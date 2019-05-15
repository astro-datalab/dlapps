/**
 *  FTABCNV -- Convert FITS Binary Tables to one or more CSV files.
 *
 *    Usage:
 *		ftabcnv [<otps>] [ <input> ..... ]
 *
 *  @file       ftabcnv.c
 *  @author     Mike Fitzpatrick
 *  @date       8/21/16
 *
 *  @brief       Convert FITS Binary Tables to one or more CSV files.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <arpa/inet.h>

#include "fitsio.h"
#include "dlApps.h"


#define MAX_CHUNK               100000
#define MAX_COLS                1024

#define	SZ_RESBUF	        8192
#define SZ_COLNAME              32
#define SZ_EXTNAME              32
#define SZ_COLVAL               1024
#define SZ_ESCBUF               1024
#define SZ_LINEBUF              10240

#define TAB_DELIMITED           0               // delimited ascii table
#define TAB_IPAC                1               // IPAC table format
#define TAB_POSTGRES            2               // SQL -- PostgreSQL
#define TAB_MYSQL               3               // SQL -- MySQL
#define TAB_SQLITE              4               // SQL -- SQLite

#define DEF_CHUNK               10000
#define DEF_ONAME               "root"

#define DEF_FORMAT              TAB_DELIMITED
#define DEF_DELIMITER           ','
#define DEF_QUOTE               '"'
#define DEF_MODE                "w+"


typedef struct {
    int       colnum;
    int       dispwidth;
    int       type;
    int       nelem;
    int       ndim;
    int       nrows;
    int       ncols;
    long      width;
    long      repeat;
    char      colname[SZ_COLNAME];
    char      coltype[SZ_COLNAME];
    char      colunits[SZ_COLNAME];
} Col, *ColPtr;

Col inColumns[MAX_COLS];
Col outColumns[MAX_COLS];

int numInCols 		= 0;		// number of input columns
int numOutCols 		= 0;		// number of output columns


char    esc_buf[SZ_ESCBUF];             // escaped value buffer
char   *obuf, *optr;                    // output buffer pointers
long    olen = 0;                       // output buffer length


char   *extname         = NULL;         // extension name
char   *iname           = NULL;         // input file name
char   *oname           = NULL;         // output file name
char   *basename        = NULL;         // base output file name
char   *rows            = NULL;         // row selection string
char   *expr            = NULL;         // selection expression string
char   *tablename       = NULL;         // database table name

char    delimiter       = DEF_DELIMITER;// default to CSV
char    quote_char      = DEF_QUOTE;    // string quote character
char   *omode           = DEF_MODE;     // output file mode

int     format          = DEF_FORMAT;   // default output format
int     mach_swap       = 0;            // is machine swapped relative to FITS?
int     do_binary       = 0;            // do binary SQL output
int     do_quote        = 1;            // quote ascii values?
int     do_escape       = 0;            // escape strings for quotes?
int     do_strip        = 0;            // strip leading/trailing whitespace?
int     do_drop         = 0;            // drop db table before creating new one
int     do_create       = 0;            // create new db table
int     do_truncate     = 0;            // truncate db table before load
int     nfiles          = 0;            // number of input files
int     noop            = 0;            // no-op ??

int     concat          = 0;            // concat input file to single output?
int     explode         = 0;            // explode arrays to new columns?
int     extnum          = -1;           // extension number
int     header          = 1;            // prepend column headers
int     number          = 0;            // number rows ?
int     chunk_size      = DEF_CHUNK;    // processing chunk size

int     debug           = 0;            // debug flag
int     verbose         = 0;            // verbose output flag

char   *pgcopy_hdr      = "PGCOPY\n\377\r\n\0\0\0\0\0";
int     len_pgcopy_hdr  = 15;

size_t  sz_char         = sizeof (char);
size_t  sz_short        = sizeof (short);
size_t  sz_int          = sizeof (int);
size_t  sz_long         = sizeof (long);
size_t  sz_longlong     = sizeof (long long);
size_t  sz_float        = sizeof (float);
size_t  sz_double       = sizeof (double);



/*  Task specific option declarations.  Task options are declared using the
 *  getopt_long(3) syntax.
static Task  self       = {  "ftabcnv",  ftabcnv,  0,  0,  0  };
 */

static char  *opts 	= "hdvnc:e:E:i:o:r:s:t:BCXHQNS012345:678";
static struct option long_opts[] = {
    { "help",         no_argument,          NULL,   'h'},
    { "debug",        no_argument,          NULL,   'd'},
    { "verbose",      no_argument,          NULL,   'v'},
    { "noop",         no_argument,          NULL,   'n'},

    { "chunk",        required_argument,    NULL,   'c'},
    { "extnum",       required_argument,    NULL,   'e'},
    { "extname",      required_argument,    NULL,   'E'},
    { "input",        required_argument,    NULL,   'i'},
    { "output",       required_argument,    NULL,   'o'},
    { "rowrange",     required_argument,    NULL,   'r'},
    { "select",       required_argument,    NULL,   's'},
    { "table",        required_argument,    NULL,   't'},

    { "binary",       no_argument,          NULL,   'B'},
    { "concat",       no_argument,          NULL,   'C'},
    { "explode",      no_argument,          NULL,   'X'},
    { "noheader",     no_argument,          NULL,   'H'},
    { "noquote",      no_argument,          NULL,   'Q'},
    { "number",       no_argument,          NULL,   'N'},
    { "singlequote",  no_argument,          NULL,   'S'},

    { "asv",          no_argument,          NULL,   '0'},
    { "bsv",          no_argument,          NULL,   '1'},
    { "csv",          no_argument,          NULL,   '2'},
    { "tsv",          no_argument,          NULL,   '3'},
    { "ipac",         no_argument,          NULL,   '4'},

    { "sql",          required_argument,    NULL,   '5'},
    { "drop",         no_argument,          NULL,   '6'},
    { "create",       no_argument,          NULL,   '7'},
    { "truncate",     no_argument,          NULL,   '8'},

    { NULL,           0,                    0,       0 }
};


/*  All tasks should declare a static Usage() method to print the help 
 *  text in response to a '-h' or '--help' flag.  The help text should 
 *  include a usage summary, a description of options, and some examples.
 */
static void Usage (void);

static void dl_escapeCSV (char* in);
static void dl_quote (char* in);
static void dl_ftabcnv (char *iname, char *oname, int filenum, int nfiles);
static void dl_printHdr (fitsfile *fptr, int firstcol, int lastcol, FILE *ofd);
static void dl_printIPACTypes (char *tablename, fitsfile *fptr, int firstcol,
                                int lastcol, FILE *ofd);
static void dl_createSQLTable (char *tablename, fitsfile *fptr, int firstcol,
                                int lastcol, FILE *ofd);
static void dl_printSQLHdr (char *tablename, fitsfile *fptr, int firstcol,
                                int lastcol, FILE *ofd);
static void dl_getColInfo (fitsfile *fptr, int firstcol, int lastcol);
static int  dl_validateColInfo (fitsfile *fptr, int firstcol, int lastcol);
static void dl_getOutputCols (fitsfile *fptr, int firstcol, int lastcol);

static unsigned char *dl_printCol (unsigned char *dp, ColPtr col, char end_ch);
static unsigned char *dl_printString (unsigned char *dp, ColPtr col);
static unsigned char *dl_printLogical (unsigned char *dp, ColPtr col);
static unsigned char *dl_printByte (unsigned char *dp, ColPtr col);
static unsigned char *dl_printShort (unsigned char *dp, ColPtr col);
static unsigned char *dl_printInt (unsigned char *dp, ColPtr col);
static unsigned char *dl_printLong (unsigned char *dp, ColPtr col);
static unsigned char *dl_printFloat (unsigned char *dp, ColPtr col);
static unsigned char *dl_printDouble (unsigned char *dp, ColPtr col);

extern int dl_atoi (char *v);
extern int dl_isFITS (char *v);
extern int dl_isGZip (char *v);
extern void dl_error (int exit_code, char *error_message, char *tag);

static char *dl_colType (ColPtr col);
static char *dl_IPACType (ColPtr col);
static char *dl_SQLType (ColPtr col);
static char *dl_makeTableName (char *fname);

extern void bswap2 (char *a, char *b, int nbytes);
extern void bswap4 (char *a, int aoff, char *b, int boff, int nbytes);
extern void bswap8 (char *a, int aoff, char *b, int boff, int nbytes);
extern char *sstrip (char *s);
extern int is_swapped (void);



/**
 *  Application entry point.  All DLApps tasks MUST contain this 
 *  method signature.
 */
int
main (int argc, char **argv)
{
    char **pargv, optval[SZ_FNAME];
    char **iflist = NULL, **ifstart = NULL;
    char  *iname = NULL, *oname = NULL;
    int    i, ch = 0, status = 0, pos = 0;


    /*  Initialize local task values.
     */
    ifstart = calloc (argc, sizeof (char *));
    iflist = ifstart;


    /*  Parse the argument list.  The use of dl_paramInit() is required to
     *  rewrite the argv[] strings in a way dl_paramNext() can be used to
     *  parse them.  The programmatic interface allows "param=value" to
     *  be passed in, but the getopt_long() interface requires these to
     *  be written as "--param=value" so they are not confused with 
     *  positional parameters (i.e. any param w/out a leading '-').
     */
    prog_name = argv[0];
    pargv = dl_paramInit (argc, argv, opts, long_opts);
    memset (optval, 0, SZ_FNAME);
    while ((ch = dl_paramNext(opts,long_opts,argc,pargv,optval,&pos)) != 0) {
        if (ch > 0) {
	    /*  If the 'ch' value is > 0 we are parsing a single letter
	     *  flag as defined in the 'opts string.
	     */
	    switch (ch) {
	    case 'h':  Usage ();			return (OK);
	    case 'd':  debug++;                         break;  // --debug
	    case 'v':  verbose++;                       break;  // --verbose
	    case 'n':  noop++;                          break;  // --noop

	    case 'c':  chunk_size = dl_atoi (optval);	break;  // --chunk_size
	    case 'e':  extnum = dl_atoi (optval);	break;  // --extnum
	    case 'E':  extname = strdup (optval);	break;  // --extname
	    case 'r':  rows = strdup (optval);		break;  // --rows
	    case 's':  expr = strdup (optval);	        break;  // --select
	    case 't':  tablename = strdup (optval);	break;  // --table

	    case 'B':  do_binary++;			break;  // --binary
	    case 'C':  concat++;			break;  // --concat
	    case 'X':  explode++;			break;  // --explode
	    case 'H':  header = 0;			break;  // --noheader
	    case 'Q':  do_quote = 0;			break;  // --noquote
	    case 'N':  do_strip = 0;			break;  // --strip
	    case 'S':  quote_char = '\'';		break;  // --quote

	    case 'i':  iname = strdup (optval);		break;  // --input
	    case 'o':  oname = strdup (optval);		break;  // --output

	    case '0':  delimiter = ' ';			break;  // ASV
	    case '1':  delimiter = '|';			break;  // BSV
	    case '2':  delimiter = ',';			break;  // CSV
	    case '3':  delimiter = '\t';		break;  // TSV
	    case '4':  delimiter = '|'; 
                       format = TAB_IPAC; 	        break;

	    case '5':  if (optval[0] == 'm')           // MySQL ouptut
                           format = TAB_MYSQL;
                       else if (optval[0] == 's')      // MySQL ouptut
                           format = TAB_SQLITE;
                       else                            // default to Postgres
                           format = TAB_POSTGRES;
                       break;
	    case '6':  do_drop++, do_create++;          break;  // --drop
	    case '7':  do_create++;                     break;  // --create
	    case '8':  do_truncate++;                   break;  // --truncate

	    default:
		fprintf (stderr, "Invalid option '%s'\n", optval);
		return (ERR);
	    }

	} else {
	    /*  All non-opt arguments are input files to process.
	     */
	    *iflist++ = strdup (optval);
            nfiles++;
	}

        memset (optval, 0, SZ_FNAME);
    }
    *iflist = NULL;


    if (debug) {
        fprintf (stderr, "do_create=%d  do_drop=%d  do_truncate=%d\n",
            do_create, do_drop, do_truncate);
        fprintf (stderr, "extnum=%d  extname='%s' rows='%s' expr='%s'\n",
            extnum, extname, rows, expr);
        fprintf (stderr, "table = '%s'\n", (tablename ? tablename : "<none>"));
        for (i=0; i < nfiles; i++)
             fprintf (stderr, "in[%d] = '%s'\n", i, ifstart[i]);
        if (noop)
            return (0);
    }

    /*  Sanity checks.  Tasks should validate input and accept stdin/stdout
     *  where it makes sense.
     */
    if (iname && *ifstart == NULL)
        *ifstart = *iflist = iname;

    if (*ifstart == NULL) {
        dl_error (2, "no input files specified", NULL);
        return (ERR);
    }

    if (format == TAB_MYSQL || format == TAB_SQLITE) {
        fprintf (stderr, "Error: database format '%s' not yet supported\n",
            (format == TAB_MYSQL ? "mysql" : "sqlite"));
        return (ERR);
    }

    if (do_binary) {
        if (format == TAB_MYSQL || format == TAB_SQLITE) {
            fprintf (stderr, "Error: binary output not supported for '%s'\n",
                (format == TAB_MYSQL ? "mysql" : "sqlite"));
            return (ERR);
        }
    }


    /*  Generate the output file lists if needed.
     */
    if (nfiles == 1 || concat) {
        /*  If we have 1 input file, output may be to stdout or to the named
         *  file only.
         */
        if (oname && strcmp (oname, "-") == 0)
            free (oname), oname = NULL;
        if (oname == NULL) 
            oname = strdup ("stdout");
    } else {
        if (oname)
            /*  For multiple files, the output arg specifies a root filename.
             *  We append the file number and a ".csv" extension.
             */
            basename = oname;
        else
            /*  If we don't specify an output name, use the input filename
             *  and replace the extension.
             */
            basename = NULL;
    }


    /* Compute and output the image metadata. */
    if (*ifstart == NULL)
        dl_error (2, "no input source specified", NULL);

    else {
        char ofname[SZ_PATH], ifname[SZ_PATH];
        int  ndigits = (int) log10 (nfiles) + 1;

        for (iflist=ifstart, i=0; *iflist; iflist++, i++) {

            memset (ifname, 0, SZ_PATH);
            memset (ofname, 0, SZ_PATH);

            strcpy (ifname, *iflist);

            if (basename) {
                if (concat) {
                    if (i == 0)
                        sprintf (ofname, "%s.csv", basename);
                    else if (i > 0)
                        sprintf (ofname, "%s%*d.csv", basename, ndigits, i);
                } else 
                    sprintf (ofname, "%s%*d.csv", basename, ndigits, i);

            } else if (oname) {
                strcpy (ofname, oname);

            } else {
                char *ip = (ifname + strlen(ifname) - 1);

                do {
                    *ip-- = '\0';
                } while (*ip != '.' && ip > ifname);
                *ip = '\0';
                sprintf (ofname, "%s.csv", ifname);
                strcpy (ifname, *iflist);
            }
            omode = ((concat && i > 0) ? "a+" : "w+");

            if (debug)
                fprintf (stderr, "ifname='%s'  ofname='%s'\n", ifname, ofname);

            /*  Do the conversion.
             */
            if (dl_isFITS (ifname) || dl_isGZip (ifname)) {
                if (verbose)
                    fprintf (stderr, "Processing file: %s\n", ifname);

                if (!noop)
                    dl_ftabcnv (ifname, ofname, i, nfiles);

            } else
                fprintf (stderr, 
                    "Error: Skipping non-FITS file '%s'.\n", ifname);
        }
    }

    if (status)
        fits_report_error (stderr, status);     // print any error message


    /*  Clean up.  Rememebr to free whatever pointers were created when
     *  parsing arguments.
     */
    for (iflist=ifstart; *iflist; iflist++)     // free the file list
        if (*iflist) 
            free ((void *) *iflist), *iflist = NULL;

    if (rows) free (rows);
    if (expr) free (expr);
    if (iname) free (iname);
    if (oname) free (oname);
    if (extname) free (extname);
    if (tablename) free (tablename);

    dl_paramFree (argc, pargv);

    return (status);	/* status must be OK or ERR (i.e. 0 or 1)     	*/
}


/**
 *  DL_FTABCNV -- Convert a FITS file to a database, i.e. actual SQL code or
 *  some ascii 'database' table like a CSV.
 */
static void
dl_ftabcnv (char *iname, char *oname, int filenum, int nfiles)
{
    fitsfile *fptr = (fitsfile *) NULL;
    int   status = 0;
    long  jj, nrows;
    int   hdunum, hdutype, ncols, ii, i, j;
    int   firstcol = 1, lastcol = 0, firstrow = 1;
    int   nelem, chunk = MAX_CHUNK;

    FILE  *ofd = (FILE *) NULL;
    ColPtr col = (ColPtr) NULL;

    unsigned char *data = NULL, *dp = NULL;
    long   naxis1, rowsize = 0, nbytes = 0, firstchar = 1, totrows = 0;


    mach_swap = is_swapped ();

    if (!fits_open_file (&fptr, iname, READONLY, &status)) {
        if ( fits_get_hdu_num (fptr, &hdunum) == 1 )
            /*  This is the primary array;  try to move to the first extension
             *  and see if it is a table.
             */
            fits_movabs_hdu (fptr, 2, &hdutype, &status);
         else 
            fits_get_hdu_type (fptr, &hdutype, &status); /* Get the HDU type */

        if (hdutype == IMAGE_HDU) 
            printf ("Error: this program only converts tables, not images\n");

        else {
            fits_get_num_rows (fptr, &nrows, &status);
            fits_get_num_cols (fptr, &ncols, &status);

            lastcol = ncols;
            chunk = (chunk > nrows ? nrows : chunk);
            nelem = chunk;
            

            /*  Get the optimal I/O row size.
             */
            fits_get_rowsize (fptr, &rowsize, &status);
            nelem = rowsize;

            /*  Open the output file.
             */
            if (strcasecmp (oname, "stdout") == 0 || oname[0] == '-')
                ofd = stdout;
            else {
                if ((ofd = fopen (oname, omode)) == (FILE *) NULL)
                    dl_error (3, "Error opening output file '%s'\n", oname);
            }

            /*  Print column names as column headers when writing a new file,
             *  skip if we're appending output.
             *
             *  FIXME -- Need to add a check that new file matches columns
             *           when we have multi-file input.
             */
            fits_read_key (fptr, TLONG, "NAXIS1", &naxis1, NULL, &status);
            if (filenum == 0) {
		dl_getColInfo (fptr, firstcol, lastcol);

                if (!tablename) 
                    tablename = dl_makeTableName (iname);

                if (format == TAB_DELIMITED)
                    dl_printHdr (fptr, firstcol, lastcol, ofd);
                else if (format == TAB_IPAC)
                    dl_printIPACTypes (iname, fptr, firstcol, lastcol, ofd);
                else {
                    // This is some sort of SQL output.
                    if (do_drop)
                        fprintf (ofd, "DROP TABLE IF EXISTS %s;\n", tablename);
                    if (do_create)
                        dl_createSQLTable (tablename, fptr, firstcol, lastcol,
                            ofd);
                    if (do_truncate)
                        fprintf (ofd, "TRUNCATE TABLE %s;\n", tablename);

                    dl_printSQLHdr (tablename, fptr, firstcol, lastcol, ofd);
                    delimiter = '\t';
                    do_quote = 0;
                }
            } else {
                // Make sure this file has the same columns.
                if (dl_validateColInfo (fptr, firstcol, lastcol)) {
                    fprintf (stderr, "Skipping unmatching table '%s'\n", 
                        oname);
                    return;
                }
            }

            /*  Allocate the I/O buffer.
             */                
            nbytes = nelem * naxis1;
            if (debug)
                fprintf (stderr, "nelem=%d  naxis1=%ld  nbytes=%ld\n",
                    nelem, naxis1, nbytes);

            data = (unsigned char *) calloc (1, nbytes * 8);
            obuf = (char *) calloc (1, nbytes * 8);

            /*  Loop over the rows in the table in optimal chunk sizes.
             */
            for (jj=firstrow; jj < nrows; jj += nelem) {
                if ( (jj + nelem) >= nrows)
                    nelem = (nrows - jj + 1);

                /*  Read a chunk of data from the file.
                 */
                nbytes = nelem * naxis1;
                fits_read_tblbytes (fptr, firstrow, firstchar, nbytes,
                    data, &status);
                if (status) {                   /* print any error message */
                    fits_report_error (stderr, status);
                    break;
                }

                /* Process the chunk by parsing the binary data and printing
                 * out according to column type.
                 */
                dp = data;
                optr = obuf;
                olen = 0;
                for (j=firstrow; j <= nelem; j++) {
                    if (format == TAB_POSTGRES && do_binary) {
                        unsigned short val = 0;
                        val = (explode ? htons ((short)numOutCols) : 
                                         htons (ncols));
                        memcpy (optr, &val, sz_short);
                        optr += sz_short;
                        olen += sz_short;
                    }

                    for (i=firstcol; i <= ncols; i++) {
                        dp = dl_printCol (dp, &inColumns[i], 
                            (i < ncols ? delimiter : '\n'));
                    }
                }
                write (fileno(ofd), obuf, olen);
                fflush (ofd);

                /*  Advance the offset counters in the file.
                 */
                firstchar += nbytes;
                totrows += nelem;
            }

            if (concat && filenum == (nfiles - 1)) {
                if (format == TAB_POSTGRES) {
                    short  eof = -1;
                    optr = obuf, olen = 0;
                    memset (optr, 0, nbytes);

                    if (do_binary) {
                        memcpy (optr, &eof, sz_short);
                        olen += sz_short;
                    } else {
                        memcpy (optr, "\\.\n", 3);
                        olen += 3;
                    }
                }
                write (fileno(ofd), obuf, olen);
                fflush (ofd);
            }


            /*  Free the column structures and data pointers.
             */
            if (data) free ((char *) data);
            if (obuf) free ((void *) obuf);
            for (ii = firstcol; ii <= lastcol; ii++) {
                col = (ColPtr) &inColumns[ii];
            }

            /*  Close the output file.
             */
            if (ofd != stdout)
                fclose (ofd);
        }
    }
    fits_close_file (fptr, &status);

    if (status)                                 /* print any error message */
        fits_report_error (stderr, status);
}


/**
 *  DL_GETCOLINFO -- Get information about the columns in teh table.
 */
static void
dl_getColInfo (fitsfile *fptr, int firstcol, int lastcol)
{
    register int i;
    char  keyword[FLEN_KEYWORD], dims[FLEN_KEYWORD];
    ColPtr icol = (ColPtr) NULL;
    int   status = 0;


    /* Gather information about the input columns.
     */
    for (i = firstcol, numInCols = 0; i <= lastcol; i++, numInCols++) {
        icol = (ColPtr) &inColumns[i];
        memset (icol, 0, sizeof(Col));

        status = 0;                             // reset CFITSIO status
        memset (keyword, 0, FLEN_KEYWORD);
        fits_make_keyn ("TTYPE", i, keyword, &status);
        fits_read_key (fptr, TSTRING, keyword, icol->colname, NULL, &status);
        fits_get_coltype (fptr, i, &icol->type, &icol->repeat,
            &icol->width, &status);
        fits_get_col_display_width (fptr, i, &icol->dispwidth, &status);
        if (icol->type == TSTRING && do_quote) 
            icol->dispwidth+= 2;
        icol->colnum = i;

        icol->ndim = 1;		                // default dimensions
        icol->nrows = 1;
        icol->ncols = icol->repeat;

        if (icol->repeat > 1 && icol->type != TSTRING && explode) {
            memset (keyword, 0, FLEN_KEYWORD);
            fits_make_keyn ("TDIM", i, keyword, &status);
            fits_read_key (fptr, TSTRING, keyword, dims, NULL, &status);
            if (status == 0)    // Dimension string usually means a 2-D array
                icol->ndim = sscanf (dims, "(%d,%d)", 
                    &icol->nrows, &icol->ncols);
            status = 0;
        }
    }

    if (debug) {
        fprintf (stderr, "Input Columns [%d]:\n", numInCols);
        for (i=1; i <= numInCols; i++) {
            icol = (ColPtr) &inColumns[i];
            fprintf (stderr, "  %d  '%s'  rep=%ld nr=%d nc=%d\n", icol->colnum, 
                icol->colname, icol->repeat, icol->nrows, icol->ncols);
        }
    }

    /* Now expand to create the output column information so we don't need 
     * to repeat this for each table type.  If we're exploding columns
     * compute all the output columns names, otherwise simply copy the input
     * names.
     */
    dl_getOutputCols (fptr, firstcol, lastcol);
}


/**
 *  DL_VALIDATECOLINFO -- Validate that this file has the same column
 *  information.
 */
static int
dl_validateColInfo (fitsfile *fptr, int firstcol, int lastcol)
{
    register int i;
    char   keyword[FLEN_KEYWORD], dims[FLEN_KEYWORD];
    Col    newColumns[MAX_COLS], *col = (ColPtr) NULL, *icol = (ColPtr) NULL;
    int    numCols, status = 0;


return (0);
    /* Gather information about the input columns.
     */
    for (i = firstcol, numCols = 0; i <= lastcol; i++, numCols++) {
        col = (ColPtr) &newColumns[i];
        memset (col, 0, sizeof(Col));
        col->colnum = i;

        status = 0;                             // reset CFITSIO status
        memset (keyword, 0, FLEN_KEYWORD);
        fits_make_keyn ("TTYPE", i, keyword, &status);
        fits_read_key (fptr, TSTRING, keyword, col->colname, NULL, &status);
        fits_get_coltype (fptr, i, &col->type, &col->repeat,
            &col->width, &status);

        col->ndim = 1;				// default dimensions
        col->nrows = 1;
        col->ncols = col->repeat;

        if (col->repeat > 1 && col->type != TSTRING && explode) {
            memset (keyword, 0, FLEN_KEYWORD);
            fits_make_keyn ("TDIM", i, keyword, &status);
            fits_read_key (fptr, TSTRING, keyword, dims, NULL, &status);
            if (status == 0)
                col->ndim = sscanf (dims, "(%d,%d)", &col->nrows, &col->ncols);
            status = 0;
        }
    }

    if (debug) {
        fprintf (stderr, "Table Columns [%d]:\n", numCols);
        for (i=1; i <= numCols; i++) {
            col = (ColPtr) &newColumns[i];
            fprintf (stderr, "  %d  '%s'  rep=%ld nr=%d nc=%d\n", col->colnum, 
                col->colname, col->repeat, col->nrows, col->ncols);
        }
    }

    /*  Check column names, dimensionality, and type for equality.
     */
    for (i=firstcol; i <= lastcol; i++) {
        col = (ColPtr) &newColumns[i];
        icol = (ColPtr) &inColumns[i];

        if (strcmp (col->colname, icol->colname))
            return (1);
        if (col->type != icol->type ||
            col->ndim != icol->ndim ||
            col->nrows != icol->nrows ||
            col->ncols != icol->ncols ||
            col->repeat != icol->repeat)
                return (1);
    }

    return (0);                                         // No error
}


/**
 *  DL_GETOUTPUTCOLS -- Get the output column information.
 */
static void
dl_getOutputCols (fitsfile *fptr, int firstcol, int lastcol)
{
    register int i, j, ii, jj;
    ColPtr icol = (ColPtr) NULL;
    ColPtr ocol = (ColPtr) NULL;


    if (explode) {
        jj = firstcol;
        for (ii = firstcol; ii <= lastcol; ii++) {
            icol = (ColPtr) &inColumns[ii];

            if (icol->repeat > 1 && icol->type != TSTRING) {
                if (icol->ndim > 1) {                           // 2-D array
                    for (i=1; i <= icol->nrows; i++) {
                        for (j=1; j <= icol->ncols; j++)  {
                            ocol = (ColPtr) &outColumns[jj++];
                            memset (ocol->colname, 0, SZ_COLNAME);
                            sprintf (ocol->colname, "%s_%d_%d",
                                icol->colname, i, j);
                            strcpy (ocol->coltype, dl_colType (icol));
                            ocol->dispwidth = icol->dispwidth;
                        }
                    }
                } else {                                        // 1-D array
                    for (i=1; i <= icol->repeat; i++) {
                        ocol = (ColPtr) &outColumns[jj++];
                        memset (ocol->colname, 0, SZ_COLNAME);
                        sprintf (ocol->colname, "%s_%d", icol->colname, i);
                        strcpy (ocol->coltype, dl_colType (icol));
                        ocol->dispwidth = icol->dispwidth;
                    }
                }
            } else {
                ocol = (ColPtr) &outColumns[jj++];
                memset (ocol->colname, 0, SZ_COLNAME);
                strcpy (ocol->colname, icol->colname);
                strcpy (ocol->coltype, dl_colType (icol));
                ocol->dispwidth = icol->dispwidth;
            }
        }
        numOutCols = jj - 1;

    } else {
        for (i = firstcol, numOutCols = 0; i <= lastcol; i++, numOutCols++) {
            icol = (ColPtr) &inColumns[i];
            ocol = (ColPtr) &outColumns[i];
            memcpy (ocol, icol, sizeof(Col));

            if (format == TAB_IPAC)
                strcpy (ocol->coltype, dl_IPACType (icol));
            else if (format == TAB_POSTGRES)
                strcpy (ocol->coltype, dl_SQLType (icol));
        }
        numOutCols = numInCols;
    }

    if (debug) {
        fprintf (stderr, "Output Columns [%d]:\n", numOutCols);
        for (i=1; i <= numOutCols; i++) {
            ocol = (ColPtr) &outColumns[i];
            fprintf (stderr, "  %d  %-24s  '%s'\n", ocol->colnum, 
                ocol->colname, ocol->coltype);
        }
    }
}


/**
 *  DL_PRINTHDR -- Print the CSV column headers.
 */
static void
dl_printHdr (fitsfile *fptr, int firstcol, int lastcol, FILE *ofd)
{
    register int i;
    ColPtr col = (ColPtr) NULL;


    if (*omode == 'a')
        return;

    if (format == TAB_IPAC)
        fprintf (ofd, "|");

    for (i=1; i <= numOutCols; i++) {           // print column types
        col = (ColPtr) &outColumns[i];

        if (format == TAB_IPAC)                 // FIXME
            fprintf (ofd, "%-*s", col->dispwidth, col->colname);
        else
            fprintf (ofd, "%-s", col->colname);
        if (i < numOutCols)
            fprintf (ofd, "%c", delimiter);
    }

    if (format == TAB_IPAC)
        fprintf (ofd, "|");

    if (format == TAB_IPAC || format == TAB_DELIMITED)
        fprintf (ofd, "\n");
    fflush (ofd);
}


/**
 *  DL_CREATESQLTABLE -- Print the SQL CREATE command.
 */
static void
dl_createSQLTable (char *tablename, fitsfile *fptr, int firstcol, int lastcol,
                    FILE *ofd)
{
    register int  i;
    ColPtr col = (ColPtr) NULL;


    fprintf (ofd, "CREATE TABLE IF NOT EXISTS %s (\n", tablename);

    for (i=1; i <= numOutCols; i++) {             // print column types
        col = (ColPtr) &outColumns[i];
        fprintf (ofd, "    %s\t%s", col->colname, col->coltype);
        if (i < numOutCols)
            fprintf (ofd, ",\n");
    }

    fprintf (ofd, "\n);\n\n");
    fflush (ofd);
}


/**
 *  DL_PRINTSQLHDR -- Print the SQL COPY headers.
 */
static void
dl_printSQLHdr (char *tablename, fitsfile *fptr, int firstcol, int lastcol,
                    FILE *ofd)
{
    if (do_binary) {
        int hdr_extn = 0;
        char copy_buf[160];

        memset (copy_buf, 0, 160);
        sprintf (copy_buf, "COPY %s FROM stdin WITH BINARY;\n", tablename);

        if (!noop)
            write (fileno(ofd), copy_buf, strlen(copy_buf));   // header string

        write (fileno(ofd), pgcopy_hdr, len_pgcopy_hdr);   // header string
        write (fileno(ofd), &hdr_extn, sz_int);            // header extn length

    } else {
        fprintf (ofd, "COPY %s (", tablename);
        dl_printHdr (fptr, firstcol, lastcol, ofd);
        fprintf (ofd, ") from stdin;\n");
        fflush (ofd);
    }
}


/**
 *  DL_PRINTIPACTYPES -- Print the IPAC column type headers.
 */
static void
dl_printIPACTypes (char *tablename, fitsfile *fptr, int firstcol, int lastcol,
                    FILE *ofd)
{
    register int i;
    ColPtr col = (ColPtr) NULL;


    if (*omode == 'a' || format != TAB_IPAC)
        return;

    dl_printHdr (fptr, firstcol, lastcol, ofd);         // print column names

    fprintf (ofd, "|");
    for (i=1; i <= numOutCols; i++) {           // print column types
        col = (ColPtr) &outColumns[i];
        fprintf (ofd, "%-*s|", col->dispwidth, col->coltype);
    }

    fprintf (ofd, "\n");
    fflush (ofd);
}


/**
 * DL_COLTYPE -- Get the type string for the column.
 */
static char *
dl_colType (ColPtr col)
{
    if (format == TAB_POSTGRES)
        return dl_SQLType (col);
    else //if (format == TAB_IPAC)
        return dl_IPACType (col);
}


/**
 * DL_SQLTYPE -- Get the SQL type string for the column.
 */
static char *
dl_SQLType (ColPtr col)
{
    static char *type = NULL;
    static char tbuf[64];


    switch (col->type) {
    case TBIT:                                          // NYI
    case TCOMPLEX:
    case TDBLCOMPLEX:      
                    break;
    case TSTRING:   type = ((col->repeat > 1) ? "character varying" : "char");
                    break;
    case TLOGICAL:  type = "smallint";
                    break;

    case TBYTE:
    case TSBYTE:    type = "smallint";
		    break;
    case TSHORT:
    case TUSHORT:   type = "smallint";
		    break;
    case TINT:
    case TUINT:
    case TINT32BIT: type = "integer";
		    break;

    case TLONGLONG: type = "bigint";
                    break;
    case TFLOAT:    type = "real";
		    break;
    case TDOUBLE:   type = "double precision";
		    break;

    default:        fprintf (stderr, "Error: unsupported type %d\n", col->type);
    }

    memset (tbuf, 0, 64);
    if (col->repeat > 1 && col->type == TSTRING)
        sprintf (tbuf, "%s(%ld)", type, col->repeat);
        //strcpy (tbuf, type);
    else if (!explode && col->repeat > 1)
        sprintf (tbuf, "%s[%ld]", type, col->repeat);
    else
        strcpy (tbuf, type);

    return (tbuf);
}


/**
 * DL_IPACTYPE -- Get the IPAC-table type string for the column.
 */
static char *
dl_IPACType (ColPtr col)
{
    char *type = NULL;

    switch (col->type) {
    case TBIT:                                          // NYI
    case TCOMPLEX:                                      // NYI
    case TDBLCOMPLEX:                                   // NYI
                    break;

    case TSTRING:   type = "char";          
                    break;
    case TLOGICAL:
    case TBYTE:
    case TSBYTE:
    case TSHORT:
    case TUSHORT:
    case TINT:
    case TUINT:
    case TINT32BIT:
    case TLONGLONG: type = "int";
                    break;
    case TFLOAT:    type = "real";
                    break;
    case TDOUBLE:   type = "double";
                    break;
    default:        type = " ";
    }

    return (type);
}


/**
 *  DL_PRINTCOL -- Print the column value(s).
 *
 *  FIXME -- We don't handled unsigned or long correctly yet.
 */

#define SZ_TXTBUF               16738

static unsigned char *
dl_printCol (unsigned char *dp, ColPtr col, char end_char)
{
    if (!explode && !do_binary && col->type != TSTRING && col->repeat > 1) {
        if (format == TAB_DELIMITED) {
            *optr++ = quote_char, *optr++ = '(';
            olen += 2;
        } else {
            *optr++ = '{';
            olen += 1;
        }
    }

                    
    if (format == TAB_IPAC && col->colnum == 1)
        *optr++ = '|', olen++;

    switch (col->type) {
    case TBIT:                          // TFORM='X'    bit
        fprintf (stderr, "Error: Unsupported column type, col[%s] = %d\n", 
            col->colname, col->type);
        break;
    case TCOMPLEX:                      // TFORM='C'    complex
    case TDBLCOMPLEX:                   // TFORM='M'    double complex
        fprintf (stderr, "Error: Unsupported column type, col[%s] = %d\n", 
            col->colname, col->type);
        break;

    case TSTRING:                       // TFORM='A'    8-bit character
        dp = dl_printString (dp, col);
        break;

    case TLOGICAL:                      // TFORM='L'    8-bit logical (boolean)
        dp = dl_printLogical (dp, col);
        break;

    case TBYTE:                         // TFORM='B'    1 unsigned byte
    case TSBYTE:                        // TFORM='S'    8-bit signed byte
        dp = dl_printByte (dp, col);
        break;

    case TSHORT:                        // TFORM='I'    16-bit integer
    case TUSHORT:                       // TFORM='U'    unsigned 16-bit integer
        dp = dl_printShort (dp, col);
        break;

    case TINT:                          // TFORM='J'    32-bit integer
    case TUINT:                         // TFORM='V'    unsigned 32-bit integer
    case TINT32BIT:                     // TFORM='K'    signed 32-bit integer
        dp = dl_printInt (dp, col);
        break;

    case TLONGLONG:                     // TFORM='K'    64-bit integer
        dp = dl_printLong (dp, col);
        break;

    case TFLOAT:                        // TFORM='E'    single precision float
        dp = dl_printFloat (dp, col);
        break;

    case TDOUBLE:                       // TFORM='D'    double precision float
        dp = dl_printDouble (dp, col);
        break;

    default:
        fprintf (stderr, "Error: Unknown column type, col[%s] = %d\n", 
            col->colname, col->type);
        break;
    }

    if (!explode && !do_binary && col->type != TSTRING && col->repeat > 1) {
        if (format == TAB_DELIMITED) {
            *optr++ = quote_char, *optr++ = ')';
            olen += 2;
        } else {
            *optr++ = '}';
            olen += 1;
        }
    }

    if (format == TAB_IPAC && end_char == '\n')
        *optr++ = '|', olen++;

    if (!do_binary && TAB_POSTGRES)
        *optr++ = end_char, olen++;     // append the comma or newline


    return (dp);
}


/**
 *  DL_PRINTSTRING -- Print the column as a string value.
 */
static unsigned char *
dl_printString (unsigned char *dp, ColPtr col)
{
    char  buf[SZ_TXTBUF], *bp;
    int   len = 0;


    if (do_binary) {
        unsigned int val = 0;
        val = htonl (col->repeat);

        //memcpy (optr, &val, sz_int);            optr += sz_int;
        //memcpy (optr, dp, col->repeat);         optr += col->repeat;

        memset (buf, 0, SZ_TXTBUF);
        memcpy (buf, dp, col->repeat);
        len = strlen ((bp = sstrip(buf)));
        val = htonl (len);

        memcpy (optr, &val, sz_int);            optr += sz_int;
        memcpy (optr, bp, len);                 optr += len;

        olen += sz_int + len;

    } else {
        memset (buf, 0, SZ_TXTBUF);
        memcpy (buf, dp, col->repeat);
        if (do_escape) {
            dl_escapeCSV ((do_strip ? sstrip (buf) : buf));
            memcpy (optr, esc_buf, (len = strlen (esc_buf)));
        } else {
            if (do_quote) {
                dl_quote ((do_strip ? sstrip (buf) : buf));
                memcpy (optr, esc_buf, (len = strlen (esc_buf)));
            } else {
                bp =  sstrip(buf);
                memcpy (optr, bp, (len = strlen(bp)));
            }
        }

        olen += len;
        optr += len;
    }
    dp += col->repeat;

    return (dp);
}


/**
 *  DL_PRINTLOGICAL -- Print the column as logical values.
 */
static unsigned char *
dl_printLogical (unsigned char *dp, ColPtr col)
{
    char ch;
    char  valbuf[8 * col->repeat];
    int   i, j, len = 0;


    if (do_binary) {
        unsigned int sz_val = 0;
        unsigned short lval = 0;
        if (explode) {
            len = sz_short;
            sz_val = htonl(sz_short);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    ch = (char) *dp++;
                    lval = ((tolower((int)ch) == 't') ? htons(1) : 0);
                    memcpy (optr, &lval, sz_short);     optr += sz_short;
                    olen += sz_int + len;
                }
            }
        } else {
            len = col->repeat * sz_short;
            sz_val = htonl(col->repeat * sz_short);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;

            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    ch = (char) *dp++;
                    lval = ((tolower((int)ch) == 't') ? htons(1) : 0);
                    memcpy (optr, &lval, sz_short);     optr += sz_short;
                    olen += sz_short + len;
                }
            }
            olen += sz_int + sz_val;
        }

    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 8 * col->repeat);

                ch = (char) *dp++;
                if (format == TAB_IPAC)
                    sprintf (valbuf, "%d", ((tolower((int)ch) == 't') ? 1 : 0));
                else
                    sprintf (valbuf, "%*d", col->dispwidth, 
                        ((tolower((int)ch) == 't') ? 1 : 0));
                memcpy (optr, valbuf, (len = strlen (valbuf)));
                olen += len;
                optr += len;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}


/**
 *  DL_PRINTBYTE -- Print the column as byte values.
 */
static unsigned char *
dl_printByte (unsigned char *dp, ColPtr col)
{
    char ch;
    unsigned char uch;
    char  valbuf[4 * col->repeat];
    int   i, j, len = 0;


    if (do_binary) {
        unsigned int sz_val = 0;
        short sval = 0;
        if (explode) {
            len = sz_short;
            sz_val = htonl(sz_short);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    sval = (short) *dp++;
                    memcpy (optr, &sval, sz_short);     optr += sz_short;
                    olen += sz_int + len;
                }
            }
        } else {
            len = col->repeat * sz_short;
            sz_val = htonl(col->repeat * sz_short);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;
            olen += sz_int;

            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    sval = (short) *dp++;
                    memcpy (optr, &sval, sz_short);     optr += sz_short;
                    olen += sz_short;
                }
            }
        }
    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 4 * col->repeat);
                if (col->type == TBYTE) {
                    uch = (unsigned char) *dp++;
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*d", col->dispwidth, uch);
                    else
                        sprintf (valbuf, "%d", uch);
                } else {
                    ch = (char) *dp++;
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*d", col->dispwidth, ch);
                    else
                        sprintf (valbuf, "%d", ch);
                }
                memcpy (optr, valbuf, (len = strlen (valbuf)));
                olen += len;
                optr += len;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}


/**
 *  DL_PRINTSHORT -- Print the column as short integer values.
 */
static unsigned char *
dl_printShort (unsigned char *dp, ColPtr col)
{
    short sval = 0.0;
    unsigned short usval = 0.0;
    char  valbuf[8 * col->repeat];
    int   i, j, len = 0;


    if (mach_swap && !do_binary)
        bswap2 ((char *)dp, (char *)dp, sz_short * col->repeat);

    if (do_binary) {
        unsigned int sz_val = 0;
        if (explode) {
            len = sz_short;
            sz_val = htonl(sz_short);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    memcpy (optr, dp, sz_short);        optr += len;
                    olen += sz_int + len;
                    dp += sz_short;
                }
            }
        } else {
            len = col->repeat * sz_short;
            sz_val = htonl(col->repeat * sz_short);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;
            memcpy (optr, dp, sz_short * col->repeat);  optr += len;
            olen += sz_int + len;
            dp += col->repeat * sz_short;
        }

    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 8 * col->repeat);
                if (col->type == TUSHORT) {
                    memcpy (&usval, dp, sz_short);
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*d", col->dispwidth, usval);
                    else
                        sprintf (valbuf, "%d", usval);
                } else {
                    memcpy (&sval, dp, sz_short);
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*d", col->dispwidth, sval);
                    else
                        sprintf (valbuf, "%d", sval);
                }
                memcpy (optr, valbuf, (len = strlen (valbuf)));
                olen += len;
                optr += len;
                dp += sz_short;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}


/**
 *  DL_PRINTINT -- Print the column as integer values.
 */
static unsigned char *
dl_printInt (unsigned char *dp, ColPtr col)
{
    int   ival = 0.0;
    unsigned int uival = 0.0;
    char  valbuf[64 * col->repeat];
    int   i, j, len = 0;


    if (mach_swap && !do_binary)
        bswap4 ((char *)dp, 1, (char *)dp, 1, sz_int * col->repeat);

    if (do_binary) {
        unsigned int sz_val = 0;
        if (explode) {
            len = sz_int;
            sz_val = htonl(sz_int);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    memcpy (optr, dp, sz_int);          optr += len;
                    olen += sz_int + len;
                    dp += sz_int;
                }
            }
        } else {
            len = col->repeat * sz_int;
            sz_val = htonl(col->repeat * sz_int);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;
            memcpy (optr, dp, sz_int * col->repeat);    optr += len;
            olen += sz_int + len;
            dp += col->repeat * sz_int;
        }

    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 64 * col->repeat);
                if (col->type == TUINT) {
                    memcpy (&uival, dp, sz_int);
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*d", col->dispwidth, uival);
                    else
                        sprintf (valbuf, "%d", uival);
                } else {
                    memcpy (&ival, dp, sz_int);
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*d", col->dispwidth, ival);
                    else
                        sprintf (valbuf, "%d", ival);
                }
                memcpy (optr, valbuf, (len = strlen (valbuf)));
                olen += len;
                optr += len;
                dp += sz_int;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}


/**
 *  DL_PRINTLONG -- Print the column as long integer values.
 */
static unsigned char *
dl_printLong (unsigned char *dp, ColPtr col)
{
    long  lval = 0.0;
    char  valbuf[64 * col->repeat];
    int   i, j, len = 0;


    if (mach_swap && !do_binary)
        //bswap8 ((char *)dp, 1, (char *)dp, 1, sizeof(long) * col->repeat);
        // FIXME -- We're in trouble if we comes across a 64-bit int column
        bswap4 ((char *)dp, 1, (char *)dp, 1, sz_long * col->repeat);

    if (do_binary) {
        unsigned int sz_val = 0;
        if (explode) {
            len = sz_long;
            sz_val = htonl(sz_long);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    memcpy (optr, dp, sz_long);         optr += len;
                    olen += sz_int + len;
                    dp += sz_long;
                }
            }
        } else {
            len = col->repeat * sz_long;
            sz_val = htonl(col->repeat * sz_long);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;
            memcpy (optr, dp, sz_long * col->repeat);   optr += len;
            olen += sz_int + len;
            dp += col->repeat * sz_long;
        }

    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 64 * col->repeat);
                memcpy (&lval, dp, sz_long);
                if (format == TAB_IPAC)
                    sprintf (valbuf, "%*ld", col->dispwidth, lval);
                else
                    sprintf (valbuf, "%ld", lval);
                memcpy (optr, valbuf, (len = strlen (valbuf)));
                olen += len;
                optr += len;
                dp += sz_long;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}


/**
 *  DL_PRINTFLOAT -- Print the column as floating-point values.
 */
static unsigned char *
dl_printFloat (unsigned char *dp, ColPtr col)
{
    float rval = 0.0;
    char  valbuf[64 * col->repeat];
    int   i, j, sign = 1, len = 0;


    if (mach_swap && !do_binary)
        bswap4 ((char *)dp, 1, (char *)dp, 1, sz_float * col->repeat);

    if (do_binary) {
        unsigned int sz_val = 0;
        if (explode) {
            len = sz_float;
            sz_val = htonl(sz_float);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    memcpy (optr, dp, sz_float);        optr += len;
                    olen += sz_int + len;
                    dp += sz_float;
                }
            }
        } else {
            len = col->repeat * sz_float;
            sz_val = htonl(col->repeat * sz_float);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;
            memcpy (optr, dp, sz_float * col->repeat);  optr += len;
            olen += sz_int + len;
            dp += col->repeat * sz_float;
        }

    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 64 * col->repeat);
                memcpy (&rval, dp, sz_float);
                if (! isnan (rval) ) {
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*f", col->dispwidth, (double) rval);
                    else
                        sprintf (valbuf, "%f", (double) rval);
                    memcpy (optr, valbuf, (len = strlen (valbuf)));
                    olen += len;
                    optr += len;

                //   FIXME -- Need to add option to filter NaN/Inf
                } else if ((sign = isinf (rval)) ) {
                    if (sign > 0)
                        memcpy (optr, "Inf", (len = strlen ("Inf")));
                    else
                        memcpy (optr, "-Inf", (len = strlen ("-Inf")));
                    olen += len;
                    optr += len;

                } else if (isnan (rval) ) {
                    memcpy (optr, "NaN", (len = strlen ("NaN")));
                    olen += len;
                    optr += len;
                }
                dp += sz_float;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}


/**
 *  DL_PRINTDOUBLE -- Print the column as double-precision values.
 */
static unsigned char *
dl_printDouble (unsigned char *dp, ColPtr col)
{
    double dval = 0.0;
    char  valbuf[64 * col->repeat];
    int   i, j, len = 0;


    if (mach_swap && !do_binary)
        bswap8 ((char *)dp, 1, (char *)dp, 1, sz_double * col->repeat);

    if (do_binary) {
        unsigned int sz_val = 0;
        if (explode) {
            len = sz_double;
            sz_val = htonl(sz_double);
            for (i=1; i <= col->nrows; i++) {
                for (j=1; j <= col->ncols; j++) {
                    memcpy (optr, &sz_val, sz_int);     optr += sz_int;
                    memcpy (optr, dp, sz_double);       optr += len;
                    olen += sz_int + len;
                    dp += sz_double;
                }
            }
        } else {
            len = col->repeat * sz_double;
            sz_val = htonl(col->repeat * sz_double);
            memcpy (optr, &sz_val, sz_int);             optr += sz_int;
            memcpy (optr, dp, sz_double * col->repeat); optr += len;
            olen += sz_int + len;
            dp += col->repeat * sz_double;
        }

    } else {
        for (i=1; i <= col->nrows; i++) {
            for (j=1; j <= col->ncols; j++) {
                memset (valbuf, 0, 64 * col->repeat);
                memcpy (&dval, dp, sizeof(double));
                if (! isnan (dval) ) {
                    if (format == TAB_IPAC)
                        sprintf (valbuf, "%*lf", col->dispwidth, dval);
                    else
                        sprintf (valbuf, "%.16lf", dval);
                    memcpy (optr, valbuf, (len = strlen (valbuf)));
                    olen += len;
                    optr += len;
                } else {
                    memcpy (optr, "NaN", (len = strlen ("NaN")));
                    olen += len;
                    optr += len;
                }
                dp += sz_double;
                if (col->repeat > 1 && j < col->ncols)
                    *optr++ = delimiter,  olen++;
            }
        }
    }

    return (dp);
}



/***********************************************************/
/****************** LOCAL UTILITY METHODS ******************/
/***********************************************************/


/**
 *  DL_MAKETABLENAME -- Name a table name from the input file name.
 */
static char *
dl_makeTableName (char *fname)
{
    char *ip, *tp, *np;

    ip = strdup (fname);                // copy the input name
    tp = strchr (ip, (int)'.');         // locate and kill the '.'
    *tp = '\0';

    for (np=ip; *np; np++) {
        if (*np == '-')
            *np = '_';
    }

    return (ip);                        // return the start of the filename
}


/**
 *  DL_ESCAPECSV -- Escape quotes for CSV printing.
 */
static void
dl_escapeCSV (char* in)
{
    int   in_len = 0;
    char *ip = in, *op = esc_buf;

    memset (esc_buf, 0, SZ_ESCBUF);
    if (in)
        in_len = strlen (in);

    *op++ = quote_char;
    for ( ; *ip; ip++) {
        *op++ = *ip;
        if (*ip == quote_char)
            *op++ = quote_char;
    }
    *op++ = quote_char;
}


/**
 *  DL_QUOTE -- Set quotes for CSV printing.
 */
static void
dl_quote (char* in)
{
    int    in_len = 0;
    char  *op = esc_buf;

    memset (esc_buf, 0, SZ_ESCBUF);
    in_len = (in ? strlen (in) : 0);

    *op++ = quote_char;
    if (in_len) {
        memcpy (op, in, in_len);
        op += in_len;
    }
    *op++ = quote_char;
}


/**
 *  USAGE -- Print task help summary.
 */
static void
Usage (void)
{
    fprintf (stderr, "\n  Usage:\n\t"
        "ftabcnv [<opts>] [ <input> ... ]\n\n"
        "  where\n"
        "      -h,--help                this message\n"
        "      -d,--debug               set debug flag\n"
        "      -v,--verbose             set verbose output flag\n"
	"\n"
        "      -e,--extnum=<N>          process table extension number <N>\n"
        "      -E,--extname=<name>      process table extension name <name>\n"
        "      -i,--input=<file>        set input filename\n"
        "      -o,--output=<file>       set output filename\n"
	"\n"
        "      -c,--chunk=<nrows>       set size of ingest chunk size\n"
        "      -C,--concat              concatenate all input files\n"
        "      -X,--explode             explode array cols to new columns\n"
        "      -H,--noheader            suppress CSV column header\n"
        "      -N,--number              number rows\n"
        "      -r,--rowrange=<range>    convert only row range\n"
        "      -s,--select=<expr>       convert rows selected by expression\n"
        "      -S--singlequote          use single quotes for strings\n"
	"\n"
 	"  Examples:\n\n"
	"    1)  Convert table to CSV on stdout\n\n"
	"	    %% ftabcnv test.fits\n"
	"\n"
        "  Example input file syntax: \n"
        "    ftabcnv tab.fits[GTI]               - list the GTI extension\n"
        "    ftabcnv tab.fits[1][#row < 101]     - list first 100 rows\n"
        "    ftabcnv tab.fits[1][col X;Y]        - list X and Y cols only\n"
        "    ftabcnv tab.fits[1][col -PI]        - list all but the PI col\n"
        "    ftabcnv tab.fits[1][col -PI][#row < 101]  - combined case\n"
        "\n"
	"\n"
    );
}
