/**
 *  FITS2CSV -- Convert FITS Binary Tables to one or more CSV files.
 *
 *    Usage:
 *		fits2csv [<otps>] [ <input> ..... ]
 *
 *  @file       fits2csv.c
 *  @author     Mike Fitzpatrick
 *  @date       8/13/16
 *
 *  @brief       Convert FITS Binary Tables to one or more CSV files.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

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

#define DEF_CHUNK               10000
#define DEF_ONAME               "root"


typedef struct {
    int       colnum;
    int       dispwidth;
    int       type;
    int       nelem;
    long      width;
    long      repeat;
    char      colname[SZ_COLNAME];
} Col, *ColPtr;

Col Columns[MAX_COLS];

char    esc_buf[SZ_ESCBUF];             // escaped value buffer
char   *obuf, *optr;                    // output buffer pointers
long    olen = 0;                       // output buffer length


char   *extname         = NULL;         // extension name
char   *iname           = NULL;         // input file name
char   *oname           = NULL;         // output file name
char   *basename        = NULL;         // base output file name
char   *rows            = NULL;         // row selection string
char   *expr            = NULL;         // selection expression string

char    quote_char      = '"';          // string quote character

int     mach_swap       = 0;            // is machine swapped relative to FITS?
int     do_quote        = 0;            // quote ascii values?
int     do_escape       = 0;            // escape strings for quotes?
int     do_strip        = 0;            // strip leading/trailing whitespace?
int     nfiles          = 0;            // number of input files

int     concat          = 0;            // concat input file to single output?
int     explode         = 0;            // explode arrays to new columns?
int     extnum          = -1;           // extension number
int     header          = 1;            // prepend column headers
int     number          = 0;            // number rows ?
int     chunk_size      = DEF_CHUNK;    // processing chunk size

int     debug           = 0;            // debug flag
int     verbose         = 0;            // verbose output flag




/*  Task specific option declarations.  Task options are declared using the
 *  getopt_long(3) syntax.
static Task  self       = {  "fits2csv",  fits2csv,  0,  0,  0  };
 */

static char  *opts 	= "hdvc:e:E:i:o:r:s:CXHNS";
static struct option long_opts[] = {
    { "help",         no_argument,          NULL,   'h'},
    { "debug",        no_argument,          NULL,   'd'},
    { "verbose",      no_argument,          NULL,   'v'},

    { "chunk",        required_argument,    NULL,   'c'},
    { "extnum",       required_argument,    NULL,   'e'},
    { "extname",      required_argument,    NULL,   'E'},
    { "input",        required_argument,    NULL,   'i'},
    { "output",       required_argument,    NULL,   'o'},
    { "rowrange",     required_argument,    NULL,   'r'},
    { "select",       required_argument,    NULL,   's'},

    { "concat",       no_argument,          NULL,   'C'},
    { "explode",      no_argument,          NULL,   'X'},
    { "noheader",     no_argument,          NULL,   'H'},
    { "number",       no_argument,          NULL,   'N'},
    { "singlequote",  no_argument,          NULL,   'S'},

    { NULL,           0,                    0,       0 }
};


/*  All tasks should declare a static Usage() method to print the help 
 *  text in response to a '-h' or '--help' flag.  The help text should 
 *  include a usage summary, a description of options, and some examples.
 */
static void Usage (void);

static void escapeCSV (char* in);
static unsigned char *printColVal (unsigned char *dp, ColPtr col,
                        char end_char);

extern int dl_atoi (char *v);
extern int dl_isFITS (char *v);
extern void dl_error (int exit_code, char *error_message, char *tag);

extern int is_swapped (void);
extern void bswap2 (char *a, char *b, int nbytes);
extern void bswap4 (char *a, int aoff, char *b, int boff, int nbytes);
extern void bswap8 (char *a, int aoff, char *b, int boff, int nbytes);
extern char *sstrip (char *s);




/**
 *  Application entry point.  All DLApps tasks MUST contain this 
 *  method signature.
 */
int
main (int argc, char **argv)
{
    fitsfile *fptr = (fitsfile *) NULL;
    char **pargv, optval[SZ_FNAME];
    char **iflist = NULL, **ifstart = NULL;
    char **oflist = NULL, **ofstart = NULL;
    char  *iname, *oname;
    int    i, ch = 0, status = 0, pos = 0;

    FILE  *ifd = (FILE *) NULL, *ofd = (FILE *) NULL;


    /*  Initialize local task values.
     */
    iname  = NULL;
    oname  = NULL;

    iflist = ifstart = calloc (argc, sizeof (char *));
    oflist = ofstart = calloc (argc, sizeof (char *));


    /*  Parse the argument list.  The use of dl_paramInit() is required to
     *  rewrite the argv[] strings in a way dl_paramNext() can be used to
     *  parse them.  The programmatic interface allows "param=value" to
     *  be passed in, but the getopt_long() interface requires these to
     *  be written as "--param=value" so they are not confused with 
     *  positional parameters (i.e. any param w/out a leading '-').
     */
    prog_name = argv[0];
    pargv = dl_paramInit (argc, argv, opts, long_opts);
    while ((ch = dl_paramNext(opts,long_opts,argc,pargv,optval,&pos)) != 0) {
        if (ch > 0) {
	    /*  If the 'ch' value is > 0 we are parsing a single letter
	     *  flag as defined in the 'opts string.
	     */
	    switch (ch) {
	    case 'h':  Usage ();			return (OK);
	    case 'd':  debug++;
	    case 'v':  verbose++;

	    case 'c':  chunk_size = dl_atoi (optval);	break;
	    case 'e':  extnum = dl_atoi (optval);	break;
	    case 'E':  extname = strdup (optval);	break;
	    case 'r':  rows = strdup (optval);		break;
	    case 's':  expr = strdup (optval);	        break;

	    case 'C':  concat++;			break;
	    case 'X':  explode++;			break;
	    case 'H':  header = 0;			break;
	    case 'N':  number++;			break;
	    case 'S':  quote_char = '\'';		break;

	    case 'i':  iname = strdup (optval);		break;
	    case 'o':  oname = strdup (optval);		break;

	    default:
		fprintf (stderr, "Invalid option '%s'\n", optval);
		return (1);
	    }

	} else if (ch == PARG_ERR) {
	    return (ERR);

	} else {
	    /*  All non-opt arguments are input files to process.
	     */
	    *iflist++ = strdup (optval);
            nfiles++;
	}
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


    /*  Generate the output file lists if needed.
     */
    if (nfiles == 1 || concat) {
        /*  If we have 1 input file, output may be to stdout or to the named
         *  file only.
         */
        if (strcmp (oname, "-") == 0)
            free (oname), oname = NULL;
        if (oname == NULL) 
            ofd = stdout;
        else if ((ofd = fopen (oname, "w+")) == (FILE *) NULL)
            dl_error (3, "Error opening output file '%s'\n", oname);
    } else {
        if (oname) {
            /*  For multiple files, the output arg specifies a root filename.
             *  We append the file number and a ".csv" extension.
             */
            basename = oname;
        } else {
            /*  If we don't specify an output name, use the input filename
             *  and replace the extension.
             */
            basename = NULL;
        }
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
                sprintf (ofname, "%s%*d.csv", basename, ndigits, i);
            } else {
                char *ip = (ifname + strlen(ifname) - 1);

                do {
                    *ip-- = '\0';
                } while ( *ip != '.' && ip > ifname);
                *ip = '\0';
                sprintf (ofname, "%s.csv", ifname);
                strcpy (ifname, *iflist);
            }

            if ((ofd = fopen (ofname, "w+")) == (FILE *) NULL)
                dl_error (3, "Error opening output file '%s'\n", oname);
printf ("ifname = '%s'   ofname = '%s'\n", ifname, ofname);

            if (dl_isFITS (ifname)) {
                if (debug)
                    fprintf (stderr, "Processing FITS = '%s'\n", ifname);
                //dl_procFITS ( *iflist, NULL );
            } else
                fprintf (stderr, 
                    "Error: Skipping non-FITS file '%s'.\n", ifname);

            if (ofd != stdout)
                fclose (ofd);
        }
    }

    if (status)
        fits_report_error (stderr, status);     // print any error message


    /*  Clean up.  Rememebr to free whatever pointers were created when
     *  parsing arguments.
     */
    for (iflist=ifstart; *iflist; iflist++)     // free the file list
        if (*iflist) 
            free ((void *) *iflist);
    for (oflist=ofstart; *oflist; oflist++)     // free the file list
        if (*oflist) 
            free ((void *) *oflist);

    if (rows) free (rows);
    if (expr) free (expr);
    if (iname) free (iname);
    if (oname) free (oname);
    if (extname) free (extname);

    dl_paramFree (argc, pargv);

    return (status);	/* status must be OK or ERR (i.e. 0 or 1)     	*/
}


/**
 *  USAGE -- Print task help summary.
 */
static void
Usage (void)
{
    fprintf (stderr, "\n  Usage:\n\t"
        "fits2csv [<opts>] [ <input> ... ]\n\n"
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
	"	    %% fits2csv test.fits\n"
	"\n"
        "  Example input file syntax: \n"
        "    fits2csv tab.fits[GTI]               - list the GTI extension\n"
        "    fits2csv tab.fits[1][#row < 101]     - list first 100 rows\n"
        "    fits2csv tab.fits[1][col X;Y]        - list X and Y cols only\n"
        "    fits2csv tab.fits[1][col -PI]        - list all but the PI col\n"
        "    fits2csv tab.fits[1][col -PI][#row < 101]  - combined case\n"
        "\n"
	"\n"
    );
}



/* *********** Local Methods ********** */


/*  ESCAPECSV -- Escape quotes for CSV printing.
 */
static void
escapeCSV (char* in)
{
    int   in_len = 0;
    char *ip = in, *op = esc_buf;

    memset (esc_buf, 0, SZ_ESCBUF);
    if (in)
        in_len = strlen (in);

    *op++ = '"';
    for ( ; *ip; ip++) {
        *op++ = *ip;
        if (*ip == '"')
            *op++ = '"';
    }
    *op++ = '"';
}


/*  PRINTCOLVAL -- Print the column value.
 */

#define SZ_TXTBUF               16738

static unsigned char *
printColVal (unsigned char *dp, ColPtr col, char end_char)
{
    char ch, buf[SZ_TXTBUF];            // unsigned char uch;
    short sval = 0;                     // unsigned short usval = 0;
    int ival = 0;                       // unsigned int uival = 0;
    long lval = 0;                      // unsigned long ulval = 0;
    float rval = 0.0;
    double dval = 0.0;

    char  valbuf[64];
    int   len = 0;


    memset (valbuf, 0, 64);
    switch (col->type) {
    case TBIT:                          // TFORM='X'    bit
        break;          // NYI
    case TBYTE:                         // TFORM='B'    1 unsigned byte
    case TSBYTE:                        // TFORM='S'    8-bit signed byte
        ch = *dp;
        break;          // NYI
    case TLOGICAL:                      // TFORM='L'    logical (boolean)
        if (mach_swap)
            bswap4 ((char *)dp, 1, (char *)&ival, 1, sizeof(int));
        else
            memcpy (&ival, dp, sizeof(int));
        break;          // NYI
    case TCOMPLEX:                      // TFORM='C'    complex
    case TDBLCOMPLEX:                   // TFORM='M'    double complex
        break;          // NYI

    case TSTRING:                       // TFORM='A'    8-bit character
        memset (buf, 0, SZ_TXTBUF);
        memcpy (buf, dp, col->repeat);
        if (do_escape) {
            escapeCSV ((do_strip ? sstrip (buf) : buf));
            memcpy (optr, esc_buf, (len = strlen (esc_buf)));
        } else
            memcpy (optr, (do_strip ? sstrip (buf) : buf), (len = col->repeat));

        olen += len;
        optr += len;
        dp += col->repeat;
        break;

    case TSHORT:                        // TFORM='I'    16-bit integer
    case TUSHORT:                       // TFORM='U'    unsigned 16-bit integer
        if (mach_swap)
            bswap2 ((char *)dp, (char *)&sval, sizeof(short));
        else
            memcpy (&sval, dp, sizeof(short));

        sprintf (valbuf, "%d", sval);
        memcpy (optr, valbuf, (len = strlen (valbuf)));
        olen += len;
        optr += len;
        dp += sizeof (short);
        break;

    case TINT:                          // TFORM='J'    32-bit integer
    case TUINT:                         // TFORM='V'    unsigned 32-bit integer
        if (mach_swap)
            bswap4 ((char *)dp, 1, (char *)&ival, 1, sizeof(int));
        else
            memcpy (&ival, dp, sizeof(int));

        sprintf (valbuf, "%d", ival);
        memcpy (optr, valbuf, (len = strlen (valbuf)));
        olen += len;
        optr += len;
        dp += sizeof (int);
        break;

    case TLONG:                         // TFORM='K'    64-bit integer
    case TULONG:                        // TFORM='X'    unsigned 64-bit integer
    case TLONGLONG:                     // TFORM='K'    64-bit integer
        if (mach_swap)
            bswap8 ((char *)dp, 1, (char *)&lval, 1, sizeof(long));
        else
            memcpy (&lval, dp, sizeof(long));

        sprintf (valbuf, "%ld", lval);
        memcpy (optr, valbuf, (len = strlen (valbuf)));
        olen += len;
        optr += len;
        dp += sizeof (long);
        break;

    case TFLOAT:                        // TFORM='E'    single precision float
        if (mach_swap)
            bswap4 ((char *)dp, 1, (char *)&rval, 1, sizeof(float));
        else
            memcpy (&rval, dp, sizeof(float));

        if (! isnan (rval) ) {
            sprintf (valbuf, "%f", rval);
            memcpy (optr, valbuf, (len = strlen (valbuf)));
            olen += len;
            optr += len;
        }
        dp += sizeof (float);
        break;

    case TDOUBLE:                       // TFORM='D'    double precision float
        if (mach_swap)
            bswap8 ((char *)dp, 1, (char *)&dval, 1, sizeof(double));
        else
            memcpy (&dval, dp, sizeof(double));

        if (! isnan (rval) ) {
            sprintf (valbuf, "%lf", dval);
            memcpy (optr, valbuf, (len = strlen (valbuf)));
            olen += len;
            optr += len;
        }
        dp += sizeof (double);
        break;

    default:
        break;
    }

    *optr++ = end_char; olen++;         // append the comma or newline

    return (dp);
}
