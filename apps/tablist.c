/**
 *
 **/

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "fitsio.h"


#define MAX_CHUNK               10000
#define MAX_COLS                1024
#define SZ_COLNAME              32
#define SZ_COLVAL               1024
#define SZ_ESCBUF               1024
#define SZ_LINEBUF              10240


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

int     mach_swap       = 0;            // is machine swapped relative to FITS?
int     do_quote        = 0;            // quote ascii values?
int     do_escape       = 0;            // escape strings for quotes?
int     do_strip        = 0;            // strip leading/trailing whitespace?


extern int is_swapped (void);
extern void bswap2 (char *a, char *b, int nbytes);
extern void bswap4 (char *a, int aoff, char *b, int boff, int nbytes);
extern void bswap8 (char *a, int aoff, char *b, int boff, int nbytes);
extern char *sstrip (char *s);





/*  ESCAPECSV -- Escape quotes for CSV printing.
 */
void
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

unsigned char *
printColVal (unsigned char *dp, ColPtr col, char end_char)
{
    char ch, buf[SZ_TXTBUF];    //unsigned char uch;
    short sval = 0;             //unsigned short usval = 0;
    int ival = 0;               //unsigned int uival = 0;
    long lval = 0;              //unsigned long ulval = 0;
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



/**
 *
 */
int
main(int argc, char *argv[])
{
    fitsfile *fptr;
    char  keyword[FLEN_KEYWORD];
    long  jj, nrows;
    int   hdunum, hdutype, ncols, ii, i, j;
    int   firstcol = 1, lastcol = 0, firstrow = 1;
    int   nelem, chunk = MAX_CHUNK;
    int   status = 0;   // CFITSIO status value MUST be initialized to zero!
                    
    ColPtr col = (ColPtr) NULL;

    unsigned char *data = NULL, *dp = NULL;
    long   naxis1, rowsize = 0, nbytes = 0, firstchar = 1, totrows = 0;


    if (argc != 2) {
        printf("Usage:  tablist filename[ext][col filter][row filter] \n");
        printf("\n");
        printf("List the contents of a FITS table \n");
        printf("\n");
        printf("Examples: \n");
        printf("  tablist tab.fits[GTI]           - list the GTI extension\n");
        printf("  tablist tab.fits[1][#row < 101] - list first 100 rows\n");
        printf("  tablist tab.fits[1][col X;Y]    - list X and Y cols only\n");
        printf("  tablist tab.fits[1][col -PI]    - list all but the PI col\n");
        printf("  tablist tab.fits[1][col -PI][#row < 101]  - combined case\n");
        printf("\n");
        printf("Display formats can be modified with the TDISPn keywords.\n");

        return(0);
    }


    mach_swap = is_swapped ();

    if (!fits_open_file (&fptr, argv[1], READONLY, &status)) {
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


            /*  Print column names as column headers
             */
            fits_read_key (fptr, TINT, "NAXIS1", &naxis1, NULL, &status);
            for (ii = firstcol; ii <= lastcol; ii++) {
                col = (ColPtr) &Columns[ii];

                fits_make_keyn ("TTYPE", ii, keyword, &status);
                fits_read_key (fptr, TSTRING, keyword, col->colname, NULL,
                    &status);
                fits_get_eqcoltype (fptr, ii, &col->type, &col->repeat,
                    &col->width, &status);
                col->colnum = ii;
                printf ("%-s%c", col->colname, (ii < lastcol ? ',' : '\n'));
            }
            fflush (stdout);


            /*  Allocate the I/O buffer.
             */                
            nbytes = nelem * naxis1;
            data = (unsigned char *) calloc (1, nbytes);
            obuf = (char *) calloc (1, nbytes * 10);

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
                    for (i=firstcol; i <= ncols; i++) {
                        //col = (ColPtr) &Columns[i];
                        dp = printColVal (dp, &Columns[i], 
                            (i < ncols ? ',' : '\n'));
                    }
                }
                write (fileno(stdout), obuf, olen);
                fflush (stdout);

                /*  Advance the offset counters in the file.
                 */
                firstchar += nbytes;
                totrows += nelem;
            }
            free ((char *) data);               // free the data pointer
            free ((void *) obuf);

            /*  Free the column structures.
             */
            for (ii = firstcol; ii <= lastcol; ii++) {
                col = (ColPtr) &Columns[ii];
            }
        }
    }
    fits_close_file(fptr, &status);

    if (status)                                 /* print any error message */
        fits_report_error(stderr, status);
    return(status);
}
