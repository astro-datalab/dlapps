/*
 *  DLPHEAD -- 
 *
 *    Usage:
 *              dlphead [<opts>] [<file> | <URI>]
 *
 *  @file       dlphead.c
 *  @author     Mike Fitzpatrick
 *  @date       1/03/16
 *
 *  @brief      Compute statistics for numeric columns of a VOTable.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

//#include "fitsio.h"
//#include "ast.h"

//#include "votParse.h"                 /* keep these in order!         */
#include "dlApps.h"


/*  All tasks must support a "--return" flag if they can return a pointer
 *  to a result as part of the programmatic invocation.  As a cmdline task
 *  'reslen' and 'result' are ignored (actually, they're thrown away), but
 *  when called from an API the 'result' is a pointer to an arbitrary memory
 *  location that is passed back to the caller.
 *
 *  The result object is defined by the task and can be anything.  The
 *  '--return' flag can be defined to take an optional argument to specify
 *  which of multiple possible objects are returned (e.g. "--return=fits"
 *  "--return=votable") but the task is responsible for creating the object.
 */
static int  do_return   = 0;            /* return result?               */

/*  Global task declarations.  These should all be defined as 'static' to
 *  avoid namespace collisions.
 */
static int foo          = 0;


/*  A result buffer should be defined to point to the result object if it is
 *  created dynamically, e.g. a list of votable columns.  The task is
 *  responsible for initially allocating this pointer and then resizing as
 *  needed.
 */
#define SZ_RESBUF       8192

static char *resbuf;


/*  Task specific option declarations.  Task options are declared using the
 *  getopt_long(3) syntax.
 */
static Task  self       = {  "dlphead",  dlphead,  0,  0,  0  };

static char  *opts      = "%hno:r";
static struct option long_opts[] = {
        { "test",         2, 0,   '%'},         /* --test is std        */
        { "help",         2, 0,   'h'},         /* --help is std        */
        { "number",       2, 0,   'n'},         /* opt w/ no arg        */
        { "output",       1, 0,   'o'},         /* opt w/ required arg  */
        { "return",       2, 0,   'r'},         /* --return is std      */
        { NULL,           0, 0,    0 }
};


/*  All tasks should declare a static Usage() method to print the help
 *  text in response to a '-h' or '--help' flag.  The help text should
 *  include a usage summary, a description of options, and some examples.
 */
static void Usage (void);
static void Tests (char *input);

static void transform (void);
static void printvec (double *vec, int n, int extraline)

    
    
/**
 *  Application entry point.  All DLApps tasks MUST contain this
 *  method signature.
 */
int
dlphead (int argc, char **argv, size_t *reslen, void **result)
{
    /*  These declarations are required for the DLApps param interface.
     */
    char **pargv, optval[SZ_FNAME];

    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    AstFitsChan *fitschan;
    AstFrameSet *wcsinfo, *simpler;
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    char *header = NULL, *fname = NULL;
    int single = 0, hdupos, nkeys, verbose=0;
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */


    if (strcmp (argv[1], "-v") == 0) {
	verbose = 1;
	fname = argv[2];
    } else
	fname = argv[1];



    if (!fits_open_file(&fptr, fname, READONLY, &status)) {
      fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

      /* List only a single header if a specific extension was given */ 
      if (hdupos != 1 || strchr(argv[1], '[')) single = 1;

      for (; !status; hdupos++)  /* Main loop through each extension */ {
        printf("Header listing for HDU #%d:\n", hdupos);

	if (fits_hdr2str (fptr, 0, NULL, 0, &header, &nkeys, &status ) )
	    printf(" Error getting header\n");
	else {
	    astBegin;

	    /* Create a FitsChan and fill it with FITS header cards. */
	    fitschan = astFitsChan( NULL, NULL, "" );
	    astPutCards( fitschan, header );

	    /* Free the memory holding the concatenated header cards. */
	    free (header);
	    header = NULL;

	    /* Read WCS information from the FitsChan. */
	    wcsinfo = astRead( fitschan );
	    if (wcsinfo == NULL)
		break;

	    /* Print the WCS info. */
	    simpler = astSimplify (wcsinfo);
	    if (simpler != NULL) {
		if (verbose)
		    astShow (simpler);
	    } else
		printf ("astSimplify failed\n");

	    /* Test coordinate transformation. */
	    (void) transform (fptr, simpler);

	    astEnd;
	}

        if (single) break;  /* quit if only listing a single header */
        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
      }

      if (status == END_OF_FILE)
	  status = 0; /* Reset after normal error */

      fits_close_file (fptr, &status);
    }

    if (status)
	fits_report_error (stderr, status); /* print any error message */



    /*  Clean up.  Rememebr to free whatever pointers were created when
     *  parsing arguments.
     */
    dl_paramFree (argc, pargv);

    return (status);    /* status must be OK or ERR (i.e. 0 or 1)       */
}


/*
 * Transform - test an N-D coordinate transform.
 */
static void
transform (fitsfile *fptr, AstFrameSet *wcs)
{
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    double in[64], out[64], pix_zero = 1.0;
    int imdim=4, npts=2, istep=2, ostep=2;
    int bitpix, naxis, ncols, i;
    char ctype[8][20], cunit[8][20], key[20];
    char comment[8][64];
    void printvec();
    long naxes[10];

    bzero ((void *)naxes, sizeof(naxes));
    bzero ((void *)ctype, sizeof(ctype));
    bzero ((void *)cunit, sizeof(cunit));
    fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);

    imdim = naxis;  /* should be naxes! */
    for (i=1;  i <= imdim;  i++) {
      sprintf (key, "CTYPE%d", i);
      if (fits_read_key (fptr, TSTRING, key, ctype[i], comment[i], &status))
        break;
      sprintf (key, "CUNIT%d", i);
      if (fits_read_key (fptr, TSTRING, key, cunit[i], comment[i], &status))
        break;
    }

    printf ("naxes=%d, naxis=[%ld,%ld,%ld,%ld]\n", imdim,
        naxes[0], naxes[1], naxes[2], naxes[3]);
    printf ("ctypes: [%s,%s,%s,%s]   --  ",
        ctype[1], ctype[2], ctype[3], ctype[4]);
    printf ("cunits: [%s,%s,%s,%s]\n",
        cunit[1], cunit[2], cunit[3], cunit[4]);

    bzero ((void *)in, sizeof(in));
    bzero ((void *)out, sizeof(out));
    astSetI (wcs, "report", 1);

    /* Forward transform. */
    astTranN (wcs, npts, imdim, istep, in, 1, imdim, ostep, out);
    printvec (out, 16, 1);

    /* Reverse transform. */
    astTranN (wcs, npts, imdim, ostep, out, 0, imdim, istep, in);
    printvec (in, 16, 1);

    /* Forward transform, using output from above. */
    astTranN (wcs, npts, imdim, istep, in, 1, imdim, ostep, out);
    printvec (out, 16, 1);
}


/* Utility to dump a vector. */
static void
printvec (double *vec, int n, int extraline)
{
    int i;

    printf ("Output:  ");
    for (i=0;  i < n;  i++)
	printf ("%g ", vec[i]);
    printf ("\n");

    if (extraline)
	printf ("\n");
}


/** ************************************************************************
*** ***********************************************************************/

/**
 *  USAGE -- Print task help summary.
 */
static void
Usage (void)
{
    fprintf (stderr, "\n  Usage:\n\t"
        "dlphead [<opts>] <input>\n\n"
        "  where\n"
        "       -%%,--test              run unit tests\n"
        "       -h,--help               this message\n"
        "       -r,--return             return result from method\n"
        "       -d,--debug              debug output\n"
        "       -v,--verbose            verbose output\n"
	"\n"
        "       -o,--output=<file>      output file\n"
        "\n"
        "  Examples:\n\n"
        "    1)  First example\n\n"
        "           %% dlphead decam.fits\n"
        "           %% cat test.fits | dlphead\n"
        "\n"
        "    2)  Second example\n\n"
        "           %% dlphead -o hdr.txt test.fits\n"
        "\n"
    );
}


/**
 *  TESTS -- Task unit tests.
 */
static void
Tests (char *input)
{
   /*  First argument must always be the 'self' variable, the last must
    *  always be a NULL to terminate the cmd args.
    */
   dl_taskTest (self, "--help", NULL);
}

