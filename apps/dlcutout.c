/**
 *  DLCUTOUT -- Create a cutout from an image.
 *
 *    Usage:
 *              dlcutout [<opts>] [ <file> | @<file>] ....
 *
 *  @file       dlcutout.c
 *  @author     Mike Fitzpatrick
 *  @date       2323/16
 *
 *  @brief      Create a cutout from an image.
 */

#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/errno.h>
#ifdef HAVE_LIBPQ
#include <libpq-fe.h>		// for postgres connection
#endif
#include "fitsio.h"
#include "ast.h"
#include "dlApps.h"
#include "dbConnect.h"


/* Axis map. Maps logical axes in the axis map to physical
 * image axes.
 */
#define	AX_SP1		0	/* First spatial axis 			*/
#define	AX_SP2		1	/* Second spatial axis 			*/
#define	AX_EM		2	/* Spectral (EM) axis 			*/
#define	AX_TIME		3	/* Time axis 				*/
#define	AX_POL		4	/* Polarization axis 			*/

/* Command line options. */
#define	RA		21	/* RA of region center, ICRS 		*/
#define	DEC		22	/* DEC of region center, ICRS 		*/
#define	WIDTH		23	/* Angular width of region, deg 	*/
#define	HEIGHT		24	/* Angular height of region, deg 	*/
#define	WAVELO		25	/* Start of spectral band, meters 	*/
#define	WAVEHI		26	/* End of spectral band, meters 	*/
#define	TIMELO		27	/* Start of time band, MJD 		*/
#define	TIMEHI		28	/* End of time band, MJD 		*/
#define	POLSTATES	29	/* List of polarization states 		*/
#define	SECTION		30	/* User-supplied image section 		*/
#define	PREVIEW		31	/* Generate PNG preview			*/


/*
#define	TEST_QUERY \
 "fileRef=c4d_130211_024118_ooi_i_v2.fits.fz&extn=12&POS=150.62,2.68&SIZE=0.005"
#define	TEST_QUERY "fileRef=c4d_130211_024118_ooi_i_v2.fits.fz&extn=12"
 */
#define	TEST_QUERY \
 "fileRef=c4d_130211_024118_ooi_i_v2.fits.fz&extn=12&POS=150.62,2.68&SIZE=0.005&PREVIEW=T"


/* Global params shared by all functions. 
 */
static int   debug 	= 0;
//const  char *prog_name  = NULL;
static char *pubdid 	= NULL;
static char *dir 	= NULL;
static char *polstates 	= NULL;
static char *imsection 	= NULL;
static char *collection = NULL;

static int spat_filter  = 0;		// are we doing spatial filtering?
static int spec_filter  = 0;		// are we doing spectral filtering?
static int time_filter 	= 0;		// are we doing temporal filtering?
static int pol_filter 	= 0;		// are we doing polarization filtering?
static int filter_term 	= 0;		// are we doing any kind of filtering?
static int pixel_term 	= 0;		// are we doing a pixel cutout?
static int extn		= 0;		// requested image extension
static int preview	= 0;		// generate a preview?
static int cutout	= 1;		// generate a cutout or use whole CCD?
static int isCGI	= 0;		// are we running as a CGI script?

static double wavelo 	= 0.0;		// low side of wavelength range
static double wavehi 	= 0.0;		// high side of wavelength range
static double timelo 	= 0.0;		// low side of time range
static double timehi 	= 0.0;		// high side of time range
static double ra 	= 361.0;	// RA of cutout (degrees)
static double dec 	= 91.0;		// Dec of cutout (degrees)
static double width 	= -1.0;		// Width of cutout (degrees)
static double height 	= -1.0;		// Height of cutout (degrees)

char  *phu_header 	= (char *) NULL;
char  *ehu_header 	= (char *) NULL;


/*  Task method declarations.
 */
static AstFrameSet *read_header (char *imagefile, 
	int *naxes, long *naxis, int *bitpix, int *axmap, 
	char ctype[][KEYLEN], char cunit[][KEYLEN],
	char comment[][COMLEN], int maxaxes);

static void  dl_cgiResponse (char *fname, char *type);
static char *dl_computeMetadata (char *imagefile);
static void  dl_readMDFile (char *mdfile, FILE *out);
static char *dl_extractImage  (char *mdfile, int force);
static char *dl_getKeyword (char *file, char *name);

static void  dl_parseQueryString (char *qs, char **fileRef, int *extn, 
	double *ra, double *dec, double *sz1, double *sz2, int *preview,
	int *cutout);
static int   dl_copyPixels (fitsfile *in, fitsfile *out);

static void task_Usage (void);

/*  Library method declarations.
 */
extern void  dl_loadPHU (fitsfile *fptr);
extern void  dl_loadEHU (fitsfile *fptr);
extern void  dl_error (int exit_code, char *error_message, char *tag);

extern int dl_FITS2PNG (char *fits_fname, char *png_fname);

extern char *dl_mergeHeaders (char *phu_hdr, char *ehu_hdr);



/* Command line arguments.
 */
static char opts[] = "Dmxp:l:d:F";
static struct option long_opts[] = {
    { "debug",		no_argument,		NULL,	'D' },
    { "metadata",	no_argument,		NULL,	'm' },
    { "extract",	no_argument,		NULL,	'x' },
    { "pubdid",		required_argument,	NULL,	'p' },
    { "mdfile",		required_argument,	NULL,	'l' },
    { "dir",		required_argument,	NULL,	'd' },
    { "force",		no_argument,		NULL,	'F' },
    { "ra",		required_argument,	NULL,	RA  },
    { "dec",		required_argument,	NULL,	DEC },
    { "width",		required_argument,	NULL,	WIDTH },
    { "height",		required_argument,	NULL,	HEIGHT },
    { "wavelo",		required_argument,	NULL,	WAVELO },
    { "wavehi",		required_argument,	NULL,	WAVEHI },
    { "timelo",		required_argument,	NULL,	TIMELO },
    { "timehi",		required_argument,	NULL,	TIMEHI },
    { "polstates",	required_argument,	NULL,	POLSTATES },
    { "section",	required_argument,	NULL,	SECTION },
    { NULL,		0,			NULL,	0 },
};




/*
 * DLCUTOUT - Cutout regions of multi-parameter space from a FITS image.
 *
 * The initial version supports only the spatial and spectral axes,
 * although support for time and polarization would be easy to add later.
 * The Starlink AST library is used for WCS processing; CFITSIO is used
 * for FITS file access.
 */
int main (int argc, char *argv[])
{
    prog_name = argv[0];
    int metadata=0, extract=0, force=0;;
    char *mdfile = NULL, *imagefile = NULL, *qs = NULL;
    int ch;


    if ((qs = getenv ("QUERY_STRING"))) {
	dir = "/dl1/temp";                              // [MACHDEP]
	//debug++;
	extract++;
    	spat_filter++;
    	filter_term = 1;
	isCGI++;

	dl_parseQueryString (qs, &imagefile, &extn, &ra, &dec, &width, &height,
	    &preview, &cutout);

	if (debug) {
	    printf ("file = '%s'  extn=%d\n", imagefile, extn);
	    printf ("POS=(%g,%g) SIZE=(%g,%g) prev=%d\n", 
	    	ra, dec, width, height, preview);
	}

    } else if ((qs = getenv ("CUTOUT_TEST"))) {
	dir = "/dl1/temp";                              // [MACHDEP]
	//debug++;
	extract++;
    	spat_filter++;
    	filter_term = 1;
	isCGI++;

	dl_parseQueryString (TEST_QUERY, &imagefile, &extn, &ra, &dec, 
	    &width, &height, &preview, &cutout);

    } else {

        /* Process command line options. */
        while ((ch = getopt_long(argc, argv, opts, long_opts, NULL)) != -1) {
    	    char *endptr;
    
    	    switch (ch) {
    	    case 'h':  task_Usage (); 	exit (0);
    
    	    case 'D':  debug++; 		break;
    	    case 'm':  metadata++; 		break;
    	    case 'x':  extract++; 		break;
    	    case 'p':  pubdid = optarg; 	break;
    	    case 'l':  mdfile = optarg; 	break;
    	    case 'd':  dir = optarg; 		break;
    	    case 'F':  force++; 		break;
    
    	    case RA:
    	        ra = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (RA, "invalid RA value", optarg);
    	        spat_filter++;
    	        filter_term = 1;
    	        break;
    	    case DEC:
    	        dec = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (DEC, "invalid DEC value", optarg);
    	        spat_filter++;
    	        filter_term = 1;
    	        break;
    	    case WIDTH:
    	        width = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (WIDTH, "invalid WIDTH value", optarg);
    	        spat_filter++;
    	        break;
    	    case HEIGHT:
    	        height = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (HEIGHT, "invalid HEIGHT value", optarg);
    	        spat_filter++;
    	        break;
    
    	    case WAVELO:
    	        wavelo = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (WAVELO, "invalid WAVELO value", optarg);
    	        spec_filter++;
    	        filter_term = 1;
    	        break;
    	    case WAVEHI:
    	        wavehi = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (WAVEHI, "invalid WAVEHI value", optarg);
    	        spec_filter++;
    	        filter_term = 1;
    	        break;
    	    case TIMELO:
    	        timelo = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (TIMELO, "invalid TIMELO value", optarg);
    	        time_filter++;
    	        filter_term = 1;
    	        break;
    	    case TIMEHI:
    	        timehi = strtod (optarg, &endptr);
    	        if (endptr == optarg)
    		    dl_error (TIMEHI, "invalid TIMEHI value", optarg);
    	        time_filter++;
    	        filter_term = 1;
    	        break;
    	    case POLSTATES:
    	        polstates = optarg;
    	        pol_filter++;
    	        filter_term = 1;
    	        break;
    	    case SECTION:
    	        imsection = optarg;
    	        pixel_term = 1;
    	        break;
    	    case PREVIEW:
    	        preview = 1;
    	        break;
    
    	    default:
    	        dl_error (1, "unknown option", optarg);
    	    }
        }
        argc -= optind;
        argv += optind;
    
        /*  Get the pathname of the image to be processed.
         */
        imagefile = argv[0];
    }
    

    if (debug) {
	printf ("prog=%s, mdfile=%s, imagefile=%s  preview=%d\n",
	    prog_name, mdfile, imagefile, preview);
	printf ("metadata=%d, extract=%d\n", metadata, extract);
    }

    /*  If we're not extracting an extension or cutout, simply return the
     *  entire file.
     */
    if (extn == 0 && preview == 0 && (ra > 360 && dec > 90.0)) {
        if (isCGI)
	    if (!debug) dl_cgiResponse (imagefile, "fits");

	/*  FIXME -- Need to copy to output file. 
         */

        if (mdfile) free ((void *) mdfile);
        if (imagefile) free ((void *) imagefile);

        return (0);
    }


    /*  If no precomputed metadata for the virtual image is referenced
     *  (by pointing to a mdfile), then we need to analyze the input image
     *  and specfified filter parameters to determine how to process the
     *  image.  The creation metadata for the virtual image will be saved
     *  to a mdfile and may be used later to instantiate the image.
     */
    if (mdfile == NULL) {
	if (imagefile == NULL)
	    dl_error (2, "no imagefile specified", NULL);
	else if ((mdfile = dl_computeMetadata (imagefile)) == NULL)
	    dl_error (3, "computation of virtual image failed", imagefile);
    }

    /*  If requested, return metadata for the virtual image to the client.
     */
    if (metadata)
	dl_readMDFile (mdfile, stdout);

    /*  If requested, compute the cutout image and return access information
     *  to the client.  Created images are saved as files in the staging area
     *  (rather than for example streaming them back dynamically), to permit
     *  caching to optimize repeated accesses, and to facilitate creation of
     *  images via asynchronous batch processing.
     */
    if (extract) {
	char *fname = (char *) NULL;

	if ((fname = dl_extractImage (mdfile, force)) == NULL)
	    dl_error (4, "image creation failed", mdfile);

        if (isCGI) {
	    if (!debug) dl_cgiResponse (fname, (preview ? "png" : "fits"));

	    /*  Delete the temp files.
	     */
            if (access (fname, F_OK) == 0)  
                unlink (fname);		
            if (access (mdfile, F_OK) == 0) 
                unlink (mdfile);
        }

	if (fname)
	    free ((void *) fname);
    }


    if (imagefile)
	free ((void *) imagefile);
    if (mdfile)
	free ((void *) mdfile);

    return (0);
}


/**
 *  CGIRESPONSE -- Return the CGI response, i.e. return the file with the
 *  appropriate HTML header.
 */
static void
dl_cgiResponse (char *fname, char *type)
{
    int   fd;
	    

    /*  Print the response header.
     */
    printf ("Content-Type:  image/%s\n\n", type);
    fflush (stdout);

    /*  Copy the file to stdout.
     */
    if ((fd = open (fname, O_RDONLY)) < 0) {
        fprintf (stderr, "Error opening cutout file '%s': (%s)\n",
            fname, strerror(errno));
    } else {
	char  buf[40960];
	int   nread = 0;
 
    	memset (buf, 0, 40960);
    	while ((nread = read (fd, buf, 40960)) > 0) {
   	    write (fileno(stdout), buf, nread);
    	    memset (buf, 0, 40960);
    	}
	close (fd);
    }
}


/**
 *  DL_COMPUTEMETADATA -- Compute the metadata for a virtual image to be
 *  generated from the referenced static image file.  The input parameters
 *  specify two terms of the image access model, the filter (world-cutout)
 *  term, and the pixel-space image section term.  The filter term specifies
 *  the bounds in multi-parameter space of the region to be extracted from
 *  the image.  The image section term applies to the virtual image resulting
 *  from application of the filter term.  Both terms are optional; the entire
 *  image is returned if both are absent.
 *
 *  Metadata is saved to a MDFILE, and the name of this file is returned as
 *  the function argument.
 */
static char *
dl_computeMetadata (char *imagefile) 
{
    long   naxis[MAXAXES], datalen;
    int    real_naxes, axmap[MAXAXES], imdim=0, naxes, bitpix, i;
    double ra_cen, dec_cen, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4;
    int    cutout1[MAXAXES], cutout2[MAXAXES], npts=2, step=2;
    char   ctype[MAXAXES][KEYLEN]; char cunit[MAXAXES][KEYLEN];
    double in[MAXAXES*2], out[MAXAXES*2];
    double  cra1, cdec1, cra2, cdec2, cra3, cdec3, cra4, cdec4;
    char   comment[MAXAXES][COMLEN];
    int    status1=0, status2=0;
    char  *mdfile;
    AstFrameSet *wcs;
    

    astBegin;

    /*  Read the image geometry and WCS.
     */
    if ((wcs = read_header (imagefile, &naxes, naxis, &bitpix, axmap,
		    ctype, cunit, comment, MAXAXES)) == NULL) {
	dl_error (5, "cannot read metadata for image", imagefile);
    }

    /*  Compute the world coordinates of the full image by taking the forward
     *  transform of the full pixel array.  Due to the way array dimensioning
     *  is used AST wants the coordinates of the two points to be interleaved,
     *  e.g. with a stepsize of 2, pt1-val1, pt2-val1, pt1-val2, pt2-val2, etc.
     */
    memset (in, 0, (MAXAXES * 2 * sizeof(double)));
    memset (out, 0, (MAXAXES * 2 * sizeof(double)));
    for (i=0, imdim=naxes;  i < imdim;  i++) {
	in[(i*step)+0] = 1.0;
	in[(i*step)+1] = naxis[i];
    }
    if (debug) {
        printf("Fwd Bounding B4: in %g %g %g %g\n", in[0], in[1], in[2], in[3]);
        printf("Fwd Bounding B4: out %g %g %g %g\n", out[0], out[1], out[2], out[3]);
    }
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (debug) {
        printf("Fwd Bounding B4: in %g %g %g %g\n", in[0], in[1], in[2], in[3]);
        printf("Fwd Bounding B4: out %g %g %g %g\n", out[0], out[1], out[2], out[3]);
    }
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);

    /*  Apply each given constraint filter.  The image is assumed to be
     *  in approximately the same region of parameter space as the requested
     *  cutout region; this should have been ensured by the first stage of
     *  selection.  This code is hardly rigorous at this point and some
     *  assumptions are made to simplify the code.
     */
    int ra_axis   = dl_findAxis (AX_SP1, axmap, naxes);
    int dec_axis  = dl_findAxis (AX_SP2, axmap, naxes);
    int em_axis   = dl_findAxis (AX_EM, axmap, naxes);
    int time_axis = dl_findAxis (AX_TIME, axmap, naxes);
    int pol_axis  = dl_findAxis (AX_POL, axmap, naxes);

    char *em_ctype = ctype[em_axis];
    if (em_ctype[0] == '\0')
	em_ctype = NULL;
    char *em_unit = cunit[em_axis];
    if (em_unit[0] == '\0')
	em_unit = NULL;

    /*  Disable filtering on an axis if we don't have enough WCS information.
     */
    if (ra_axis < 0 || dec_axis < 0)
	spat_filter = 0;
    if (em_axis < 0)
	spec_filter = 0;
    if (time_axis < 0)
	time_filter = 0;
    if (pol_axis < 0)
	pol_filter = 0;

    /*
    if (debug) {
	// RA, DEC at image edge.
	double rs = out[(ra_axis*step)+0] * RAD2DEG;
	double re = out[(ra_axis*step)+1] * RAD2DEG;
	double ds = out[(dec_axis*step)+0] * RAD2DEG;
	double de = out[(dec_axis*step)+1] * RAD2DEG;

	rs = ((rs < 0) ? (rs + 360) : rs);
	re = ((re < 0) ? (re + 360) : re);

	printf ("Spatial edge:  ra_start %.*g ra_end %.*g\n",
	    DBL_DIG, rs, DBL_DIG, re);
	printf ("          dec_start %.*g dec_end %.*g\n",
	    DBL_DIG, ds, DBL_DIG, de);
    }
    */

retry_:
    if (spat_filter) {
	/*  AST expresses spatial coords in radians.  ICRS is assumed here,
	 *  and is enforced in read_header().  If either axis of the image
	 *  is flipped the axis will be flipped to the standard orientation
	 *  in the generated cutout.  If the image is rotated the cutout
	 *  will also be rotated, and will only approximate the requested
	 *  region.
	 */
	if (height <= 0.0)
	    height = width;

	double ra1  = ra  - (width  / 2.0);
	double ra2  = ra  + (width  / 2.0);
	double dec1 = dec - (height / 2.0);
	double dec2 = dec + (height / 2.0);

	if (debug) {
	    printf ("Spatial:  ra %.*g dec %.*g size %.*g %.*g\n",
		DBL_DIG, ra, DBL_DIG, dec, DBL_DIG, width, DBL_DIG, height);
	    printf ("       :  r1 %.*g r2 %.*g d1 %.*g d2 %.*g\n",
		DBL_DIG, ra1, DBL_DIG, ra2, DBL_DIG, dec1, DBL_DIG, dec2);
	}

	out[(ra_axis*step)+0]  = ra1 / RAD2DEG;
	out[(ra_axis*step)+1]  = ra2 / RAD2DEG;
	out[(dec_axis*step)+0] = dec1 / RAD2DEG;
	out[(dec_axis*step)+1] = dec2 / RAD2DEG;
    }

    if (spec_filter) {
	/*  Only WAVE and FREQ image coords are currently supported, in a
	 *  variety of units.  If the computation fails, due to a sufficiently
	 *  wierd image or whatever, we disable filtering of the axis and
	 *  merely return the entire spectral axis.
	 */
	double em1 = dl_wave2img (wavelo, em_ctype, em_unit, &status1);
	double em2 = dl_wave2img (wavehi, em_ctype, em_unit, &status2);
	if (status1 == 0 && status2 == 0) {
	    out[(em_axis*step)+0] = em1;
	    out[(em_axis*step)+1] = em2;
	}

	if (debug) {
	    printf ("Spectral:  wavelo %.*g em1 %.*g wavehi %.*g em2 %.*g\n",
		DBL_DIG, wavelo, DBL_DIG, em1, DBL_DIG, wavehi, DBL_DIG, em2);
	}
    }

    if (time_filter) {
	/* ignore; not yet implemented. */
    }
    if (pol_filter) {
	/* ignore; not yet implemented. */
    }

    /*  We now have the bounding world coordinates of the (unclipped) cutout
     *  region.  Axes which were not filtered are unchanged.  We perform the
     *  reverse transformation to get back to image pixel coordinates.
     */
    if (debug) {
        printf("Bounding B4: in %g %g %g %g\n", in[0], in[1], in[2], in[3]);
        printf("Bounding B4: out %g %g %g %g\n", out[0], out[1], out[2], out[3]);
    }
    astTranN (wcs, npts, imdim, step, out, 0, imdim, step, in);
    if (debug) {
        printf("Bounding after: in %g %g %g %g\n", in[0], in[1], in[2], in[3]);
        printf("Bounding after: out %g %g %g %g\n", out[0], out[1], out[2], out[3]);
    }

    /*  Clip the cutout region, which may extend beyond the bounds of the
     *  target image, to the physical image pixel matrix of the target image.
     *  If any projected axis lies completely outside the target image the
     *  cutout has no coverage in the ROI (it must have coverage in all dim-
     *  ensions), and is rejected.
     *
     *  The final region in pixel coordinates is specified by cutout1,2[n] 
     *  where N is the axis number (0-4) and cutout1-2 are the starting and
     *  ending pixel values for that axis.
     */
    for (i=0;  i < naxes;  i++) {
        double v1 = in[(i*step)+0];
        double v2 = in[(i*step)+1];

        if (v1 < 0.0) 
            v1 += naxis[i], v2 += naxis[i];

        // Check that we have some coverage.
        if (v1 > naxis[i] || v2 < 1.0) {
            char region[20];  sprintf(region, "axis%d: %s", i+1, ctype[i]);
            if (debug)
                printf ("no coverage ax %d: v1=%g axlim=%ld v2=%g, axlim=1\n",
		    i, v1, naxis[i], v2);
        }

        if (debug)
            printf("Cutout: axis %d: %f %f\n", i, v1, v2);

        in[(i*step)+0] = cutout1[i] = max(1, min(naxis[i], nint(v1)));
        in[(i*step)+1] = cutout2[i] = max(1, min(naxis[i], nint(v2)));

        if (debug)
            printf("Cutout: clip %d: %d %d\n", i, cutout1[i], cutout2[i]);
    }

    /*  Do a final forward transform back to World coordinates to define
     *  the physical coverage of the clipped cutout region.
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
     */ 

    /*  Compute RA,DEC in degrees of image corner 1,1.
     */
    in[(ra_axis*step)+0] = cutout1[ra_axis];
    in[(dec_axis*step)+0] = cutout1[dec_axis];
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    ra1 = out[(ra_axis*step)+0] * RAD2DEG;
    dec1 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Compute RA,DEC in degrees of image corner 2,2.
     */
    in[(ra_axis*step)+0] = cutout2[ra_axis];
    in[(dec_axis*step)+0] = cutout1[dec_axis];
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    ra2 = out[(ra_axis*step)+0] * RAD2DEG;
    dec2 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Compute RA,DEC in degrees of image corner 3,3.
     */
    in[(ra_axis*step)+0] = cutout2[ra_axis];
    in[(dec_axis*step)+0] = cutout2[dec_axis];
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    ra3 = out[(ra_axis*step)+0] * RAD2DEG;
    dec3 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Compute RA,DEC in degrees of image corner 4,4.
     */
    in[(ra_axis*step)+0] = cutout1[ra_axis];
    in[(dec_axis*step)+0] = cutout2[dec_axis];
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    ra4 = out[(ra_axis*step)+0] * RAD2DEG;
    dec4 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Make sure that RA values are in the range 0-360 (not +/-180).
     *  Can also be done this way:  double pt1[2]; astNorm(wcs, pt1);
     */
    ra1 = (ra1 < 0.0) ? (360.0 + ra1) : ra1;
    ra2 = (ra2 < 0.0) ? (360.0 + ra2) : ra2;
    ra3 = (ra3 < 0.0) ? (360.0 + ra3) : ra3;
    ra4 = (ra4 < 0.0) ? (360.0 + ra4) : ra4;

    /*  Image center. */
    if (ra3 > ra1)
	ra_cen = (ra3 - ra1) / 2.0 + ra1;
    else
	ra_cen = (ra1 - ra3) / 2.0 + ra3;
    if (dec3 > dec1)
	dec_cen = (dec3 - dec1) / 2.0 + dec1;
    else
	dec_cen = (dec1 - dec3) / 2.0 + dec3;


//if (debug) {
    /*  Compute RA,DEC in degrees of image corner 1,1.
     */
    in[(ra_axis*step)+0] = 1.0;
    in[(dec_axis*step)+0] = 1.0;
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    cra1 = out[(ra_axis*step)+0] * RAD2DEG;
    cdec1 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Compute RA,DEC in degrees of image corner 2,2.
     */
    in[(ra_axis*step)+0] = naxis[ra_axis];
    in[(dec_axis*step)+0] = 1.0;
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    cra2 = out[(ra_axis*step)+0] * RAD2DEG;
    cdec2 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Compute RA,DEC in degrees of image corner 3,3.
     */
    in[(ra_axis*step)+0] = naxis[ra_axis];
    in[(dec_axis*step)+0] = naxis[dec_axis];
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    cra3 = out[(ra_axis*step)+0] * RAD2DEG;
    cdec3 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Compute RA,DEC in degrees of image corner 4,4.
     */
    in[(ra_axis*step)+0] = 1.0;
    in[(dec_axis*step)+0] = naxis[dec_axis];
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    if (!astOK)
	dl_error (6, "cannot read metadata for image", imagefile);
    cra4 = out[(ra_axis*step)+0] * RAD2DEG;
    cdec4 = out[(dec_axis*step)+0] * RAD2DEG;

    /*  Make sure that RA values are in the range 0-360 (not +/-180).
     *  Can also be done this way:  double pt1[2]; astNorm(wcs, pt1);
     */
    cra1 = (cra1 < 0.0) ? (360.0 + cra1) : cra1;
    cra2 = (cra2 < 0.0) ? (360.0 + cra2) : cra2;
    cra3 = (cra3 < 0.0) ? (360.0 + cra3) : cra3;
    cra4 = (cra4 < 0.0) ? (360.0 + cra4) : cra4;


    /*  If we're not doing a cutout, set the query position to 
     *  the image/ccd center.
     */
    if (ra > 360.0 && dec > 90.0) {
	if (debug)
	    printf ("No query POS set, retrying.\n");

        if (cra3 > cra1) {
	    ra_cen = (cra3 - cra1) / 2.0 + cra1;
	    width = (cra3 - cra1) + 0.01;
        } else {
	    ra_cen = (cra1 - cra3) / 2.0 + cra3;
	    width = (cra1 - cra3) + 0.01;
	}
        if (cdec3 > cdec1) {
	    dec_cen = (cdec3 - cdec1) / 2.0 + cdec1;
	    height = (cdec3 - cdec1) + 0.01;
        } else {
	    dec_cen = (cdec1 - cdec3) / 2.0 + cdec3;
	    height = (cdec1 - cdec3) + 0.01;
	}
	ra = ra_cen;
	dec = dec_cen;

	ra1 = cra1;   dec1 = cdec1;
	ra2 = cra2;   dec2 = cdec2;
	ra3 = cra3;   dec3 = cdec3;
	ra4 = cra4;   dec4 = cdec4;

	goto retry_;
    }


    if (debug) {
        printf("CCD:\n");
        printf("     Corn1:  %lg %lg\n", cra1, cdec1);
        printf("     Corn2:  %lg %lg\n", cra2, cdec2);
        printf("     Corn3:  %lg %lg\n", cra3, cdec3);
        printf("     Corn4:  %lg %lg\n", cra4, cdec4);
        printf("Cutout:\n");
        printf("    Center:  %lg %lg\n", ra_cen, dec_cen);
        printf("     Corn1:  %lg %lg\n", ra1, dec1);
        printf("     Corn2:  %lg %lg\n", ra2, dec2);
        printf("     Corn3:  %lg %lg\n", ra3, dec3);
        printf("     Corn4:  %lg %lg\n", ra4, dec4);
    }


    /*  Image scale.  We assume that pixels are square.  DEC is used to
     *  avoid the cos(dec) term required to convert delta-RA to angular
     *  extent on the sky.
    double im_scale = ((dec3 - dec1) / dec_len);
    im_scale = im_scale * (60 * 60);
     */

    /*  Metadata for a virtual image is saved to a file in the staging
     *  area, and may be referenced later to create the described virtual
     *  image (as when metadata computed for a virtual image in a query
     *  response is later used to create and retrieve the virtual image).
     *  Saving the metadata in a file is a simple way to do this for now.
     *  A more scalable and persistent solution would be to create a DBMS
     *  record describing the virtual image, and reference it externally
     *  with a PubDID.
     */
    long cutaxis[MAXAXES];
    FILE *fp = NULL;
    int fd = -1;
    long npix;

    /*  Create the MDFILE (metadata link storage file).
     */
    if (dir == NULL)
	dl_error (8, "root directory of staging area not specified", NULL);

    mdfile = (char *) malloc ((size_t)256);
    sprintf (mdfile, "%s/image-XXXXXX", dir);
    if ((fd = mkstemp (mdfile)) < 0)
	dl_error (9, "cannot create metadata link file", mdfile);
    else
	fp = fdopen (fd, "r+");

    /*  Save the MDFILE filename, minus the staging directory prefix.
     */
    fprintf (fp, "# Virtual image definition file.\n\n");
    fprintf (fp, "MDFILE = %s\n", mdfile+strlen(dir)+1);

    /*  Save the filename of the archive image.
     */
    if (extn > 0)
        fprintf (fp, "image = %s[%d]\n", imagefile, extn);
    else
        fprintf (fp, "image = %s\n", imagefile);


    /* ---- FILTER TERM ---- */
    fprintf (fp, "\n[filter]\n\n");

    if (filter_term) {
	fprintf (fp, "filter_term = true\n");

	/*  Record the filter term required to extract the cutout region.
	 */
	switch (naxes) {
	    case 2:
		fprintf (fp, "cutout = [%d:%d,%d:%d]\n",
		    cutout1[0], cutout2[0], cutout1[1], cutout2[1]);
		break;
	    case 3:
		fprintf (fp, "cutout = [%d:%d,%d:%d,%d:%d]\n",
		    cutout1[0], cutout2[0], cutout1[1], cutout2[1],
		    cutout1[2], cutout2[2]);
		break;
	    case 4:
		fprintf (fp, "cutout = [%d:%d,%d:%d,%d:%d,%d:%d]\n",
		    cutout1[0], cutout2[0], cutout1[1], cutout2[1],
		    cutout1[2], cutout2[2], cutout1[3], cutout2[3]);
		break;
	    default:
		fprintf (fp, "cutout = [*]\n");
		break;
	}

	/*  Count the number of axes and compute their lengths.
	 */
	for (i=0, real_naxes=0, datalen=1;  i < naxes;  i++) {
	    int v1 = cutout1[i]; int v2 = cutout2[i];
	    if (v2 > v1)
		npix = v2 - v1 + 1;
	    else
		npix = v1 - v2 + 1;

	    if (npix > 1) {
		cutaxis[real_naxes] = npix;
		real_naxes++;
	    }
	    datalen *= npix;
	}

	/*  A single pixel is an image with a single axis of length 1.
	 */
	fprintf (fp, "im_naxes = %d\n", max(1, real_naxes));
	switch (real_naxes) {
	    case 1:
		fprintf (fp, "im_naxis = %ld\n", cutaxis[0]);
		break;
	    case 2:
		fprintf (fp, "im_naxis = %ld %ld\n", cutaxis[0], cutaxis[1]);
		break;
	    case 3:
		fprintf (fp, "im_naxis = %ld %ld %ld\n",
		    cutaxis[0], cutaxis[1], cutaxis[2]);
		break;
	    case 4:
		fprintf (fp, "im_naxis = %ld %ld %ld %ld\n",
		    cutaxis[0], cutaxis[1], cutaxis[2], cutaxis[3]);
		break;
	    default:
		fprintf (fp, "im_naxis = ****\n");
	}

	fprintf (fp, "access_estsize = %ld\n", (datalen*(abs(bitpix)/8)/1024));
	fprintf (fp, "dataset_length = %ld\n", datalen);
	fprintf (fp, "obs_creation_type = cutout\n");

	if (spat_filter) {
	    /*  Image corners.
	     *  [Does not change in a cutout, omitted.]
	     *
	    fprintf (fp, "im_ra1 = %.*g\n", DBL_DIG, ra1);
	    fprintf (fp, "im_dec1 = %.*g\n", DBL_DIG, dec1);
	    fprintf (fp, "im_ra2 = %.*g\n", DBL_DIG, ra2);
	    fprintf (fp, "im_dec2 = %.*g\n", DBL_DIG, dec2);
	    fprintf (fp, "im_ra3 = %.*g\n", DBL_DIG, ra3);
	    fprintf (fp, "im_dec3 = %.*g\n", DBL_DIG, dec3);
	    fprintf (fp, "im_ra4 = %.*g\n", DBL_DIG, ra4);
	    fprintf (fp, "im_dec4 = %.*g\n", DBL_DIG, dec4);
	     */

	    /*  Image scale.
	     *  [Does not change in a cutout, omitted.]
	     *
	    fprintf (fp, "im_scale = %g\n", im_scale);
	     */

	    /*  Image center.
	     */
	    fprintf (fp, "s_ra = %.*g\n", DBL_DIG, ra_cen);
	    fprintf (fp, "s_dec = %.*g\n", DBL_DIG, dec_cen);

	    /*  Image field of view.  Let's use the smaller extent.
	     */
	    double ra_width = abs(ra3 - ra1) * cos(dec3 / RAD2DEG);
	    double dec_width = abs(dec3 - dec1);
	    fprintf(fp, "s_fov = %g\n", min(ra_width, dec_width));

	    /*  Image footprint (STC AstroCoordArea).
	     */
	    fprintf (fp, "s_region = polygon icrs");
	    fprintf (fp, " %.*g %.*g", DBL_DIG, ra1, DBL_DIG, dec1);
	    fprintf (fp, " %.*g %.*g", DBL_DIG, ra2, DBL_DIG, dec2);
	    fprintf (fp, " %.*g %.*g", DBL_DIG, ra3, DBL_DIG, dec3);
	    fprintf (fp, " %.*g %.*g", DBL_DIG, ra4, DBL_DIG, dec4);
	    fprintf (fp, "\n");
	}

	if (spec_filter) {
	    /*  Spectral axis, if given.  This will be skipped if both the CTYPE
	     *  and CUNIT are not defined, or have invalid values.
	     */
	    if (em_ctype && em_unit) {
		double em1 = out[(em_axis*step)+0];
		double em2 = out[(em_axis*step)+1];

		em1 = dl_img2wave (em1, em_ctype, em_unit, &status1);
		em2 = dl_img2wave (em2, em_ctype, em_unit, &status2);
		if (em2 > em1) {
		    double temp = em1;
		    em1 = em2;
		    em2 = temp;
		}

		if (status1 == 0 && status2 == 0) {
		    /* Spectral coverage. */
		    fprintf (fp, "em_min = %.*g\n", DBL_DIG, em1);
		    fprintf (fp, "em_max = %.*g\n", DBL_DIG, em2);

		    if ((em2-em1) > EPSILOND) {
			/* Spectral resolution. */
			double em_res = (em2 - em1) / naxis[em_axis];
			double em_loc = (em2 - em1) / 2.0 + em1;
			fprintf(fp, "em_resolution = %g\n", em_res);

			/* Spectral resolving power. */
			fprintf (fp, "em_res_power = %g\n", em_loc / em_res);
		    }
		}
	    }
	}
    } else
	fprintf (fp, "filter_term = false\n");

    /*  ---- WCS TERM ---- */
    fprintf (fp, "\n[wcs]\n\n");
    fprintf (fp, "wcs_term = false\n");

    /*  ---- PIXEL TERM ---- */
    fprintf (fp, "\n[pixel]\n\n");

    if (pixel_term) {
	fprintf (fp, "pixel_term = true\n");
	fprintf (fp, "section = %s\n", imsection);
    } else
	fprintf (fp, "pixel_term = false\n");

    /*  ---- FUNCTION TERM ---- */
    fprintf (fp, "\n[function]\n\n");
    fprintf (fp, "function_term = false\n");

    fprintf (fp, "END\n");
    fclose (fp);
    astEnd;

    return (mdfile);
}


/* Read the image header and extract essential metadata.
 */
static AstFrameSet *
read_header (char *imagefile, int *naxes, long *naxis, int *bitpix, int *axmap,
    char ctype[][KEYLEN], char cunit[][KEYLEN], char comment[][COMLEN], int maxaxes) 
{
    AstFrameSet *wcsinfo = NULL;
    AstFitsChan *fitschan = NULL;
    fitsfile *fptr = NULL;
    int   status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    //int   hdupos = 0, nkeys = 0;
    int   i;
    char *header = (char *) NULL;


    if (fits_open_file (&fptr, imagefile, READONLY, &status) != 0) {
fprintf (stderr, "fits_open_file() failed w/ status = %d\n", status);
fits_report_error (stderr, status);
	dl_error (14, "error reading FITS header", imagefile);
	exit (2);

    } else {
	char key[20];

	dl_loadPHU (fptr); 			// load the primary header

	if (extn > 0) {
	    fitsfile *eptr;
	    char  epath[1024];
	    int   status = 0, nkeys;

	    memset (epath, 0, 1024);
	    sprintf (epath, "%s[%d]", imagefile, extn);
    	    if (fits_open_file (&eptr, epath, READONLY, &status) != 0) {
	    	ehu_header = strdup (" ");
	    } else {
		if (ehu_header == NULL)
		    ehu_header = calloc (1, 2880000);

		status = 0;
		if (fits_hdr2str (eptr, 0, NULL, 0, &ehu_header, &nkeys,
		    &status ) != 0) {
			strcpy (ehu_header, " ");
		}
	    	dl_loadEHU (eptr);
	    }

            /*  Merge the PHU and EHU headers so we can make use of
             *  inherited keywords to find the WCS.
             */
            header = dl_mergeHeaders (phu_header, ehu_header);

	    status = 0;
	    fits_close_file (fptr, &status);
	    fptr = eptr;
    	} else {
	    header = calloc (1, strlen (phu_header));
	    strcpy (header, phu_header);
	}


	/*  Create a FitsChan and fill it with FITS header cards.
  	 */
	fitschan = astFitsChan (NULL, NULL, "");
	astPutCards (fitschan, header);

	/* Read WCS information from the FitsChan. */
	wcsinfo = astRead (fitschan);
	int wcsdim = astGetI (wcsinfo, "Naxes");

	/* Initialize output arrays. */
	memset(naxis, 0, sizeof(*naxis) * maxaxes);
	for (i=0;  i < maxaxes;  i++) {
	    ctype[i][0] = '\0';
	    cunit[i][0] = '\0';
	}

	/* Get the primary FITS header metadata. */
	fits_get_img_param (fptr, maxaxes, bitpix, naxes, naxis, &status);

	/* Since we are working with the WCS here, we want WCSDIM not NAXES */
	*naxes = max(*naxes, wcsdim);

	for (i=0;  i < *naxes;  i++) {
	    sprintf (key, "CTYPE%d", i + 1);
	    if (fits_read_key (fptr,TSTRING,key,ctype[i],comment[i],&status))
	        break;
	    sprintf (key, "CUNIT%d", i + 1);
	    if (fits_read_key (fptr,TSTRING,key,cunit[i],NULL,&status)) {
	        cunit[i][0] = '\0';
		status = 0;
	    }
	}

	/* Compute the axis map.  This tells which image axis
	 * corresponds to the 1st or 2nd spatial axis, time axis,
	 * spectral axis, or polarization axis.  For the moment we
	 * use simple string comparison, and support for esoteric
	 * time scales is limited.
	 */
	for (i=AX_SP1;  i <= AX_POL;  i++)
	    axmap[i] = -1;

	for (i=0;  i < *naxes;  i++) {
	    if (debug)
		printf("axis %d, CTYPE = %s\n", i, ctype[i]);

	    if (strncasecmp(ctype[i],      "RA", 2) == 0)
		axmap[i] = AX_SP1;
	    else if (strncasecmp(ctype[i], "DEC", 3) == 0)
		axmap[i] = AX_SP2;
	    else if (strncasecmp(ctype[i], "GLON", 4) == 0)
		axmap[i] = AX_SP1;
	    else if (strncasecmp(ctype[i], "GLAT", 4) == 0)
		axmap[i] = AX_SP2;
	    else if (strncasecmp(ctype[i], "STOKES", 6) == 0)
		axmap[i] = AX_POL;
	    else if (strncasecmp(ctype[i], "FREQ", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "WAVE", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "ENER", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "WAVN", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "VRAD", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "VOPT", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "ZOPT", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "AWAV", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "FELO", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "VELO", 4) == 0)
		axmap[i] = AX_EM;
	    else if (strncasecmp(ctype[i], "UTC", 3) == 0)
		axmap[i] = AX_TIME;
	    else if (strncasecmp(ctype[i], "TT", 2) == 0)
		axmap[i] = AX_TIME;
	    else if (strncasecmp(ctype[i], "TAI", 3) == 0)
		axmap[i] = AX_TIME;
	    else if (strncasecmp(ctype[i], "GMT", 3) == 0)
		axmap[i] = AX_TIME;
	    else
		dl_error (15, "unknown CTYPE value", ctype[i]);
	}

	/* Free the memory holding the concatenated header cards. */
	if (header) 
	    free (header);
	header = NULL;
    }

    if (status == END_OF_FILE)
	status = 0; /* Reset after normal error */
    fits_close_file(fptr, &status);

    if (status || wcsinfo == NULL || !astOK) {
	fits_report_error (stderr, status); /* print any error message */
fprintf (stderr, "final fail w/ status = %d\n", status);
fprintf (stderr, "final fail w/ wcsinfo = %d\n", wcsinfo);
fprintf (stderr, "final fail w/ astOK = %d\n", astOK);
	dl_error (16, "error reading FITS header", imagefile);
    }

    return (wcsinfo);
}


/**
 *  DL_READMDFILE -- Copy the MDFILE content to the indicated output stream.
 */
static void
dl_readMDFile (char *mdfile, FILE *out) 
{
    char buf[4096];
    FILE *fp;

    if ((fp = fopen (mdfile, "r")) == NULL)
	dl_error (17, "cannot open metadata storage file", mdfile);
    while (fgets (buf, 4096, fp) != NULL)
	fputs(buf, out);
}


/**
 *  DL_EXTRACTIMAGE --  Extract the image cutout specified in MDFILE.  We
 *  don't need most of the metadata, only the image file name and the image
 *  section, specified in pixel coordinates, required to extract the filter
 *  term or cutout region.
 *
 *  The MDFILE records the information needed to generate the image filename.
 *  If the image does not exist it is created and the filename returned.  If
 *  the image exists we merely return the filename.  If the flag "force" is
 *  set the image is regenerated.
 *
 *  The filename of the image associated with an MDFILE is created by merely
 *  appending ".fits" to the MDFILE filename.
 */
static char *
dl_extractImage (char *mdfile, int force) 
{
    char  tempimage[256], newimage[256], pngimage[256], mdpath[256];
    char *cutout=NULL, *section=NULL, *sval=NULL, *ip=NULL;
    char  imageref[256], imageroot[256], card[81], *imagefile=NULL;
    int   bitpix, naxis = 0, nkeys, stat;
    long  naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};



    /*  Generate the imagefile pathname.
     */
    if (dir == NULL)
	dl_error (18, "root directory of staging area not specified", NULL);
    //sprintf (mdpath, "%s/%s", dir, mdfile);
    //sprintf (newimage, "%s/%s.fits", dir, mdfile);
    //sprintf (tempimage, "%s/%s-temp.fits", dir, mdfile);

    memset (mdpath, 0, 256);
    memset (newimage, 0, 256);
    memset (pngimage, 0, 256);
    memset (tempimage, 0, 256);

    sprintf (mdpath, "%s", mdfile);
    sprintf (newimage, "%s.fits", mdfile);
    sprintf (pngimage, "%s.png", mdfile);
    sprintf (tempimage, "%s-temp.fits", mdfile);

    if (debug) {
        printf ("mdpath = %s\n", mdfile);
        printf ("newimage = %s.fits\n", mdfile);
        printf ("pngimage = %s.png\n", mdfile);
        printf ("tempimage = %s-temp.fits\n", mdfile);
    }

    /*  Create the image if it does not already exist, or if the
     *  "force" flag is set.
     */
    if ((access (newimage, F_OK) < 0) || force) {
        if (debug) printf ("unlinking %s\n", newimage);
	unlink (newimage);

	fitsfile  *infptr, *outfptr, *phufptr;
	int status = 0, ii = 1;

	/* Get the image pathname.
	 */
	if ((sval = dl_getKeyword (mdpath, "image")) != NULL)
	    imagefile = strdup(sval);

	/* --- Filter Term ---- */

	sval = dl_getKeyword (mdpath, "filter_term");
	if (sval != NULL && strcmp(sval,"true") == 0) {

	    /*  Apply the filter term.  This amounts to a simple cutout, e.g.,
	     *  "image[cutout]".
	     */
	    if ((sval = dl_getKeyword (mdpath, "cutout")) != NULL)
		cutout = strdup(sval);
	    if (imagefile == NULL || cutout == NULL)
		dl_error (19, "cannot identify imagefile", mdfile);


	    /*  Open the input image PHU for the filter term. 
	     */
	    memset (imageref, 0, 256);
	    if (imagefile[0] == '/')
	        sprintf (imageref, "%s%s", imagefile, cutout);
	    else
	        sprintf (imageref, "%s/%s%s", dir, imagefile, cutout);

	    memset (imageroot, 0, 256);
	    strcpy (imageroot, imageref);
    	    if ((ip = strchr (imageroot, '[')))         // get PHU pointer
	        *ip = '\0';
	    strcat (imageroot, "[0]");

	    status = 0;
	    fits_open_file (&phufptr, imageroot, READONLY, &status);
	    if (status) {
	        dl_error (20, "cannot open image root file", imageroot);
	        fits_report_error (stderr, status);
	    }


	    /*  Open the input image section for the filter term. 
	     */
	    status = 0;
	    if (fits_open_file (&infptr, imageref, READONLY, &status) != 0) {
		dl_error (20, "cannot open imageref file", imageref);
	        if (status)
		    fits_report_error (stderr, status);
		status = 0;
	    }


	    /*  Create the output file and copy the data from the specified
	     *  image section for the filter term.  This should automatically
	     *  update the WCS to reflect the smaller image coverage (copying
	     *  the entire image is a special case).
	     */
	    if (!fits_create_file (&outfptr, newimage, &status)) {

          	/*  Explicitly create new image, to support compression.
		 */
                fits_get_img_param (infptr, 9, &bitpix, &naxis, naxes, &status);
          	fits_create_img (outfptr, bitpix, naxis, naxes, &status);

          	/*  Copy all the PHU user keywords (except structural keywords)
		 */
          	fits_get_hdrspace (phufptr, &nkeys, NULL, &status);
          	for (ii = 1; ii <= nkeys; ii++) {
                    fits_read_record (phufptr, ii, card, &status);
                    if (fits_get_keyclass (card) > TYP_CMPRS_KEY)
                        fits_write_record (outfptr, card, &status);
          	}

          	/*  Copy all the user keywords (except structural keywords).
		 */
          	fits_get_hdrspace (infptr, &nkeys, NULL, &status);
          	for (ii = 1; ii <= nkeys; ii++) {
                    fits_read_record (infptr, ii, card, &status);
                    if (fits_get_keyclass (card) > TYP_CMPRS_KEY)
                        fits_write_record (outfptr, card, &status);
          	}

		/*  Copy the cutout pixels to the output image.
	 	 */
		status = dl_copyPixels (infptr, outfptr);

		/*  Reset status after normal error.
		 */
		if (status == END_OF_FILE)
		    status = 0;

		fits_close_file (outfptr,  &status);
	    }

	    /*  Close the PHU and input pointer.
	     */
	    status = 0; fits_close_file (infptr, &status);
	    status = 0; fits_close_file (phufptr, &status);

	    /*  If error occured, print out error message.
	     */
	    if (status) {
		fits_report_error(stderr, status);
		dl_error (21, "cannot create image cutout", imageref);
	    }

	    /*  Apply the Pixel term if any to "newimage".
	     */
	    imagefile = newimage;
	}


	/* --- Pixel Term ---- */

	/*  FIXME -- Add the cutout code from spatial filter above.
	 */

	sval = dl_getKeyword (mdpath, "pixel_term");
	if (sval != NULL && strcmp(sval,"true") == 0) {

	    /*  Apply the pixel term.  This is a pixel space operation manip-
	     *  ulating the pixel matrix.  Currently this is just another image
	     *  section, but this time it is a client-defined section.  It is
	     *  applied after the filter and WCS terms if any, otherwise it is
	     *  applied to the original input image.
	     */

	    if ((sval = dl_getKeyword (mdpath, "section")) != NULL)
		section = strdup(sval);
	    if (newimage[0] == '\0' || section == (char *) NULL)
		dl_error (19, "cannot apply pixel term", mdfile);

	    if (imagefile[0] == '/')
		sprintf (imageref, "%s%s", imagefile, section);
	    else
		sprintf (imageref, "%s/%s%s", dir, imagefile, section);

	    /*  Open the input image section for the filter term.
	     */
	    status = 0;  ii = 1;
	    if (fits_open_file (&infptr, imageref, READONLY, &status) != 0)
		dl_error (20, "cannot open imagefile", imageref);

	    /*  Create the output file and copy the data from the specified
	     *  image section for the filter term.  This should automatically
	     *  update the WCS to reflect the smaller image coverage (copying
	     *  the entire image is a special case).
	     */
	    if (!fits_create_file (&outfptr, tempimage, &status)) {
		/*  Copy every HDU until we get an error.
		while (!fits_movabs_hdu (infptr, ii++, NULL, &status))
		 */
		    fits_copy_hdu (infptr, outfptr, 0, &status);
	 
		/* Reset status after normal error */
		if (status == END_OF_FILE)
		    status = 0;

		fits_close_file (outfptr,  &status);
	    }

	    fits_close_file (infptr, &status);

	    /*  If error occured, print out error message.
	     */
	    if (status) {
		fits_report_error (stderr, status);
		dl_error (21, "cannot create image cutout", imageref);
	    }

	    stat = unlink (newimage);
            if (debug) printf ("extImg unlinked %s = %d\n", newimage, stat);
	    stat = rename (tempimage, newimage);
            if (debug) printf ("extImg renamed %s -> %s %d\n", tempimage, newimage, stat);
	    if (stat != 0)
		dl_error (22, "error renaming temporary image", imageref);
	}
    }

    if (preview) {
	dl_FITS2PNG (newimage, pngimage);
	stat = unlink (newimage);
        if (debug) printf ("preview unlinked %s = %d\n", newimage, stat);
        return ( strdup (pngimage) );
    } else
        return ( strdup (newimage) );
}


/*  DL_COPYPIXELS -- Copy the cutout pixels to the output image.
 */
static int
dl_copyPixels (fitsfile *in, fitsfile *out)
{
    int  bitpix, bytepix, naxis=0, datatype=0, anynul, iteration, status=0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    long first, totpix = 0, npix;
    double *array, *a, bscale = 1.0, bzero = 0.0, nulval = 0.;


    fits_get_img_param (in, 9, &bitpix, &naxis, naxes, &status);
    totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] *
        naxes[5] * naxes[6] * naxes[7] * naxes[8];

    switch(bitpix) {
    case BYTE_IMG: 	datatype = TBYTE; 	break;
    case SHORT_IMG: 	datatype = TSHORT; 	break;
    case LONG_IMG: 	datatype = TLONG; 	break;
    case FLOAT_IMG: 	datatype = TFLOAT; 	break;
    case DOUBLE_IMG: 	datatype = TDOUBLE; 	break;
    }
    bytepix = abs(bitpix) / 8;

    npix = totpix;
    iteration = 0;

    /*  Try to allocate memory for the entire image, use double type to force
     *  memory alignment.
     */
    array = a = (double *) calloc (npix, 2 * bytepix);

    /* if allocation failed, divide size by 2 and try again */
    while (!array && iteration < 10)  {
        iteration++;
        npix = npix / 2;
        array = (double *) calloc(npix, bytepix);
    }

    if (!array)  {
        printf ("dl_copyPixels: Memory allocation error\n");
        return (0);
    }

    /* turn off any scaling so that we copy the raw pixel values */
    fits_set_bscale (in,  bscale, bzero, &status);
    fits_set_bscale (out, bscale, bzero, &status);

    first = 1;
    while (totpix > 0 && !status) {
       /*  Read all or part of image then write it back to output file.
	*/
       fits_read_img (in, datatype, first, npix,
               &nulval, array, &anynul, &status);
       fits_write_img (out, datatype, first, npix, array, &status);

       totpix = totpix - npix;
       first  = first  + npix;
    }
    free ((double *) a);

    return (status);
}


/**
 *  DL_GETKEYWORD --  Small utility to read a keyword from the MDFILE (or any
 *  other keyword=value formatted file.  The returned string if overwritten
 *  in each call.
 */
static char *
dl_getKeyword (char *file, char *name) 
{
    char token[64], *ip, *op;
    static char buf[256];
    FILE *fp;

    if ((fp = fopen(file, "r")) == NULL)
	dl_error (22, "cannot open file", file);

    while (fgets(buf, 256, fp) != NULL) {
	for (ip=buf, op=token;  (ip-buf) < 63 && !isspace(*ip);  )
	    *op++ = *ip++;
	*op++ = '\0';

	if (strcasecmp (token, name) == 0) {
	    while (*ip && (isspace(*ip) || *ip == '='))
		ip++;
	    for (op=ip;  *op && *op != '\n';  op++)
		;
	    *op = '\0';

	    fclose (fp);
	    return (ip);
	}
    }

    fclose (fp);
    return (NULL);
}


/**
 * TASK_USAGE - Print task usage.
 */

static void
task_Usage (void)
{
    fprintf (stderr, "\n  Usage:\n\t"
        "dlmeta [<opts>] [ <file> | @<file> | <dir> ] .... \n\n"
        "  where <opts> includes:\n\n"
        "       -h,--help               this message\n"
        "       -d,--debug              set debug flag\n"
        "       -v,--verbose            set verbose output flag\n"
        "\n"
        "       :\n"
        "\n"
    );
    exit (0);
}


/**
 *  DL_FILEPATH -- Convert the 'fileRef' identifier to a local path in the 
 *  mass store.
 */
char *
dl_filePath (char *fname)
{
#ifdef HAVE_LIBPQ
    PGconn      *conn;
    PGresult    *res;
    int         row = 0;
    static char path[256];
    static char query[256];
    char        *value = NULL;


    memset (path, 0, 256);                      // initialize
    memset (query, 0, 256);

    /*  Format the query string, substituting the returned path for the one
     *  used on the machine.
     */
#ifdef USE_NSA
    sprintf (query,
        "select replace(data_path,'Volumes','nsa') from r_data_main where data_name='%s'", fname);
#else
    if (collection)
        sprintf (query, "select path from ivoa_%s.exposure where fileRef='%s'", 
            collection, fname);
    else
        sprintf (query, "select path from ivoa_nsa.exposure where fileRef='%s'",
            fname);
#endif


    conn = PQconnectdb (conn_info);             // create the DB connection
    if (PQstatus(conn) == CONNECTION_BAD) {
        fprintf (stderr, "ERROR: Bad database connection\n");
        return (NULL);
    }

//fprintf (stderr, "collection='%s'\n", (collection ? 'none' : collection));
//printf ("conn = '%s'\n", conn_info);
//printf ("query = '%s'\n", query);
    res = PQexec(conn, query);                  // execute the query
    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        fprintf (stderr, "ERROR: Empty query result\n");
        return (NULL);
    }

    value = PQgetvalue(res, row, 0);     // fetch result
//printf ("value = '%s'\n", value);
    memset (path, 0, sizeof(path));
    if (value)
        strncpy (path, value, strlen(value));

    PQclear(res);
    PQfinish(conn);
    return (path);

#else
    return ("/nsa/archive/pipeline/Q20141110/DEC13A/20130209/c4d_130211_024118_ooi_i_v2.fits.fz");
#endif
}



/**
 *  DL_PARSEQUERYSTRING -- Parse the QUERY_STRING for supplied options.
 */
#define	SZ_QUERY_ARG		255

static void
dl_parseQueryString (char *qs, char **fileRef, int *extn, 
	double *ra, double *dec, double *sz1, double *sz2, int *preview,
	int *cutout)
{
    char  *ip = qs, *kp, *vp;
    char   keyw[SZ_QUERY_ARG], val[SZ_QUERY_ARG], path[1024];


    if (debug)
	printf ("query = '%s'\n", qs);

    memset (path, 0, 1024);
    for (ip=qs; *ip; ) {
	memset (keyw, 0, SZ_QUERY_ARG);
	memset (val, 0, SZ_QUERY_ARG);

	for (kp=keyw; *ip && *ip != '='; ip++)
	    *kp++ = *ip;
	if (*ip == '=') ip++; else break;

	for (vp=val; *ip && *ip != '&'; ip++)
	    *vp++ = *ip;
	if (*ip == '&') ip++;

	if (debug)
	    printf ("key = '%s'   val = '%s'\n", keyw, val);


	if (strncasecmp (keyw, "fileRef", 3) == 0) {
	    strcpy (path, dl_filePath (val));

	} else if (strncasecmp (keyw, "siaRef", 3) == 0) {
	    strcpy (path, dl_filePath (val));

	} else if (strncasecmp (keyw, "col", 3) == 0) {
            if (val[0])
	        collection = strdup (val);

	} else if (strncasecmp (keyw, "extn", 3) == 0) {
	    *extn = atoi (val);

	} else if (strncasecmp (keyw, "pos", 3) == 0) {
	    if (strchr (val, (int)',')) {
		sscanf (val, "%lg,%lg", ra, dec);
	    } else {
		/*  Assume an object name and try to resolve it.  NYI */
		;
	    }

	} else if (strncasecmp (keyw, "size", 3) == 0) {
	    if (strchr (val, (int)',')) {
		sscanf (val, "%lg,%lg", sz1, sz2);
	    } else {
		sscanf (val, "%lg", sz1);
		*sz2 = *sz1;
	    }

	} else if (strncasecmp (keyw, "cutout", 3) == 0) {
	    *cutout = ( (val[0] == 't' || val[0] == 'T') ? 1 : 0 );

	} else if (strncasecmp (keyw, "preview", 3) == 0) {
	    *preview = ( (val[0] == 't' || val[0] == 'T') ? 1 : 0 );

	} else 
	    ;				// ignore unknown arguments
    }

    /*  Append the extension to the filename if only processing one extension.
     */
    if (*fileRef == NULL)
	*fileRef = calloc (1, 1024);

    strcpy (*fileRef, path);

    /*  If we're not doing cutouts, force the size to be 10.0 degrees, i.e.
     *  large enough to ensure we clip at the whole ccd/image.
     */
    if (!cutout)			
	*sz1 = *sz2 = 10.0;
}
