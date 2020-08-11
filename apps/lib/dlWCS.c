/** 
 *  DLWCS -- WCS Utility procedures.
 *
 *    Usage:
 *              dlmeta [<opts>] [ <file> | @<file>] ....
 *
 *  @file       dlmeta.c
 *  @author     Mike Fitzpatrick
 *  @date       2/23/16
 *
 *  @brief      WCS Utility procedures.
 */

#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "fitsio.h"
#include "ast.h"
#include "dlApps.h"


/* Axis map. Maps logical axes in the axis map to physical
 * image axes.
 */
#define	AX_SP1		0		/* First spatial axis */
#define	AX_SP2		1		/* Second spatial axis */
#define	AX_EM		2		/* Spectral (EM) axis */
#define	AX_TIME		3		/* Time axis */
#define	AX_POL		4		/* Polarization axis */


#define DEBUG		0


AstFrameSet *dl_readWCSHeader (fitsfile *fptr, char *imagefile,
	int *naxes, long *naxis, int *bitpix, int *axmap, char ctype[][KEYLEN],
	char cunit[][KEYLEN], char bunit[KEYLEN], char comment[][COMLEN],
	int maxaxes, int is_phu);

int dl_findAxis (int axis_type, int *axmap, int naxes);

double dl_img2wave (double imval, char *ctype, char *unit, int *status);
double dl_wave2img (double wave, char *ctype, char *unit, int *status);

extern char  *ehu_header, *phu_header;

extern int     dl_strCard (fitsfile *fptr, char *keyw, char *value);
extern int     dl_intCard (fitsfile *fptr, char *keyw, int *value);
extern int     dl_dblCard (fitsfile *fptr, char *keyw, double *value);
extern int     dl_floatCard (fitsfile *fptr, char *keyw, float *value);
extern int     dl_dmsCard (fitsfile *fptr, char *keyw, double *value);
extern int     dl_hmsCard (fitsfile *fptr, char *keyw, double *value);

extern char  *dl_mergeHeaders (char *phu, char *ehu);

extern double dl_sexa (char *s);
extern int dl_atoi (char *v);

extern void dl_error (int exit_code, char *error_message, char *tag);




/**
 *  DL_READWCSHEADER -- Read the image header and extract essential metadata.
 */
AstFrameSet *
dl_readWCSHeader (fitsfile *fptr, char *imagefile, int *naxes, long *naxis,
    int *bitpix, int *axmap, char ctype[][KEYLEN], char cunit[][KEYLEN], 
    char bunit[KEYLEN], char comment[][COMLEN], int maxaxes, int is_phu)
{
    AstFrameSet *wcsinfo = NULL;
    AstFitsChan *fitschan = NULL;
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int nkeys = 0;
    char *header = (char *) NULL, *ehdr = (char *) NULL;


    if (ehu_header == NULL)
	ehu_header = calloc (1, 28800000);
    else
	memset (ehu_header, 0, 28800000);

    if (fits_hdr2str (fptr, 0, NULL, 0, &ehdr, &nkeys, &status ) != 0) {
	dl_error (14, "readWCSHeader B4: error reading FITS header", imagefile);
	return ( (AstFrameSet *) NULL );

    } else {
	int i;

	/*  Merge the PHU and EHU headers so we can make use of inherited 
    	 *  keywords to find the WCS.
    	 */
//printf ("------------- PHU HDR ----------------\n%s\n", phu_header);
//printf ("------------- EHU HDR ----------------\n%s\n", ehdr);
	strcpy (ehu_header, ehdr);
	header = dl_mergeHeaders (phu_header, ehdr);
	if (header) 
            free ((char *) ehdr);


	/* Initialize output arrays.
   	 */
	memset(naxis, 0, sizeof(*naxis) * maxaxes);
	for (i=0;  i < maxaxes;  i++) {
	    ctype[i][0] = '\0';
	    cunit[i][0] = '\0';
	    bunit[0] = '\0';
	}

#ifdef READ_FITS_HEADER
	/* Get the primary FITS header metadata. */
	fits_get_img_param(fptr, maxaxes, bitpix, naxes, naxis, &status);

	for (i=0;  i < *naxes;  i++) {
            char key[81];

            memset (key, 0, 81);
	    sprintf (key, "CTYPE%d", i + 1);
	    if (fits_read_key (fptr, TSTRING, key, ctype[i], comment[i],
		&status))
	            break;

            memset (key, 0, 81);
	    sprintf (key, "CUNIT%d", i + 1);
	    if (fits_read_key (fptr, TSTRING, key, cunit[i], NULL, &status)) {
	        cunit[i][0] = '\0';
		status = 0;
	    }
	}
	
	/* Get BUNIT.
	 */
	if (fits_read_key (fptr, TSTRING, "BUNIT", bunit, NULL, &status)) {
	    bunit[0] = '\0';
	    status = 0;
	}
#else
	/*  We assume BITPIX/NAXES/NAXIS are present
	 *  Get the primary FITS header metadata for the extension.
	 */
	status = 0;
	fits_get_img_param(fptr, maxaxes, bitpix, naxes, naxis, &status);

	if (fits_read_key (fptr, TSTRING, "CTYPE1", ctype[0], comment[0],
	    &status) > 0) {
	        strcpy (ctype[0], phu.ctype1);
		status = 0;
	}
	if (fits_read_key (fptr, TSTRING, "CTYPE2", ctype[1], comment[1],
	    &status) > 0) {
	        strcpy (ctype[1], phu.ctype2);
		status = 0;
	}
	if (fits_read_key (fptr, TSTRING, "CUNIT1", cunit[0], comment[0],
	    &status) > 0) {
	        strcpy (cunit[0], phu.cunit1);
		status = 0;
	}
	if (fits_read_key (fptr, TSTRING, "CUNIT2", cunit[1], comment[1],
	    &status) > 0) {
	        strcpy (cunit[1], phu.cunit2);
		status = 0;
	}

	if (fits_read_key (fptr, TSTRING, "BUNIT", bunit, comment[0],
	    &status) > 0) {
	        strcpy (bunit, phu.bunit);
		status = 0;
	}
#endif

	/* Compute the axis map.  This tells which image axis
	 * corresponds to the 1st or 2nd spatial axis, time axis,
	 * spectral axis, or polarization axis.  For the moment we
	 * use simple string comparison, and support for esoteric
	 * time scales is limited.
	 */
	for (i=AX_SP1;  i <= AX_POL;  i++)
	    axmap[i] = -1;

	for (i=0;  i < *naxes;  i++) {
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
	    else {
		//dl_error (15, "unknown CTYPE value '%s'", ctype[i]);
                if (ctype[i][0])
		    fprintf (stderr, "unknown CTYPE value '%s'\n", ctype[i]);
                else
		    fprintf (stderr, "NULL CTYPE value\n");
            }
	}

	/* Create a FitsChan and fill it with FITS header cards.
	 */
	fitschan = astFitsChan (NULL, NULL, "", &status);
        if (status != 0) {
            fprintf (stderr, "readWCSHeader: err in astFitsChan = %d\n",status);
	    fitschan = astDelete (fitschan);
	    astClearStatus;
        }

 
        /* Replace INDEF values to avoid error in astPutCards().
         */
        char *p;
        while ((p = strstr(header, "INDEF")))
            memcpy (p, "  0.0", 5);

	astPutCards (fitschan, header);

	/* Free the memory holding the concatenated header cards.
	 */
	if (header) free (header);
	header = NULL;

	/* Read WCS information from the FitsChan.
	 */
	wcsinfo = NULL;
	wcsinfo = astRead (fitschan);
	astClearStatus;

	if (wcsinfo != NULL)
	    fitschan = astDelete (fitschan);
        //else
        //    fprintf (stderr, "readWCSHeader:  wcsinfo = NULL\n");
    }

    if (status == END_OF_FILE)
	status = 0; 		/* Reset after normal error */

    if (!is_phu && (status || wcsinfo == NULL || !astOK)) {
	fits_report_error (stderr, status); /* print any error message */
	dl_error (16, "readWCSHeader: error reading FITS WCS", imagefile);
	if (DEBUG)
	    fprintf (stderr, "status = %d  wcsinfo = %ld  astOK = %d\n",
	        status, (long)wcsinfo, astOK);
	astClearStatus;
	return ( (AstFrameSet *) NULL );
    }

    return (wcsinfo ? astSimplify (wcsinfo) : (AstFrameSet *) NULL);
}



/************************************************
 * UTILITY METHODS
 ************************************************/

/**
 *  DL_FINDAXIS -- Find the given world axis in the axis map.  The axis map
 *  tells, for each physical image axis, what type of world coordinate is
 *  on the axis.  -1 is returned if the axis type is not found.
 */
int
dl_findAxis (int axis_type, int *axmap, int naxes)
{
    int i;
    for (i=0;  i < naxes;  i++)
	if (axmap[i] == axis_type)
	    return (i);
    return (-1);
}


/**
 *  DL_WAVE2IMAGE --  Convert a wavelength value in meters into the units
 *  used in the image WCS.  (Skip energy units for now).
 */
double
dl_wave2img (double wave, char *ctype, char *unit, int *status)
{
    double val=wave, scale=1.0;

    if (!ctype) {
	*status = -1;

    } else if (strncasecmp(ctype, "FREQ", 4) == 0) {
        if (!ctype) {
	    scale = 1.0;
	    *status = -1;

	} else if (strncasecmp(unit, "Hz", 2) == 0)
	    scale = 1.0;
	else if (strncasecmp(unit, "KHz", 3) == 0)
	    scale = 1.0E3;
	else if (strncasecmp(unit, "MHz", 3) == 0)
	    scale = 1.0E6;
	else if (strncasecmp(unit, "GHz", 3) == 0)
	    scale = 1.0E9;
	else
	    dl_error (10, "unrecognized frequency unit", unit);
	val = (299792458.0 / wave) * scale;

    } else if (strncasecmp(ctype, "WAVE", 4) == 0) {
        if (!ctype) {
	    scale = 1.0;
	    *status = -1;

	} else if (strncasecmp(unit, "m", 1) == 0)
	    scale = 1.0;
	else if (strncasecmp(unit, "cm", 2) == 0)
	    scale = 1.0E2;
	else if (strncasecmp(unit, "mm", 2) == 0)
	    scale = 1.0E3;
	else if (strncasecmp(unit, "um", 2) == 0)
	    scale = 1.0E6;
	else if (strncasecmp(unit, "nm", 2) == 0)
	    scale = 1.0E9;
	else
	    dl_error (11, "unrecognized wavelength unit", unit);
	val = wave * scale;

    } else
	*status = -1;

    return (val);
}


/**
 *  DL_IMG2WAVE -- Convert a wavelength value in image WCS units to meters.
 *  (Skip energy units for now).
 */
double
dl_img2wave (double imval, char *ctype, char *unit, int *status)
{
    double val=imval, scale=1.0;


    if (!ctype) {
	*status = -1;

    } else if (strncasecmp(ctype, "FREQ", 4) == 0) {
        if (!ctype) {
	    scale = 1.0;
	    *status = -1;

	} else if (strncasecmp(unit, "Hz", 2) == 0)
	    scale = 1.0;
	else if (strncasecmp(unit, "KHz", 3) == 0)
	    scale = 1.0E3;
	else if (strncasecmp(unit, "MHz", 3) == 0)
	    scale = 1.0E6;
	else if (strncasecmp(unit, "GHz", 3) == 0)
	    scale = 1.0E9;
	else
	    dl_error (12, "unrecognized frequency unit", unit);
	val = 299792458.0 / (imval / scale);

    } else if (strncasecmp(ctype, "WAVE", 4) == 0) {
        if (!ctype) {
	    scale = 1.0;
	    *status = -1;

	} else if (strncasecmp(unit, "m", 1) == 0)
	    scale = 1.0;
	else if (strncasecmp(unit, "cm", 2) == 0)
	    scale = 1.0E2;
	else if (strncasecmp(unit, "mm", 2) == 0)
	    scale = 1.0E3;
	else if (strncasecmp(unit, "um", 2) == 0)
	    scale = 1.0E6;
	else if (strncasecmp(unit, "nm", 2) == 0)
	    scale = 1.0E9;
	else
	    dl_error (13, "unrecognized wavelength unit", unit);
	val = imval / scale;

    } else
	*status = -1;

    return (val);
}
