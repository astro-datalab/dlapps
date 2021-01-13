/** 
 *  DLMETA -- Generate metadata from a FITS image for DAL services.
 *
 *    Usage:
 *              dlmeta [<opts>] [ <file> | @<file>] ....
 *
 *  @file       dlmeta.c
 *  @author     Mike Fitzpatrick
 *  @date       2/23/16
 *
 *  @brief      Generate metadata from a FITS image for DAL services.
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

#define	KEYLEN		FLEN_KEYWORD
#define	COMLEN		FLEN_COMMENT
#define RAD2DEG		57.2957795
#define	EPSILONR	(1.19e-7)
#define	EPSILOND	(2.22e-16)


#define	BASE_SIA_URL	"https://datalab.noao.edu/svc/cutout"
#define	BASE_FILE_REF	"fileRef"
#define	SIA_COLLECTION	"nsa"
#define	CREATION_TYPE	"Stack"
#define	CREATOR	        "NOAO Data Lab"


/* Axis map. Maps logical axes in the axis map to physical
 * image axes.
 */
#define	AX_SP1		0		/* First spatial axis           */
#define	AX_SP2		1		/* Second spatial axis          */
#define	AX_EM		2		/* Spectral (EM) axis           */
#define	AX_TIME		3		/* Time axis                    */
#define	AX_POL		4		/* Polarization axis            */

#define	HDR_PHU		0		/* Primary header               */
#define	HDR_EHU		1		/* Extension header             */


/* Global params shared by all functions.
*/
static int debug 	= 0;
static int verbose 	= 0;
static int extn 	= -1;

static char *cur_name  	= NULL;		// current file name
static char *cur_path  	= NULL;		// current file's full path
static int   cur_extn  	= 0;		// current MEF extension
static long  cur_fsize  = 0L;		// current file size
static int   use_parent = 0;		// use parent dir in fileref

static char *keyw_file  = NULL;
static char *file_ref   = NULL;
static char *instrument = NULL;
static char *telescope  = NULL;
static char *propid     = NULL;
static char *ctype      = NULL;
static char *creator      = NULL;

static char *sia_tbl  	= NULL;
static char *exp_tbl  	= NULL;
static char *obs_tbl  	= NULL;
static char *collection = NULL;

static FILE *efd	= (FILE *) NULL;
static FILE *ofd	= (FILE *) NULL;
static FILE *s1fd	= (FILE *) NULL;
static FILE *s2fd	= (FILE *) NULL;

char *phu_header 	= (char *) NULL;
char *ehu_header 	= (char *) NULL;
int   wcs_in_phu        = 0;

AstFrameSet *phu_wcs    = (AstFrameSet *) NULL;
AstFrameSet *ehu_wcs    = (AstFrameSet *) NULL;


#define EXP_HDR "obs_collection,obs_id,obs_pub_did,product_type,calib_level,access_url,access_format,access_estsize,ra_j2000,dec_j2000,raj2000_,decj2000_,date_obs,mjd_obs,expnum,exptime,filter,filt_str,fileref,path,proctype,prodtype,obstype,telescope,instrument,obsid,object,proposer,propid,ha,zd,airmass,seeing,plver,photflag,moonangle,fwhm,elliptic,magzero"

#define OBS_HDR "obs_collection,obs_id,obs_pub_did,product_type,calib_level,access_url,access_format,access_estsize,target_name,s_ra,s_dec,s_fov,s_region,s_resolution,s_xel1,s_xel2,t_min,t_max,t_exptime,t_resolution,t_xel,em_min,em_max,em_res_power,em_xel,o_ucd,pol_states,pol_xel,facility_name,instrument_name"

#define SIAV1_HDR "obs_collection,obs_id,obs_pub_did,access_url,access_format,access_estsize,filter,calib_level,datalength,nsubarrays,naxes,naxis1,naxis2,pixtype,wcsaxes1,wcsaxes2,title,creator,collection,pubdid,creationtype,spat_loc1,spat_loc2,spat_lolimit1,spat_lolimit2,spat_hilimit1,spat_hilimit2,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,spat_res1,spat_res2,pix_res1,pix_res2,spec_location,spec_start,spec_stop,spec_res,spec_respower,time_location,time_start,time_stop,fluxaxis_ucd,fluxaxis_unit,el_naxes,el_naxis1,el_naxis2,el_pixtype,el_length,el_size,extname,ccdnum,plver,photflag,fileref"

#define SIAV2_HDR "obs_collection,obs_id,obs_pub_did,access_url,access_format,access_estsize,filter,calib_level,datalength,nsubarrays,naxes,naxis1,naxis2,pixtype,wcsaxes1,wcsaxes2,title,creator,collection,pubdid,creationtype,spat_loc1,spat_loc2,spat_lolimit1,spat_lolimit2,spat_hilimit1,spat_hilimit2,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,spat_res1,spat_res2,pix_res1,pix_res2,spec_location,spec_start,spec_stop,spec_res,spec_respower,time_location,time_start,time_stop,fluxaxis_ucd,fluxaxis_unit,el_naxes,el_naxis1,el_naxis2,el_pixtype,el_length,el_size,extname,ccdnum,plver,photflag,fileref,proctype,prodtype"



/*  Task method declarations.
 */
static FILE *dl_csvOpen (char *tname, char *defname, char *hdr);

static int dl_procDir (char *dirname, int recurse);
static int dl_procFile (char *fname);
static int dl_procFITS (char *fname, char *dir);
static int dl_getMetadata (fitsfile *fptr, int hdu, AstFrameSet *wcs, FILE *fp, char *imagefile);

static void dl_writeObsCore (void);
static void dl_writeExposure (void);
static void dl_writeSIAV1 (void);
static void dl_writeSIAV2 (void);

static char *dl_fixProctype (char *proc);

static void task_Usage (void);


/*  Library method declarations.
 */
extern AstFrameSet *dl_readWCSHeader (fitsfile *fptr, char *imagefile,
	int *naxes, long *naxis, int *bitpix, int *axmap, char ctype[][KEYLEN],
	char cunit[][KEYLEN], char bunit[KEYLEN], char comment[][COMLEN],
	int maxaxes, int is_phu);
extern void dl_error (int exit_code, char *error_message, char *tag);

extern int dl_findAxis (int axis_type, int *axmap, int naxes);
extern double dl_img2wave (double imval, char *ctype, char *unit, int *status);
extern double dl_wave2img (double wave, char *ctype, char *unit, int *status);

extern void dl_loadKeywConfig (char *fname);
extern int dl_getCard (fitsfile *fptr, char *keyw, char *card, int *status);
extern int dl_strCard (fitsfile *fptr, char *card, char *value);
extern int dl_intCard (fitsfile *fptr, char *card, int *value);
extern int dl_dblCard (fitsfile *fptr, char *card, double *value);
extern int dl_floatCard (fitsfile *fptr, char *card, float *value);
extern int dl_dmsCard (fitsfile *fptr, char *keyw, double *value);
extern int dl_hmsCard (fitsfile *fptr, char *keyw, double *value);
extern void dl_loadPHU (fitsfile *fptr);
extern void dl_loadEHU (fitsfile *fptr, int extnum);
extern void dl_printKeywConfig (void);

extern int dl_atoi (char *v);
extern int dl_isFITS (char *v);
extern int dl_isDir (char *v);
extern double dl_sexa (char *s);



/* Command line arguments. */
static char opts[] = "hdvc:f:e:i:o:t:C:E:I:K:O:P:S:T:";
static struct option long_opts[] = {
    { "help",		no_argument,		NULL,	'h' },
    { "debug",		no_argument,		NULL,	'd' },
    { "verbose",	no_argument,		NULL,	'v' },

    { "collection",	required_argument,	NULL,	'c' },
    { "file_ref",	required_argument,	NULL,	'f' },
    { "extn",		required_argument,	NULL,	'e' },
    { "image",		required_argument,	NULL,	'i' },
    { "outfile",	required_argument,	NULL,	'o' },
    { "creation-type",	required_argument,	NULL,	't' },
    { "creator",	required_argument,	NULL,	'C' },
    { "parent",	        no_argument,	        NULL,	'p' },

    { "instrument",	required_argument,	NULL,	'I' },
    { "keywords",	required_argument,	NULL,	'K' },
    { "propid",		required_argument,	NULL,	'P' },
    { "telescope",	required_argument,	NULL,	'T' },

    { "obscore",	required_argument,	NULL,	'O' },
    { "exposure",	required_argument,	NULL,	'E' },
    { "siav2",		required_argument,	NULL,	'S' },

    { NULL,		0,			NULL,	0 },
};



/**
 *  DLMETA - Compute standard ObsTap/SIAV3 metadata for a FITS image.
 *
 *  This uses the generic FITS image and WCS metadata to compute standard VO
 *  (SIAV2) metadata values.  No collection-specific metadata is computed.
 *  If insufficient information is available to compute a given metadata value
 *  then it is omitted.
 */

int main (int argc, char *argv[])
{
    //char *imagefile = NULL;
    //char *ifname = NULL;
    char **flist = NULL, **fstart = NULL;
    int  status = 0;			/* status must be initialized	*/
    int  pos;
    char **pargv, ch, optval[SZ_LINE];


    collection = strdup (SIA_COLLECTION);
    file_ref = strdup (BASE_FILE_REF);

    /* Process command line options. */
    prog_name = argv[0];
    flist = fstart = calloc (argc, sizeof (char *));
    pargv = dl_paramInit (argc, argv, opts, long_opts);
    while ((ch = dl_paramNext (opts,long_opts,argc,pargv,optval,&pos)) != 0) {
	if (ch > 0) {
	    switch (ch) {
	    case 'h':	task_Usage (); 			return (0);
	    case 'd':	debug++; 			break;
	    case 'v':	verbose++; 			break;

	    case 'c':	collection = strdup (optval); 	break;
	    case 'f':	file_ref = strdup (optval); 	break;
	    case 'e':	extn = dl_atoi (optval); 	break;
	    case 'i': 	*flist++ = strdup (optval); 	break;
	    case 't': 	ctype = strdup (optval); 	break;
	    case 'p':	use_parent++; 			break;

	    case 'C':	creator = strdup (optval);	break;
	    case 'K':	keyw_file = strdup (optval);	break;
	    case 'E':	exp_tbl = strdup (optval);	break;
	    case 'O':	obs_tbl = strdup (optval);	break;
	    case 'S':	sia_tbl = strdup (optval);	break;

	    case 'I':	instrument = strdup (optval);	break;
	    case 'T':	telescope = strdup (optval);	break;
	    case 'P':	propid = strdup (optval);	break;

	    default:
		fprintf (stderr, "Invalid option '%s'\n", optval);
		return (ERR);
	    }
	} else if (ch == PARG_ERR) {
            return (ERR);

        } else {
	    /*  Add the argument to the list of files/dirs to process.
	     */
            *flist++ = strdup (optval);
        }
    }


    /*  Sanity checks.
     */
    if (*fstart == NULL) {
	dl_error (2, "no input files specified", NULL);
	return (ERR);
    }

    /*  Open the output CSV file containing the metadata.
     */
    efd = dl_csvOpen (exp_tbl, "exposure.csv", EXP_HDR);
    ofd = dl_csvOpen (obs_tbl, "obscore.csv", OBS_HDR);
    s1fd = dl_csvOpen (sia_tbl, "siav1.csv", SIAV1_HDR);
    s2fd = dl_csvOpen (sia_tbl, "siav2.csv", SIAV2_HDR);

    if (keyw_file) {
        dl_loadKeywConfig (keyw_file);
        if (debug > 1)
            dl_printKeywConfig ();
    }

    cur_path = calloc (1, SZ_PATH);

    /* Compute and output the image metadata. */
    if (*fstart == NULL)
	dl_error (2, "no input source specified", NULL);

    else {
	for (flist=fstart; *flist; flist++) {

	    if (*flist && *flist[0] == '@') {
                /* Process an @-file containing file/directory names.
                 */
                char  *fname = &(flist[0][1]);

	        if (debug || verbose)
		    fprintf (stderr, "Processing FILE = '%s'\n", fname);
	        dl_procFile ( fname );

	    } else if (dl_isDir (*flist)) {
	        if (debug || verbose)
		    fprintf (stderr, "Processing DIR = '%s'\n", *flist);
	        dl_procDir ( *flist, TRUE );

	    } else if (dl_isFITS (*flist) || dl_isGZip (*flist) || dl_isBZip2 (*flist)) {
	        if (debug || verbose)
		    fprintf (stderr, "Processing FITS = '%s'\n", *flist);
	        dl_procFITS ( *flist, NULL );
	    }
	}
    }

    if (status)
        fits_report_error (stderr, status); 	// print any error message


    /*  Clean up.
     */
    for (flist=fstart; *flist; flist++) {	// free the file list
        if (*flist) 
    	    free ((void *) *flist);
    }
    if (fstart) free ((void *) fstart);
    if (cur_path) free ((void *) cur_path);

    if (exp_tbl) free ((void *) exp_tbl);
    if (obs_tbl) free ((void *) obs_tbl);
    if (sia_tbl) free ((void *) sia_tbl);
    if (collection) free ((void *) collection);

    // close the output files
    if (efd) fclose (efd);
    if (ofd) fclose (ofd);
    if (s1fd) fclose (s1fd);
    if (s2fd) fclose (s2fd);

    return (status);
}


/**
 *  DL_CSVOPEN -- Open and initialize a CSV output file.
 */
static FILE *
dl_csvOpen (char *tname, char *defname, char *hdr)
{
    FILE *fd = (FILE *) NULL;
    char *fname = (tname ? tname : defname);
    int   do_init = 0;


    /* If the file doesn't already exist, initialize with the CSV header, 
     * otherwise we are appending the file.
     */
    if (access (fname, W_OK) != 0) 
	do_init = 1;

    fd = fopen (fname, "a+");
    if (do_init)
	fprintf (fd, "%s\n", hdr);

    return (fd);
}


/**
 *  DL_PROCFILE -- Process metadata for the named @-file.
 */
static int
dl_procFile (char *atfile)
{
    register FILE *fp;
    char fname[SZ_FNAME];


    if ((fp = fopen (atfile, "r"))) {
	while (1) {
	    memset (fname, 0, SZ_FNAME);
            if (fscanf (fp, "%s", fname) == EOF)
		break;

	    if (verbose)
		fprintf (stderr, "Processing file: '%s'\n", fname);
	    if (access (fname, R_OK) == 0) {
	        if (dl_isDir (fname)) {
	            if (debug)
		        fprintf (stderr, "Processing DIR = '%s'\n", fname);
   	            if (verbose == 1)
		        fprintf (stderr, "Processing dir: '%s'\n", fname);
	            dl_procDir (fname, TRUE);

		} else if (dl_isFITS (fname) || dl_isGZip (fname) || dl_isBZip2 (fname))
	            dl_procFITS (fname, NULL);

	        else
	            fprintf (stderr, "File '%s' is not a FITS file.\n", fname);
	    } else
	        fprintf (stderr, "File '%s' does not exist.\n", fname);
	}
        fclose (fp);
    }

    return (OK);
}


/**
 *  DL_PROCDIR -- Process metadata for the named directory.
 */
static int
dl_procDir (char *dirname, int recurse)
{
    DIR    *dp;
    struct  dirent *entry;
    struct  stat st;
    char    *cwd, *ip, *op, buf[SZ_PATH], rdir[SZ_PATH], cdir[SZ_PATH];
    char    newpath[SZ_PATH];


    /* Scan through the directory.
     */
    if ((dp = opendir (dirname)) == (DIR *) NULL)
        return (ERR);

    memset (cdir, 0, SZ_PATH);
    memset (rdir, 0, SZ_PATH);
    cwd = getcwd (buf, (size_t) SZ_PATH);

    if (dirname[0] == '/') {				// absolute path
	sprintf (rdir, "%s/", dirname);
    } else {
	if (dirname[0] == '.' && ! dirname[1])
	    sprintf (rdir, "%s/", cwd);
	else
	    sprintf (rdir, "%s/%s/", cwd, dirname);
    }

    for (ip=rdir, op=cdir; *ip; ip++,op++) {		// clean up path
	if (*ip == '.' && (*(ip+1) == '.' || *(ip+1) == '/')) {
	    ip += 2;
	    continue;
	} else if (*ip == '/' && *(ip+1) == '/')
	    ip++;
	*op = *ip;
    }


    while ( (entry = readdir (dp)) != (struct dirent *) NULL ) {
	/* Skip the 'dot' directories.
   	 */
        if (strcmp (entry->d_name, ".") == 0 ||         
            strcmp (entry->d_name, "..") == 0)
                continue;

        memset (newpath, 0, SZ_PATH);
        sprintf (newpath, "%s%s", rdir, entry->d_name);
        (void) stat (newpath, &st);

        if (S_ISDIR(st.st_mode) && recurse) {
            /* Recursively check the subdirs.
            */
            memset (newpath, 0, SZ_PATH);
            sprintf (newpath, "%s%s", rdir, entry->d_name);

	    strcpy (cur_path, newpath);

	    if (verbose>=1)
		fprintf (stderr, "Processing dir '%s'\n", newpath);
            (void) dl_procDir (newpath, recurse);

        } else {
	    if (dl_isFITS (newpath) || dl_isGZip (newpath) || dl_isBZip2 (newpath)) {
	        if (verbose>1)
		    fprintf (stderr, "Processing file '%s'\n", newpath);

	        dl_procFITS (newpath, cdir);
	    }
        }
    }
    closedir (dp);                      /* clean up     */

    return (OK);
}


/**
 *  DL_PROCFITS -- Process metadata for the named FITS file.
 */
static int
dl_procFITS (char *fname, char *dir)
{
    fitsfile *fptr = (fitsfile *) NULL;
    int   hdupos = 0, single = 0, status = 0;
    char  imgfile[SZ_PATH];


    memset (&phu, 0, sizeof (phu));
    memset (&ehu, 0, sizeof (phu));

    cur_name = fname;
    if (fname[0] == '/') {
	/*  Absolute path in filename.
 	 */
	char *ip;

	strcpy (cur_path, fname);

	/*  Extract filename from path.
         */
	for (ip=&fname[strlen(fname)-1]; *ip != '/' && ip >= fname; ip--)
	    ;
	if (*ip == '/')
            if (use_parent) {
	        for (; *ip != '/' && ip >= fname; ip--)
	            ;
            }
	    ip++;
	cur_name = ip;

    } else {
	/*  Relative (or no) path in filename.
 	 */
        static char  *cwd, buf[SZ_PATH];

        memset (buf, 0, SZ_PATH);
	if (dir == NULL) {
            cwd = getcwd (buf, (size_t) SZ_PATH);
	    sprintf (cur_path, "%s/%s", cwd, fname);
	} else {
	    if (dir[strlen(dir)-1] == '/')
	        dir[strlen(dir)-1] = '\0';
	    sprintf (cur_path, "%s/%s", dir, fname);
	}
    }

    long fsize = 0;  struct stat statbuf;
    if (stat(fname, &statbuf) == 0)
	fsize = statbuf.st_size;
    cur_fsize = fsize / 1024;


    /*  Append the extension to the filename if only processing one extension
     *  specified by the user and the "--extn" flag.
     */
    memset (imgfile, 0, SZ_PATH);
    if (extn >= 0) {
	sprintf (imgfile, "%s[%d]", fname, extn);
    } else
        strcpy (imgfile, fname);


    if (!fits_open_file (&fptr, imgfile, READONLY, &status)) {

//fits_get_hdu_num (fptr, &hdupos);  	// get current HDU position
//printf ("after openeing FITS file: HDUPOS = %d\n", hdupos);
	/*  Load the primary header information.
	 */
	dl_loadPHU (fptr);
//printf ("after 1st loadPHU, phu.exptime = %g....\n", phu.exptime);
	if (dl_getMetadata (fptr, HDR_PHU, phu_wcs, stderr, imgfile) == OK)
            ; //wcs_in_phu++;
        else {
            phu_wcs = (AstFrameSet *) NULL;
            memset (&phu, 0, sizeof (phu));
	    dl_loadPHU (fptr);
        }
/*
 */

	/*  Apply the constraint filters before proceeding.  Require we have
	 *  an RA/DEC keyword in the PHU as well
	 */
        if (debug > 1) {
            fprintf (stderr, "inst = '%s' '%s'  |  ", 
                (instrument ? instrument : ""), 
                (phu.instrument ? phu.instrument : ""));
            fprintf (stderr, "tele = '%s' '%s'  |  ", 
                (telescope ? telescope : ""), 
                (phu.telescope ? phu.telescope : ""));
            fprintf (stderr, "ra = '%s'  dec = '%s'\n", 
                (phu.ra ? phu.ra : ""), (phu.dec ? phu.dec : ""));
        }

	//if ( (instrument && strcasecmp (instrument, phu.instrument) != 0) ||
	//     (telescope && strcasecmp (telescope, phu.telescope) != 0) ||
	//     (!phu.ra[0] && !phu.dec[0]) ) {
	if ( (instrument && strcasecmp (instrument, phu.instrument) != 0) ||
	     (telescope && strcasecmp (telescope, phu.telescope) != 0) ) {
                // No required keyword, skip this file.
                status = 0;  /* reset status after error */
	        fits_close_file (fptr, &status);	// close the FITS file
                if (debug || verbose)
	            fprintf (stderr, "Skipping '%s', keyword mismatch\n",
                        imgfile);

		if (ehu_header) {
	    	    free ((char *) ehu_header); ehu_header = NULL;
		}
		if (phu_header) {
	    	    //free ((char *) phu_header); phu_header = NULL;
	    	    fits_free_memory ((char *) phu_header, &status); phu_header = NULL;
		}
	        return (OK);
	}

        fits_get_hdu_num (fptr, &hdupos);  	// get current HDU position
//printf ("HDUPOS = %d\n", hdupos);

        /*  List only a single header if a specific extension was given.
         */
        if (hdupos != 1 || extn >= 0) {
            if (debug)
                fprintf (stderr, "Doing SIF image or extn = %d\n", extn);
            single = 1;
        } else {
            /*  Skip PHU for an MEF image.  If this throws an EOF then try
	     *  to process the primary HDU as a SIF image.
             */
            if (debug > 1)
                fprintf (stderr, "Doing MEF image, skipping to first extn\n");

            if (fits_movrel_hdu (fptr, 1, NULL, &status) == END_OF_FILE) {
fits_get_hdu_num (fptr, &hdupos);
//printf ("B4 loadEHU: HDUPOS = %d\n", hdupos);
	        dl_loadEHU (fptr, 1);			// load 1st EHU header
//printf ("after 1st loadEHU, phu.exptime = %g....\n", phu.exptime);
	        if (dl_getMetadata (fptr, HDR_EHU, ehu_wcs, stderr, imgfile) == OK) {
		    dl_writeObsCore ();
		    dl_writeExposure ();
	            dl_writeSIAV1 ();
	            dl_writeSIAV2 ();
		} else {
                    if (debug || verbose)
	                fprintf (stderr, "Skipping '%s', EHU loading error\n",
                            imgfile);
		}

		if (ehu_header) {
	    	    free ((char *) ehu_header); ehu_header = NULL;
		}
		if (phu_header) {
	    	    //free ((char *) phu_header); phu_header = NULL;
	    	    fits_free_memory ((char *) phu_header, &status); phu_header = NULL;
		}

		status = 0;				// must be reset 
		fits_close_file (fptr, &status);	// close the FITS file
		return (OK);
	    }
        }

	cur_extn = 0;           // CFITSIO is one-indexed, i.e. PHU = 1

        /*  Main loop through each extension.
         */
        int exp_wrote = 0, hdutype;
        for (; !status; hdupos++)  {
	    if (debug > 1)
	    	fprintf (stderr, "Header listing for HDU #%d:\n", hdupos);

	    cur_extn++;

            /* Process only the image extensions.
             */
            fits_get_hdu_type (fptr, &hdutype, &status); /* Get the HDU type */
            if (hdutype != IMAGE_HDU) {
		astClearStatus;
	        fits_movrel_hdu (fptr, 1, NULL, &status); 
                continue;
            }

	    dl_loadEHU (fptr, cur_extn);		// load EHU header
//printf ("after loop loadEHU, phu.exptime = %g....\n", phu.exptime);
	    if (dl_getMetadata (fptr, HDR_EHU,ehu_wcs,stderr,imgfile) == ERR) {
		fprintf (stderr, "Error reading Header #%d\n", hdupos);
		astClearStatus;
	    } else {
                if (! exp_wrote++) {
	            dl_writeObsCore ();
	            dl_writeExposure ();
                }
	        dl_writeSIAV1 ();
	        dl_writeSIAV2 ();
	    }

//	    if (ehu_header) {				// free the header str
//printf ("loop procFITS %d:  ehu_header = 0x%x\n", hdupos, ehu_header);
//	        free ((char *) ehu_header); ehu_header = NULL;
//	    }

	    if (single)
		break;
	    else {
		// try to move to next HDU
	        fits_movrel_hdu (fptr, 1, NULL, &status); 	
	    }

	    if (ehu_header) {
	        free ((char *) ehu_header); ehu_header = NULL;
	    }
	    if (phu_header) {
	        //free ((char *) phu_header); phu_header = NULL;
	    	fits_free_memory ((char *) phu_header, &status); phu_header = NULL;
	    }
	}

	if (status == END_OF_FILE)
	    status = 0; 			// reset after normal error
	fits_close_file (fptr, &status);	// close the FITS file

    } else {
	fprintf (stderr, "Cannot open FITS file '%s'\n", imgfile);
	return (ERR);
    }



    if (ehu_header) {				// free the header str
	free ((char *) ehu_header); ehu_header = NULL;
    }
    if (phu_header) {
	//free ((char *) phu_header); phu_header = NULL;
	fits_free_memory ((char *) phu_header, &status); phu_header = NULL;
    }

    return (OK);
}


/* DL_FIXPROCTYPE -- Enforce type names.
 */
static char *
dl_fixProctype (char *proc)
{
    if (strncasecmp (proc, "raw", 3) == 0) 
        return "Raw";
    if (strncasecmp (proc, "stack", 5) == 0) 
        return "Stack";
    return (proc);
}


/**
 *  DL_WRITEOBSCORE -- Write the entry into the ObsCore table.
 */
static void
dl_writeObsCore (void)
{
    if (!ofd)
	return;


    fprintf (ofd, "\"%s\",\"%s\",\"ivo://datalab.noao/%s/%s\",\"image\",%d,", 
	phu.propid, phu.obsid, collection, cur_name,
	((strncasecmp (phu.proctype,"raw",3) == 0) ? 1 :
	    ((strncasecmp (phu.proctype,"stack",5) == 0) ? 3 : 2))
        );
    fprintf (ofd, "\"%s?col=%s&%s=%s\",", BASE_SIA_URL, 
        (collection ? collection : phu.obsid), file_ref, cur_name);
    fprintf (ofd, "\"image/fits\",%ld,\"%s\",%.16lg,%.16lg,%.3g,", cur_fsize,
	phu.object, phu.ra_j2000, phu.dec_j2000, ehu.fov);
    fprintf (ofd, "\"CIRCLE 'icrs geocenter' %g %g %.3g\",%.5g,%d,%d,", 
	phu.ra_j2000, phu.dec_j2000, ehu.fov, ehu.scale, 
        ehu.naxis1, ehu.naxis2);
    fprintf (ofd, "%s,%lg,%g,1.0,1,", phu.mjd_obs_str,
	(phu.mjd_obs + (phu.exptime / 86400.)), phu.exptime);

    fprintf (ofd, "%g,%g,", phu.filt_start, phu.filt_end);

    fprintf (ofd, "1.0,1,\"phot.count\",,,\"%s\",\"%s\"\n", 
	phu.telescope, phu.instrument);

    fflush (ofd);
}


/**
 *  DL_WRITEEXPOSURE -- Write the entry into the Exposure table.
 */
static void
dl_writeExposure (void)
{
    if (!efd)
	return;


    fprintf (efd, "\"%s\",\"%s\",\"ivo://datalab.noao/%s/%s\",\"%s\",%d,", 
	phu.propid, phu.obsid, collection, cur_name, ehu.prodtype,
	((strncasecmp (phu.proctype,"raw",3) == 0) ? 1 :
	    ((strncasecmp (phu.proctype,"stack",5) == 0) ? 3 : 2))
        );
    fprintf (efd, "\"%s?col=%s&%s=%s\",", BASE_SIA_URL, 
        (collection ? collection : phu.obsid), file_ref, cur_name);
    fprintf (efd, "\"image/fits\",%ld,", cur_fsize);
    fprintf (efd, "%.16lg,%.16lg,\"%s\",\"%s\",", 
	phu.ra_j2000, phu.dec_j2000, phu.ra,phu.dec);
    if (phu.filter[0])
        fprintf (efd, "\"%s\",%s,%d,%g,\"%s\",\"%s\",", phu.date_obs, 
	    phu.mjd_obs_str, phu.expnum,
	    phu.exptime, phu.filter, phu.filt_card);
    else if (ehu.filter[0])
        fprintf (efd, "\"%s\",%s,%d,%g,\"%s\",\"%s\",", phu.date_obs, 
	    phu.mjd_obs_str, phu.expnum,
	    phu.exptime, ehu.filter, ehu.filt_card);
    else
        fprintf (efd, "\"%s\",%s,%d,%g,\"\",\"%s\",", phu.date_obs, 
	    phu.mjd_obs_str, phu.expnum, phu.exptime, phu.filt_card);
    fprintf (efd, "\"%s\",\"%s\",", cur_name, cur_path);

    fprintf (efd, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",", 
	dl_fixProctype(phu.proctype), ehu.prodtype, phu.obstype,
	phu.telescope, phu.instrument);
    fprintf (efd, "\"%s\",\"%s\",\"%s\",\"%s\",%g,", 
	phu.obsid, phu.object, phu.proposer,
	phu.propid, phu.ha);
    fprintf (efd, "%g,%g,%g,", phu.zd, phu.airmass, phu.dimsee);
    fprintf (efd, "\"%s\",%d,", phu.plver, phu.photflag);
    fprintf (efd, "%g,%g,%g,%g", phu.moonangle, phu.fwhm, 
	phu.elliptic, phu.magzero);

    fprintf (efd, "\n");
    fflush (efd);
}


/**
 *  DL_WRITESIAV2 -- Write the entry into the SIAV2 table.
 */
static void
dl_writeSIAV2 (void)
{
    int  i = 0;


    if (!s2fd)
	return;

    fprintf (s2fd, "\"%s\",\"%s\",\"ivo://datalab.noao/%s/%s#%d\",", 
	phu.propid, phu.obsid, collection, cur_name, cur_extn);
    fprintf (s2fd, "\"%s?col=%s&%s=%s&extn=%d\",\"image/fits\",%ld,", 
	BASE_SIA_URL, 
        (collection ? collection : phu.obsid), file_ref, cur_name, cur_extn, 
        cur_fsize);

    if (phu.filter[0])
        fprintf (s2fd, "\"%s\",", phu.filter);
    else if (ehu.filter[0])
        fprintf (s2fd, "\"%s\",", ehu.filt_card);
    else
        fprintf (s2fd, ",");

    fprintf (s2fd, "%d,%d,0,2,%d,%d,\"%s\",", 
	((strncasecmp (phu.proctype,"raw",3) == 0) ? 1 :
	    ((strncasecmp (phu.proctype,"stack",5) == 0) ? 3 : 2)),
	(ehu.naxis1 * ehu.naxis2 * (abs (ehu.bitpix / 8))),
	ehu.naxis1, ehu.naxis2, ehu.pixtype);
    fprintf (s2fd, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",",
        ehu.ctype1,ehu.ctype2,phu.program,creator,collection,phu.propid,ctype);

    fprintf (s2fd, "%lg,%lg,%lg,%lg,%lg,%lg,", ehu.ra_cen, ehu.dec_cen,
	ehu.ra_limit[0], ehu.dec_limit[0], ehu.ra_limit[1], ehu.dec_limit[1]);
    for (i=0; i < 4; i++)
        fprintf (s2fd, "%lg,%lg,", ehu.ra_corner[i], ehu.dec_corner[i]);
    fprintf (s2fd, "%.5g,%.5g,1,1,", ehu.scale, ehu.scale);
    fprintf (s2fd, "%g,%g,%g,1.0,,", phu.filt_cen,phu.filt_start,phu.filt_end);

    fprintf (s2fd, "\"%s\",%s,%.8lg,\"phot.count\",,",
	phu.date_obs, phu.mjd_obs_str,
	(phu.mjd_obs + (phu.exptime / 86400.)) );

    fprintf (s2fd, "2,%d,%d,\"%s\",%d,%d,",
	ehu.naxis1, ehu.naxis2, ehu.pixtype,
	(ehu.naxis1 * ehu.naxis2),
	(ehu.naxis1 * ehu.naxis2 * abs(ehu.bitpix / 8)) );

    fprintf (s2fd, "\"%s\",%d,", ehu.extname, ehu.ccdnum);
    fprintf (s2fd, "\"%s\",%d,\"%s\",\"%s\",\"%s\"\n", 
        phu.plver, phu.photflag, cur_name, dl_fixProctype(ehu.proctype),
        ehu.prodtype);
    fflush (s2fd);
}


/**
 *  DL_WRITESIAV1 -- Write the entry into the SIAV1 table.
 */
static void
dl_writeSIAV1 (void)
{
    int  i = 0;


    if (!s1fd)
	return;

    fprintf (s1fd, "\"%s\",\"%s\",\"ivo://datalab.noao/%s/%s#%d\",", 
	phu.propid, phu.obsid, collection, cur_name, cur_extn);
    fprintf (s1fd, "\"%s?col=%s&%s=%s&extn=%d\",\"image/fits\",%ld,", 
	BASE_SIA_URL, 
        (collection ? collection : phu.obsid), file_ref, cur_name, cur_extn, 
        cur_fsize);

    if (phu.filter[0])
        fprintf (s1fd, "\"%s\",", phu.filter);
    else if (ehu.filter[0])
        fprintf (s1fd, "\"%s\",", ehu.filt_card);
    else
        fprintf (s1fd, ",");

    fprintf (s1fd, "%d,%d,0,2,%d,%d,\"%s\",", 
	((strncasecmp (phu.proctype,"raw",3) == 0) ? 1 :
	    ((strncasecmp (phu.proctype,"stack",5) == 0) ? 3 : 2)),
	(ehu.naxis1 * ehu.naxis2 * (abs (ehu.bitpix / 8))),
	ehu.naxis1, ehu.naxis2, ehu.pixtype);
    fprintf (s1fd, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",",
        ehu.ctype1,ehu.ctype2,phu.program,creator,collection,phu.propid,ctype);
    
    fprintf (s1fd, "%lg,%lg,%lg,%lg,%lg,%lg,", ehu.ra_cen, ehu.dec_cen,
	ehu.ra_limit[0], ehu.dec_limit[0], ehu.ra_limit[1], ehu.dec_limit[1]);
    for (i=0; i < 4; i++)
        fprintf (s1fd, "%lg,%lg,", ehu.ra_corner[i], ehu.dec_corner[i]);
    fprintf (s1fd, "%.5g,%.5g,1,1,", ehu.scale, ehu.scale);
    fprintf (s1fd, "%g,%g,%g,1.0,,", phu.filt_cen,phu.filt_start,phu.filt_end);

    fprintf (s1fd, "\"%s\",%s,%.8lg,\"phot.count\",,",
	phu.date_obs, phu.mjd_obs_str,
	(phu.mjd_obs + (phu.exptime / 86400.)) );

    fprintf (s1fd, "2,%d,%d,\"%s\",%d,%d,",
	ehu.naxis1, ehu.naxis2, ehu.pixtype,
	(ehu.naxis1 * ehu.naxis2),
	(ehu.naxis1 * ehu.naxis2 * abs(ehu.bitpix / 8)) );

    fprintf (s1fd, "\"%s\",%d,", ehu.extname, ehu.ccdnum);
    fprintf (s1fd, "\"%s\",%d,\"%s\"\n", phu.plver, phu.photflag, cur_name);
    fflush (s1fd);
}


/**
 *  DL_GETMETADATA -- Compute the SIAV2 image table metadata for an image,
 *  and write it to the given output streams.
 */

#define	DBG	if(debug>1)fprintf

static int
dl_getMetadata (fitsfile *fptr, int hdu, AstFrameSet *wcs, FILE *fp, char *imagefile)
{
    long   naxis[MAXAXES], datalen;
    int    real_naxes=0, npts=2, step=2;
    int    axmap[MAXAXES], imdim, naxes, bitpix, i;
    double ra_cen, dec_cen;
    double ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4;
    char   ctype[MAXAXES][KEYLEN]; char cunit[MAXAXES][KEYLEN];
    double in[MAXAXES*2], out[MAXAXES*2];
    char   comment[MAXAXES][COMLEN], bunit[KEYLEN];
    int    status1=0, status2=0;


//printf ("hdu = %d  wcs_in_phu = %d\n", hdu, wcs_in_phu);
    if (hdu == HDR_EHU && wcs_in_phu) {
        /*  Do a deep copy of the WCS values found in the PHU.
         */
//printf ("doing ehu copy\n");
        memcpy (ehu.ra_corner, phu.ra_corner, sizeof(phu.ra_corner));
        memcpy (ehu.dec_corner, phu.dec_corner, sizeof(phu.dec_corner));
        memcpy (ehu.ra_limit, phu.ra_limit, sizeof(phu.ra_limit));
        memcpy (ehu.dec_limit, phu.dec_limit, sizeof(phu.dec_limit));
        ehu.ra_cen = phu.ra_cen;
        ehu.dec_cen = phu.dec_cen;
        ehu.fov = phu.fov;
        ehu.scale = phu.scale;
        ehu.em_res = phu.em_res;
        ehu.em_loc = phu.em_loc;
        ehu.em_res_power = phu.em_res_power;
        strcpy (ehu.pixtype, phu.pixtype);
        return (OK);
    }

    astBegin;

    /*
     * Collect standard information from the FITS header.
     * ---------------------------------------------------
     */

    /* Read the image geometry and WCS.
     */
    if ((wcs = dl_readWCSHeader (fptr, imagefile, &naxes, naxis,
	    &bitpix, axmap, ctype, cunit, bunit, comment, MAXAXES,
            (int)(hdu == HDR_PHU))) == (AstFrameSet *) NULL) {
                if (hdu == HDR_EHU) {
	            dl_error (5, "cannot read metadata for image", imagefile);
	            astEnd;
	            return (ERR);
                }
    } else if (hdu == HDR_PHU && wcs)
        wcs_in_phu++;

    if (naxes == 0) 
        return (OK);

    /* Identify the image axes.
     */
    int ra_axis = dl_findAxis (AX_SP1, axmap, naxes);
    int dec_axis = dl_findAxis (AX_SP2, axmap, naxes);
    int em_axis = dl_findAxis (AX_EM, axmap, naxes);

    int ra_len = naxis[ra_axis];
    int dec_len = naxis[dec_axis];

    char *em_ctype = ctype[em_axis];
    if (em_ctype[0] == '\0')
	em_ctype = NULL;
    char *em_unit = cunit[em_axis];
    if (em_unit[0] == '\0')
	em_unit = NULL;

    /* Compute the world coordinates of the full image by taking the forward
     * transform of the full pixel array.  Due to the way array dimensioning
     * is used AST wants the coordinates of the two points to be interleaved,
     * e.g. with a stepsize of 2, pt1-val1, pt2-val1, pt1-val2, pt2-val2, etc.
     */
    for (i=0, imdim=naxes;  i < imdim;  i++) {
	in[(i*step)+0] = 1.0;
	in[(i*step)+1] = naxis[i];
    }

    /* Compute RA,DEC in degrees of image corner 1,1.			LL
     */
    in[(ra_axis*step)+0] = 1.0;
    in[(dec_axis*step)+0] = 1.0;
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    //if (!astOK)
    //	dl_error(6, "cannot read metadata for corner 1,1", imagefile);
    ra1 = out[(ra_axis*step)+0] * RAD2DEG;
    dec1 = out[(dec_axis*step)+0] * RAD2DEG;

    /* Compute RA,DEC in degrees of image corner 2,2.			LR
     */
    in[(ra_axis*step)+0] = ra_len;
    in[(dec_axis*step)+0] = 1.0;
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    //if (!astOK)
    //	dl_error(6, "cannot read metadata for corner 2,2", imagefile);
    ra2 = out[(ra_axis*step)+0] * RAD2DEG;
    dec2 = out[(dec_axis*step)+0] * RAD2DEG;

    /* Compute RA,DEC in degrees of image corner 3,3.			UR
     */
    in[(ra_axis*step)+0] = ra_len;
    in[(dec_axis*step)+0] = dec_len;
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    //if (!astOK)
    //	dl_error(6, "cannot read metadata for corner 3,3", imagefile);
    ra3 = out[(ra_axis*step)+0] * RAD2DEG;
    dec3 = out[(dec_axis*step)+0] * RAD2DEG;

    /* Compute RA,DEC in degrees of image corner 4,4.			UL
     */
    in[(ra_axis*step)+0] = 1.0;
    in[(dec_axis*step)+0] = dec_len;
    astTranN (wcs, npts, imdim, step, in, 1, imdim, step, out);
    //if (!astOK)
    //	dl_error(6, "cannot read metadata for corner 4,4", imagefile);
    ra4 = out[(ra_axis*step)+0] * RAD2DEG;
    dec4 = out[(dec_axis*step)+0] * RAD2DEG;


    /* Make sure that RA values are in the range 0-360 (not +/-180).
     */
    ra1 = (ra1 < 0.0) ? (360.0 + ra1) : ra1;
    ra2 = (ra2 < 0.0) ? (360.0 + ra2) : ra2;
    ra3 = (ra3 < 0.0) ? (360.0 + ra3) : ra3;
    ra4 = (ra4 < 0.0) ? (360.0 + ra4) : ra4;

    if (hdu == HDR_PHU) {
        phu.ra_corner[0] = ra1;	phu.dec_corner[0] = dec1;	
        phu.ra_corner[1] = ra2;	phu.dec_corner[1] = dec2;	
        phu.ra_corner[2] = ra3;	phu.dec_corner[2] = dec3;	
        phu.ra_corner[3] = ra4;	phu.dec_corner[3] = dec4;	

        phu.ra_limit[0]  = min( min (ra1,ra2), min(ra3,ra4) );
        phu.ra_limit[1]  = max( max (ra1,ra2), max(ra3,ra4) );
        phu.dec_limit[0] = min( min (dec1,dec2), min(dec3,dec4) );
        phu.dec_limit[1] = max( max (dec1,dec2), max(dec3,dec4) );
    } else {
        ehu.ra_corner[0] = ra1;	ehu.dec_corner[0] = dec1;	
        ehu.ra_corner[1] = ra2;	ehu.dec_corner[1] = dec2;	
        ehu.ra_corner[2] = ra3;	ehu.dec_corner[2] = dec3;	
        ehu.ra_corner[3] = ra4;	ehu.dec_corner[3] = dec4;	

        ehu.ra_limit[0]  = min( min (ra1,ra2), min(ra3,ra4) );
        ehu.ra_limit[1]  = max( max (ra1,ra2), max(ra3,ra4) );
        ehu.dec_limit[0] = min( min (dec1,dec2), min(dec3,dec4) );
        ehu.dec_limit[1] = max( max (dec1,dec2), max(dec3,dec4) );
    }


    /* Count the number of axes and compute their lengths, and the total
     * number of pixels.
     */
    for (i=0, datalen=1;  i < naxes;  i++) {
	long npix = abs(naxis[i]);
	if (npix > 1)
	    real_naxes++;
	datalen *= npix;
    }

    /* Image center.
     */
    if (ra3 > ra1)
	ra_cen = (ra3 - ra1) / 2.0 + ra1;
    else
	ra_cen = (ra1 - ra3) / 2.0 + ra3;
    if (dec3 > dec1)
	dec_cen = (dec3 - dec1) / 2.0 + dec1;
    else
	dec_cen = (dec1 - dec3) / 2.0 + dec3;

    if (hdu == HDR_PHU) {
        phu.ra_cen = ra_cen;
        phu.dec_cen = dec_cen;
    } else {
        ehu.ra_cen = ra_cen;
        ehu.dec_cen = dec_cen;
    }


    /* Image scale.  We assume that pixels are square.  DEC is used to
     * avoid the cos(dec) term required to convert delta-RA to angular
     * extent on the sky.
     */
    double im_scale = ((dec3 - dec1) / dec_len);
    im_scale = im_scale * (60 * 60);     /* arcsec */

    /*
     * Output the standard image metadata.
     * ---------------------------------------------------
     */

    DBG(fp, "\n[%s]\n", imagefile);    /* Start a new context. */

    DBG(fp, "access_format = %s\n", "image/fits");
    long fsize = 0;  struct stat statbuf;
    if (stat(imagefile, &statbuf) != 0)
	fsize = datalen * (abs(bitpix) / 8);
    else
	fsize = statbuf.st_size;
    DBG(fp, "access_estsize = %ld\n", fsize / 1024);
    cur_fsize = fsize / 1024;

    DBG(fp, "dataproduct_type = %s\n", real_naxes > 2 ? "cube" : "image");
    DBG(fp, "dataset_length = %ld\n", datalen);

    /* Currently we support only a single subarray.
     */
    DBG(fp, "im_nsubarrays = 1\n");

    /* Image geometry.
     */
    DBG(fp, "im_naxes = %d\n", real_naxes);
    for (i=0;  i < naxes;  i++)
	DBG(fp, "im_naxis%d = %ld\n", i+1, naxis[i]);

    /* Pixel datatype as in VOTable.
     */
    char pixtype[SZ_VAL];
    memset (pixtype, 0, SZ_VAL);
    switch (bitpix) {
    case 8:
	DBG(fp, "im_pixtype = %s\n", "unsignedByte");
	strcpy (pixtype, "unsignedByte");
	break;
    case 16:
	DBG(fp, "im_pixtype = %s\n", "short");
	strcpy (pixtype, "short");
	break;
    case 32:
	DBG(fp, "im_pixtype = %s\n", "int");
	strcpy (pixtype, "int");
	break;
    case -32:
	DBG(fp, "im_pixtype = %s\n", "float");
	strcpy (pixtype, "float");
	break;
    case -64:
	DBG(fp, "im_pixtype = %s\n", "double");
	strcpy (pixtype, "double");
	break;
    default:
	dl_error(6, "invalid value for bitpix", NULL);
    }
    strcpy ((hdu == HDR_PHU ? phu.pixtype : ehu.pixtype), pixtype);


    /* WCS axes.
     */
    for (i=0;  i < naxes;  i++)
	DBG(fp, "im_wcsaxes%d = %s\n", i+1, ctype[i]);

    /* Image corners.
     */
    DBG(fp, "im_ra1 = %.*g\n", DBL_DIG, ra1);
    DBG(fp, "im_dec1 = %.*g\n", DBL_DIG, dec1);
    DBG(fp, "im_ra2 = %.*g\n", DBL_DIG, ra2);
    DBG(fp, "im_dec2 = %.*g\n", DBL_DIG, dec2);
    DBG(fp, "im_ra3 = %.*g\n", DBL_DIG, ra3);
    DBG(fp, "im_dec3 = %.*g\n", DBL_DIG, dec3);
    DBG(fp, "im_ra4 = %.*g\n", DBL_DIG, ra4);
    DBG(fp, "im_dec4 = %.*g\n", DBL_DIG, dec4);

    /* Image scale.
     */
    DBG(fp, "im_scale = %g\n", im_scale);
    phu.scale = im_scale;
    ehu.scale = im_scale;

    /* Image center.
     */
    DBG(fp, "s_ra = %.*g\n", DBL_DIG, ra_cen);
    DBG(fp, "s_dec = %.*g\n", DBL_DIG, dec_cen);

    /* Image field of view.  Let's use the smaller extent.
     */
    double ra_width = abs(ra3 - ra1) * cos(dec3 / RAD2DEG);
    double dec_width = abs(dec3 - dec1);
    DBG(fp, "s_fov = %g\n", min(ra_width, dec_width));
    if (hdu == HDR_PHU)
        phu.fov = min(ra_width, dec_width);
    else
        ehu.fov = min(ra_width, dec_width);

    /* Image footprint (STC AstroCoordArea).
     */
    DBG(fp, "s_region = polygon icrs");
    DBG(fp, " %g %g", ra1, dec1); DBG(fp, " %g %g", ra2, dec2);
    DBG(fp, " %g %g", ra3, dec3); DBG(fp, " %g %g", ra4, dec4);
    DBG(fp, "\n");

    /* Spectral axis, if given.  This will be skipped if both the CTYPE 
     * and CUNIT are not defined, or have invalid values.
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
            double em_res, em_loc, em_res_power;

	    /* Spectral coverage. */
	    DBG(fp, "em_min = %.*g\n", DBL_DIG, em1);
	    DBG(fp, "em_max = %.*g\n", DBL_DIG, em2);

	    if (em2-em1 > EPSILOND) {
		/* Spectral resolution. */
		em_res = (em2 - em1) / naxis[em_axis];
		em_loc = (em2 - em1) / 2.0 + em1;
		DBG(fp, "em_resolution = %g\n", em_res);

		/* Spectral resolving power. */
		em_res_power = em_loc / em_res;
		DBG(fp, "em_res_power = %g\n", em_res_power);

                if (hdu == HDR_PHU) {
                    phu.em_res = em_res;
                    phu.em_loc = em_loc;
                    phu.em_res_power = em_res_power;
                } else {
                    ehu.em_res = em_res;
                    ehu.em_loc = em_loc;
                    ehu.em_res_power = em_res_power;
                }
	    }
	}
    }

    /* Polarization axis, if given.
     */
    /* to be added. */

    /* Observable axis.
     */
    if (bunit[0] != '\0')
	DBG(fp, "o_unit = %s\n", bunit);

    astEnd;

    return (OK);
}


/************************************************
 * UTILITY METHODS
 ************************************************/

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
        "	:\n"
        "\n"
    );
    exit (0);
}

