/**
 *  DLAPPS.h -- Task declarations for the Data Lab Package applications.
 *
 *  @file       dlApps.h
 *  @author     Mike Fitzpatrick
 *  @date       1/03/16
 *
 *  @brief      Task declarations for the Data Lab Package applications.
 */

#include <getopt.h>

#ifdef  SZ_FORMAT
#undef  SZ_FORMAT
#endif
#define	SZ_FORMAT	32

#ifdef  SZ_FNAME
#undef  SZ_FNAME
#endif
#define SZ_FNAME	256

#ifdef  SZ_PATH
#undef  SZ_PATH
#endif
#define SZ_PATH		512

#ifdef  SZ_LINE
#undef  SZ_LINE
#endif
#define	SZ_LINE		4096

#define PARG_ERR	-512000000

#ifndef OK                              /* Utility values               */
#define OK              0
#endif
#ifndef ERR
#define ERR             1
#endif

#ifndef TRUE
#define TRUE            1
#endif
#ifndef FALSE
#define FALSE           0
#endif


#define MAXAXES         7
#define KEYLEN          FLEN_KEYWORD
#define COMLEN          FLEN_COMMENT
#define RAD2DEG         57.29577951308232
#define DEG2RAD         0.017453292519943
#define EPSILONR        (1.19e-7)
#define EPSILOND        (2.22e-16)

#define abs(x)          ((x)<0 ? -(x):(x))
#define min(x,y)        ((x)<(y)?x:y)
#define max(x,y)        ((x)>(y)?x:y)
#define nint(x)         ((int)((x)+0.5))



#define SZ_VAL          69

/*  Primary header metedata.
 */
struct PHU {
    char proctype[SZ_VAL];
    char prodtype[SZ_VAL];
    char obsid[SZ_VAL];
    char object[SZ_VAL];
    char ra[SZ_VAL];
    char dec[SZ_VAL];
    char date_obs[SZ_VAL];
    char mjd_obs_str[SZ_VAL];
    char filt_card[SZ_VAL];
    char proposer[SZ_VAL];
    char program[SZ_VAL];
    char propid[SZ_VAL];
    char obstype[SZ_VAL];
    char telescope[SZ_VAL];
    char instrument[SZ_VAL];
    char pixtype[SZ_VAL];

    char ctype1[SZ_VAL];
    char ctype2[SZ_VAL];
    char cunit1[SZ_VAL];
    char cunit2[SZ_VAL];
    char bunit[SZ_VAL];
    char plver[SZ_VAL];

    double ra_cen;
    double dec_cen;
    double ra_corner[4];                // Counter-CW from LL
    double dec_corner[4];
    double ra_limit[2];                 // RA min/max
    double dec_limit[2];                // DEC min/max
    double fov;                         // FOV
    double scale;                       // plate scale (""/pix)

    double em_res;
    double em_loc;
    double em_res_power;

    int expnum;
    int photflag;

    float moonangle;
    float fwhm;
    float elliptic;
    float magzero;

    char filter[SZ_VAL];
    float filt_cen;
    float filt_start;
    float filt_end;

    double mjd_obs;
    double ra_j2000;
    double dec_j2000;
    float ha;
    float zd;
    float exptime;
    float airmass;
    float dimsee;
} phu;


/*  Extension header metedata.
 */
struct EHU {
    int extnum;
    int ccdnum;
    int bitpix;
    int naxis1;
    int naxis2;
    char ctype1[SZ_VAL];
    char ctype2[SZ_VAL];
    char cunit1[SZ_VAL];
    char cunit2[SZ_VAL];
    char pixtype[SZ_VAL];
    char bunit[SZ_VAL];
    char extname[SZ_VAL];

    char prodtype[SZ_VAL];
    char proctype[SZ_VAL];
    char filt_card[SZ_VAL];
    char filter[SZ_VAL];

    double ra_cen;
    double dec_cen;
    double ra_corner[4];                // Counter-CW from LL
    double dec_corner[4];
    double ra_limit[2];                 // RA min/max
    double dec_limit[2];                // DEC min/max
    double fov;                         // FOV
    double scale;                       // plate scale (""/pix)

    double em_res;
    double em_loc;
    double em_res_power;
} ehu;


// keep the full PHU header for merging
//char *phu_header = (char *) NULL;       
//char *ehu_header = (char *) NULL;

static char *prog_name	= (char *) NULL;


#define	dabs(x)		((x<0.0?-x:x))

/*  Debug and verbose flags.
 */
#define DLAPP_DEBUG  (getenv("DLAPP_DBG")||access("/tmp/DLAPP_DBG",F_OK)==0)
#define DLAPP_VERB   (getenv("DLAPP_VERB")||access("/tmp/DLAPP_VERB",F_OK)==0)

/**
 *  Output formats.
 */
#define FORMATS "|vot|asv|bsv|csv|tsv|html|shtml|fits|ascii|xml|raw|"

#define VOT     0                       /* A VOTable                    */
#define ASV     1                       /* ascii separated values       */
#define BSV     2                       /* bar separated values         */
#define CSV     3                       /* comma separated values       */
#define TSV     4                       /* tab separated values         */
#define HTML    5                       /* standalone HTML document     */
#define SHTML   6                       /* single HTML <table>          */
#define FITS    7                       /* FITS binary table            */
#define ASCII   8                       /* ASV alias                    */
#define XML     9                       /* VOTable alias                */
#define RAW     10                      /*    "      "                  */


/*  File Formats.
 */
#define	DL_FITS			0
#define	DL_VOTABLE		1
#define	DL_CSV			2
#define	DL_FITS_SPEC		3
#define	DL_VOTABLE_SPEC		4

/*  URL Formats.
 */
#define	DL_LOCALURL		0	/* e.g. http://127.0.0.1/foo	*/
#define	DL_LOCALURI		1	/* e.g. file:///path/foo	*/
#define	DL_LOCALFILE		2	/* e.g. /path/foo		*/
#define	DL_REMOTE 		3	/* e.g. http://foo.bar/junk	*/



/******************************************************************************
 *  Image information structure.  Note that although we can report 3-D
 *  images, we're really only setup to deal with equatorial sky coord
 *  systems.
 *****************************************************************************/
typedef struct {
    char   *imname;                             /* image name 		    */
    int     is_image;                           /* is it an image?  	    */
    int     is_table;                           /* is it a table?  	    */
    int     has_wcs;                            /* image has wcs 	    */

    int     extnum;                             /* extension number 	    */
    int     naxis;                              /* number of axes 	    */
    int     naxes[3];                           /* axis dimensions 	    */
    int     bitpix;                             /* pixel size 		    */
    int     axflip;                             /* are axes flipped?	    */

    double  xc[4], yc[4];                       /* corner positions (wcs)   */
    double  cx, cy;                             /* center position (wcs)    */
    double  lx, ly;                             /* LL corner (wcs) 	    */
    double  ux, uy;                             /* UR corner (wcs) 	    */
    double  xrval, yrval;                       /* CRVAL values 	    */
    double  xrpix, yrpix;                       /* CRPIX values 	    */

    double  width, height;                      /* image width/height (deg) */
    double  radius;                             /* cone radius (deg) 	    */
    double  rotang;                             /* rotation angle (deg)     */
    double  scale;                              /* plate scale (""/pix)     */
    char    ctype[16];                          /* coordinate type 	    */
} frameInfo, *frameInfoP;

typedef struct {
    char        imname[SZ_PATH];                /* full image name 	    */
    int         nextend;                        /* number of extensions     */

    frameInfo   frame;                          /* full frame information   */
    frameInfo  *extns;                          /* image extn information   */
} ImInfo, *ImInfoP;


ImInfo *vot_imageInfo (char *name, int do_all);
void    vot_printImageInfo (FILE *fd, ImInfo *im);
int     vot_imageNExtns (char *image);
void    vot_freeImageInfo (ImInfo *img);



/*  Task structure.
 */
typedef struct {
   char	 *name;				/* task name		      	*/
   int  (*func)(int argc, char **argv, size_t *len, void **result);

   int   ntests;			/* number of unit tests		*/
   int   npass;				/* number of passed tests	*/
   int   nfail;				/* number of failed tests	*/
} Task;


/*  Tasking execution procedure.
 */
int  dl_runTask (char *method, Task *apps, int argc, char **argv, size_t *len, 
		    void **result);
int  dl_taskTest (Task *self, char *arg, ...);
void dl_taskTestFile (char *str, char *fname);
void dl_taskTestReport (Task self);
void dl_taskDbg (void);

int  dl_setResultFromFile (char *fname, size_t *len, void **data);
int  dl_setResultFromString (char *str, size_t *len, void **data);
int  dl_setResultFromInt (int value, size_t *len, void **data);
int  dl_setResultFromReal (float value, size_t *len, void **data);
int  dl_appendResultFromString (char *str, size_t *len, void **data, 
		size_t *maxlen);



/*  Tasking parameter procedures.
 */
char **dl_paramInit (int argc, char *argv[],
                char *opts, struct option long_opts[]);
int    dl_paramNext (char *opts, struct option long_opts[], 
		int argc, char *argv[], char *optval, int *posindex);
void   dl_paramFree (int argc, char *argv[]);


/* Utility Method declarations.
 */
int 	dl_findAxis (int axis_type, int *axmap, int naxes);

double 	dl_wave2img (double wave, char *ctype, char *unit, int *status);
double 	dl_img2wave (double imval, char *ctype, char *unit, int *status);
