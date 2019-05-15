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

       
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
/*
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
*/
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <dirent.h>

#include "fitsio.h"
#include "ast.h"
#include "dlApps.h"


/*  Keyword/constant mapping
 */
#define	SZ_KWVAL	        12
#define	SZ_CONSTANT	        64

struct keywMap {
    char field[SZ_KWVAL];
    char keyw[SZ_KWVAL];
    char constant[SZ_CONSTANT];
} keywCfg[] = {
    { "airmass",  "AIRMASS",  "\0" },
    { "bitpix",   "BITPIX",   "\0" },
    { "bunit",    "BUNIT",    "\0" },
    { "ccdnum",   "CCDNUM",   "\0" },
    { "ctype1",   "CTYPE1",   "\0" },
    { "ctype2",   "CTYPE2",   "\0" },
    { "cunit1",   "CUNIT1",   "\0" },
    { "cunit2",   "CUNIT2",   "\0" },
    { "date-obs", "DATE-OBS", "\0" },
    { "dateobs1", "DATEOBS1", "\0" },
    { "dec",      "DEC",      "\0" },
    { "dimmsee",  "DIMMSEE",  "\0" },
    { "elliptic", "ELLIPTIC", "\0" },
    { "expnum",   "EXPNUM",   "\0" },
    { "exptime",  "EXPTIME",  "\0" },
    { "extname",  "EXTNAME",  "\0" },
    { "extnum",   "EXTNUM",   "\0" },
    { "filter",   "FILTER",   "\0" },
    { "fwhm",     "FWHM",     "\0" },
    { "ha",       "HA",       "\0" },
    { "instrume", "INSTRUME", "\0" },
    { "magzero",  "MAGZERO",  "\0" },
    { "mjd-obs",  "MJD-OBS",  "\0" },
    { "moonangl", "MOONANGL", "\0" },
    { "naxis1",   "NAXIS1",   "\0" },
    { "naxis2",   "NAXIS2",   "\0" },
    { "object",   "OBJECT",   "\0" },
    { "obsid",    "OBSID",    "\0" },
    { "obstype",  "OBSTYPE",  "\0" },
    { "photflag", "PHOTFLAG", "\0" },
    { "plver",    "PLVER",    "\0" },
    { "proctype", "PROCTYPE", "\0" },
    { "prodtype", "PRODTYPE", "\0" },
    { "program",  "PROGRAM",  "\0" },
    { "propid",   "PROPID",   "\0" },
    { "proposer", "PROPOSER", "\0" },
    { "ra",       "RA",       "\0" },
    { "telescop", "TELESCOP", "\0" },
    { "zbitpix",  "ZBITPIX",  "\0" },
    { "zd",       "ZD",       "\0" },
    { "znaxis1",  "ZNAXIS1",  "\0" },
    { "znaxis2",  "ZNAXIS2",  "\0" }
};


/* Axis map. Maps logical axes in the axis map to physical
 * image axes.
 */
#define	AX_SP1		0		/* First spatial axis */
#define	AX_SP2		1		/* Second spatial axis */
#define	AX_EM		2		/* Spectral (EM) axis */
#define	AX_TIME		3		/* Time axis */
#define	AX_POL		4		/* Polarization axis */


#define	DEBUG		0
#define	VERBOSE		0


void   dl_loadPHU (fitsfile *fptr);
void   dl_loadEHU (fitsfile *fptr);


int 	dl_getCard (fitsfile *fptr, char *keyw, char *card, int *status);
int 	dl_strCard (fitsfile *fptr, char *keyw, char *value);
int 	dl_intCard (fitsfile *fptr, char *keyw, int *value);
int 	dl_dblCard (fitsfile *fptr, char *keyw, double *value);
int 	dl_floatCard (fitsfile *fptr, char *keyw, float *value);
int 	dl_dmsCard (fitsfile *fptr, char *keyw, double *value);
int 	dl_hmsCard (fitsfile *fptr, char *keyw, double *value);

char  *dl_mergeHeaders (char *phu, char *ehu);

extern	char *ehu_header, *phu_header;

extern double dl_sexa (char *s);
extern int dl_atoi (char *v);

extern void dl_error (int exit_code, char *error_message, char *tag);




/**
 *  DL_LOADPHU -- Load the PHU header structure.
 */
void
dl_loadPHU (fitsfile *fptr)
{
    double  dval = 0.0;
    int	    nkeys = 0, status = 0;


    /*
    if (phu_header == NULL)
	phu_header = calloc (1, 2880000);
    else
        memset (&phu_header, 0, 2880000);
    */

    if (fits_hdr2str (fptr, 0, NULL, 0, &phu_header, &nkeys, &status ) != 0) {
	fprintf (stderr, "Error reading PHU FITS header\n");
	return;
    }

//printf ("loadPHU:  header str:\n%s\n======================\n", phu_header);
    //memset (&phu, 0, sizeof (struct PHU));

    dl_strCard (fptr, "PROCTYPE", phu.proctype);
    dl_strCard (fptr, "PRODTYPE", phu.prodtype);
    dl_strCard (fptr, "OBSID", phu.obsid);
    dl_strCard (fptr, "OBJECT", phu.object);
    dl_strCard (fptr, "FILTER", phu.filt_card); 
    dl_strCard (fptr, "PROPOSER", phu.proposer);
    dl_strCard (fptr, "PROGRAM", phu.program);
    dl_strCard (fptr, "PROPID", phu.propid);
    dl_strCard (fptr, "OBSTYPE", phu.obstype);
    dl_strCard (fptr, "TELESCOP", phu.telescope);
    dl_strCard (fptr, "INSTRUME", phu.instrument);
    dl_strCard (fptr, "PLVER", phu.plver);
    if (dl_strCard (fptr, "DATE-OBS", phu.date_obs) != OK)
        dl_strCard (fptr, "DATEOBS1", phu.date_obs);

    dl_intCard (fptr, "EXPNUM", &phu.expnum);
    dl_intCard (fptr, "PHOTFLAG", &phu.photflag);

    dl_floatCard (fptr, "MAGZERO", &phu.magzero);
    dl_floatCard (fptr, "MOONANGL", &phu.moonangle);
    dl_floatCard (fptr, "FWHM", &phu.fwhm);
    dl_floatCard (fptr, "ELLIPTIC", &phu.elliptic);

    if (dl_hmsCard (fptr, "HA", &dval) == OK)
        phu.ha = (float) dval;

    dl_strCard (fptr, "MJD-OBS", &phu.mjd_obs_str[0]);
    dl_dblCard (fptr, "MJD-OBS", &phu.mjd_obs);
    if (phu.mjd_obs == 0.0)				// No MJD-OBS keyword
	phu.mjd_obs_str[0] = '\0';

    dl_floatCard (fptr, "ZD", &phu.zd);
    dl_floatCard (fptr, "EXPTIME", &phu.exptime);
    dl_floatCard (fptr, "AIRMASS", &phu.airmass);
    dl_floatCard (fptr, "DIMMSEE", &phu.dimsee);

    memset (phu.filter, 0, SZ_VAL);
    phu.filter[0] = phu.filt_card[0];			// DECam filters
    if (strchr ("ugrizyUGRIZY", (int)phu.filter[0]) == NULL)
        strcpy (phu.filter, phu.filt_card);

    /* Positional information.
     */
    dl_strCard (fptr, "RA", phu.ra);	
    dl_strCard (fptr, "DEC", phu.dec);	
    if ( strchr (phu.ra, (int)':'))  {
        phu.ra_j2000 = dl_sexa (phu.ra) * 15.;		// Sexigesimal value
        phu.dec_j2000 = dl_sexa (phu.dec);
    } else {
        phu.ra_j2000 = atof (phu.ra);		        // Decimage degrees
        phu.dec_j2000 = atof (phu.dec);
    }

    char  f = '\0', j1[12], j2[12], j3[12];
    float wcen, wwidth;

    /*  Parse a filter string such as 'i DECam SDSS c0003 7835.0 1470.0'
     */
    sscanf (phu.filt_card, "%c %s %s %s %g %g", &f, j1, j2, j3, &wcen, &wwidth);

    if (f != 'E' && f != 'N') {			// Empty or None values
        phu.filt_cen = wcen;
        phu.filt_start = wcen - wwidth / 2.0;
        phu.filt_end = wcen + wwidth / 2.0;
    } else {
	memset (phu.filt_card, 0, SZ_VAL);
	memset (phu.filter, 0, SZ_VAL);
    }

    /*  Inherit WCS keywords from PHU for older image files/
     */
    dl_strCard (fptr, "CTYPE1", phu.ctype1);
    dl_strCard (fptr, "CTYPE2", phu.ctype2);
    dl_strCard (fptr, "CUNIT1", phu.cunit1);
    dl_strCard (fptr, "CUNIT2", phu.cunit2);
    dl_strCard (fptr, "BUNIT",  phu.bunit);

    if (DEBUG) 
	fprintf (stderr, "PHU ctype = '%s','%s' cunit='%s','%s' bunit = '%s'\n",
	    phu.ctype1, phu.ctype2, phu.cunit1, phu.cunit2, phu.bunit);
}


/**
 *  DL_LOADEHU -- Load the EHU header structure.
 */
void
dl_loadEHU (fitsfile *fptr)
{
    int  stat = 0;

    memset (&ehu, 0, sizeof (struct EHU));
    if (phu_header == (char *) NULL)
        dl_loadPHU (fptr);

    if (((stat = dl_intCard (fptr, "BITPIX", &ehu.bitpix)) != OK) &&
        ((stat = dl_intCard (fptr, "ZBITPIX", &ehu.bitpix)) != OK)) {
	  if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get BITPIX keyw from EHU\n", stat);
	  ;
    }
    if (((stat = dl_intCard (fptr, "NAXIS1", &ehu.naxis1)) != OK) &&
        ((stat = dl_intCard (fptr, "ZNAXIS1", &ehu.naxis1)) != OK)) {
	  if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get NAXIS1 keyw from EHU\n", stat);
	  ;
    }
    if (((stat = dl_intCard (fptr, "NAXIS2", &ehu.naxis2)) != OK) &&
        ((stat = dl_intCard (fptr, "ZNAXIS2", &ehu.naxis2)) != OK)) {
	  if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get NAXIS2 keyw from EHU\n", stat);
	  ;
    }

    if ((stat = dl_strCard (fptr, "CTYPE1", ehu.ctype1)) != OK) {
	if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get CTYPE1 keyw from EHU\n", stat);
	strcpy (ehu.ctype1, phu.ctype1);
    }
    if ((stat = dl_strCard (fptr, "CTYPE2", ehu.ctype2)) != OK) {
	if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get CTYPE2 keyw from EHU\n", stat);
	strcpy (ehu.ctype2, phu.ctype2);
    }

    dl_strCard (fptr, "EXTNAME", ehu.extname);
    dl_intCard (fptr, "CCDNUM", &ehu.ccdnum);
   
    dl_strCard (fptr, "FILTER", ehu.filt_card); 
    memset (ehu.filter, 0, SZ_VAL);
    ehu.filter[0] = ehu.filt_card[0];			// DECam filters
    if (strchr ("ugrizyUGRIZY", (int)ehu.filter[0]) == NULL)
        strcpy (ehu.filter, ehu.filt_card);

    if ((stat = dl_strCard (fptr, "PRODTYPE", ehu.prodtype)) != OK) {
	if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get PRODTYPE keyw from EHU\n", stat);
	    
	strcpy (ehu.prodtype, phu.prodtype);
    }
    if ((stat = dl_strCard (fptr, "PROCTYPE", ehu.proctype)) != OK) {
	if (DEBUG > 2 || VERBOSE > 2)
	    fprintf (stderr, "Err %d: Cannot get PROCTYPE keyw from EHU\n", stat);
	strcpy (ehu.proctype, phu.proctype);
    }
}


/** 
 *  DL_MERGEHEADERS -- Merge the PHU and EHU headers to implement inheritence
 *  for the AST library to find the WCS.  The caller must free the returned
 *  pointer.
 */
char *
dl_mergeHeaders (char *phu_hdr, char *ehu_hdr)
{
    char *header = (char *) NULL;
    char *ehu = ehu_hdr;
    char *phu = phu_hdr;
    char *hc, *pc, *ec, keyw[8];
    int  i, j, pnkeys, enkeys;
    int  found = 0;


    /*  Trim the END keyword from the PHU.
     */
    if ((ec = strstr (phu, "END")))
        *ec = '\0';

    /*  Skip the XTENSION keyword in a MEF file.
     */
    if (strncmp (ehu, "XTENSION", 8) == 0)
	ehu += 80;

    header = calloc (1, 2 * (strlen (phu) + strlen (ehu) + 80));
    pnkeys = strlen (phu) / 80;
    enkeys = strlen (ehu) / 80;

    if (DEBUG) {
	fprintf (stderr, "pnkeys = %d   enkeys = %d  ", pnkeys, enkeys);
	fflush (stderr);
    }

    /*  Use the EHU as the default, copy everything but the END card.
     */
    if (ehu && strlen(ehu) > 80)
        memcpy (header, ehu, (int)(strlen(ehu)-80));

    /*  FIXME -- Should check for INHERIT value here....
     */
    /*  FIXME -- Need to remove INDEF values from the EHU....
     */

    /*  Not loop through PHU, copying over those keyword not present in the EHU.
     */
    pc = phu;
    hc = &header[strlen(header)];		// end of 'header' pointer
    for (i=0; i < pnkeys; i++) {
	memcpy (keyw, pc, 8);

	if (strncmp (keyw, "END", 3) == 0)
	    break;

	/* foreach keyw in EHU .... */
	found = 0;
	for (j=0, ec=ehu; j < enkeys; j++, ec+=80) {
	    if (ec[0] == keyw[0]) {		    // match first letter
		if (strncmp (ec, keyw, 8) == 0)	{   // PHU keyword not in EHU
		    found = 1;
		    break;
		}
	    }
	}
	if (!found) {
	    if (strstr (pc, "INDEF") == NULL)
	        memcpy (hc, pc, 80);
	    hc += 80;
	}
	pc += 80;
    }

    if (DEBUG) {
	fprintf (stderr, "nkeys = %d\n", (int) strlen (header) / 80);
	fprintf (stderr, "header = '%s'\n", header);
	fflush (stderr);
    }
    return (header);
}


/************************************************
 * UTILITY METHODS
 ************************************************/

/**
 *  DL_xxxCARD - Utilities to return card values into supplied variables. 
 *  Return OK if the card was found, the 'status' value on error.
 */

int
dl_strCard (fitsfile *fptr, char *keyw, char *value)
{
    char  card[FLEN_CARD], comment[FLEN_COMMENT], val[FLEN_VALUE];
    char  str[FLEN_VALUE], *ip, *op;
    int   status = 0;


    memset (card, 0, FLEN_CARD);
    memset (val, 0, FLEN_VALUE);
    memset (str, 0, FLEN_VALUE);
    memset (comment, 0, FLEN_COMMENT);

    //if (fits_read_card (fptr, keyw, card, &status))
    if (dl_getCard (fptr, keyw, card, &status))
	return (status);

    if (card[0]) {
	fits_parse_value (card, val, comment, &status);

        for (ip=val, op=str; *ip; ip++) {
	    if (*ip != '"' && *ip != '\'') {
	        *op = *ip;
	        op++;
	    }
        }

        /*  Remove quotes and trailing whitespace.
         */
        strncpy (value, str, strlen (str));
        for (ip=&value[strlen(value)-1]; isspace (*ip) && ip > value; ip--) {
	    *ip = '\0';
        }

        return (OK);
    }
    return (ERR);
}

int
dl_intCard (fitsfile *fptr, char *keyw, int *value)
{
    char  card[SZ_LINE], comment[FLEN_COMMENT], val[FLEN_VALUE];
    int   status = 0;

    memset (val, 0, FLEN_VALUE);
    //*value = 0;

    //if (fits_read_card (fptr, keyw, card, &status))
    if (dl_getCard (fptr, keyw, card, &status))
	return (status);

    if (card[0]) {
	fits_parse_value (card, val, comment, &status);
        *value = (int) dl_atoi (val);
        return (OK);
    }
    return (ERR);
}

int
dl_dblCard (fitsfile *fptr, char *keyw, double *value)
{
    char  card[SZ_LINE], comment[FLEN_COMMENT], val[FLEN_VALUE];
    int   status = 0;
    extern double atof(const char *s);

    memset (val, 0, FLEN_VALUE);
    //if (fits_read_card (fptr, keyw, card, &status))
    if (dl_getCard (fptr, keyw, card, &status))
	return (status);

    if (card[0]) {
	fits_parse_value (card, val, comment, &status);
        if (strstr (val, "INDEF") != (char *) NULL) 
            *value = (double) 0.0;
        else
            *value = (double) atof (val);
        return (OK);
    }
    fprintf(stderr, "dblCard returning error: %s\n", card);
    return (ERR);
}

int
dl_floatCard (fitsfile *fptr, char *keyw, float *value)
{
    double  dval = -0.999;

    if (dl_dblCard (fptr, keyw, &dval) == OK) {
        *value = (float) dval;
        return (OK);
    }
    return (ERR);
}

int
dl_dmsCard (fitsfile *fptr, char *keyw, double *value)
{
    char  val[FLEN_VALUE];

    //*value = 0.0;
    memset (val, 0, FLEN_VALUE);
    if (dl_strCard (fptr, keyw, val) == OK) {
        *value = dl_sexa (val);
        return (OK);
    } 
    return (ERR);
}

int
dl_hmsCard (fitsfile *fptr, char *keyw, double *value)
{
    //*value = 0.0;
    if (dl_dmsCard (fptr, keyw, value) == OK) {
        *value *= 15.0;
        return (OK);
    } 
    return (ERR);
}

/* ----------------------------------------------------------------------- */


int
dl_getCard (fitsfile *fptr, char *field, char *card, int *status)
{
    register int i = 0;
    int nmap = sizeof(keywCfg) / sizeof (struct keywMap);
    char val[80], line[80], *ip = NULL, *op = NULL;

			
    memset (val,  0, 80);
    memset (card, 0, 80);
    memset (line, 0, 80);

    // Loop over each of the keyword fields.
    for (i=0; i < nmap; i++) {

	if (strcasecmp (field, keywCfg[i].field) == 0) {
	    if (keywCfg[i].constant[0]) {
		// Config defines a constant value for the keyword.  Format
		// a FITS header card with the value.
		if (isdigit(keywCfg[i].constant[0])) {
		    sprintf (line, "%-8.8s= %s", field, keywCfg[i].constant);

		} else if (keywCfg[i].constant[0] == '"' ||
		    keywCfg[i].constant[0] == '\'') {
			for (ip = &keywCfg[i].constant[1], op = val;
			    *ip && *ip != '"' && *ip != '\''; ip++, op++) 
			        *op = *ip;
		        sprintf (line, "%-8.8s= '%s'", field, val);
		} else
		     sprintf (line, "%-8.8s= %s", field, keywCfg[i].constant);

		memset (card, ' ', 80);
		memcpy (card, line, strlen (line));
		return ( (*status = 0) );

	    } else {
		// Config has a keyword mapping so use the mapped keyword
		// and get the card from the header.
    		if (fits_read_card (fptr, keywCfg[i].keyw, card, status))
                    memset (card, 0, 80);
        	    return (*status);
	    }
	}
    }
    return (0);
}

char *
dl_getKeyw (char *field)
{
    register int i = 0;
    int nmap = sizeof(keywCfg) / sizeof (struct keywMap);

    for (i=0; i < nmap; i++) {
	if (strcasecmp (field, keywCfg[i].field) == 0)
	    return (keywCfg[i].constant[0] ? 
		keywCfg[i].constant : keywCfg[i].keyw);
    }
    return (NULL);
}


void
dl_keywSetCfg (char *line)
{
    char keyw[80], val[80], *ip = line, *op;
    int  i, nmap = sizeof(keywCfg) / sizeof (struct keywMap);
    int  is_str = 0, is_constant = 0;


    memset (keyw, 0, 80);
    memset (val, 0, 80);

    for (op=keyw; !isspace(*ip); ip++, op++)		// get keyword
	*op = *ip;
    while (isspace (*ip))				// skip whitespace
	ip++;

    is_str = (*ip == '"' || *ip == '\'');

    /*  Extract the value as a string, numeric, or keyword value.
    */
    if (is_str) {
        op = val;
        do {
            *op++ = *ip++;
        } while ( (*ip != '\'' && *ip != '"' && *ip != '\n') );
        *op = *ip;
        is_constant++;
    } else if (isdigit(*ip)) {
        for (op=val; *ip && *ip != '\n'; ip++, op++)
            *op = *ip;
        is_constant++;
    } else {
        for (op=val; *ip && *ip != '\n'; ip++, op++)
            *op = *ip;
    }

    for (i=0; i < nmap; i++) {
	if (strcasecmp (keyw, keywCfg[i].field) == 0) {
	    if (is_constant) {
		keywCfg[i].keyw[0] = '\0';
		strcpy (keywCfg[i].constant, val);
	    } else if (isalpha (val[0]))		// keyword mapping
		strcpy (keywCfg[i].keyw, val);
	    else 
		printf ("Unknown value '%s'\n", val);
	}
    }
}

void
dl_loadKeywConfig (char *fname)
{
    FILE *fd = (FILE *) NULL;
    char  *line = NULL, *ip = NULL;
    size_t len = 80;


    if ( (fd = fopen (fname, "r")) != (FILE *) NULL) {
	while (1) {
            line = NULL;
	    if ((int)getline (&line, &len, fd) < 0)  {
                if (line) free (line);
		break;
            }
	    if (*line == '#')				// skip comments
		continue;
            for (ip=line; *ip; ip++) {                  // inline comments too
                if (*ip == '#') {
                    *(--ip) = '\0';
                    break;
                }
            }
	    dl_keywSetCfg (line);			// set config value
            if (line) free (line);
	}
        fclose (fd);

    } else
	fprintf (stderr, "Error: cannot open keyword config '%s', errno=%d\n",
	    fname, errno);
}


/**  Debug routine
 **/
void
dl_printKeywConfig ()
{
    int  i, nmap = sizeof(keywCfg) / sizeof (struct keywMap);

    printf ("\nnmap = %d\n\n", (int)(sizeof(keywCfg) / sizeof(struct keywMap)));

    for (i=0; i < nmap; i++)
        printf ("%12s | %-12s | %s\n", 
            keywCfg[i].field, 
            keywCfg[i].keyw, 
            (keywCfg[i].constant[0] ?  keywCfg[i].constant : ""));
    printf ("\n===========================\n\n");
}
