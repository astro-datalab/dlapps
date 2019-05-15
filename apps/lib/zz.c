/** 
 *  DLKEYW -- Keyword Mapping Utility procedures.
 *
 *  @file       dlKeyw.c
 *  @author     Mike Fitzpatrick
 *  @date       10/3/17
 *
 *  @brief      Keyword Mapping Utility procedures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
/*
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
*/

#include "fitsio.h"
#include "dlApps.h"

#define	SZ_VAL		12
#define	SZ_CONSTANT	64

struct keywMap {
    char field[SZ_VAL];
    char keyw[SZ_VAL];
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


char *
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
		    sprintf (line, "%-8.8s = %s", field, keywCfg[i].constant);

		} else if (keywCfg[i].constant[0] == '"' ||
		    keywCfg[i].constant[0] == '\'') {
			for (ip = &keywCfg[i].constant[1], op = val;
			    *ip && *ip != '"' && *ip != '\''; ip++, op++) 
			        *op = *ip;
		        sprintf (line, "%-8.8s = '%s'", field, val);
		} else
		     sprintf (line, "%-8.8s = %s", field, keywCfg[i].constant);

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
    return (NULL);
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
    char  *line = NULL;
    int    len = 80;


    if ( (fd = fopen (fname, "r")) != (FILE *) NULL) {
	while (1) {
            line = NULL;
	    if (getline (&line, &len, fd) < 0) 
		break;
	    if (*line == '#')				// skip comments
		continue;
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

    printf ("\nnmap = %d\n\n", sizeof(keywCfg) / sizeof (struct keywMap));

    for (i=0; i < nmap; i++)
        printf ("%12s | %-12s | %s\n", 
            keywCfg[i].field, 
            keywCfg[i].keyw, 
            (keywCfg[i].constant[0] ?  keywCfg[i].constant : ""));
    printf ("\n===========================\n\n");
}


int
main (int argc, char *argv[])
{
    fitsfile *fptr;
    int status = 0;
    char card[80];


    memset (card, ' ', 80);

    dl_loadKeywConfig ("keyw.dat");
    dl_printKeywConfig ();

    if (fits_open_file (&fptr, "zztest.fits", READONLY, &status) == 0) {

        if (! dl_getCard (fptr, "ra", card, &status))
            printf ("%8s : %-32.32s : %d\n", "ra", card, status);
        else
            printf ("%8s : %-32.32s : %d  NO KEYWORD\n", "ra", card, status);

        if (! dl_getCard (fptr, "dec", card, &status))
            printf ("%8s : %-32.32s : %d\n", "dec", card, status);
        else
            printf ("%8s : %-32.32s : %d  NO KEYWORD\n", "dec", card, status);

        if (! dl_getCard (fptr, "zd", card, &status))
            printf ("%8s : %-32.32s : %d\n", "zd", card, status);
        else
            printf ("%8s : %-32.32s : %d  NO KEYWORD\n", "zd", card, status);

        if (! dl_getCard (fptr, "mjd-obs", card, &status))
            printf ("%8s : %-32.32s : %d\n", "mjd-obs", card, status);
        else
            printf ("%8s : %-32.32s : %d  NO KEYWORD\n", "mjd-obs", card, status);

        if (! dl_getCard (fptr, "expnum", card, &status))
            printf ("%8s : %-32.32s : %d\n", "expnum", card, status);
        else
            printf ("%8s : %-32.32s : %d  NO KEYWORD\n", "expnum", card, status);

        if (! dl_getCard (fptr, "telescop", card, &status))
            printf ("%8s : %-32.32s : %d\n", "telescop", card, status);
        else
            printf ("%8s : %-32.32s : %d  NO KEYWORD\n", "telescop", card, status);
    }
}
