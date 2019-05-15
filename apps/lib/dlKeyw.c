/** 
 *  DLKEYW -- Keyword Mapping Utility procedures.
 *
 *  @file       dlKeyw.c
 *  @author     Mike Fitzpatrick
 *  @date       10/3/17
 *
 *  @brief      Keyword Mapping Utility procedures.
 */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
/*
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "dlApps.h"
*/

#define	SZ_VAL		12
#define	SZ_CONSTANT	64

typedef struct keywMap {
    char field[SZ_VAL];
    char keyw[SZ_VAL];
    char constant[SZ_CONSTANT];
} keywCfg[] = {
    { "airmass",  "AIRMASS",  "" },
    { "bitpix",   "BITPIX",   "" },
    { "bunit",    "BUNIT",    "" },
    { "ccdnum",   "CCDNUM",   "" },
    { "ctype1",   "CTYPE1",   "" },
    { "ctype2",   "CTYPE2",   "" },
    { "cunit1",   "CUNIT1",   "" },
    { "cunit2",   "CUNIT2",   "" },
    { "date-obs", "DATE-OBS", "" },
    { "dateobs1", "DATEOBS1", "" },
    { "dec",      "DEC",      "" },
    { "dimmsee",  "DIMMSEE",  "" },
    { "elliptic", "ELLIPTIC", "" },
    { "expnum",   "EXPNUM",   "" },
    { "exptime",  "EXPTIME",  "" },
    { "extname",  "EXTNAME",  "" },
    { "extnum",   "EXTNUM",   "" },
    { "filter",   "FILTER",   "" },
    { "fwhm",     "FWHM",     "" },
    { "ha",       "HA",       "" },
    { "instrume", "INSTRUME", "" },
    { "magzero",  "MAGZERO",  "" },
    { "mjd-obs",  "MJD-OBS",  "" },
    { "moonangl", "MOONANGL", "" },
    { "naxis1",   "NAXIS1",   "" },
    { "naxis2",   "NAXIS2",   "" },
    { "object",   "OBJECT",   "" },
    { "obsid",    "OBSID",    "" },
    { "obstype",  "OBSTYPE",  "" },
    { "photflag", "PHOTFLAG", "" },
    { "plver",    "PLVER",    "" },
    { "proctype", "PROCTYPE", "" },
    { "prodtype", "PRODTYPE", "" },
    { "program",  "PROGRAM",  "" },
    { "propid",   "PROPID",   "" },
    { "proposer", "PROPOSER", "" },
    { "ra",       "RA",       "" },
    { "telescop", "TELESCOP", "" },
    { "zbitpix",  "ZBITPIX",  "" },
    { "zd",       "ZD",       "" },
    { "znaxis1",  "ZNAXIS1",  "" },
    { "znaxis2",  "ZNAXIS2",  "" }
};


char *
dl_getKeyw (char *field)
{
    register int i = 0;
    keywMap *map;

    for (&keywCfg[0]; *map; map++) {
	if (strcasecmp (field, map.field) == 0)
	    return (*map.constant ? map.constant : map.keyw);
    }
    return (NULL)
}



