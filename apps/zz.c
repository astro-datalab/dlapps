#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

void parse_query_string (char *qs, char *fileRef, int *extn, 
	double *ra, double *dec, double *sz1, double *sz2, int *preview);


int main (int argc, char *argv[])
{
    char   *s = NULL;
    char    fileRef[180];
    int     extn = -1, preview = -1;
    double  ra = 0.0, dec = 0.0, sz1 = 0.0, sz2 = 0.0;


    memset (fileRef, 0, 180);

    if ((s = getenv ("QUERY_STRING"))) 
	parse_query_string (s, fileRef, &extn, &ra, &dec, &sz1, &sz2, &preview);
    printf ("0: fileRef='%s'  extn=%d  pos=(%g,%g)  size=(%g,%g)  preview=%d\n",
	fileRef, extn, ra, dec, sz1, sz2, preview);

    parse_query_string ("fileRef=/tmp/foo.fits&extn=3", 
	fileRef, &extn, &ra, &dec, &sz1, &sz2, &preview);
    printf ("1: fileRef='%s'  extn=%d  pos=(%g,%g)  size=(%g,%g)  preview=%d\n",
	fileRef, extn, ra, dec, sz1, sz2, preview);

    parse_query_string ("fileRef=/tmp/foo.fits&POS=12.34,56.78", 
	fileRef, &extn, &ra, &dec, &sz1, &sz2, &preview);
    printf ("2: fileRef='%s'  extn=%d  pos=(%g,%g)  size=(%g,%g)  preview=%d\n",
	fileRef, extn, ra, dec, sz1, sz2, preview);

    parse_query_string ("fileRef=/tmp/foo.fits&POS=12.34,56.78&SIZE=0.2", 
	fileRef, &extn, &ra, &dec, &sz1, &sz2, &preview);
    printf ("3: fileRef='%s'  extn=%d  pos=(%g,%g)  size=(%g,%g)  preview=%d\n",
	fileRef, extn, ra, dec, sz1, sz2, preview);

    parse_query_string ("fileRef=/tmp/foo.fits&POS=12.34,56.78&SIZE=0.2,0.1", 
	fileRef, &extn, &ra, &dec, &sz1, &sz2, &preview);
    printf ("4: fileRef='%s'  extn=%d  pos=(%g,%g)  size=(%g,%g)  preview=%d\n",
	fileRef, extn, ra, dec, sz1, sz2, preview);

    parse_query_string ("fileRef=/tmp/foo.fits&extn=2&preview=true", 
	fileRef, &extn, &ra, &dec, &sz1, &sz2, &preview);
    printf ("5: fileRef='%s'  extn=%d  pos=(%g,%g)  size=(%g,%g)  preview=%d\n",
	fileRef, extn, ra, dec, sz1, sz2, preview);

}


#define	SZ_ARG		32

void
parse_query_string (char *qs, char *fileRef, int *extn, 
	double *ra, double *dec, double *sz1, double *sz2, int *preview)
{
    char *ip = qs, *kp, *vp;
    char  keyw[SZ_ARG], val[SZ_ARG];

    //if (debug)
    //	printf ("\nqstring = '%s'\n", qs);

    for (ip=qs; *ip; ) {
	memset (keyw, 0, SZ_ARG);
	memset (val, 0, SZ_ARG);

	for (kp=keyw; *ip && *ip != '='; ip++)
	    *kp++ = *ip;
	if (*ip == '=') ip++; else break;

	for (vp=val; *ip && *ip != '&'; ip++)
	    *vp++ = *ip;
	if (*ip == '&') ip++;

	//if (debug)
	//  printf ("  keyword = '%s' val = '%s'\n", keyw, val);

	if (strncasecmp (keyw, "fileRef", 3) == 0) {
	    strcpy (fileRef, val);
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
	} else if (strncasecmp (keyw, "preview", 3) == 0) {
	    *preview = ( (val[0] == 't' || val[0] == 'T') ? 1 : 0 );
	} else 
	    ;				// ignore unknown arguments

    }
}
