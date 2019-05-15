/**
 *  VOSUTIL - Utility routines for the VOSAMP tools.
 *
 *  @file       vosUtil.c
 *  @author     Mike Fitzpatrick
 *  @date       6/03/11
 *
 *  @brief      Utility routines for the VOSAMP tools.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <time.h>

#include <netdb.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <curl/curl.h>
#ifdef OLD_CURL
#include <curl/types.h>
#endif
#include <curl/easy.h>


#include "dlApps.h"				/* dlApps interface	    */
#include "dlAppsP.h"				/* dlApps interface	    */


#define	SZ_MTYPE	64
#define	SZ_BUF		128
#define	SZ_HOSTIP	16
#define	SZ_APPNAME	16

#define	MAX_ROWS	256
#define	MAX_CLIENTS	256

#define DL_DEBUG  (getenv("DL_DBG")!=NULL||access("/tmp/DL_DBG",F_OK)==0)


char *dl_toURL (char *arg);
char *dl_optArg (char *arg);
char *dl_getLocalIP (void);
char *dl_getFName (char *path);

int   dl_isDir (char *path);
int   dl_isFITS (char *fname);
int   dl_isGZip (char *fname);
int  *dl_toIntArray (char *arg, int *nrows);
int   dl_urlType (char *url);
int   dl_getURL (char *url, char *fname);

int    dl_atoi (char *val);
long   dl_atol (char *val);
double dl_atof (char *val);





/****************************************************************************
 */


/**
 * DL_ERROR - Process an error and exit the program.
 */

extern	char *prog_name;

void
dl_error (int exit_code, char *error_message, char *tag)
{
    if (tag != NULL && strlen (tag) > 0)
	fprintf (stderr, "ERROR %s: %s (%s)\n", prog_name, error_message,
tag);
    else
	fprintf (stderr, "ERROR %s: %s\n", prog_name, error_message);
    fflush (stdout);

    //exit (exit_code);
}


/****************************************************************************
 */


/**
 *  DL_URLTYPE -- Determine the type of a URL parameter
 */
int
dl_urlType (char *url)
{
    if (strncasecmp (url, "http://127.0.0.1", 16) == 0)
	return (DL_LOCALURL);
    else if (strncasecmp (url, "file://", 7) == 0)
	return (DL_LOCALURI);
    else if (strncasecmp (url, "http://", 7) == 0)
	return (DL_REMOTE);
    else if (access (url, F_OK) == 0)
	return (DL_LOCALFILE);

    return (-1);
}


/**
 *  DL_GETFNAME -- Get a filename from a path or URL.
 */
char *
dl_getFName (char *path)
{
    static char fname[SZ_FNAME];
    static int  filenum = 0;

    memset (fname, 0, SZ_FNAME);
    if (access (path, R_OK) == 0) {
        int i, len = strlen (path);

        for (i=len-1; i >=0 && path[i] != '/'; i--) ;  /* get filename */
        strcpy (fname, &path[i+1]);
    }

    if (!fname[0])
	sprintf (fname, "vos%d_%03d", (int)getpid(), filenum++);

    return (fname);
}


/** 
 *  DL_GETURL -- Utility routine to do a simple URL download to the file.
 */
int 
dl_getURL (char *url, char *fname)
{
    int  stat = 0;
    char errBuf[SZ_LINE];
    FILE *fd;
    CURL *curl_handle;

	
    if (access (fname, F_OK) == 0)	/* see if file already exists	*/
	unlink (fname);


    /*  For the CURL operation to download the file.
     */
    curl_global_init (CURL_GLOBAL_ALL);     	/* init curl session	*/
    curl_handle = curl_easy_init ();

    /*  Open the output file.
     */
    if ((fd = fopen (fname, "wb")) == NULL) { 	
	fprintf (stderr, "Error: cannot open output file '%s'\n", fname);
        curl_easy_cleanup (curl_handle);
        return -1;
    }

    /*  Set cURL options
     */
    curl_easy_setopt (curl_handle, CURLOPT_URL, url);
    curl_easy_setopt (curl_handle, CURLOPT_NOPROGRESS, 1L);
    curl_easy_setopt (curl_handle, CURLOPT_WRITEDATA, fd);
    curl_easy_setopt (curl_handle, CURLOPT_ERRORBUFFER, errBuf);
    curl_easy_setopt (curl_handle, CURLOPT_FOLLOWLOCATION, 1);
    curl_easy_setopt (curl_handle, CURLOPT_FAILONERROR, 1);

    /*  Do the download.
     */
    if ((stat = curl_easy_perform (curl_handle)) != 0) {
	/*  Error in download, clean up.
	 */
	fprintf (stderr, "Error: can't download '%s' : %s\n", url, errBuf);
	unlink (fname);
        fclose (fd); 			    	/* close the file 	*/
        curl_easy_cleanup (curl_handle);    	/* cleanup curl stuff 	*/
	return (-1);
    }

    fflush (fd);
    fclose (fd); 			    	/* close the file 	*/
    curl_easy_cleanup (curl_handle); 	    	/* cleanup curl stuff 	*/

    return (1);
}


/**
 *  DL_OPTARG -- Input command arguments are allowed to be of the form
 *  'param=value', but the SAMP interface only wants the value string. 
 *  Skip past the '=' and return the value, or just return the value if
 *  there is no parameter name.
 */
char *
dl_optArg (char *arg)
{
    char *ip, first = (arg ? *arg : 0);

    if (!arg || !first) return ("");
    return ( ((ip = strchr (arg, (int) '=')) ? ++ip : arg) );
}


/**
 *  DL_TOURL -- Convert the argument to a URL suitable for a message.
 */
char *
dl_toURL (char *arg)
{
    /*  If we have an existing protocol simply return the argument.
     */
    if ((strncmp (arg, "http:", 5) == 0) ||
        (strncmp (arg, "file:", 5) == 0) ||
        (strncmp (arg, "ftp:", 4) == 0))
    	    return (arg);

    if (access (arg, F_OK) == 0) {
	static char buf[SZ_FNAME];

 	memset (buf, 0, SZ_FNAME);
	if (arg[0] != '/') {
	    char  cwd[SZ_FNAME];

 	    memset (cwd, 0, SZ_FNAME);
	    getcwd (cwd, (unsigned long) SZ_FNAME);
	    sprintf (buf, "file://%s/%s", cwd, arg);
	} else
	    sprintf (buf, "file://%s", arg);

	return (buf);
    }

    return (arg);
}


/**
 *  DL_TOINTARRAY -- Convert a range string to an unpacked array of ints.
#define MAX_RANGES	 256
 */

int *
dl_toIntArray (char *arg, int *nrows)
{
    int  i, val, nvalues;
    static int  ranges[MAX_RANGES], values[MAX_ROWS];
    extern int  dl_decodeRanges(), get_next_number();


    memset (values, 0, (sizeof(int) * MAX_ROWS));
    memset (ranges, 0, (sizeof(int) * MAX_RANGES));

    if (dl_decodeRanges (arg, ranges, MAX_RANGES, &nvalues) < 0)
        fprintf (stderr, "Error decoding range string.\n");

    for (i=0, val=0; (val = get_next_number (ranges, val)) > 0; i++ )
	values[i] = val;

    *nrows = nvalues;
    return (values);
}


/**
 *  DL_STRSUB -- Do a string subsitution.
 */
int
dl_strsub (char *in, char *from, char *to, char *outstr, int maxch)
{
    int   flen = strlen (from);
    int   nsub = 0;
    char  *ip, *op;

    if (!from || !to)
        return (0);

    for (ip=in, op=outstr; *ip; ip++) {
        if (! *ip || (ip - in) > maxch)
            break;
        if (*ip == '$') {
            /* Start of a macro.
             */
            if (strncasecmp (ip, from, flen) == 0) {
                /* Our macro, do the substitution.
                 */
                char *tp = to;

                ip += flen - 1;         /* skip the input macro         */
                while (*tp)             /* copy replacement string      */
                    *op++ = *tp++;
                nsub++;
            } else {
                /* Not our macro, just pass it through.
                 */
                *op++ = *ip;
            }
        } else {
            *op++ = *ip;
        }
    }
    *op = '\0';

    return (nsub);
}


/**
 *  DL_ISDIR -- Check whether a path represents a directory.
 *
 *  @brief  Check whether a path represents a directory.
 *  @fn     int dl_isDir (char *path)
 *
 *  @param  path        pathname to be checked
 *  @return             1 (one) if path is a directory, 0 (zero) otherise
 */
int
dl_isDir (char *path)
{
    struct stat stb;
    int     res = -1;

    if ((res = stat(path, &stb)) == 0)
        return ((int) S_ISDIR (stb.st_mode));
    return (0);
}


/**
 *  DL_ISFITS -- Test a file to see if it is a simple FITS file.
 */
int 
dl_isFITS (char *fname)
{
    register FILE *fp;
    int value = 0;
    char keyw[8], val;

    memset (keyw, 0, 8);
    if ((fp = fopen (fname, "r"))) {
        fscanf (fp, "%6s = %c", keyw, &val);
        if (strcmp ("SIMPLE", keyw) == 0 && val == 'T')
            value = 1;
        fclose (fp);
    }
    return value;
}


/**
 *  DL_ISGZIP -- Test a file to see if it is GZip compressed.
 */
int 
dl_isGZip (char *fname)
{
    int   fp, value = 0;
    unsigned short sval = 0;

    if ((fp = open (fname, O_RDONLY))) {
        read (fp, &sval, 2);
        if (sval == 35615)              // ushort value of "\037\213"
            value = 1;
        close (fp);
    }
    return value;
}


/**
 *  DL_ATOI -- System atoi() with lexical argument checking.
 */
int
dl_atoi (char *val)
{
    char *ip;

    for (ip = val; *ip; ip++) {
        if (isalpha ((int) *ip)) {
            fprintf (stderr, "Warning: value '%s' is not an integer\n", val);
            break;
        }
    }
    return (atoi (val));
}


/**
 *  DL_ATOL -- System atol() with lexical argument checking.
 */
long
dl_atol (char *val)
{
    char *ip;

    for (ip = val; *ip; ip++) {
        if (isalpha ((int) *ip) && *ip != '-' && *ip != '+') {
            fprintf (stderr, "Warning: value '%s' is not an integer\n", val);
            break;
        }
    }
    return (atol (val));
}


/**
 *  DL_ATOF -- System atoi() with lexical argument checking.
 */
double
dl_atof (char *val)
{
    char *ip, c;

    for (ip = val; *ip; ip++) {
        c = *ip;
        if (! (tolower(c) == 'e' || tolower(c) == 'd' || isspace(c) ||
              (c == '-' || c == '+' || c == '.') || isdigit(c))) {
                fprintf (stderr, 
                    "Warning: value '%s' is not a floating point value\n", val);
            break;
        }
    }
    return (atof (val));
}




/**
 * DL_DEXA - Convert a sexagesimal string to decimal.
 */
double
dl_sexa (char *str)
{
    int     n, sign;
    int     hr, minutes;
    double  sec, val;
    char   *s = str;


    if (!str || (*str == '\0'))
        return (0.0);

    /*  Skip leading whitespace and garbage.
     */
    while (isspace (*s) || ((!isdigit (*s)) && (*s != '-' && *s != '+')))
        s++;
    sign = (*s == '-') ? (s++, -1) : 1; /* get the sign                 */

    minutes = 0;
    sec = 0.;
    n = sscanf (s, "%d:%d:%lg", &hr, &minutes, &sec);
    if (n < 1 || minutes < 0 || sec < 0)
        val = -999.0;
    else
        /*  Beware: Evaluation here can produce roundoff errors!
        */
        val = sign * (hr + ((double)minutes)/60. + sec/3600.);

    return (val);
}

