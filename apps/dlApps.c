/**
 *  DLAPPS.C -- Generic task main() for the Data Lab Package applications.
 *
 *  @file       dlpps.c
 *  @author     Mike Fitzpatrick
 *  @date       1/03/16
 *
 *  @brief      Generic task main() for the Data Lab Package applications.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*  Run as a detached child process?
*/
#define	RUN_DETACHED		0


#include "dlApps.h"

Task  *app    = (Task *) NULL;


/*  Task entry-point declarations.
 */
extern int  dlcutout (int argc, char **argv, size_t *len, void **result);
extern int  dlmeta (int argc, char **argv, size_t *len, void **result);
extern int  dlphead (int argc, char **argv, size_t *len, void **result);

extern int  fits2csv (int argc, char **argv, size_t *len, void **result);



/*  Task application table.  Applications must be declared here to be found
 *  and used with the generic main().
 */

Task dlApps[] = {
// { "dlcutout",        dlcutout    },          /* Data Lab Apps      	      */
// { "dlmeta",          dlmeta      },
// { "dlphead",         dlphead     },
   { "fits2csv",        fits2csv    },
   { NULL,              NULL        }
};



/**
 *  Program main().
 */
int
main (int argc, char *argv[])
{
    char  *task   = (char *) NULL;
    int    status = 0;
    size_t reslen = 0;
    void  *result = (void *) NULL;


    /*  Get the task name minus any leading path prefix.
     */
    task = argv[0] + strlen (argv[0]) - 1;
    while (task > argv[0]) {
	if ( *(task - 1) != '/')
	    task--;
	else
	    break;
    }

    /*  Loop through the application table.  If the names match, call
     *  the entry-point function.
     */
    if (getenv ("DLAPP_CONNECTED") || (! RUN_DETACHED) ) {
	for (app = dlApps; app->name; app++) {
	    if (strcmp (app->name, task) == 0) {
	        status = (*app->func) (argc, argv, &reslen, &result);
		break;
	    }
	}
    } else {
	/*  Run as a detached sub-process.
	 */
        status = dl_runTask (task, dlApps, argc, argv, &reslen, &result);
    }


    if (status && result)		/* print the error message */
	fprintf (stderr, "%s\n", (char *) result);
    if (reslen && result)		/* free result pointer	   */
	free (result);

    return (status);
}
