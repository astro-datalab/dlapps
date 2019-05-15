/*
 *   Compile as:
 *
 *      cc -o zzarch -I/usr/pgsql-9.3/include zzarch.c -L/usr/pgsql-9.3/lib -lpq
 */

#include <stdio.h>
#include <libpq-fe.h>
#include <string.h>
    
char *conn_info = 
	"dbname=ICAT host=dsan3.tuc.noao.edu user=pipeline password=Mosaic_DHS";

char *dl_getFilePath (char *fname);


int main() 
{
    char *fname = "c4d_130211_024118_ooi_i_v2.fits.fz", *s;

    if ((s = dl_getFilePath (fname)))
	printf ("%s\n", s);
    else
	printf ("Cannot resolve path.\n");
}


char *
dl_getFilePath (char *fname)
{
    PGconn      *conn;
    PGresult    *res;
    int         rec_count;
    int         row;
    static char path[256];
    static char query[256];


    memset (path, 0, 256);			// initialize
    memset (query, 0, 256);

    /*  Format the query string, substituting the returned path for the one
     *  used on the machine.
     */
    sprintf (query, 
        "select replace(data_path,'Volumes','nsa') from r_data_main where data_name='%s'", fname);

    conn = PQconnectdb (conn_info);		// create the DB connection
    if (PQstatus(conn) == CONNECTION_BAD)
        return (NULL);

    res = PQexec(conn, query);			// execute the query
    if (PQresultStatus(res) != PGRES_TUPLES_OK)
        return (NULL);

    strcpy (path, PQgetvalue(res, row, 0));	// fetch result

    PQclear(res);
    PQfinish(conn);

    return (path);
}
