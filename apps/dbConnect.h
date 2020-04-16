/**
 *  DBCONNECT.H -- Connection string for the NSA database to resolve
 *  'fileRef' to paths in the mass store.
 *
 *  NOTE -- This file SHOULD NOT be included in any public repository
 */

#ifdef USE_NSA

static char *conn_info =
    "dbname=ICAT host=dsan3.tuc.noao.edu user=pipeline password=Mosaic_DHS";

#else

static char *conn_info =
    "dbname=tapdb host=db02.datalab.noao.edu user=dlquery password=datalab";

#endif
