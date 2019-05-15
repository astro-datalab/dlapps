/**
 *  VOAPPPROTO.H -- VOApps prototype headers.
 */


/**
 *  VOACLIST.C -- Procedures for handling the AccessList of images/data
 */
void dl_addToAclist (char *url, char *outfile);
void dl_freeAclist (void);
void dl_procAclist (void);


/**
 *  VODALUTIL.C  -- Utility procedures for the DAL interface worker procedures.
*/
int   dl_extractResults (char *result, char delim, svcParams *pars);
char *dl_openExFile (svcParams *pars, int nrows, char *extn, FILE **ofd);
char *dl_getOFName (svcParams *pars, char *extn, int pid);
char *dl_getOFIndex (svcParams *pars, char *extn, int pid);
int   dl_countResults (char *result);
void  dl_dalExit (int code, int count);
void  dl_printHdr (int fd, svcParams *pars);
void  dl_printCountHdr (void);
int   dl_printCount (Query query, svcParams *pars, int *res_count);
void  dl_printCountLine (int nrec, svcParams *pars);

char  dl_svcTypeCode (int type);
char *dl_getExtn (void);

char *dl_procTimestamp (void);
void  dl_concat (void);


/**
 *  VODALUTIL.C  -- Utility procedures for the DAL interface worker procedures.
 */
void  dl_initHTML (FILE *fd, svcParams *pars);
void  dl_printHTMLRow (FILE *fd, char *line, int isHdr, int rownum);
void  dl_closeHTML (FILE *fd);


/**
 *  VOINV.C -- VOInventory service routines.
 */
char *dl_doInventory (void);
char *dl_execInv (double ra, double dec, double radius, char *sources,
    		    char *resources, char *id, char *rettype, FILE *outfile);


/**
 *  VOKML.C  -- Utility procedures for writing Google KML files.
 */
void  dl_initKML (FILE *fd, svcParams *pars);
void  dl_printKMLPlacemark (FILE *fd, char *id, double ra, double dec, 
	char *line, char *acref, svcParams *pars);
void  dl_mkPlaceDescr (FILE *fd, char *line, char *acref, svcParams *pars);
void  dl_closeKML (FILE *fd);
void  dl_concatKML (char *fname);
void  dl_concatKMLByObject (FILE *fd);
void  dl_concatKMLByService (FILE *fd);
char *dl_getSName (char *root);
char *dl_getOName (char *root);
int   dl_copyKMLFile (char *root, char *name, FILE *fd);
void  dl_cleanKML (void);


/**
 *  VOLOG.C -- VOApps logging interface.  
 */
void  dl_appLog (FILE *fd, char *format, ...);
void  dl_encodeString (char *buf, char *format, va_list *argp);
char *dl_doarg (va_list **argp, int dtype);
char *dl_logtime (void);


/**
 *  VOARGS.C -- Procedures for commandline argument handling.  We also do
 */
int   dl_parseObjectList (char *list, int isCmdLine);
void  dl_freeObjectList (void);
int   dl_countObjectList (void);
int   dl_printObjectList (FILE *fd);
void  dl_readObjFile (char *fname);


/**
 *  VOPARAMS.C -- Interface to manage cmdline options or library parameters.
 */
char **dl_paramInit (int argc, char *argv[]);
int    dl_paramNext (char *opts, struct option long_opts[], int argc, 
		char *argv[], char *optval, int *posindex);
void   dl_paramFree (int argc, char *argv[]);

 
/**
 *  VORANGES -- Simple range-specification package to decode lists of numbers
 *  or ranges of the form:
 */
int dl_decodeRanges (range_string, ranges, max_ranges, nvalues);
int get_next_number (int ranges[], int number);
int is_in_range (int ranges[], int number);


/**
 *  VOSCS.C -- Worker procedure to query a Simple Cone Search service.
 */
int dl_callConeSvc (svcParams *pars);


/**
 *  VOSIAP.C -- Worker procedure to make a query to an SIAP service.
 */
int dl_callSiapSvc (svcParams *pars);
char *dl_validateFile (char *fname);


/**
 *  VOSSAP.C -- Worker procedure to make a query to an SSAP service.
 */
int dl_callSsapSvc (svcParams *pars);


/**
 *  VOSVC.C -- Procedures for commandline argument and DAL service handling.
 */
int   dl_parseServiceList (char *list, int dalOnly);
void  dl_freeServiceList (void);
void  dl_resetServiceCounters (void);
void  dl_addToSvcList (char *name, char *ident, char *url, char *type, 
				char *title);
int   dl_countServiceList (void);
int   dl_printServiceList (FILE *fd);
int   dl_printServiceVOTable (FILE *fd);
void  dl_readSvcFile (char *fname, int dalOnly);


/**
 *  VOTASK.C -- Utilities to run a VOApps task as a connected subprocess.
 */
int dl_runTask (char *method, Task *apps, int argc, char **argv, size_t *len, 
		void **result);
int  dl_taskTest (Task self, char *arg, ...);

int dl_setResultFromFile (char *fname, size_t *len, void **data);
int dl_setResultFromInt (int value, size_t *len, void **data);
int dl_setResultFromReal (float value, size_t *len, void **data);
int dl_setResultFromString (char *str, size_t *len, void **data);


/**
 *  VOUTIL.C -- Utility procedures for the VO-CLI tasks.
 */
int   dl_regResolver (char *term, char *svctype, char *bpass, char *subject, 
		 char *fields, int index, int exact, int dalOnly, char **res);
int   dl_regSearch (char **ids, int nids, char *svctype, char *bpass, 
		   char *subject, int orValues, int votable, FILE *dl_fd,
		   int dalOnly, int sortRes, int terse);
void  pretty_print (char *result, int nresults);

void  ppResSummary (char *result, int nresults);
void  ppMultiLine (char *result, int poffset, int pwidth, int maxchars);

void  pretty_print_table (char *result, int nresults, char *fields);
char *dl_parseSvcType (char *svctype, int exact);

char *dl_parseBandpass (char *bpass);
char *dl_parseSubject (char *subject);

char *dl_urlFname (char *url);
void  dl_printAttrs (char *fname, Query query, char *ident);
void  dl_printRegVOTableHdr (FILE *fd);

void  dl_printRegVOTableRec (FILE *fd, RegResult resource, int recnum);
void  dl_printRegVOTableTail (FILE *fd);

char *xmlEncode (char *in);
char *dl_getline (FILE *fd);
char *dl_normalizeCoord (char *coord);
char *dl_normalize (char *str);
char *dl_toURL (char *arg);
void  dl_setArg (char **argv, int *argc, char *value);

int   isVOTable (char *fname); 			/* utility functions  */
int   isSexagesimal (char *str);
int   isDecimal (char *str);
float sexa (char *s);
char *toSexa (double pos);
char *toSexaTime (int nsec);
char *dl_mktemp (char *root);
char *dl_copyStdin (void);
void  dl_skipHdr (FILE *fd);

char *dl_getTableCol (char *line, int col, int span);
int   dl_isNumericField (handle_t field);
int   dl_fileType (char *fname);
int   dl_sum32 (char *str);
int   strdic (char *in_str, char *out_str, int maxchars, char *dict);


/**
**  VOXML.C  -- Utility procedures for writing XML files, i.e. the raw
*/

void  dl_concatXML (char *fname);
int   dl_copyXMLFile (char *root, char *name, FILE *fd);
void  dl_cleanXML (void);


/**
 *  VOSUTIL.C - Utility routines for the VOSAMP tools.
 */
int   vos_urlType (char *url);
char *vos_getFName (char *path);
char *vos_typeName (int type);
int   vos_getURL (char *url, char *fname);
char *vos_optArg (char *arg);
char *vos_toURL (char *arg);
int  *vos_toIntArray (char *arg, int *nrows);

int   vos_openServerSocket (int port);
int   vos_openClientSocket (char *host, int port, int retry);
int   vos_testClientSocket (char *host, int port);
int   vos_sockReadHdr (int fd, int *len, char *name, int *type, int *mode);
int   vos_sockWriteHdr (int fd, int len, char *name, int type, int mode, 
		char *to);
void  vos_sockPrintHdr (char *msg, int fd);
int   vos_sockRead (int fd, void *vptr, int nbytes);
int   vos_sockWrite (int fd, void *vptr, int nbytes);
int   vos_fileRead (int fd, void *vptr, int nbytes);
int   vos_fileWrite (int fd, void *vptr, int nbytes);
void  vos_setNonBlock (int sock);
char *vos_getLocalIP (void);

struct hostent *vos_getHostByName (char *name);
struct hostent *vos_dupHostent (struct hostent *hentry);

int   vos_strsub (char *in, char *from, char *to, char *outstr, int maxch);
