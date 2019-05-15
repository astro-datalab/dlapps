/**
 *  DLAPPS.I -- SWIG Interface definition file.
 *
 *  @file       dlApps.i
 *  @author     Mike Fitzpatrick
 *  @date       4/03/16
 *
 *  @brief       SWIG Interface definition file.
 */

%module voapps
%{

/*  Cutout Tool
 */
extern int  dl_cutout (int argc, char **argv, size_t *len, void **result);

/*  Metadata Collection Tool
 */
extern int  dl_meta (int argc, char **argv, size_t *len, void **result);

%}


/*  Cutout Tool
 */
extern int  dl_cutout (int argc, char **argv, size_t *len, void **result);

/*  Metadata Collection Tool
 */
extern int  dl_meta (int argc, char **argv, size_t *len, void **result);


