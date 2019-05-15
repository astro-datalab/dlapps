Generic computational tasks.

    dlmeta.c
	Extract Image table metadata from a FITS image.

    dltaskd.c
	Tasking daemon.  Used to asynchronously subit, run, and manage
	multiple concurrent tasks on a remote server (or the localhost).

    dlcutout.c
	Used by the SIAV2 service to both compute (plan) and generate
	cutouts of images.
        SimpleLink:  cc dlcutout.c -lcfitsio `ast_link` -o dlcutout
	Requires: cfitsio, Starlink AST library

    dlphead.c
    	A simple utility used to print image header information, to
	examine images used as input to vocutout.c

