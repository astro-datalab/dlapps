#///////////////////////////////////////////////////////////////////////////////
#//
#//  Makefile for the Data Lab Applications
#//
#///////////////////////////////////////////////////////////////////////////////

# primary dependencies

NAME       	 = DataLab
VERSION    	 = 1.0
PLATFORM  	:= $(shell uname -s)
PLMACH  	:= $(shell uname -m)
HERE      	:= $(shell /bin/pwd)

BINDIR  	:= $(shell ./install_env --bindir)
LIBDIR  	:= $(shell ./install_env --libdir)
INCDIR  	:= $(shell ./install_env --incdir)
MANDIR  	:= $(shell ./install_env --mandir)



# includes, flags and libraries
CC 	  	= gcc
CINCS  	  	= -I$(INCDIR) -I./
ifeq  ($(PLATFORM), "Darwin")
    ifeq  ($(PLATFORM), "x86_64")
        CARCH	= -m64 -mmacosx-version-min=10.7
    else
        CARCH	= -arch i386 -arch ppc -m32 -mmacosx-version-min=10.7
    endif
else
    CARCH	= 
endif

CFLAGS 		= -g -Wall $(CARCH) -D$(PLATFORM) $(CINCS) -L./
LIBS		= -lm -lc -lpthread


all:
	(cd ast-8.6.2   ; \
	     ./configure --prefix=$(HERE) --without-stardocs ;\
	     make all  ; make install ; make clean ; \
	     /bin/rm -rf help man manifests news share )
	(cd cfitsio     ; \
	     ./configure --prefix=$(HERE) ; \
	     make all  ; make install ; make clean)
	(cd apps        ; make all  ; make install)
	/bin/rm -rf bin/curl* lib/pkgconfig
	/bin/rm -rf lib/*.dylib lib/*.la lib/*.so

libs:
	(cd ast-8.6.2   ; make all  ; make install)
	(cd cfitsio     ; make all  ; make install)
	(cd apps        ; make libs)

apps:
	(cd apps        ; make all)

examples:

install:
	(cd ast-8.6.2   ; make all  ; make install)
	(cd cfitsio     ; make all  ; make install)
	(cd apps        ; make all  ; make install)
	/bin/rm -rf bin/curl* lib/pkgconfig
	/bin/rm -rf lib/*.dylib lib/*.la lib/*.so
	#cp -rp bin/*      $(BINDIR)
	#cp -rp lib/*      $(LIBDIR)
	#cp -rp include/*  $(INCDIR)
	#cp -rp doc/*.man  $(MANDIR)

clean:
	(cd ast-8.6.2   ; make clean)
	(cd cfitsio     ; make clean)
	(cd apps        ; make clean)
	/bin/rm -rf bin/* lib/* include/* *spool* */*spool* */*/*spool*
	/bin/rm -rf */config.log */*/config.log */*/*/config.log
	/bin/rm -rf */pkgconfig



###############################################################################
# Leave this stuff alone.
###############################################################################

$(STATICLIB): $(SRCS:%.c=Static/%.o)
	/usr/bin/ar rv $@ $?
Static/%.o: %.c $(INCS)
	/usr/bin/gcc $(CINCS) $(CFLAGS) -c $< -o $@
Static:
	/bin/mkdir $@
	chmod 777 $@

$(SHAREDLIB): $(SRCS:%.c=Shared/%.o)
	/usr/bin/ld -shared -o $@ $? -lc -ldl
Shared/%.o: %.c $(INCS)
	/usr/bin/gcc $(CINCS) $(CFLAGS) -fpic -shared -c $< -o $@
Shared:
	/bin/mkdir $@
	chmod 777 $@
