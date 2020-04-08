#///////////////////////////////////////////////////////////////////////////////
#//
#//  Makefile for the Data Lab Package Applications
#//
#///////////////////////////////////////////////////////////////////////////////

# primary dependencies

NAME       	= DLApps
VERSION    	= 1.0
PLATFORM        := $(shell uname -s)
PLMACH          := $(shell uname -m)
HERE            := $(shell /bin/pwd)
BINDIR          := ../bin/
LIBDIR          := ../lib/
INCDIR          := ../include/


# secondary dependencies
LIBBASE     = lib$(NAME)
STATICLIB   = $(HERE)/$(LIBBASE).a
SHAREDLIB   = $(HERE)/$(LIBBASE).so.$(VERSION)


# stuff that's precious to keep
.PRECIOUS:	$(STATICLIB) $(SHAREDLIB)
.KEEP_STATE:


# includes, flags and libraries
CC              = gcc
CINCS           = -I$(INCDIR)  -I.


LIBAST		= -last_grf3d -last_pgplot -last
LIBAST		= -last
LIBAST		= /usr/local/lib/libast.a /usr/local/lib/libast_pal.a /usr/local/lib/libast_err.a /usr/local/lib/libast_grf3d.a /usr/local/lib/libast_pgplot3d.a /usr/local/lib/libast_grf_5.6.a /usr/local/lib/libast_grf_2.0.a /usr/local/lib/libast_grf_3.2.a
LIBAST		= ../lib/libast.a ../lib/libast_pal.a ../lib/libast_err.a ../lib/libast_grf3d.a ../lib/libast_pgplot3d.a ../lib/libast_grf_5.6.a ../lib/libast_grf_2.0.a ../lib/libast_grf_3.2.a
LIBCURL		= -lcurl
LIBCFITSIO	= -lcfitsio
LIBPNG		= -lpng
#LIBCURL	= ../lib/libcurl.a
#LIBCFITSIO	= ../lib/libcfitsio.a

        
ifeq ($(PLATFORM),Darwin)
    ifeq  ($(PLMACH),x86_64)
        #CARCH   = -m64 -mmacosx-version-min=10.5
        CARCH   = 
    else
        #CARCH   = -arch i386 -arch ppc -m32 -mmacosx-version-min=10.4
        #CARCH   = -mmacosx-version-min=10.5
        CARCH   = 
    endif
else
    CLIBS       = -lm -lc $(LIBCURL) -lpthread -L/usr/pgsql-9.3/lib -lpq
    CARCH       = -DHAVE_LIBPQ -I/usr/pgsql-9.3/include
    LIBCURL	= -lcurl -lrt
endif
        
DEPLIBS           = $(LIBCURL) $(LIBAST) $(LIBCFITSIO) -lpthread -lm
CLIBS           = -L/usr/local/lib -L../lib/ -lm -lc $(DEPLIBS)
PSQL_LIBS   	= -L/usr/pgsql-9.3/lib -lpq
CFLAGS          = -g -Wall $(CARCH) -D$(PLATFORM) $(CINCS)
#CFLAGS          = -g -Wall $(CARCH) -D$(PLATFORM) $(CINCS) -DDL_DIRECT


# includes, flags and libraries
CC 	    = gcc
CINCS  	    = -I$(HERE) -I./ -I../include -I/usr/local/include

#F77 	    = g77
F77 	    = gfortran
FFLAGS 	    = -g -Wall


# list of source and include files

C_SRCS 	    = dlcutout.c dlmeta.c dlphead.c
C_OBJS 	    = dlcutout.o dlmeta.o dlphead.o
C_INCS 	    = dbConnect.h

F77_SRCS    =
F77_OBJS    =
F77_INCS    =  

SPP_SRCS    =
SPP_OBJS    =
SPP_INCS    =

SPP_TASKS   =
F77_TASKS   =
C_TASKS	    = dlcutout dlmeta

	      
TARGETS	    = $(F77_TASKS) $(SPP_TASKS) $(C_TASKS)

SRCS	    = $(C_SRCS) $(F77_SRCS) $(SPP_SRCS) 
OBJS	    = $(C_OBJS) $(F77_OBJS) $(SPP_OBJS) 
HOST_LIBS   = $(LIBCURL) $(LIBAST) $(LIBCFITSIO) $(LIBPNG) -lpthread -lm
#LIBS        = -L../lib/ $(HOST_LIBS)
LIBS        = -L/usr/local/lib -L../lib/ lib$(NAME).a $(HOST_LIBS)



all: lib lib$(NAME).a
	(make c_progs)
	cp lib$(NAME).a $(LIBDIR)
	cp $(TARGETS) $(BINDIR)

World:  lib

install: all
	cp lib$(NAME).a $(LIBDIR)
	mv $(TARGETS) $(BINDIR)

objs:   $(OBJS)



# Targets

c_progs:    $(C_TASKS)
spp_progs:  $(SPP_TASKS)
f77_progs:  $(F77_TASKS)

distclean:
	make clean
	/bin/rm -rf *.fits *.xml

clean:
	/bin/rm -rf .make.state .nse_depinfo *.[aeo] *.dSYM
	/bin/rm -rf $(TARGETS)

everything:
	make clean
	make all
	make install

help: HELP

install: all 


####################################
#  LIBDLAPPS dependency libraries.
####################################

lib:
	(cd lib ; make all)

lib$(NAME).a:
	(cd lib ; make all)


###############################################################################
# Unit test programs to be built.
###############################################################################

demo:	$(TARGETS)

zztest:  zztest.c
	$(CC) $(CFLAGS) -o zztest zztest.c $(LIBS)


###########################
#  C Test programs
###########################

dlcutout: dlcutout.o lib$(NAME).a
	$(CC) $(CFLAGS) -o dlcutout dlcutout.o $(LIBS)  $(PSQL_LIBS) 
	/bin/rm -rf dlcutout.dSYM

dlmeta: dlmeta.o lib$(NAME).a
	$(CC) $(CFLAGS) -o dlmeta dlmeta.o $(LIBS)
	/bin/rm -rf dlmeta.dSYM

dlphead: dlphead.o lib$(NAME).a
	$(CC) $(CFLAGS) -o dlphead dlphead.o $(LIBS)
	/bin/rm -rf dlphead.dSYM

ftabcnv: ftabcnv.o lib$(NAME).a
	#$(CC) $(CFLAGS) -pg -o ftabcnv ftabcnv.o $(LIBS)
	$(CC) $(CFLAGS) -o ftabcnv ftabcnv.o $(LIBS)
	/bin/rm -rf ftabcnv.dSYM

fits2db: fits2db.o lib$(NAME).a
	$(CC) $(CFLAGS) -o fits2db fits2db.o $(LIBS)
	/bin/rm -rf fits2db.dSYM

fits2sql: fits2sql.o lib$(NAME).a
	$(CC) $(CFLAGS) -o fits2sql fits2sql.o $(LIBS)
	/bin/rm -rf fits2sql.dSYM



tablist: tablist.o lib$(NAME).a
	$(CC) $(CFLAGS) -o tablist tablist.o $(LIBS)
	/bin/rm -rf tablist.dSYM

liststruct: liststruct.o lib$(NAME).a
	$(CC) $(CFLAGS) -o liststruct liststruct.o $(LIBS)
	/bin/rm -rf liststruct.dSYM

zzpng: zzpng.o
	$(CC) $(CFLAGS) -o zzpng zzpng.c $(LIBS) -lpng



###########################
#  SPP Test programs
###########################



###########################
#  Fortran Test programs.
###########################





###############################################################################
# Leave this stuff alone.
###############################################################################

%.o: %.c $(C_INCS)
	/usr/bin/gcc $(CINCS) $(CFLAGS) -g -c $< -o $@

$(STATICLIB): $(C_SRCS:%.c=Static/%.o)
	/usr/bin/ar rv $@ $?
Static/%.o: %.c $(C_INCS)
	/usr/bin/gcc $(CINCS) $(CFLAGS) -g -c $< -o $@
Static:
	/bin/mkdir $@
	chmod 777 $@

$(SHAREDLIB): $(C_SRCS:%.c=Shared/%.o)
	/usr/bin/ld -shared -o $@ $? -lc -ldl
Shared/%.o: %.c $(C_INCS)
	/usr/bin/gcc $(CINCS) $(CFLAGS) -fpic -shared -c $< -o $@
Shared:
	/bin/mkdir $@
	chmod 777 $@
