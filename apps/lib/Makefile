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
BINDIR          := ../../bin/
LIBDIR          := ../../lib/
INCDIR          := ../../include/


# secondary dependencies
LIBBASE     = ../lib$(NAME)
STATICLIB   = $(HERE)/$(LIBBASE).a
SHAREDLIB   = $(HERE)/$(LIBBASE).so.$(VERSION)


# stuff that's precious to keep
.PRECIOUS:	$(STATICLIB) $(SHAREDLIB)
.KEEP_STATE:


# includes, flags and libraries
CC              = gcc
CINCS           = -I$(INCDIR)  -I.

ifeq ($(PLATFORM), "Darwin")
    ifeq  ($(PLMACH), "x86_64")
        CARCH   = -m64 -mmacosx-version-min=10.5
    else
        CARCH   = -arch i386 -arch ppc -m32 -mmacosx-version-min=10.4
    endif
else
    CLIBS       = -lm -lc -lcurl -lpthread
    CARCH       =
endif

CLIBS           = -lm -lc -lcurl -lpthread -lm
CFLAGS          = -g -Wall $(CARCH) -D$(PLATFORM) $(CINCS)
CFLAGS          = -g -Wall $(CARCH) -D$(PLATFORM) $(CINCS) -DDL_DIRECT




# includes, flags and libraries
CC 	    = gcc
CINCS  	    = -I$(HERE) -I../ -I../../ -I../../include 


# list of source and include files

SRCS 	    = dlAclist.c dlRanges.c dlLog.c \
              dlTask.c dlParams.c dlUtil.c dlWCS.c dlHdr.c \
	      dlZScale.c dlPNG.c util.c
OBJS 	    = dlAclist.o dlRanges.o dlLog.o \
              dlTask.o dlParams.o dlUtil.o dlWCS.o dlHdr.o \
	      dlZScale.o dlPNG.o util.o
INCS 	    = ../dlApps.h ../dlAppsP.h


all:    lib

install: all

objs:   $(OBJS)


# Targets

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

lib: objs
	ar rv ../lib$(NAME).a $(OBJS)
	$(RM) $(OBJS)



###############################################################################
# Unit test programs to be built.
###############################################################################

demo:	$(TARGETS)

zz:  zz.c
	$(CC) $(CFLAGS) -o zz zz.c $(LIBS) $(CLIBS) -L../../lib -lcfitsio

zztest:  zztest.c
	$(CC) $(CFLAGS) -o zztest zztest.c $(LIBS)



###############################################################################
# Leave this stuff alone.
###############################################################################

$(STATICLIB): $(C_SRCS:%.c=Static/%.o)
	/usr/bin/ar rv $@ $?
Static/%.o: %.c $(C_INCS)
	/usr/bin/gcc $(CINCS) $(CFLAGS) -c $< -o $@
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
