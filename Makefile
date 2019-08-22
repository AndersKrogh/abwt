# This file is part of the abwt package.
# Copyright 2016-2019 by Anders Krogh.
# The abwt package is licensed under the GPLv3, see the file LICENSE.

# for profilling use "make PROF=-pg"

#OFLAGS = -O3
OFLAGS = -g

# The FMITYPE's are 1:compact (standard) 2: simple 3: funny

CFLAGS  = $(OFLAGS) $(PROF) -DFMITYPE=1 -I ../aklib -Wall
LDFLAGS  = -lpthread -lm -L../aklib/ $(PROF)
LDLIBS  = -laklib

PERL = /usr/bin/perl

SRC = src
TOOLDIR = ../tools
TOOLS = $(CURDIR)/$(TOOLDIR)
BACKUP = backup
VPATH = $(SRC)


EXEFILES = makeabwt readfmi searchbwt exactRepeats

ALL: $(EXEFILES)

makeabwt: makeabwt.o bwt.o fmi.o readLongFasta.o suffixArray.o multikeyqsort.o interval.o

exactRepeats: exactRepeats.o bwt.o interval.o suffixArray.o fmi.o

readfmi: readfmi.o bwt.o suffixArray.o fmi.o interval.o

searchbwt: searchbwt.o bwt.o suffixArray.o fmi.o interval.o



# o files

makeabwt.o: makeabwt.c makeabwt_vars.inc common.h multikeyqsort.h fmi.h version.h
	$(CC) $(CFLAGS) -c $<

searchbwt.o: searchbwt.c searchbwt_vars.inc bwt.h fmi.h common.h version.h
	$(CC) $(CFLAGS) -c $<

readfmi.o: readfmi.c readfmi_vars.inc bwt.h fmi.h common.h version.h
	$(CC) $(CFLAGS) -c $<

exactRepeats.o: exactRepeats.c exactRepeats_vars.inc bwt.h fmi.h common.h version.h
	$(CC) $(CFLAGS) -c $<

readLongFasta.o: readLongFasta.c readLongFasta.h common.h

fmi.o: fmi.c fmi.h common.h

suffixArray.o: suffixArray.c suffixArray.h common.h

bwt.o: bwt.c bwt.h fmi.h common.h interval.h

suffixArray.o: suffixArray.c suffixArray.h bwt.h fmi.h common.h

multikeyqsort.o: multikeyqsort.c multikeyqsort.h 

interval.o: interval.c interval.h

# h files

makeabwt_vars.inc: makeabwt.variables
	cd $(SRC); $(TOOLS)/OptionsAndArguments makeabwt.variables > $@

searchbwt_vars.inc: searchbwt.variables
	cd $(SRC); $(TOOLS)/OptionsAndArguments searchbwt.variables > $@

readfmi_vars.inc: readfmi.variables
	cd $(SRC); $(TOOLS)/OptionsAndArguments readfmi.variables > $@

exactRepeats_vars.inc: exactRepeats.variables
	cd $(SRC); $(TOOLS)/OptionsAndArguments exactRepeats.variables > $@


cleansrc:
	- cd $(SRC); rm -f *.old *~ *_vars.inc

clean: cleansrc
	- rm -f *.o *.old

cleanall: clean
	- rm -f *~ $(EXEFILES)
