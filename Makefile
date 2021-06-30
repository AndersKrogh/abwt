# This file is part of the abwt package.
# Copyright 2016-2021 by Anders Krogh.
# The abwt package is licensed under the GPLv3, see the file LICENSE.


# for profilling use "make PROF=-pg"

OFLAGS = -O3
#OFLAGS = -g

# The FMITYPE's are 1:compact (standard) 2: simple 3: funny

CFLAGS  = $(OFLAGS) $(PROF) -DFMITYPE=1 -I ../aklib -Wall
LDFLAGS  = -L../aklib/ $(PROF)
LDLIBS  = -laklib -lm -lpthread

BINDIR = ${HOME}/usr/bin/

PERL = /usr/bin/perl

SRC = src
TOOLDIR = ../aklib/tools
TOOLS = $(CURDIR)/$(TOOLDIR)
BACKUP = backup
VPATH = $(SRC)


EXEFILES = makeabwt readfmi searchbwt exactRepeats predictDNA editSequence

ALL: $(EXEFILES)

INSTALL: ALL
	cp -p $(EXEFILES) ${BINDIR}


makeabwt: makeabwt.o bwt.o fmi.o readLongFasta.o suffixArray.o multikeyqsort.o interval.o

exactRepeats: exactRepeats.o bwt.o interval.o suffixArray.o fmi.o

readfmi: readfmi.o bwt.o suffixArray.o fmi.o interval.o

searchbwt: searchbwt.o bwt.o suffixArray.o fmi.o interval.o

predictDNA: predictDNA.o bwt.o suffixArray.o fmi.o interval.o

editSequence: editSequence.o

# o files

makeabwt.o: makeabwt.c makeabwt_vars.inc common.h multikeyqsort.h fmi.h version.h
	$(CC) $(CFLAGS) -c $<

searchbwt.o: searchbwt.c searchbwt_vars.inc bwt.h fmi.h common.h version.h
	$(CC) $(CFLAGS) -c $<

readfmi.o: readfmi.c readfmi_vars.inc bwt.h fmi.h common.h version.h
	$(CC) $(CFLAGS) -c $<

exactRepeats.o: exactRepeats.c exactRepeats_vars.inc bwt.h fmi.h common.h version.h
	$(CC) $(CFLAGS) -c $<

predictDNA.o: predictDNA.c predictDNA_vars.inc bwt.h fmi.h common.h version.h

editSequence.o: editSequence.c editSequence_vars.inc

readLongFasta.o: readLongFasta.c readLongFasta.h common.h

fmi.o: fmi.c fmi.h common.h

suffixArray.o: suffixArray.c suffixArray.h common.h

bwt.o: bwt.c bwt.h fmi.h common.h interval.h

suffixArray.o: suffixArray.c suffixArray.h bwt.h fmi.h common.h

multikeyqsort.o: multikeyqsort.c multikeyqsort.h 

interval.o: interval.c interval.h


# h files

makeabwt_vars.inc: makeabwt.variables
	$(TOOLS)/OptionsAndArguments $< > $(SRC)/$@

searchbwt_vars.inc: searchbwt.variables
	$(TOOLS)/OptionsAndArguments $< > $(SRC)/$@

readfmi_vars.inc: readfmi.variables
	$(TOOLS)/OptionsAndArguments $< > $(SRC)/$@

exactRepeats_vars.inc: exactRepeats.variables
	$(TOOLS)/OptionsAndArguments $< > $(SRC)/$@

predictDNA_vars.inc: predictDNA.variables
	$(TOOLS)/OptionsAndArguments $< > $(SRC)/$@

editSequence_vars.inc: editSequence.variables
	$(TOOLS)/OptionsAndArguments $< > $(SRC)/$@

cleansrc:
	- cd $(SRC); rm -f *.old *~ *_vars.inc

clean: cleansrc
	- rm -f *.o *.old

cleanall: clean
	- rm -f *~ $(EXEFILES)
