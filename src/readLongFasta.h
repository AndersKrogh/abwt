/*
This file is part of the abwt package.
Copyright 2016-2020 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#ifndef READLONGFASTA_h
#define READLONGFASTA_h

#include "common.h"

/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
Sequence *readLongFasta(FILE *fp, long length, AlphabetStruct *astruct, int revcomp, char term, int padding);
/* FUNCTION PROTOTYPES END */

#endif
