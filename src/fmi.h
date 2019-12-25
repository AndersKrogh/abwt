/*
This file is part of the abwt package.
Copyright 2016-2020 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/


/*

Three different codings of the FM index identified by FMITYPE:
1: "Campact fmi", which is the default and the one most tested
2: "Simple fmi" - the straight-forward counting FMI
3: "Funny fmi"

You can define it at compile time with -DFMITYPE=3

 */


#ifndef FMITYPE
#define FMITYPE 1
#endif

#ifndef FMINDEX_h
#define FMINDEX_h


#include "common.h"
#include "aklib.h"
#include "interval.h"


/* FM index with a double checkpoint structure */
typedef struct {
  int alen;           // Length of alphabet
  IndexType bwtlen;   // Total length of BWT
  uchar *bwt;         // BWT string
  int N1;             // Total number of entries in index 1 (bwtlen>>ex1 +1);
  int N2;             // Total number of entries in index 2 (bwtlen>>ex2 +1);
  IndexType **index1; // FM index1 (index,letter)
  ushort **index2;    // Counts relative to index1 checkpoints (assuming 16 bit int)
  // Only used in campactfmi
  int *startLcode;    // start numbers for byte encoding of letter and number
  int bwt_page;       // The start offset of the bwt string from page boundary (mmap)
} FMI;

typedef struct {
  int minlen; // Min length of matches in core
  int north;  // Max length of expansion to smaller suffixes
  int south;  // Max length of expansion to larger suffixes
  int greedy; // Flag for greedy mode
  double frac; // Min fraction of majority letter
} expandParam;


FMI *alloc_FMI(uchar *bwt, IndexType bwtlen, int alen);
void free_FMI(FMI *f);
FMI *read_fmi(FILE *fp);
void write_fmi(const FMI *f, FILE *fp);
IndexType FMindex(FMI *f, uchar ct, IndexType k);
IndexType FMindexCurrent(FMI *f, uchar *c, IndexType k);
void FMindexAll(FMI *f, IndexType k, IndexType *fmia);
FMI *makeIndex(uchar *bwt, long bwtlen, int alen);
void UpdateInterval(FMI *f, interval *si);
IndexType FMIforward(FMI *f, uchar c, IndexType k);
int SAchar(const FMI *f, IndexType k);
int ExpandSI(FMI *fmi, interval *parent, expandParam *param);
void trackSIpositions(FMI *fmi, interval *parent);
uchar *get_bwt(const FMI *f, const IndexType start, const IndexType length);
uchar get_bwt_letter(uchar *bwt, IndexType k);


#endif
