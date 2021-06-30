/*
This file is part of the abwt package.
Copyright 2016-2021 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#ifndef BWT_h
#define BWT_h

#include "common.h"
#include "aklib.h"
#include "fmi.h"
#include "suffixArray.h"
#include "interval.h"

typedef struct {
  IndexType len;      // Length of bwt (not counting initial zeros)
  int nseq;
  char reverse;       // Flag where bit one is reverse and bit 2 is complement
  uchar *bwt;

  AlphabetStruct *astruct;
  // Alphabet
  // int alen;
  // char *alphabet;

  FMI *f;
  suffixArray *s;

} BWT;


typedef struct _SI_ {
  IndexType start;  // Start of suffix interval
  int len;          // Interval length
  int qi;           // Position in query
  int ql;           // Length in query (if relevant)
  int count;        // Used to count matches below current
  int score;
  struct _SI_ *next;
  struct _SI_ *samelen;
} SI;



/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
BWT *alloc_BWT(IndexType len, int nseq, uchar reverse, AlphabetStruct *astruct);
void write_BWT_header(BWT *b, FILE *bwtfile);
BWT *read_BWT_header(FILE *bwtfile);
BWT *read_BWT(FILE *bwtfile);
BWT *readIndexes(FILE *fp);
void free_BWT(BWT *b);
void get_suffix(FMI *fmi, suffixArray *s, IndexType i, int *iseq, IndexType *pos);
uchar *retrieve_seq(int snum, BWT *b);
IndexType InitialSI(FMI *f, uchar ct, IndexType *si);
IndexType UpdateSI(FMI *f, uchar ct, IndexType *si, IndexType *newsi);
void recursive_free_SI(SI *si);
IndexType bwtLookupRaw(FMI *f, char *str, int len, IndexType *si);
SI *bwtLookup(FMI *f, char *str, int len);
SI *maxMatches(FMI *f, char *str, int len, int L, int max_matches);
SI *greedyExact(FMI *f, char *str, int len, int L, int jump);
void construct_fmi_from_files(char *filename, int removefiles, int msg);
int depthFirstOptions(BWT *bwt, interval *si, int depth, void *keeptrack, 
       int (*evalFunc)(interval*, int, void*), int track_si, expandParam *expand, int use_term);
void countNwords (int len, BWT *bwt);
/* FUNCTION PROTOTYPES END */

#endif

static inline int depthFirst(BWT *bwt, interval *si, int depth, void *keeptrack, 
	      int (*evalFunc)(interval*, int, void*) ) {
  return depthFirstOptions(bwt, si, depth, keeptrack, evalFunc, 0, NULL, 0);
}

