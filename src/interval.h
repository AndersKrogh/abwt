/*
This file is part of the abwt package.
Copyright 2016-2020 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#ifndef INTERVAL_H
#define INTERVAL_H

#include "common.h"

typedef struct __interval__ {
  IndexType i[2];       // Actual SI interval
  IndexType r;          // Start of the reverse interval
  int n;                // Number of paths using it
  unsigned char c;      // First letter of suffix
  unsigned char flag;   // Used to flag if this SI is identical to another
  struct __interval__ *next;
  struct __interval__ **kids;  // Array of asize pointers to SI extensions
  void *ptr;            // Additional pointer used for e.g. storing paths
} interval;


typedef struct {
  int asize;
  int Nident;    // Number of identical SIs found (function alphasi_identical)
  interval **si; // Array of asize pointers to linked lists of intervals
  interval **c;  // Pointers to current (typically last si for each letter)
} alphasi;

#define SI_IDENT 255


/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
interval *alloc_interval(IndexType start, IndexType end, int c);
interval *clone_interval(interval *si);
void interval_make_kids(interval *si, int asize, IndexType *si_low, IndexType *si_high, int protect);
interval *free_interval(interval *si);
void interval_free_kids(interval *si, int asize);
void interval_update_reverse(interval *si, char *compTrans, int asize);
alphasi *alloc_alphasi(int asize);
void free_alphasi(alphasi *asi);
void free_alphasi_all(alphasi *asi);
void add_si_to_alphasi(interval *new, alphasi *asi);
interval *iterate_alphasi(alphasi *asi, interval *si);
void alphasi_promote_kids(alphasi *asi);
void alphasi_identical(alphasi *asi);
/* FUNCTION PROTOTYPES END */


/* Releasing an interval means to decrease n by 1
   Called when freing a path
   If n<=0, it is freed
*/
static inline void release_interval(interval *si) {
  si->n -= 1;
  if (si->n<=0 && !si->kids && !si->next) free_interval(si);
}

static inline void interval_alloc_kids_array(interval *si, int asize) {
  si->kids=(interval **)calloc(asize,sizeof(interval*));
}

/* Return length of interval */
static inline IndexType interval_length(interval *si) { if (si) return si->i[1]-si->i[0]; return 0; }

// interval_width is synonymous with above
#define interval_width interval_length
//static inline IndexType interval_width(interval *si) { return si->i[1]-si->i[0]; }

// Duplicate interval, but don't copy n, etc
static inline interval *dup_interval(interval *si) {
  return alloc_interval(si->i[0], si->i[1], si->c);
}


#endif
