/*
This file is part of the abwt package.
Copyright 2016-2020 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

/*
  Structure of suffix intervals (SIs) sorted alphabetically and on interval start
*/

#include <stdlib.h>
#include <string.h>

#include "interval.h"

interval *alloc_interval(IndexType start, IndexType end, int c) {
  interval *r = (interval *)malloc(sizeof(interval));
  r->n=0;
  r->i[0] = start;
  r->i[1] = end;
  r->r=0;
  r->c = c;
  r->flag = 0;
  r->kids = NULL;
  r->next = NULL;
  r->ptr = NULL;
  return r;
}


interval *clone_interval(interval *si) {
  interval *r = (interval *)malloc(sizeof(interval));
  r->n    = si->n;
  r->i[0] = si->i[0];
  r->i[1] = si->i[1];
  r->r    = si->r;
  r->c    = si->c;
  r->flag = si->flag;
  r->kids = NULL;
  r->next = NULL;
  r->ptr = NULL;
  return r;
}


// Takes arrays of FMI values and creates all asize kids for an interval
void interval_make_kids(interval *si, int asize, IndexType *si_low, IndexType *si_high, int protect) {
  int i;
  interval_alloc_kids_array(si,asize);
  for (i=0; i<asize; ++i) {
    if (si_high[i]-si_low[i]<=0) si->kids[i] = NULL;
    else {
      si->kids[i] = alloc_interval(si_low[i], si_high[i], i);
      if (protect) si->kids[i]->next = si;  // Pointer is set to flag that this SI cannot be deleted
    }
  }
}


// Free interval and return the one it points to next
interval *free_interval(interval *si) {
  interval *next=si->next;
  if (si->kids) free(si->kids);
  free(si);
  return next;
}


// Free all kids with n==0 and the array itself
void interval_free_kids(interval *si, int asize) {
  int i;
  // Kids have no ->next set
  for (i=0; i<asize; ++i) if (si->kids[i] && si->kids[i]->n<=0) free_interval(si->kids[i]);
  free(si->kids);
  si->kids=NULL;
}

/*
  If you have updated all kids (by e.g. calling "UpdateInterval", you can
  calculate the start of the reverse intervals of the kids
  ->r must be correct in parent.

  revcomp_order contains the letter ordering in revcomp alphabet.
  if alphabet is *ACGT the array contains 0,4,3,2,1 (*TGCA)
  In the alphabet structure, it is contained in compTrans

  Note that the intervals are done from behind, because we sometimes ignore
  the terminating (first) interval
 */
void interval_update_reverse(interval *si, char *compTrans, int asize) {
  int i,j;
  IndexType k = si->r+si->i[1] - si->i[0]; // end of interval

  // If 0th interval exists
  if (si->kids[0]) si->kids[0]->r = si->r;

  if (compTrans) {
    for (i=asize-1; i>0; --i) {
      j = compTrans[i];
      if (si->kids[j]) {
	k -= si->kids[j]->i[1] - si->kids[j]->i[0];
	si->kids[j]->r=k;
      }
    }
  }
  else {
    for (i=asize-1; i>0; --i) {
      if (si->kids[i]) {
	k -= si->kids[i]->i[1] - si->kids[i]->i[0];
	si->kids[i]->r=k;
      }
    }
  }
}




alphasi *alloc_alphasi(int asize) {
  alphasi *r = (alphasi *)malloc(sizeof(alphasi));
  r->asize = asize;
  r->si = (interval **)calloc(asize,sizeof(interval*));
  r->c  = (interval **)calloc(asize,sizeof(interval*));
  return r;
}

/* Do not free the intervals in the structure (these are pointed to by paths) */
void free_alphasi(alphasi *asi) {
  free(asi->si);
  free(asi->c);
  free(asi);
}

/* Also free intervals (assuming no kids) */
void free_alphasi_all(alphasi *asi) {
  int i;
  interval *si, *tmp;
  for (i=0; i<asi->asize; ++i) {
    si = asi->si[i];
    while (si) { tmp=si; si=si->next; free(tmp); }
  }
  free_alphasi(asi);
}


void add_si_to_alphasi(interval *new, alphasi *asi) {
  int c = new->c;
  if (asi->si[c]==NULL) asi->si[c] = new;
  else asi->c[c]->next = new;
  asi->c[c] = new;
}


interval *iterate_alphasi(alphasi *asi, interval *si) {
  int i;

  // If first call
  if (si==NULL) i = 0;
  else {
    i = si->c+1;
    si=si->next;
  }

  while ( si==NULL && i<asi->asize ) si = asi->si[i++];

  return si;
}



/*
  Go through alphasi struct and 
  - Move kids to new alphasi
  - Delete unused
  - Free old alphasi if they are unused
 */
void alphasi_promote_kids(alphasi *asi) {
  int i,j;
  interval *si, *next, *kid, **old = asi->si;

  asi->si = (interval **)calloc(asi->asize,sizeof(interval*));

  for (i=0; i<asi->asize; ++i) {
    // Go through alphabet
    si = old[i];
    while (si) {
      if (si->kids) {
	for (j=0; j<asi->asize; ++j) {
	  // For each si there are asize kids to process
	  kid = si->kids[j];
	  if (kid) {
	    kid->next=NULL;
	    if (kid->n) add_si_to_alphasi(kid,asi);
	    else free_interval(kid);
	  }
	}
      }
      // Free old interval if unused and move to next
      next = si->next;
      if (!si->n) free_interval(si);
      else {
	free(si->kids);
	si->kids=NULL;
	si->next=NULL;
      }
      si = next;
    }
  }

  free(old);
}




/*
  Go through alphasi struct and identify kids identical to other parents
  If a child SI is identical to a parent, the flag is set to SI_IDENT
  and it's parent is set to c (the initial char of the child).
  for the parent it is identical to, the flag is set to SI_IDENT.
 */
void alphasi_identical(alphasi *asi) {
  int i,j;
  interval *si, *kid, **parents;

  parents = (interval **)malloc(asi->asize*sizeof(interval*));
  memcpy((void*)parents, (void*)(asi->si),asi->asize*sizeof(interval*));

  asi->Nident=0;
  for (i=0; i<asi->asize; ++i) {
    // Go through alphabet
    si = asi->si[i];
    while (si) {
      if (si->kids) {
	for (j=0; j<asi->asize; ++j) {
	  // For each si there are asize kids to process
	  kid = si->kids[j];
	  if (kid) {
	    while (parents[j] && parents[j]->i[0] < kid->i[0]) parents[j]=parents[j]->next;
	    if (parents[j] && parents[j]->i[0] == kid->i[0] && parents[j]->i[1] == kid->i[1] ) {
	      si->flag = j;
	      kid->flag = SI_IDENT;
	      parents[j]->flag = SI_IDENT;
	      asi->Nident += 1;
	    }
	  }
	}
      }
      si=si->next;
    }
  }

  free(parents);
}





#ifdef TEST

#include <stdio.h>

/* */ int main() {
  const int asize=5;
  alphasi *asi = alloc_alphasi(asize);
  int i, j, n;
  interval *si;

  for (i=0; i<asize; ++i) {
    for (n=0; n<100; ++n) {
      si = alloc_interval(n+i*200,10,i);
      add_si_to_alphasi(si,asi);
    }
  }

  si=NULL;
  while ( (si = iterate_alphasi(asi, si)) ) {
    si->kids = (interval **)malloc(asize*sizeof(interval *));
    for (j=0;j<asize;++j) {
      si->kids[j] = alloc_interval(si->i[0]+j,20,j);
      si->kids[j]->n = rand()%2;
    }
  }
  alphasi_promote_kids(asi);
  free_alphasi(asi);
}

#endif
