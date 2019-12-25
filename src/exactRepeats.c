/*
This file is part of the abwt package.
Copyright 2016-2020 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "aklib.h"
#include "common.h"
#include "fmi.h"
#include "bwt.h"
#include "OptionsAndArguments.h"


/**LOCALSTRUCT Repeat interval *si, ->slen, ->frac, ->n
  int slen;
  int n;    // Non-terminal count
  int flag;
  double frac;  // Fraction of most abundant neighbor letter
  char *seq; !!= (char*)malloc((slen+1)*sizeof(char)) !!-FREE
  IndexType si[2]; !!=NONE
  char rev; // is it reverse?
  struct Repeat *rptr;
ALLOC END SELF->si[0]=si->i[0]; SELF->si[1]=si->i[1];
DUP
**/
//BEGIN AUTO GENERATED

typedef struct Repeat {
  int slen;
  int n; //  Non-terminal count
  int flag;
  double frac; //  Fraction of most abundant neighbor letter
  char *seq;
  IndexType si[2];
  char rev; //  is it reverse?
  struct Repeat *rptr;
} Repeat;
Repeat *alloc_Repeat(interval *si, int slen, double frac, int n) {
  Repeat *r;
  r = (Repeat *)malloc(sizeof(Repeat));
  r->slen = slen;
  r->n = n;
  r->flag = (int)0;
  r->frac = frac;
  r->seq = (char*)malloc((slen+1)*sizeof(char));
  r->rev = (char)0;
  r->rptr = NULL;
  r->si[0]=si->i[0]; r->si[1]=si->i[1];
  return r;
}
void free_Repeat(Repeat *r) {
  if (!r) return;
  if (r->seq) free(r->seq);
  free(r);
}
Repeat *duplicate_Repeat(Repeat *s) {
  Repeat *r = (Repeat *)malloc(sizeof(Repeat));
  memcpy((void*)r,(void*)s,sizeof(Repeat));
  return r;
}
//END AUTO GENERATED


/*
  Return reversed Repeat
 */
static Repeat *reverseRepreat(Repeat *rep, interval *si, char *compTrans) {
  int i;
  Repeat *rev = duplicate_Repeat(rep);
  rev->si[0] = si->r;
  rev->si[1] = si->i[1] - si->i[0] + si->r;
  rev->seq = malloc((rev->slen+1)*sizeof(char));
  for (i=0; i<rep->slen; ++i) rev->seq[rep->slen-1-i] = compTrans[(int)(rep->seq[i])];
  rev->rev=1;
  rev->rptr = rep;
  rep->rptr = rev;
  return rev;
}




static int _compare_Repeat_(const Repeat *a, const Repeat *b) {
  if (a->si[0] > b->si[0]) return -1;
  if (a->si[0] < b->si[0]) return 1;
  // Otherwise the top boundary is the same, so then sort on bottom
  if (a->si[1] < b->si[1]) return -1;
  if (a->si[1] > b->si[1]) return 1;
  // Two SIs are identical, let length determine order
  return -(int)(a->slen - b->slen);
}


/*
  Heap containing repeats
 */
#define HEAPNAME repHeap
#define ITEM Repeat*
#define COMPARE(x,y) _compare_Repeat_( (x),(y) )
#include "heap.template"


/**LOCALSTRUCT repeatParams  ->minlen, ->maxlen, ->min_interval, ->maxfrac, BWT *bwt
   int minlen;
   int maxlen;
   int min_interval;
   double maxfrac;
   char *word; !!= (char*)malloc((->maxlen+1)*sizeof(char)) !!-FREE
   AlphabetStruct *alphabet; !!= bwt->astruct
   FMI *fmi; !!= bwt->f
   repHeap *heap; !!=repHeap_alloc(256) !!-repHeap_free(->heap)
**/
//BEGIN AUTO GENERATED

typedef struct repeatParams {
  int minlen;
  int maxlen;
  int min_interval;
  double maxfrac;
  char *word;
  AlphabetStruct *alphabet;
  FMI *fmi;
  repHeap *heap;
} repeatParams;
repeatParams *alloc_repeatParams(int minlen, int maxlen, int min_interval, double maxfrac, BWT *bwt) {
  repeatParams *r;
  r = (repeatParams *)malloc(sizeof(repeatParams));
  r->minlen = minlen;
  r->maxlen = maxlen;
  r->min_interval = min_interval;
  r->maxfrac = maxfrac;
  r->word = (char*)malloc((r->maxlen+1)*sizeof(char));
  r->alphabet = bwt->astruct;
  r->fmi = bwt->f;
  r->heap = repHeap_alloc(256);
  return r;
}
void free_repeatParams(repeatParams *r) {
  if (!r) return;
  if (r->word) free(r->word);
  if (r->heap) repHeap_free(r->heap);
  free(r);
}
//END AUTO GENERATED



/*
  To be called by the depthFirst function.
 */
static int findRepeats(interval *si, int depth, void *keepTrack) {
  repeatParams *cp=(repeatParams*)keepTrack;
  int c, n, k, tot, max;
  double fraction;
  Repeat *rep;

  // Stop if depth is above max
  if (depth>cp->maxlen) { return 0;}

  // Record the letter
  if (depth>0) cp->word[depth-1]=si->c;

  // Do nothing if word is less than ->minlen and return 1, so search continues
  if ( depth<cp->minlen) return 1;

  // Now we have a word of length ->minlen or larger with a non-zero suffix interval

  // If interval is too small, do nothing
  if ( (int)interval_length(si) < cp->min_interval ) { return 0; }

  // Otherwise do stats on next letter in BWT 
  UpdateInterval(cp->fmi, si);
  // We're not using terminating symbol (0)
  if (si->kids[0]) { free_interval(si->kids[0]); si->kids[0]=NULL; }

  max = tot = 0;
  // Do not count terminations symbols (c=0)
  for (c=1; c<cp->fmi->alen; ++c) {
    if (si->kids[c]) {
      n=(int)interval_length(si->kids[c]);
      tot += n;
      if (n>max) max=n;
      free_interval(si->kids[c]);
      si->kids[c]=NULL;
    }
  }
  free(si->kids);
  si->kids=NULL;

  // If count of non-terminators is less than required, give up
  if (tot < cp->min_interval) { return 0; }

  // If fraction larger than cut-off, the repeat continues, so continue
  fraction = (double)max/(double)tot;
  if ( fraction >= cp->maxfrac ) return 1;

  // Otherwise record the word.
  rep = alloc_Repeat(si, depth, fraction, tot);
  for (k=0; k<depth; ++k) rep->seq[k]=cp->word[depth-k-1];
  repHeap_add(cp->heap,rep);
  repHeap_add(cp->heap,reverseRepreat(rep, si, cp->alphabet->compTrans));

  // Stop here and return 0
  return 0;
}



int main (int argc, char **argv) {
  FILE *fp;
  repeatParams *cp;
  interval *si;
  BWT *bwt;
  int i, n, sortsize;
  Repeat *rep, *prev, **sorted;


  /* Parsing options and arguments */
#include "exactRepeats_vars.inc"
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  // for (i=0; i<MAXLEN; ++i) wl[i]=1;

  // Check args
  if ( minlen<10 || minlen>50 || maxlen>100 || maxlen<minlen ) {
    fprintf(stderr, "minlen=%d maxlen=%d ", minlen, maxlen);
    ERROR("doesn't seem right",1);
  }

  // Open index file
  fp = open_file_read(indexfile,"fmi",NULL);
  fprintf(stderr,"Reading index from file %s\n",indexfile);
  bwt = readIndexes(fp);
  fclose(fp);
  fprintf(stderr,"BWT of length %ld has been read with %d sequences\n",
	  bwt->len,bwt->s->nseq);

  // Alloc struct to hold info
  cp = alloc_repeatParams(minlen, maxlen+1, min_interval, maxfrac, bwt);
 
  /* Depth first traversal of BWT */
  si = alloc_interval(0, bwt->f->bwtlen, 0);
  depthFirst(bwt, si, 0, (void *)cp, findRepeats);

  // Go through intervals and discard overlapping sequences from same intervals
  // Make sorted array from heap
  // Note that 1 was added to maxlen at some point
  sorted = repHeap_sort(cp->heap, &sortsize);
  for ( i=1; i<sortsize; ++i) {
    rep = sorted[i];
    prev = sorted[i-1];
    // The new SI is contained in (or equal to) previous <=> sequence is longer (or identical),
    // so prev is marked for deletion
    if (rep->si[1]<=prev->si[1] || prev->slen>=cp->maxlen) prev->flag=1;
  }
  prev=rep;
  if (rep->slen>=cp->maxlen) prev->flag=1;

  for ( i=0, n=0; i<sortsize; ++i) {
    rep = sorted[i];
    if (printall || !rep->flag) {
      ++n;
      fprintf(stdout,">%s%d",prefix,n);
      if (rep->flag) fprintf(stdout," SUBSTRING ");
      fprintf(stdout," strand=%c SI=%ld-%ld length=%d nseq=%d frac=%lf\n",
	      (rep->rev?'-':'+'),rep->si[0],rep->si[1],rep->slen, rep->n, rep->frac);
      printSeqRaw(stdout,rep->seq,rep->slen,cp->alphabet->a,0,rep->slen);
      fprintf(stdout,"\n");
    }
    free_Repeat(rep);
  }
  free(sorted);

  free_repeatParams(cp);
  free_BWT(bwt);
}
