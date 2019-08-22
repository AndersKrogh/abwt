/*
This file is part of the abwt package.
Copyright 2016-2019 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>

#include "common.h"
#include "aklib.h"
#include "bwt.h"
#include "fmi.h"
#include "suffixArray.h"
#include "interval.h"



/***********************************************************************
 *
 * The BWT struct collects everything and points to both the suffix array
 * and the FMI
 *
 ***********************************************************************/





/***********************************************
 *
 * I/O stuff
 *
 ***********************************************/



// BWT *alloc_BWT(IndexType len, int nseq, int reverse, int alen, char *alphabet) {

BWT *alloc_BWT(IndexType len, int nseq, uchar reverse, AlphabetStruct *astruct) {
  BWT *b=(BWT *)malloc(sizeof(BWT));
  b->len=len;
  b->nseq=nseq;
  b->reverse=reverse;
  b->astruct=astruct;
  // b->alen=alen;
  // b->alphabet=strdup(alphabet);
  b->f=NULL;
  b->s=NULL;
  b->bwt=NULL;
  return b;
}



/*
  Write BWT header
*/
void write_BWT_header(BWT *b, FILE *bwtfile) {
  fwrite(&(b->len),sizeof(IndexType),1,bwtfile);
  fwrite(&(b->nseq),sizeof(int),1,bwtfile);
  fwrite(&(b->reverse),sizeof(char),1,bwtfile);
  write_AlphabetStruct(b->astruct,bwtfile);
  //fwrite(&(b->alen),sizeof(int),1,bwtfile);
  //fwrite(b->alphabet,sizeof(char),b->alen,bwtfile);
}


/*
  Read BWT from file (made by mkbwt) and resturn in BWT struct
*/
BWT *read_BWT_header(FILE *bwtfile) {
  BWT *b=(BWT *)malloc(sizeof(BWT));

  fread(&(b->len),sizeof(IndexType),1,bwtfile);
  fread(&(b->nseq),sizeof(int),1,bwtfile);
  fread(&(b->reverse),sizeof(char),1,bwtfile);
  b->astruct = read_AlphabetStruct(bwtfile);
  // fread(&(b->alen),sizeof(int),1,bwtfile);
  // b->alphabet=(char *)calloc(sizeof(char),b->alen+1);
  // fread(b->alphabet,sizeof(char),b->alen,bwtfile);

  return b;
}



/*
  Read BWT from file (made by mkbwt) and return in BWT struct
*/
BWT *read_BWT(FILE *bwtfile) {
  BWT *b=read_BWT_header(bwtfile);
  b->bwt = (uchar *)malloc(b->len*sizeof(uchar));
  fread(b->bwt,sizeof(uchar),b->len,bwtfile);
  return b;
}






/*
  Read indexes from one file (made by mkfmi) and resturn in BWT struct
*/
BWT *readIndexes(FILE *fp) {
  BWT *b=read_BWT_header(fp);

  b->bwt=NULL;

  b->s = read_suffixArray_header(fp);
  read_suffixArray_body(b->s, fp);
  b->f = read_fmi(fp);

  return b;
}



void free_BWT(BWT *b) {
  if (b->s) free_suffixArray(b->s);
  if (b->f) free_FMI(b->f);
  if (b->bwt) free(b->bwt);
  if (b->astruct) free_AlphabetStruct(b->astruct);
  // if (b->alphabet) free(b->alphabet);
  free(b);
}




/***********************************************
 *
 * Query SA && FMI
 *
 ***********************************************/



/* Find suffix for suffix number i
   Return sequence number in *iseq and position in *pos
*/
void get_suffix(FMI *fmi, suffixArray *s, IndexType i, int *iseq, IndexType *pos) {
  IndexType k=0;
  uchar c=1;

  // CHECKORDER
  // Special case: char==0 && you are at a checkpoint
  if (i<s->nseq) {
    *iseq = reverse_seqTermOrder((int)i,s);
    *pos = s->seqlengths[*iseq];
    if (*pos<0) *pos = -*pos;              // Reverse complement seqs have negative length
    return;
  }

  while ( c && (i & s->check) ) {
    i = FMindexCurrent(fmi,&c,i);
    ++k;
  }


  if (c) {
    suffixArray_decode_number(iseq, pos,
                              (i>>s->chpt_exp)-((s->nseq-1)>>s->chpt_exp)-1, s);
    *pos += k;
  }
  else { *iseq = i; *pos=k-1; }

}


// CHECKORDER
/*
  Reconstruct sequence number snum
  snum is the number of the sorted sequences (the order in the BWT) - which was returned
  from the FMI.
  seqTermOrder will give the original ordering used for retrieval

  If you want to retrieve seq i from the read-order, use reverse look-up to get
  snum, e.g.
  snum = reverse_seqTermOrder(i, sa);
*/
uchar *retrieve_seq(int snum, BWT *b) {
  IndexType l, k;
  uchar *seq, *cur;

  l = b->s->seqlengths[snum];
  if (l<0) l=-l;                              // Reverse complement seqs have negative length
  k = (IndexType)(b->s->seqTermOrder[snum]);
  seq = (uchar *)malloc((l+2)*sizeof(uchar));

  cur = seq+l+1;
  *cur-- = 0; *cur-- = 0;
  while (cur>=seq) {
    k = FMindexCurrent(b->f, cur, k);
    --cur;
  }
  return seq;
}


/* Find whole (initial) suffix interval for letter ct */
IndexType InitialSI(FMI *f, uchar ct, IndexType *si) {
  IndexType r = f->N1-1;
  si[0]=f->index1[r][ct];
  if (ct<f->alen-1) si[1]=f->index1[r][ct+1];
  else si[1]=f->bwtlen;
  return si[1]-si[0]+1;
}



/* Finds the suffix interval for letter ct
   If newsi==NULL and an SI is found, it OVERWRITES values in si
   Note that the actual SI is from si[0] to si[1]-1
*/
IndexType UpdateSI(FMI *f, uchar ct, IndexType *si, IndexType *newsi) {
  IndexType nsi[2];

  nsi[0] = FMindex(f, ct, si[0]);
  nsi[1] = FMindex(f, ct, si[1]);

  if (nsi[0]>=nsi[1]) return 0;

  if (!newsi) newsi=si;
  newsi[0] = nsi[0];
  newsi[1] = nsi[1];

  return nsi[1] - nsi[0];
}




static SI *alloc_SI(IndexType *si, int query_pos, int query_len){
  SI *r = (SI *)malloc(sizeof(SI));
  r->start = si[0];
  r->len=(int)(si[1]-si[0]);
  r->qi = query_pos;
  r->ql = query_len;
  r->count = 0;
  r->next = NULL;
  r->samelen = NULL;
  return r;
}



void recursive_free_SI(SI *si) {
  if (!si) return;
  if (si->next) recursive_free_SI(si->next);
  if (si->samelen) recursive_free_SI(si->samelen);
  free(si);
}



/* Lookup the SI for sequence str */
SI *bwtLookup(FMI *f, char *str, int len) {
  IndexType si[2], l;

  // Go through the sequence from the back
  InitialSI(f, str[--len], si);
  while ( len-- > 0 ) {
    l = UpdateSI(f, str[len], si, NULL);
    if (l <= 0) break;
  }
  if (l>0) return alloc_SI(si, 0 , 0);
  return NULL;
}




/*
  Free matches of each length until there are at least max (we would get <max
  if more were freed).
  Returns min length of retained matches
*/
static inline int free_until_max_SI(SI *si, int max) {
  int n;
  SI *cur;
  if (!si || si->count<=max ) return 0;
  n = si->count;
  cur = si;
  while (cur->next && n - cur->next->count < max) cur=cur->next;
  // Now cur->next==NULL OR totalcount>=max for cur
  if (cur->next) {
    n = cur->next->count;
    recursive_free_SI(cur->next);
    cur->next=NULL;
    while (si) {si->count -= n; si=si->next; }
  }
  return cur->ql;
}



// Brute force stupid sorting (assuming short lists)
static SI *insert_SI_sorted(SI *base, SI *new) {
  SI *tmp;
  new->count=new->len;               // ->len is the si length = number of matches
  if (base==NULL) { return new; }
  if (base->ql<new->ql) {
    new->next=base;
    new->count += base->count;
    return new;
  }
  tmp=base;
  while (tmp->next && tmp->next->ql >= new->ql) {
    tmp->count += new->len;
    tmp=tmp->next;
  }
  tmp->count +=new->len;
  // Now tmp is >= new AND (tmp->next<new OR tmp->next==NULL)
  if (tmp->ql==new->ql) {
    new->samelen=tmp->samelen;
    if (tmp->samelen) new->count += tmp->samelen->count;
    tmp->samelen=new;
  }
  else {
    new->next=tmp->next;
    if (tmp->next) new->count += tmp->next->count;
    tmp->next=new;
  }
  return base;
}




/* Find maximal matches of str of length L in a sorted linked list
   Returns null if there are no matchesBWT *b
   If max_matches==0, not limit imposed
*/
SI *maxMatches(FMI *f, char *str, int len, int L, int max_matches) {
  SI *first=NULL, *cur=NULL;
  IndexType si[2], l;
  int i, j, k;

  // Go through the sequence from the back
  for (j=len-1; j>=L-1; --j) {
    i=j;
    InitialSI(f, str[i], si);
    // Extend backward
    while ( i-- > 0 ) {
      if ( UpdateSI(f, str[i], si, NULL) == 0) break;
    }
    i+=1;
    l = j-i+1;
    if (l>=L) {
      // If the begin of the match (i) equals the the previous, it is within previous match
      if ( !cur || i < cur->qi ) {
	// fprintf(stderr,"MATCH %d %d-%d %ld\n",(int)(si[1]-si[0]), i,j, l);
	cur = alloc_SI(si, i, l);
	first = insert_SI_sorted(first, cur);
	// If max_matches is set, check to see if max is reached and reset L
	if (max_matches>0) {
	  k = free_until_max_SI(first, max_matches);
	  if (k>L) L=k;
	  // The latest si may be freed if too short - then set it to NULL
	  if (l<k) cur=NULL;
	}
      }
    }
    // If the last match reached beginning of sequence, no need to continue
    if (i<=1) break;
  }

  return first;
}


/* Find maximal matches (longer than L) of str of length len in a linked list.
   Returns all matches of maximal length.
   Returns null if there are no matches

   if jump is positive, it jumps by L-jump after having a match of length L
   Note that L is dynamic (max length found)
*/
SI *greedyExact(FMI *f, char *str, int len, int L, int jump) {
  SI *first=NULL, *cur=NULL;
  IndexType si[2], l;
  int i, j, delta=1;

  if (jump>=0) delta=L-jump;

  // Go through the sequence from the back
  for (j=len-1; j>=L-1; j-=delta) {
    i=j;
    InitialSI(f, str[i], si);
    // Extend backward
    while ( i-- > 0 ) {
      if ( UpdateSI(f, str[i], si, NULL) == 0) break;
    }
    i+=1; // Start of match
    l = j-i+1;

    if (l>=L) {
      if (l>L) {
	recursive_free_SI(first);  // Free the shorter ones
	first = NULL;
	L=l;
	if (jump>=0) delta=L-jump;
      }
      cur = first;
      first = alloc_SI(si, i, l);
      first->samelen=cur;
    }
  }

  return first;
}





/* This function reads bwt and SA and, constructs the FMI, writes it and
   deletes the bwt and sa files.
*/
void construct_fmi_from_files(char *filename, int removefiles, int msg) {
  FILE *fp;
  BWT *b;
  char *bwtfilename, *safilename;

  bwtfilename = strconcat2(filename,".bwt");
  safilename = strconcat2(filename,".sa");

  /* Read BWT */
  if (msg) fprintf(stderr,"Reading BWT from file %s\n",bwtfilename);
  fp = open_file_read(bwtfilename,NULL,"makeabwt BWT file");
  b = read_BWT(fp);
  fclose(fp);

  /* Read SA */
  if (msg) fprintf(stderr,"Reading suffix array from file %s\n",safilename);
  fp = open_file_read(safilename,NULL,"makeabwt SA file");
  b->s = read_suffixArray_header(fp);
  /* ?????? If the whole SA is saved, don't read it! */
  if (b->s->chpt_exp > 0) read_suffixArray_body(b->s,fp);
  fclose(fp);

  /* Concatenate stuff in fmi file */
  fp = open_file_write(filename,"fmi","mkfmi FMI file");

  /* Make new file */
  if (msg) fprintf(stderr,"Writing BWT header and SA to file  %s.fmi\n",filename);
  write_BWT_header(b, fp);
  write_suffixArray(b->s,fp);

  if (msg) fprintf(stderr,"Constructing FM index\n");
  b->f = makeIndex(b->bwt, b->len, b->astruct->len);

  if (msg) fprintf(stderr,"\nWriting FM index to file\n");
  write_fmi(b->f,fp);
  fclose(fp);

  if (removefiles) {
    if (msg) fprintf(stderr,"Removing temporary files %s %s\n",bwtfilename,safilename);
    remove(bwtfilename); free(bwtfilename);
    remove(safilename); free(safilename);
  }

  if (msg) fprintf(stderr,"DONE\n");

  // return b;
}



/**********************************************************************
  Depth-first
  Recursively traverse BWT depth first

  evalFunc must return 0 when search should be terminated (e.g. no match)

  keepTrack is a pointer (void*) to an object that is passed to evalFunc

  evalFunc needs to do whatever needs to be done.

  HACK: If track_si>1, SIs are also EXPANDED
 **********************************************************************/


// This should eventually replace the recursive version of the function
/*
  Use ->next and ->n to navigate tree, so do not use for other stuff!!
  If track_si or expand are active, ->ptr is used
  To expand, this struct is used (defined in fmi.h)
struct {
  int minlen; // Min length of matches in core
  int north;  // Max length of expansion to smaller suffixes
  int south;  // Max length of expansion to larger suffixes
  int greedy; // Flag for greedy mode
  double frac; // Min fraction of majority letter
} expandParam;

 */
int depthFirstOptions_nr(BWT *bwt, interval *si, int depth, void *keeptrack, 
		      int (*evalFunc)(interval*, int, void*), int track_si,
		      expandParam *expand, int use_term) {
  int c;
  int alen = bwt->astruct->len;

  c=0;
  while (si) {
    if (si->kids) {  // We have backed up the tree
      // Select next kid
      si->n += 1;
      while ( si->n<alen && si->kids[si->n]==NULL ) si->n += 1;
      // If si is exhausted, back up further
      if (si->n>=alen) {
	// Free all kids
	for (c=0; c<alen; ++c) if (si->kids[c]) {
	    if (si->kids[c]->ptr) free(si->kids[c]->ptr);
	    free_interval(si->kids[c]);
	  }
	free(si->kids);
	si->kids=NULL;
	si = si->next;
	depth -= 1;
      }
      else si = si->kids[si->n];
      continue;
    }

    // back up if not accepted by eval func
    if ( evalFunc(si, depth, keeptrack)==0 ) { si = si->next; depth -= 1; continue; }

    // First visit to this si. calculate si->kids
    UpdateInterval(bwt->f, si);
    if (!use_term && si->kids[0]) { free_interval(si->kids[0]); si->kids[0]=NULL; }

    if (bwt->reverse) interval_update_reverse(si, bwt->astruct->compTrans, alen);

    if (track_si) trackSIpositions(bwt->f, si);

    // Only expand if depth>minlen
    //if (expand && depth>expand->minlen) ExpandSI(bwt->f, si, expand);
    if (expand) ExpandSI(bwt->f, si, expand);

    // Choose first exisitng kid
    si->n = 0;
    while ( si->n<alen && si->kids[si->n]==NULL) { si->n += 1; }
    if (si->n < alen) { si = si->kids[si->n]; depth += 1; }
  }

  return 1;
}


int depthFirstOptions(BWT *bwt, interval *si, int depth, void *keeptrack, 
		      int (*evalFunc)(interval*, int, void*), int track_si,
		      expandParam *expand, int use_term) {
  int c;
  int alen = bwt->astruct->len;

  if ( evalFunc(si, depth, keeptrack)==0 ) return 0;

  // calculate si->kids
  UpdateInterval(bwt->f, si);
  if (!use_term && si->kids[0]) { free_interval(si->kids[0]); si->kids[0]=NULL; }

  if (bwt->reverse) interval_update_reverse(si, bwt->astruct->compTrans, alen);

  if (track_si) trackSIpositions(bwt->f, si);

  // Only expand if depth>minlen
  //if (expand && depth>expand->minlen) ExpandSI(bwt->f, si, expand);
  if (expand) ExpandSI(bwt->f, si, expand);

  // recursively call depthFirst for all letters
  for (c=0; c<alen; ++c) if (si->kids[c]) {
      depthFirstOptions(bwt,si->kids[c],depth+1,keeptrack,evalFunc,track_si,expand,use_term);
    }
  for (c=0; c<alen; ++c) if (si->kids[c]) {
      if (si->kids[c]->ptr) free(si->kids[c]->ptr);
      free_interval(si->kids[c]);
    }
  free(si->kids);
  si->kids=NULL;

  return 1;
}



/*
  Count occurrences of all words of length len using the function depthFirst
*/
typedef struct {
  char *word;
  int len;
  int clen;
  AlphabetStruct *alphabet;
  FILE *ofile;
} countWords;

static countWords *alloc_countWords(int len, AlphabetStruct *alphabet, FILE *ofile) {
  countWords *r = malloc(sizeof(countWords));
  r->len=len;
  r->clen=0;
  r->alphabet=alphabet;
  r->ofile = ofile;
  r->word = (char*)malloc(sizeof(uchar)*(len+1));
  return r;
}

/*
  Count occurrences of all words of length N
  The function is called from depthFirst
  *keepTrack points to a char array, which is the current word

  SI must be non-empty
  When depth==len, print word and its count
*/
static int _countNwords_(interval *si, int depth, void *keepTrack) {
  countWords *cw=(countWords*)keepTrack;
  // int i;
  // char *test = "IMSCTNPTYY";

  // Record the letter
  if (depth>0) cw->word[depth-1]=si->c;

  /*
  // IMSCTNPTYY -> YYTPNTCSMI     174011
  if (depth>2) {
    for (i=depth-1; i>=0; --i) if ( cw->alphabet->a[(int)(cw->word[i])] != test[9-i] ) break;
    if ( i<0 ) {
      fprintf(stderr,"IMSCTNPTYY ");
      i=depth;
      while (--i>=0) fprintf(stderr,"%c",cw->alphabet->a[(int)(cw->word[i])]);
      fprintf(stderr," %d\n",(int)interval_length(si));
    }
  }
  */

  // Do nothing if word is less than ->len and return 1, so search continues
  if ( depth<cw->len) return 1;

  // Now we have a word of length ->len with a non-zero suffix interval
  // Print word (backwards)
  while (--depth>=0) fprintf(cw->ofile,"%c",cw->alphabet->a[(int)(cw->word[depth])]);
  fprintf(cw->ofile," %d\n",(int)interval_length(si));
  // Since maxdepth is reached, return 0
  return 0;
}

void countNwords (int len, BWT *bwt) {
  interval *si = alloc_interval(0, bwt->f->bwtlen, 0);
  countWords *cw = alloc_countWords(len, bwt->astruct, stdout);
  depthFirst(bwt, si, 0, (void *)cw, _countNwords_);
}
