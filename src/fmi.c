/*
This file is part of the abwt package.
Copyright 2016-2021 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#define USE_MMAP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_MMAP
#include <unistd.h>     // To get the page size
#include <sys/mman.h>
#endif

#include "fmi.h"

/*

Three different codings of the FM index identified by FMITYPE:
1: "Campact fmi", which is the default and the one most tested
2: "Simple fmi" - the straight-forward counting FMI
3: "Funny fmi"

You can define it at compile time with -DFMITYPE=3

Note that in the bottom of this file one of three different files
are included depending on the FMITYPE

 */

// Exponent for checkpoints index 2 (ex2<ex1!!)
#if FMITYPE==1
#define ex2 8
#elif FMITYPE==2
#define ex2 8
#else
#define ex2 8     // It must be 8 for "funnyfmi"
#endif

// Exponent for checkpoint index 1 (dist 2^e1 between them)
#define ex1 16   // Exponent for checkpoint index 1 (dist 2^e1 between them)


/* If ex2=3 the binary numbers are:
   size2 =   ...001000   (2^ex2=8)
   check2 =  ...000111 = (1<<ex2)-1
   round2 = ~...000111 = ...111000 (used to round -> set the lowest ex2 bits to zero)
*/
const IndexType size1 = (IndexType)1 << ex1;
//const IndexType check1 = size1-1;
const IndexType check1 = ((IndexType)1 << ex1)-1;
const IndexType size2 = (IndexType)1 << ex2;
// const IndexType check2 = size2-1;
const IndexType check2 = ((IndexType)1 << ex2)-1;
// const IndexType round2 = ~(size2-1);          /* 1 at every bit above bit ex2-1 */
const IndexType round2 = ~(((IndexType)1 << ex2)-1);          /* 1 at every bit above bit ex2-1 */


#if FMITYPE!=1
/*
  Find length to next index2 or end of BWT from position k.
*/
static inline int fmi_end_length(IndexType k, IndexType bwtlen) {
  IndexType k2 = k&round2;
  if ( k2+size2>bwtlen ) return bwtlen-k;
  else return k2+size2-k;
}
#endif


/* Where is the nearest chkpt2? Forward (+1) or backward (-1)
   If the bit ex2-1 is set, the number is half way or more from the previous
   checkpoint. So with mask = 1 << (ex2-1) the check is k&mask
*/
static inline int fmi_direction(IndexType k) {
  const IndexType mask = (IndexType)1 << (ex2-1);
  if ( k&mask ) return 1;
  else return -1;
}



/* Get the checkpointed FMI value for k and direction
   If direction=-1, the checkpoints are k>>ex2 and k>>ex1
   if          = 1, the checkpoints are chpt2=k>>ex2+1 and chpt1 = chpt2>>(ex1-ex2)
*/
static inline IndexType fmi_chpt_value_with_dir(const FMI *f, const IndexType k, const uchar c, const int direction) {
  IndexType chpt1, chpt2;

  chpt2 = k>>ex2;
  if (direction>0) chpt2 += 1;

  /* Note that if direction is +1 and we're just below a checkpoint1, then
     the checkpoint 1 to use is (k>>ex1)+1 rather than k>>ex1
     We can always use chpt1=chpt2>>(ex1-ex2)
  */
  chpt1 = chpt2>>(ex1-ex2);

  return f->index1[chpt1][c] + f->index2[chpt2][c];
}





/* Methods differ in allocation for index2
   + some additionals
*/
static FMI *alloc_FMI_common(uchar *bwt, IndexType bwtlen, int alen, size_t index2_size) {
  int i;
  FMI *f = (FMI*)malloc(sizeof(FMI));
  f->alen = alen;
  f->bwt = bwt;
  f->bwtlen = bwtlen;
  f->N1 = ((bwtlen-1)>>ex1)+2;
  if (f->N1<<ex1 == bwtlen) f->N1 -= 1;
  f->N2 = ((bwtlen-1)>>ex2)+2;
  if (f->N1<<ex2 == bwtlen) f->N2 -= 1;
  f->index1 = (IndexType **)malloc(f->N1*sizeof(IndexType *));
  for (i=0;i<f->N1;++i) f->index1[i]=(IndexType *)malloc(alen*sizeof(IndexType));
  f->index2 = (ushort**)malloc(f->N2*sizeof(ushort*));
  for (i=0;i<f->N2;++i) f->index2[i]=(ushort *)malloc(index2_size*alen);
  f->startLcode = NULL;  // Only used in campactfmi
  return f;
}

static void free_FMI_common(FMI *f) {
  int i;

#ifndef USE_MMAP
  if (f->bwt) free(f->bwt);
#endif

  if (f->index1) {
    for (i=0; i<f->N1; ++i) free(f->index1[i]);
    free(f->index1);
  }
  if (f->index2) {
    for (i=0; i<f->N2; ++i) free(f->index2[i]);
    free(f->index2);
  }
  free(f);
}


/* This function sets the values at index1 and index2.
   Each FMI method may have to additionally recode the BWT
*/
static FMI *makeIndex_common(uchar *bwt, long bwtlen, int alen) {
  IndexType ii, R1, R2, *total;
  int a, *current;
  // uchar *sbwt;
  FMI *fmi = alloc_FMI(bwt,bwtlen,alen);

  current = (int *)calloc(alen,sizeof(int));
  total = (IndexType *)calloc(alen,sizeof(IndexType));
  for (a=0;a<alen;++a) fmi->index2[0][a]=0;

  /*
  for (a=0;a<alen;++a) fmi->index1[0][a]=0;
  for (a=0;a<alen;++a) current[a]=total[a]=0;
  */

  // Note that current char is NOT counted
  // i=0;
  R1=R2=0;
  // sbwt=bwt;
  IndexType ten_percent = bwtlen/10;
  for (ii=0; ii<bwtlen; ++ii) {
    if ( (ii+1)%ten_percent==0 ) fprintf(stderr,"%d%% ... ",10*(int)((ii+1)/ten_percent));
    /* Check if we are at a checkpoint 1 */
    if ( !(ii&check1) ) {
      R1 = ii>>ex1;
      for (a=0;a<alen;++a) fmi->index1[R1][a]=total[a];
    }
    /* Check if we are at a checkpoint 2 */
    if ( ii>0 && !(ii&check2) ) {
      R2 = ii>>ex2;
      /* Checkpoint values */
      for (a=0; a<alen; ++a) fmi->index2[R2][a]=(ushort)(total[a]-fmi->index1[R1][a]);
      /* Reset counter */
      for (a=0;a<alen;++a) current[a]=0;
      //      i=0;
      // sbwt+=size2;
    }
    // a = sbwt[i];
    a = bwt[ii];
    if (a<0 || a>=alen) {
      fprintf(stderr,"makeIndex_common: letter %d not in range at %ld, alen=%d\n",a,ii,alen);
      exit(199);
    }
    total[a] += 1;
    current[a] += 1;
    // ++i;
  }

  // The last entry of index2 (doesn't matter if last was already a power of ex2)
  // R2 = 1+((ii-1)>>ex2);
  // R2 = 1+(ii>>ex2);
  R2 = fmi->N2-1;
  /* Insert values in checkpoint 2 */
  for (a=0;a<alen;++a) fmi->index2[R2][a]=(ushort)(total[a]-fmi->index1[R1][a]);

  fprintf(stderr,"index2 done ... ");

  // Save the letter starts in the last index1
  // Add to all of index1
  fmi->index1[fmi->N1-1][0]=0;
  for (a=1;a<alen;++a) fmi->index1[fmi->N1-1][a]=fmi->index1[fmi->N1-1][a-1]+total[a-1];
  for (R1=0; R1 < fmi->N1-1; ++R1) for (a=1;a<alen;++a) fmi->index1[R1][a] += fmi->index1[fmi->N1-1][a];

  free(current);
  free(total);

  return fmi;
}


/* Write the FMI in file (binary) */
static void write_fmi_common(const FMI *f, int index2_size, FILE *fp) {
  int i;
  fwrite(&(f->alen),sizeof(int),1,fp);
  fwrite(&(f->bwtlen),sizeof(IndexType),1,fp);
  fwrite(&(f->N1),sizeof(int),1,fp);
  fwrite(&(f->N2),sizeof(int),1,fp);
  fwrite(f->bwt,sizeof(uchar),f->bwtlen,fp);
  for (i=0; i<f->N1; ++i) fwrite(f->index1[i],sizeof(IndexType),f->alen,fp);
  for (i=0; i<f->N2; ++i) fwrite(f->index2[i],1,f->alen*index2_size,fp);
}


/* Read the FMI in file (binary)
*/
static FMI *read_fmi_common(int index2_size, FILE *fp) {
  int i;
  FMI *f = (FMI *)malloc(sizeof(FMI));
#ifdef USE_MMAP
  long pos, offset;
  int pagesize=getpagesize();
#endif

  f->bwt=NULL;

  fread(&(f->alen),sizeof(int),1,fp);
  fread(&(f->bwtlen),sizeof(IndexType),1,fp);
  fread(&(f->N1),sizeof(int),1,fp);
  fread(&(f->N2),sizeof(int),1,fp);

#ifdef USE_MMAP
  pos = ftell(fp); // Current postition where reading of bwt starts
  fseek(fp,f->bwtlen,SEEK_CUR); // Move to end of bwt
  offset = (pos/pagesize)*pagesize; // offset must be a multiple of the page size
  f->bwt_page = pos-offset;
  // fprintf(stderr,">>>>>>>>>>>>>>>> %d %ld %ld %ld\n",pagesize,pos,offset,f->bwt_page);
  f->bwt = ((uchar *)mmap(0,f->bwtlen+f->bwt_page, PROT_READ, MAP_SHARED, fileno(fp), offset)) + f->bwt_page;
#else
  f->bwt=(uchar *)malloc(f->bwtlen*sizeof(uchar));
  fread(f->bwt,sizeof(uchar),f->bwtlen,fp);
#endif

  f->index1 = (IndexType **)malloc(f->N1*sizeof(IndexType *));
  for (i=0;i<f->N1;++i) {
    f->index1[i]=(IndexType *)malloc(f->alen*sizeof(IndexType));
    fread(f->index1[i],sizeof(IndexType),f->alen,fp);
  }

  f->index2 = (ushort**)malloc(f->N2*sizeof(uchar*));
  for (i=0;i<f->N2;++i) {
    f->index2[i]=(ushort *)malloc(f->alen*index2_size);
    fread(f->index2[i],1,f->alen*index2_size,fp);
  }

  return f;
}


#if FMITYPE==1
#include "fmi_compact_inc.c"
#elif FMITYPE==2
#include "fmi_simple_inc.c"
#else
#include "fmi_funny_inc.c"
#endif



/*
  Expand intervals for majority child (for assembly)
  Assumes that children are already calculated and that trackSIpositions has been run
  Core suffixes are those that either
    - correpond to a numbered starting suffix (from trackSIpositions)
    - was included more than min_len steps ago
  The function deletes kids that contain no core suffixes
  It checks if there is a dominant letter in the core BWTfile
  If so, it expands the whole suffix as long as
    - the BWT corresponds to the dominant letter
    - or it meets the max (north & south)

  If param->greedy is set, ONLY the dominant child interval is kept
  In this case, search is discontinued if there is no core interval

  THIS DOES NOT WORK WITH FUNNY FMI (there is no fmi_get_letter)
*/
int ExpandSI(FMI *fmi, interval *parent, expandParam *param) {
  //int ExpandSI(FMI *fmi, interval *parent, double frac, int north, int south, int min_len) {
  int *index, *new;
  int c, cmax;
  int i, j, l;
  int count[128];
  interval *si;
  int north, south, nowinner=0;

  // First find maximal child inteval and remove kids that have no part in core interval
  count[0]=0;  // Count 0 is sum  of counts for 1..alen
  cmax=1;
  for (c=1; c<fmi->alen; ++c) {
    count[c] = 0;
    si = parent->kids[c];
    if (si) {
      l = interval_width(si);
      index = (int*)(si->ptr);
      // count[c] is number of "trustable" suffixes
      for (i=0; i<l; ++i) if ( index[i]>=0 || -index[i] >= param->minlen ) ++count[c];
      count[0] += count[c];
      if (count[c]==0) {
	free(index);
	free_interval(si);
	parent->kids[c] = NULL;
      }
      if (count[c]>count[cmax]) cmax=c;
    }
  }

  // Return if there is not a clear majority
  if ( count[cmax]==0 || count[cmax]<param->frac*count[0] ) nowinner=1;

  // If greedy, delete all intervals that are not majority
  if (param->greedy) {
    for (c=0; c<fmi->alen; ++c) {
      if ( c==cmax && nowinner==0 ) continue;
      si = parent->kids[c];
      if (si) {
	free(index);
	free_interval(si);
	parent->kids[c] = NULL;
      }
    }
  }

  // We now have the maximum child interval
  if ( nowinner ) {
    return -1;
  }

  // There is a majority letter, so we expand
  // First look north in parent BWT
  north=0;
  while ( north<param->north && cmax == fmi_get_letter(fmi->bwt[parent->i[0]-north-1])) ++north;
  // Then south
  south=0;
  while (south<param->south && cmax == fmi_get_letter(fmi->bwt[parent->i[1]+south])) ++south;

  // Now expand the appropriate child
  if (north>0 || south>0) {
    si = parent->kids[cmax];
    l = interval_width(si);
    index = (int *)(si->ptr);
    si->ptr = malloc((l+north+south)*sizeof(int));
    new = (int*)(si->ptr);

    for (j=0; j<north; ++j) new[j]=-1;
    for (j=0;j<l;++j) new[north+j]=index[j];
    for (j=0;j<south;++j) new[north+l+j]=-1;
    parent->kids[cmax]->i[0] -= north;
    parent->kids[cmax]->i[1] += south;
    free(index);
  }
  return north+south;
}


/*

  THIS SHOULD BE REPLACED BY FUNCTION UpdateInterval_track

  Function to record letters from the parent to each child
  For each child, an int array is created (si->ptr) that will contain either
    - the negative number of backward extensions for that sequence, or
    - the number of the parent (non-negative)

  Assumes that children are already calculated

  If parent->ptr==NULL, it is assumed positions are 0,1,2,...

  THIS DOES NOT WORK WITH FUNNY FMI (there is no fmi_get_letter)
*/
void trackSIpositions(FMI *fmi, interval *parent) {
  int *parent_track = (int *)(parent->ptr);
  int i, c, wparent, k;
  int kidindex[128], wkids[128], *kids_track[128];
  interval *si;

  wparent = interval_width(parent);

  for (c=0; c<fmi->alen; ++c) {
    si = parent->kids[c];
    kids_track[c] = NULL;
    kidindex[c] = wkids[c] = 0;
    if (si) {
      wkids[c] = interval_width(si);
      si->ptr = calloc(wkids[c],sizeof(int));
      kids_track[c] = (int*)si->ptr;
    }
  }
  // Go through bwt

  // parent_track is the parent index
  if (parent_track) {
    for (i=0; i<wparent; ++i) {
      c = fmi_get_letter(fmi->bwt[parent->i[0]+i]);
      if (wkids[c]) { // If a kid exists
	k = parent_track[i];
	if (k<0) k -= 1;     // Negative numbers count the number of backward extensions
	kids_track[c][kidindex[c]++] = k;
      }
    }
  }
  else {
    for (i=0; i<wparent; ++i) {
      c = fmi_get_letter(fmi->bwt[parent->i[0]+i]);
      // If the parent array does not exist, the position in the parent is recorded as 0...wparent-1
      if (wkids[c]) kids_track[c][kidindex[c]++] = i;    
    }
  }
}


/*
  Get letter at position k in BWT.
  THIS DOES NOT WORK WITH FUNNY FMI (there is no fmi_get_letter)
 */
uchar get_bwt_letter(uchar *bwt, IndexType k) {
  return fmi_get_letter(bwt[k]);
}


/*
  Get all BWT letters starting at start with length length
  Returns an allocated array
  THIS DOES NOT WORK WITH FUNNY FMI (there is no fmi_get_letter)
*/
uchar *get_bwt(const FMI *f, const IndexType start, const IndexType length) {
  IndexType i=0, l;
  uchar *bwt;
  if (start>f->bwtlen) return NULL;
  l = f->bwtlen-start;
  if (l>length) l=length;
  bwt=(uchar *)malloc(l*sizeof(uchar));
#if FMITYPE==1
  for (i=0; i<l;++i) bwt[i] = fmi_get_letter(f->bwt[start+i]);
#elif FMITYPE==2
  memcpy(bwt,f->bwt+start,l);
#endif
  return bwt;
}



/*
  Get all FMI values for letters in an interval
  The si->kids[c] are assumed to be NULL

  THIS DOES NOT WORK WITH FUNNY FMI (there is no fmi_get_letter)
*/
static void FMIcountSI(FMI *f, interval *si) {
  uchar c;
  uchar *bwt;
  IndexType k, fmival, end;

  k = si->i[0];
  bwt = f->bwt+k;
  end = si->i[1];
			  
  while ( k<end ) {
    c = fmi_get_letter(*bwt);
    // fmival = FMindexHere(f,bwt,c,k);
    if ( !si->kids[(int)c] ) {
      fmival = FMindex(f,c,k);
      si->kids[(int)c] = alloc_interval(fmival, fmival+1, (int)c);
      si->kids[(int)c]->next = si;
    }
    else {
      si->kids[c]->i[1] += 1;
    }
    ++k;
    ++bwt;
  }
}



/*
  Calculate child SIs from si
  child->next will point to parent
*/
void UpdateInterval(FMI *f, interval *si) {
  int c;
  const IndexType limit=200;   // Set arbitrarily - should be optimized...
  IndexType fmia[128];

  if (!si->kids) interval_alloc_kids_array(si, f->alen);

  if ( (si->i[1] - si->i[0]) < limit ) {
    FMIcountSI(f,si);
  }
  else {
    FMindexAll(f, si->i[0], fmia);
    for (c=0; c<f->alen; ++c) si->kids[c] = alloc_interval(fmia[c], 0, c);
    FMindexAll(f, si->i[1], fmia);
    for (c=0; c<f->alen; ++c) {
      if (fmia[c]<=si->kids[c]->i[1]) { free_interval(si->kids[c]); si->kids[c]=NULL; }
      else {
	si->kids[c]->i[1]=fmia[c];
	si->kids[c]->next = si;
      }
    }
  }
}





/*
  Forward search in BWT/FM index
  Find SA index for the next letter for a given letter c and its
  SI position k
 */
IndexType FMIforward(FMI *f, uchar c, IndexType k) {
  IndexType k1, k2, b_c, sai;
  int i, mini, maxi, debug=0;
  double a, nLettersInverse;

  //char emessage[10000]; //DEBUGGING
  //char *err = emessage; //DEBUGGING

  // Start position for letter c
  b_c = f->index1[f->N1-1][(int)c];
  // Inverse number of letters
  nLettersInverse = 1.0/(f->index1[f->N1-1][(int)c+1]-b_c);
  
  // FIRST LEVEL CHECK POINTS
  // Maxi is NOT included, so i<maxi
  mini=0;
  maxi = f->N1-1;
  // Fraction of interval per letter c
  a = ((double)(f->N1))*nLettersInverse;

  // Initial guess
  i = (1 + k - b_c)*a;
  if (i>=maxi) i = maxi-1;

  //err+=sprintf(err,"k=%ld ",k);  //DEBUGGING

  while (1) {
    if (k<f->index1[i][(int)c]) {
      // err+=sprintf(err,"|A1 k<%ld ",f->index1[i][(int)c]);  //DEBUGGING
      maxi = i;
      i = maxi-1-(f->index1[i][(int)c]-k)*a;
      if (i<mini) i = mini;
    }
    else if ( i+1<maxi && k>=f->index1[i+1][(int)c] ) {
      // err+=sprintf(err,"|A2 k>%ld ",f->index1[i+1][(int)c]);  //DEBUGGING
      mini = i+1;
      i = mini+(k-f->index1[i+1][(int)c])*a;
      if (i>=maxi) i = maxi-1;
    }
    else {
      // err+=sprintf(err,"|A3 %ld<=k<%ld ",f->index1[i][(int)c],f->index1[i+1][(int)c]);  //DEBUGGING
      break;
    }
  }


  k1 = f->index1[i][(int)c];
  sai = ((IndexType)i)<<ex1;

  //err+=sprintf(err,"|%ld<=k<%ld ",k1,f->index1[i+1][(int)c]);  //DEBUGGING

  // SECOND LEVEL CHECK POINTS
  mini = sai>>ex2;
  maxi = mini + (1<<(ex1-ex2));
  if (maxi>f->N2-1) maxi=f->N2-1;
  a = ((double)(f->N2))*nLettersInverse;

  k2 = k-k1;
  i = mini+k2*a;
  if (i>=maxi) i = maxi-1;
  
  //err+=sprintf(err,"|k2=%ld i=%d mini=%d maxi=%d ",k2,i,mini,maxi);  //DEBUGGING

  while (1) {
    if (k2<f->index2[i][(int)c]) {
      //err+=sprintf(err,"|B1 k1+k2=%ld ", k1+f->index2[i][(int)c]);  //DEBUGGING
      maxi = i;
      i = maxi-1-(f->index2[i][(int)c]-k2)*a;
      if (i<mini) i = mini;
    }
    else if (i+1<maxi && k2>=f->index2[i+1][(int)c]) {
      // err+=sprintf(err,"|B2 k1+k2=%ld %ld ", k1+f->index2[i][(int)c], k1+f->index2[i+1][(int)c]);  //DEBUGGING
      mini = i+1;
      i = mini+(k2-f->index2[i+1][(int)c])*a;
      if (i>=maxi) i = maxi-1;
    }
    else {
      // err+=sprintf(err,"|B3 k1+k2=%ld ", k1+f->index2[i][(int)c]);  //DEBUGGING
      break;
    }
  }
  sai = ((IndexType)i)<<ex2;
  k2 = f->index2[i][(int)c];

  
  // err+=sprintf(err,"|C1 k1+k2=%ld ", k1+f->index2[i][(int)c]);  //DEBUGGING

  // Straight-forward search through BWT for the right position
  i = match_fmi_and_char(f->bwt+sai,k1+k2,k,c);

  if (i<0) ERROR("FMIforward: Something is wrong",1);

  return sai+i;

}

// Return start character for SA position k
int SAchar(const FMI *f, IndexType k) {
  int i;
  for (i=0; i<f->alen-1; ++i) {
    //fprintf(stderr,"%ld %ld %d\n",k,f->index1[f->N1-1][i],i);
    if ( k < f->index1[f->N1-1][i+1] ) break;
  }
  return i;
}
