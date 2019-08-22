/*
This file is part of the abwt package.
Copyright 2016-2019 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/


/*
  The FM index (fmi) is used to identify the number of times the letter c occurred BEFORE
  the current location k. (This is how it is in this implementation, others may include
  current location in count).

  In this implementation, the BWT is recoded. For each segment of 256, for each letter
  (alphabetically) it stores the letter its position k modulus 256. For interval
  k in [i*256,(i+1)*256[.
  Assume the BWT is "224310330335421122" it looks like this (k0=i*256)
     position    letter    value
    rel. to k0           rel. to k0
        0           0        5
        1           0        8
        :           :        :
       n0-1         0        x

       n0           1        4
       n0+1         1       14
       n0+2         1       15
        :           :        :
     n0+n1-1        1        x

     n0+n1          2        0
     n0+n1+1        2        1
     n0+n1+2        2       13
     n0+n1+3        2       16
     n0+n1+4        2       17
        :           :        :
        :           :        :

  For each 256 the count so far and the number of each letter in the interval are
  stored.

  Additionally, we use the actual suffix array (SA)
  position of the corresponding suffix rather than the number from the
  first suffix starting with that letter (which is how the FMI is
  normally described).

  So we supply these functions:
  IndexType FMindex(FMI *fmi, uchar c, IndexType k):
           returns the FMI value at position k for letter c
  IndexType FMindexCurrent(FMI *fmi, IndexType k):
           returns the FMI value at position k for the letter in bwt[k]
           (for efficiency - could have used FMindex with the proper c)
  void FMindexAll(FMI *f, IndexType k, IndexType *fmia) {
           returns the FMI value at position k for all letters in alphabet
*/


FMI *alloc_FMI(uchar *bwt, IndexType bwtlen, int alen) {
  FMI *f = alloc_FMI_common(bwt, bwtlen, alen, sizeof(ushort)+sizeof(uchar));
  return f;
}


void free_FMI(FMI *f) {
  free_FMI_common(f);
}


FMI *read_fmi(FILE *fp) {
  FMI *f = read_fmi_common(sizeof(ushort)+sizeof(uchar),fp);
  return f;
}




void write_fmi(const FMI *f, FILE *fp) {
  write_fmi_common(f, sizeof(ushort)+sizeof(uchar), fp);
  // fwrite(f->startLcode,sizeof(int),f->alen+1,fp);
}


/***********************************************
 *
 * Querying FMI
 *
 *
 ***********************************************/


/*
  At each checkpoint 2, an array of bytes (length alen) holds information about
  the beginning of position array for each letter.
  The begining of the first array is always zero.
  The first byte gives the number of the last array with non-zero length.
  Arrays after the last have (nonsense) begining set to 0

  The array starts after the array of shorts in the checkpoint

  This function takes the short array of the checkpoint, identifies the array
  and returns the beginning of the array and the length (in *l).

  If l==0, returns arbitrary beginning.

  For the very last entry at the end of the bwt, the actual begin is recorded for
  all letters (since it is not 256 even if length is zero).
*/
static inline int get_start_and_length(uchar c, FMI *f, IndexType k, int *l) {
  uchar *a;
  int b, e, L;
  IndexType c2 = (k>>ex2);

  a = (uchar *)(f->index2[c2]+f->alen);
  L = fmi_end_length(k&round2, f->bwtlen);

  if (c>a[0]) { *l=0; return 0; }
  if (c==0) b=0;
  else b=a[c];
  if (c==a[0]) *l = L-b;
  else *l=a[c+1]-b;
  return b;
}


/* Get the checkpointed FMI value for k  */
static inline IndexType fmi_chpt_value(const FMI *f, const IndexType k, const uchar c) {
  return f->index1[k>>ex1][c] + f->index2[k>>ex2][c];
}



#define SEARCH gsearch


/*
  Binary search for integer closest to k and < k in sorted array *a of
  length l all numbers are positive <256

  Returns index+1 (=the number of letters before k)
*/
static inline int ubsearch(uchar *a, int l, uchar k) {
  int i,min,max;

  if (l==0 || a[0]>=k) return 0;
  if (a[l-1] < k) return l;

  min=0; max=l-1;
  while (min<=max) {
    i = (max+min)>>1;  // midpoint
    if (a[i] >= k) max=i-1;
    else {    // a[i]<k
      if (a[i+1]>=k) break;
      else min=i+1;
    }
  }
  return i+1;
}



/*
  Very simple linear search for integer closest to k and < k in sorted array *a of
  length l. All numbers are positive <256

  This function assumes that the numbers in the array are pretty uniformly distributed
  between 0 and 256.
  Initial guess of position in array is l* k/256

  Returns index+1 (=the number of letters before k)
*/
static inline int gsearch(uchar *a, int l, uchar k) {
  int i,min,max;
  const double c=1.0/256;

  if (l==0 || a[0]>=k) return 0;
  if (a[l-1] < k) return l;

  // Now a[0] < k <= a[l-1]

  // Initial guess
  i = (int) c*(k*l);

  while ( a[i]< k ) ++i;   // Makes a[i]>=k
  while ( a[i]>=k ) --i;   // Makes a[i]<k && a[i+1]>k

  return i+1;
}




/* Return the FMI value for target letter ct at position k */
IndexType FMindex(FMI *f, uchar ct, IndexType k) {
  IndexType fmi;
  IndexType k2;
  int b, l;


  // FMI values from check points
  fmi = fmi_chpt_value(f, k, ct);

  // Begin and length of array for letter ct
  b = get_start_and_length(ct, f, k, &l);

  // Find number of letter ct before k
  k2 = k&round2; // Suffix position for checkpoint 2 closest to k 
  if (l>=0) fmi += SEARCH(f->bwt+k2+b, l, (uchar)(k-k2));

  return fmi;
}




/* Return the letter (in *c) and the FMI value for the BWT letter at
   position k

   This is the most problematic task with this index structure:
   the BWT letter at pos. k is not immediately accessible.

   We have to search through all letters!
*/
IndexType FMindexCurrent(FMI *f, uchar *c, IndexType k) {
  IndexType fmi;
  IndexType k2, ch2;
  int b, l, i;
  uchar ct, *bwt, klocal;

  // Suffix position for checkpoint 2 closest to k 
  k2 = k&round2;
  // Index of checkpoint 2
  ch2 = k>>ex2;

  bwt = f->bwt+k2;
  klocal = (uchar)(k-k2);

  for (ct=0; ct<f->alen; ++ct) {
    b = get_start_and_length(ct, f, k, &l);
    if (l==0) continue;
    i = SEARCH(bwt+b, l, klocal);
    // Did we find the letter?
    if ( i<l && *(bwt+b+i)==klocal ) break;
  }

  // Return the letter in *c
  *c = ct;

  // FMI values from check points
  fmi = i + f->index1[k>>ex1][ct] + f->index2[ch2][ct];

  return fmi;
}




/*
  Return the FMI value for all letters at position k
  A result (fmia) array of length alen must be supplied (not checked!)
  Not much optimization can be done, so it is just using the FMindex function.
*/
void FMindexAll(FMI *f, IndexType k, IndexType *fmia) {
  uchar c;
  for (c=0; c<f->alen; ++c) fmia[c] = FMindex(f, c, k);
}





/***********************************************
 *
 * Building FMI
 *
 *
 ***********************************************/



/*
  
*/

FMI *makeIndex(uchar *bwt, long bwtlen, int alen) {
  IndexType i, j, ii, R1, R2;
  int a, l, s, *current;
  uchar **delta;
  uchar *sbwt;
  uchar *carray;
  FMI *fmi;

  fmi = makeIndex_common(bwt, bwtlen, alen);

  current = (int *)malloc(alen*sizeof(int));
  delta = (uchar **)malloc(alen*sizeof(uchar *));
  for (a=0;a<alen;++a) {
    current[a]=0;
    delta[a] = (uchar *)malloc(size2*sizeof(uchar));
  }


  // Note that current char is NOT counted
  i=0;
  R1=R2=0;
  sbwt=bwt;
  for (ii=0; ii<bwtlen; ++ii) {
    /* Check if we are at a checkpoint 2 */
    if ( ii>0 && !(ii&check2) ) {
      R2 = ii>>ex2;
      carray = (uchar*)(fmi->index2[R2-1]+alen);
      /* Values replacing bwt */
      for (a=0, s=0; a<alen; ++a) {
	/* Copy new values for bwt */
	l=current[a];
	memcpy(sbwt+s,delta[a],l*sizeof(uchar));
	/* Start value for letter */
	carray[a]=s;
	if (l>0) carray[0]=a;   // First byte in array gives last letter with non-zero length
	s+=l;
      }
      /* Reset counter */
      for (a=0;a<alen;++a) current[a]=0;
      i=0;
      sbwt+=size2;
    }
    a = sbwt[i];
    delta[a][current[a]] = i;
    current[a] += 1;
    ++i;
  }

  // The last entry of index2 (doesn't matter if last was already a power of ex2)
  R2 = ii>>ex2;
  carray = (uchar*)(fmi->index2[R2]+alen);
  for (a=0, s=0; a<alen; ++a) {
    /* Copy new values for bwt */
    l=current[a];
    memcpy(sbwt+s,delta[a],l*sizeof(uchar));
    /* Start value for letter */
    carray[a]=s;
    if (l>0) carray[0]=a;   // First byte in array gives last letter with non-zero length
    s+=l;
  }

  for (a=1;a<alen;++a) free(delta[a]);
  free(delta);
  free(current);

  return fmi;
}

