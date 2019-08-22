/*
This file is part of the abwt package.
Copyright 2016-2019 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

static inline uchar fmi_get_letter(uchar x) { return x; }

/*
  The FM index (fmi) is used to identify the number of times the letter c occurred BEFORE
  the current location k. (This is how it is in this implementation, others may include
  current location in count).
*/

FMI *alloc_FMI(uchar *bwt, IndexType bwtlen, int alen) {
  FMI *f = alloc_FMI_common(bwt, bwtlen, alen, sizeof(ushort));
  return f;
}


void free_FMI(FMI *f) {
  free_FMI_common(f);
}


FMI *read_fmi(FILE *fp) {
  FMI *f = read_fmi_common(sizeof(ushort),fp);
  return f;
}



void write_fmi(const FMI *f, FILE *fp) {
  write_fmi_common(f, sizeof(ushort), fp);
  // fwrite(f->startLcode,sizeof(int),f->alen+1,fp);
}



/***********************************************
 *
 * Querying FMI
 *
 *
 ***********************************************/




static int lettercount(const uchar letter, const uchar *str, const int len) {
  int i, c=0;
  for (i=0; i<len; ++i) if (str[i]==letter) ++c;
  return c;
}



/* Get the checkpointed FMI value for k  */
static inline IndexType fmi_chpt_value(const FMI *f, const IndexType k, const uchar c) {
  return f->index1[k>>ex1][c] + f->index2[k>>ex2][c];
}





/* Return the FMI value for target letter ct at position k */
static inline IndexType FMindex_undirectional(FMI *f, uchar ct, IndexType k) {
  IndexType k2 = k&round2;

  return fmi_chpt_value(f,k,ct)+lettercount(ct,f->bwt+k2,k-k2);
}



/* Return the FMI value for target letter ct at position k */
static inline IndexType FMindex_directional(FMI *f, uchar ct, IndexType k) {
  int dir = fmi_direction(k);
  IndexType k2 = k&round2;

  if (dir<0) return fmi_chpt_value_with_dir(f,k,ct,dir)+lettercount(ct,f->bwt+k2,k-k2);
  else return fmi_chpt_value_with_dir(f,k,ct,dir)-lettercount(ct,f->bwt+k,fmi_end_length(k, f->bwtlen));
}



IndexType FMindex(FMI *f, uchar ct, IndexType k) {
  return FMindex_directional(f, ct, k);
}


/* Return the letter (in *c) and the FMI value for the BWT letter at
   position k
*/
IndexType FMindexCurrent(FMI *f, uchar *c, IndexType k) {
  *c=f->bwt[k];
  return FMindex(f, *c, k);
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


FMI *makeIndex(uchar *bwt, long bwtlen, int alen) {
  return makeIndex_common(bwt, bwtlen, alen);
}

