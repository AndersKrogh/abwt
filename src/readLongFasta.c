/*
This file is part of the abwt package.
Copyright 2016-2021 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

/*
  Functions to read a fasta file into one long sequence
  (including reverse complement if specified). Using a linked list of
  Sequence structs (defined in sequence.h)

  The first struct points to the long sequence (->s), has the
  alphabet in the ->id and the total length of the whole sequence
  (->len). The following contains the individual sequences.

  Letters in in the sequence are translated according to an array
  (translation) with an entry for each of the 128 ASCII codes.
  Symbols that are not letters (spaces, newlines etc) should be
  negative.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "aklib.h"
#include "readLongFasta.h"
#include "common.h"

/* Return the size of a file in bytes */
static long file_size(FILE *fp) {
  long size;
  fseek(fp, 0L, SEEK_END);
  size = ftell(fp);
  rewind(fp);
  return size;
}



static inline void test_seq_alloc(char *seq, char *seqend) {
  if (seq>=seqend) ERROR("Not enough memory allocated for sequence",1);
}




// Skip until stopchar  (copied from sequence.c)
static int skip_until_char(FILE *fp, int stopchar) {
  int c;
  while ( ( c=fgetc(fp)) ) if ( c==EOF || c==stopchar ) break;
  return c;
}


/* Assume that \n> has been read and read to end of line
   ignore anything afer space if save_descr==0
   returns 1 on success

   This could be replaced by a function in sequence.c
*/
int readFastaID(FILE *fp, Sequence *seq, int save_descr) {
  const int id_read_size=256;
  iString *is;
  int c, l, lastc;

  // Read ID (until space)
  is = alloc_iString(id_read_size);
  c = read_line_iString(fp,is,' ',NULL); // Returns 0 on EOF
  lastc = is->lastread;
  l=0;
  if (c && is->len) {
    seq->id = convert_iString(is,&l,1);
    toggleBit(seq->flag,seq_flag_id);
  }
  else free_iString(is);

  if (c && lastc != '\n') { // Read rest of line, if not at end
    if (save_descr) {
      // Read rest of line into description
      is = alloc_iString(id_read_size);
      c = read_line_iString(fp,is,0,NULL);
      l=0;
      if (is->len) seq->descr = convert_iString(is,&l,1);
      else free_iString(is);
      toggleBit(seq->flag,seq_flag_descr);
    }
    else lastc = skip_until_char(fp,'\n'); // Skip to EOL
  }
  if (lastc==EOF) c=0;
  return c;
}



/*
  Read all sequences into one long array and ids into a linked list of IDs.
  string translate holds the char used for each ascii symbol

  Sequences are each terminated by padding term symbols

  Note: seq must point to enough allocated memory!!

  Returns Sequence that starts a linked list. In first element:
     ->s            points to beginning of sequence allocation  
     ->len          (long) holds the total length of concattenated sequence
     ->sort_order   (int) holds the total number of sequences read

*/
static Sequence *read_fasta_forward(FILE *fp, char *seq, long length, char *translate, char term, int padding) {
  const int MaxIDlen = 10000;
  int c, i;
  Sequence *ret=NULL, *ss;
  char *seqend = seq+length;
  char *idline = (char *)malloc((MaxIDlen+2)*sizeof(char));


  /* Absorb until "\n>" (in beginning of file) */
  do c=getc(fp); while (c!='>' && c!=EOF);
  if (c==EOF) return NULL;

  /* Alloc sequence structure to hold anchor and hold total length of
     allocation. */
  ret=ss=alloc_Sequence();
  ret->s = seq;
  ret->sort_order = 0;

  /* Add padding in beginning */
  for (i=0; i<padding; ++i) *seq++ = term;

  while (c!=EOF) {
    /* Now at ">" in beginning of line. Read id line */
    ss->next = alloc_Sequence();
    if (!readFastaID(fp,ss->next,0)) ERROR("Failed to read ID",3);
    ss=ss->next;

    ss->s = seq;
    ss->pos = (long)(seq - ret->s);
    ret->sort_order += 1;

    /* Read sequence */
    int lastc=0;
    c=getc(fp);
    while ( c!=EOF ) {
      /* New fasta entry is starting */
      if (c=='>' && lastc=='\n') {
	ss->len = seq - ss->s;
	// Add padding
	for (i=0; i<padding; ++i) { *seq++ = term; test_seq_alloc(seq,seqend); }
	break;
      }
      if (translate[c]>=0) {
        *seq++ = translate[c];
	test_seq_alloc(seq,seqend);
	ss->len += 1;
      }
      lastc=c;
      c=getc(fp);
    }
  }

  if (seq != ret->s) for (i=0; i<padding; ++i) { *seq++ = term; test_seq_alloc(seq,seqend); }
  ret->len = seq - ret->s;

  free(idline);

  //fprintf(stderr,"read_fasta_forward returning\n");

  return ret;
}


/*
  The ->id points to the original sequence!
 */
static Sequence *revcompLongSequence(Sequence *ss, char *s, AlphabetStruct *a) {
  Sequence *rr = alloc_Sequence();
  long i;
  char *r, *f;

  rr->id = (char*)ss;
  setBit(rr->flag,seq_flag_rev);
  setBit(rr->flag,seq_flag_comp);

  rr->len = ss->len;
  rr->s = s;
  f = ss->s;

  r=rr->s+rr->len-1;
  for (i=0; i<ss->len; ++i) *r-- = a->compTrans[(int)(*f++)];

  return rr;
}



static void revcomp_all(Sequence *fbase, AlphabetStruct *a, char term, int padding) {
  int i;
  Sequence *rbase=NULL, *ss, *rr;
  char *s=fbase->s+fbase->len;

  /* Reverse complement each */
  ss = fbase;
  while (ss->next) {
    ss=ss->next;
    if (!rbase) rr=rbase=revcompLongSequence(ss,s, a);
    else { rr->next=revcompLongSequence(ss,s,a); rr=rr->next; }
    rr->pos = (long)(s-fbase->s);
    s+=rr->len;
    for (i=0; i<padding; ++i) *s++ = term;
    fbase->sort_order += 1;  // Holds number of seqs!
  }
  ss->next = rbase;
  fbase->len = (long)(s-fbase->s);
}




/*
  If length is not given, the file size is used. Otherwise length MUST be at least as long
  as the concatenated string!!

  Sequence is padded with "padding" term chars
*/

Sequence *readLongFasta(FILE *fp, long length, AlphabetStruct *astruct, int revcomp, char term, int padding) {

  char *seq, *tmp;
  long totlen;
  Sequence *ss, *sst;

  /* alloc space for long sequence */
  if (length) totlen = length;
  else totlen = file_size(fp);
  if (revcomp) totlen *= 2;
  /* Add 1% to make room for padding (ad hoc!) */
  totlen += 0.01*totlen;
  seq = (char *)malloc(totlen*sizeof(char));

  if (!seq) {
    fprintf(stderr,"readLongFasta: Failed to alloc seq of length %ld\n",totlen);
    return NULL;
  }

  ss = read_fasta_forward(fp, seq, totlen, astruct->trans, term, padding);

  if (!ss) ERROR("readLongFasta: No sequences read",5);

  if (revcomp) {
    if (totlen<2*ss->len) ERROR("readLongFasta: Not enough space allocated for reverse complement",2);
    revcomp_all(ss, astruct, term, padding);
  }

  /* Down-size to actual size */
  tmp = ss->s;
  ss->s = (char *)realloc((void*)tmp,ss->len);

  /* If down-sizing by realloc actually copies to new array: */
  if (tmp != ss->s) {
    sst=ss;
    while (sst->next) { sst=sst->next; sst->s=ss->s+sst->pos; }
  }

  // fprintf(stderr,"readLongFasta returning\n");

  return ss;
}









