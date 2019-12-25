/*
This file is part of the abwt package.
Copyright 2016-2020 by Anders Krogh.
The abwt package is licensed under the GPLv3, see the file LICENSE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "version.h"
#include "common.h"
#include "bwt.h"
#include "fmi.h"
#include "OptionsAndArguments.h"

/*
Program to reconstruct suffix array from fmi (made by mkbwt and mkfmi)
Used only for debugging
*/

#if FMITYPE==1
#define PRGNAME "readfmi (using compact fmi)"
#elif FMITYPE==2
#define PRGNAME "readfmi (using simple fmi)"
#else
#define PRGNAME "readfmi (using funny fmi)"
#endif



static inline int readnumber(char **string_ptr) {
  char x[10];
  int i=0;
  // printf("\n%s -- ",*string_ptr);
  ++(*string_ptr);
  while (isdigit(**string_ptr)) { x[i++]=**string_ptr; ++(*string_ptr); }
  // printf("%s\n",*string_ptr);
  if (i>0) {
    x[i++]='\0';
    return atoi(x);
  }
  else return -10;
}


static inline void printformat(char **string_ptr, char *format, int def, char type) {
  int n;
  n = readnumber(string_ptr);
  if (n<0) n=def;
  sprintf(format,"%%%d%c",n,type);
}


int main (int argc, char **argv) {
  int i, iseq, z, m;
  IndexType pos, k;
  FILE *fp;
  BWT *b;
  uchar c;
  char *seq, *fptr, pformat[20];
  int si_begin, si_length;

  /* Parsing options and arguments */
#include "readfmi_vars.inc"
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  if (!bwtfile) { fprintf(stderr,"You have to specify index name\n"); exit(1); }
  fp = open_file_read(bwtfile,"fmi","readfmi");
  b = readIndexes(fp);
  fclose(fp);

  fprintf(stderr,"BWT of length %ld has been read with %d sequencs, alphabet=%s\n",
	  b->len,b->nseq, b->astruct->a); 


  if (forward) {
    pos = forward;
    c=1;
    i=0;
    while (c && i++ < 500) {
      // Get the character at the position
      c = SAchar(b->f,pos);
      printf("%c",b->astruct->a[c]);
      // Go to next char
      pos = FMIforward(b->f, c, pos);
    }
    printf("\n");
  }

  if (lookup) {
    i = strlen(lookup);
    printf("Suffix interval (start and length) for %s:",lookup);
    translate2numbers(lookup, i, b->astruct);
    SI *si = bwtLookup(b->f, lookup, i);
    if (si) printf("% ld %d\n",si->start, si->len);
    else printf(" %d %d\n",0, 0);
    if (si) {
      if (!sarray) sarray = (int*)malloc(2*sizeof(int));
      sarray[0] = si->start;
      sarray[1] = si->len;
    }
    else {
      if (sarray) free(sarray);
      sarray=NULL;
    }
  }

  if (Lookup) {
    i = strlen(lookup);
    printf("Suffix interval (start and length) for %s:",lookup);
    translate2numbers(lookup, i, b->astruct);
    SI *si = bwtLookup(b->f, lookup, i);
    if (si) printf("% ld %d\n",si->start, si->len);
    else printf(" %d %d\n",0, 0);
  }

  // Print seq info in the order of the SA
  if (PrintId) {
    printf("# Number of seqs: %d\n",b->s->nseq);
    printf("#  length order ID\n");
    for (i=0; i<b->s->nseq; ++i) printf("%9ld %5d %s\n",b->s->seqlengths[i],b->s->seqTermOrder[i],b->s->ids[i]);
  }

  // Print seqs in fasta format
  if (PrintSeq) {
    for (i=0; i<b->s->nseq; ++i) {
      seq = (char*)retrieve_seq(i, b);
      if (PrintSeq>0) {
	printf(">%s len=%ld termOrder=%d\n",b->s->ids[i],b->s->seqlengths[i],b->s->seqTermOrder[i]);
	for (pos=0; pos<b->s->seqlengths[i]; ++pos) {
	  printf("%c",b->astruct->a[(int)seq[pos]]);
	  if ((pos+1)%PrintSeq==0 || pos == b->s->seqlengths[i]-1) printf("\n");
	}
      }
      else {
	printf("%s ",b->s->ids[i]);
	for (pos=0; pos<b->s->seqlengths[i]; ++pos) printf("%c",b->astruct->a[(int)seq[pos]]);
	printf("\n");
      }
    }
  }
  
  if (sarray) {
    si_begin=sarray[0];
    si_length=sarray[1];
    if (marray) {
      si_begin -= marray[0];
      si_length += marray[0] + marray[1];
    }
  }
  else si_length=0;

  // Print part of SA if si_length>0
  for (i=0; i<si_length; ++i) {
    k=FMindexCurrent(b->f,&c,(IndexType)(si_begin+i));
    get_suffix(b->f, b->s, (IndexType)(si_begin+i), &iseq, &pos);
    seq=NULL;


    fptr = format;
    while (*fptr) {
      switch (*fptr) {
      case 'b':
	printformat(&fptr,pformat,10,'d');
	printf(pformat,si_begin+i);
	break;
      case 'i':
	printformat(&fptr,pformat,10,'d');
	printf(pformat,(int)pos);
	break;
      case 'l':
	printformat(&fptr,pformat,10,'d');
	printf(pformat,(int)(b->s->seqlengths[iseq]));
	break;
      case 'f':
	printformat(&fptr,pformat,10,'d');
	printf(pformat,(int)k);
	break;
      case 'p':    // Print the prefix
	if (!seq) seq = (char*)retrieve_seq(iseq, b);
	k = readnumber(&fptr);
	if (k<0) k=6;
	for ( ; k>1; --k) {
	  if (pos-k<0) printf("*");
	  else printf("%c",b->astruct->a[(int)seq[pos-k]]);
	}
	break;
      case 's':    // Print the suffix
	z=0;
	m = readnumber(&fptr);
	if (m<0) m=30;
	if (!seq) seq = (char*)retrieve_seq(iseq, b);
	for (k=0; k<m; ++k) {
	  if (!z && !seq[pos+k]) z=1;
	  if (z) printf("*");
	  else printf("%c",b->astruct->a[(int)seq[pos+k]]);
	}
	break;
      case 'B':    // Print char in BWT
	printf("%c",b->astruct->a[c]);
	++fptr;
	break;
      case 'n':    // Sequence id
	printf("%s",b->s->ids[iseq]);
	++fptr;
	break;
      default:
	printf("%c",*fptr++);
	break;
      }

    }
    printf("\n");

    /*
    printf("%10d%6d%10ld%10ld%10ld  ",si_begin+i,iseq,pos,b->s->seqlengths[iseq],k);
    for (k=5; k>0; --k) {
      if (pos-k<0) printf("*");
      else printf("%c",b->alphabet[(int)seq[pos-k]]);
    }
    printf(" -%c- ",b->alphabet[c]);
    int z=0;
    for (k=0; k<30; ++k) {
      if (!z && !seq[pos+k]) z=1;
      if (z) printf("*");
      else printf("%c",b->alphabet[(int)seq[pos+k]]);
    }
    printf("   %s\n",b->s->ids[iseq]);
    */

    free(seq);
  }

  return 0;
}
