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
#include "aklib.h"
#include "common.h"
#include "bwt.h"
#include "fmi.h"
#include "OptionsAndArguments.h"

#if FMITYPE==1
#define PRGNAME "searchabwt (using compact fmi)"
#elif FMITYPE==2
#define PRGNAME "searchabwt (using simple fmi)"
#else
#define PRGNAME "searchabwt (using funny fmi)"
#endif

void print_si(SI *si, BWT *b, FILE *fp) {
  IndexType k, pos;
  int iseq;
  for (k=si->start; k<si->start+si->len; ++k) {
    get_suffix(b->f, b->s, k, &iseq, &pos);
    //fprintf(fp,"%s:%d ",b->s->ids[iseq],si->ql);
    fprintf(fp,"%s pos:%ld start:%ld qi:%d ql:%d\n",b->s->ids[iseq],pos,si->start,si->qi,si->ql);
  }
}


//LEIKLDIMQYLDTRAEEIKAIKNLDMHPRW
void recursive_print_SI(SI *si, BWT *b, FILE *fp) {
  if (!si) return;
  recursive_print_SI(si->samelen, b, fp);
  print_si(si, b, fp);
  recursive_print_SI(si->next, b, fp);
}


int main (int argc, char **argv) {
  int l;
  char *seq;
  FILE *fp;
  BWT *b;
  SI *si;
  AlphabetStruct *astruct;

  /* Parsing options and arguments */
#include "searchbwt_vars.inc"
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  /*
  char *filenm;
  filenm = (char*)malloc(sizeof(char)*strlen(indexfile)+10);
  sprintf(filenm,"%s.fmi",indexfile);
  fp = open_file_read(filenm,"searchbwt");
  free(filenm);
  */

  // reading index
  if (!indexfile) { fprintf(stderr,"You have to specify index name\n"); exit(1); }
  fp = open_file_read(indexfile,"fmi","searchbwt");
  b = readIndexes(fp);
  fclose(fp);

  fprintf(stderr,"BWT of length %ld has been read with %d sequencs, alphabet=%s\n",
	  b->len,b->s->nseq, b->astruct->a); 

  astruct = alloc_AlphabetStruct(b->astruct->a, caseSens, 0, 0, 0);

  seq = malloc(10000);

  while ( fgets(seq,999,stdin) ) {
    l = strlen(seq)-1;
    translate2numbers((char*)seq, l, astruct);

    si = maxMatches(b->f, seq, l, minLen, maxMatch);
    recursive_print_SI(si, b, stdout);
    fprintf(stdout,"\n");
    recursive_free_SI(si);
  }

  return 0;
}
