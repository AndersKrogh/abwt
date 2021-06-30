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

#include "version.h"
#include "aklib.h"
#include "OptionsAndArguments.h"


/*
  This follows the bed convention that numbering starts from 0 and does not
  include the end base is not included. So start=0 and end=1 means the first
  base of the sequence only.
*/

/**LOCALSTRUCT seqFeature char *id,->start,->end
  char *id; !!= strdup(id) !!- FREE
  long start;
  long end;
  long len; !!= SELF->end - ->start
**/
//BEGIN AUTO GENERATED

typedef struct seqFeature {
  char *id;
  long start;
  long end;
  long len;
} seqFeature;
seqFeature *alloc_seqFeature(char *id,long start,long end) {
  seqFeature *r;
  r = (seqFeature *)malloc(sizeof(seqFeature));
  r->id = strdup(id);
  r->start = start;
  r->end = end;
  r->len = r->end - r->start;
  return r;
}
void free_seqFeature(seqFeature *r) {
  if (!r) return;
  if (r->id) free(r->id);
  free(r);
}
//END AUTO GENERATED



/* Read bed-like entry. Ignores everything but first 3 columns
   You have to specify number of first base
   Bed format starts at 0
 */
seqFeature *read_bed_like(FILE *fp, int firstbasenumber, int includeLast) {
  const int linelen = 1<<16;
  char line[linelen], *id;
  long start, end;
  List *strlist;

  while (fgets(line,linelen-1,fp)) {
    if (line[0]!='#' && !isspace(line[0])) {
      strlist = splitString(line,'\t',0);
      id = (char*)popList(strlist);
      start = atol((char*)popList(strlist))-firstbasenumber;
      end = atol((char*)popList(strlist))-firstbasenumber+includeLast;
      freeList(strlist);
      return alloc_seqFeature(line,start,end);
    }
  }
  return NULL;
}

/*
  Continue printing a fasta entry from where it left off
  Call with line_pos=0 on first call

  If line_pos is 0, it behaves like printFasta

  returns new linepos (n+line_pos)%linelen+linelen

  assume linelen=70
  if it ended just before a line break with (n+line_pos)%linelen=69,
  n is incremented after print and it returns 70 and will give a line
  break on next call

  if it ended with a line break with (n+line_pos)%linelen=0,
  it returns 71.
*/
static int continue_printFasta(FILE *file, Sequence *seq, char *alphabet, int linelen, int line_pos, int final) {
  int n=0;
  char *s = seq->s;

  if (seq->len<=0) return line_pos;

  // fprintf(stderr,"CPRINT: %s l=%ld n=%d linepos=%d\n",seq->id,seq->len,n,line_pos);

  //printf("*%d*",line_pos);

  if (linelen<=0) linelen=70;

  // First call
  if (line_pos==0) {
    if (seq->id != NULL) fprintf(file,">%s",seq->id);
    if (seq->descr != NULL) fprintf(file," %s",seq->descr);
    fprintf(file,"\n");
  }

  while (n<seq->len) {
    if ((n+line_pos)>0 && (n+line_pos)%linelen==0) {
      fputc('\n',file);
    }
    fputc(alphabet[(int)s[n++]],file);
  }

  if (final) fputc('\n',file);

  line_pos = linelen + ((n+line_pos)%linelen);

  //printf("!%d!",line_pos);

  return line_pos;
}



int main(int argc, char **argv) {
  AlphabetStruct *astruct, *protein_alph=NULL;
  char *ialph, *balph;
  FILE *fp;
  char eof = 0;
  seqFeature *f=NULL;
  long last_end=-1;
  int k=0, linepos=0, forwardstrand=1;
  char newid[1024];
  Sequence *seq=NULL, *newseq=alloc_Sequence(), *replace_seq=NULL;

  /* Parsing options and arguments */
#include "editSequence_vars.inc"
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  // OPT_print_vars(stderr, opt_struct, "# ", 0);

  // Interpret alphabet and make a version for each
  astruct = bio_AlphabetStruct(Alphabet);
  ialph = strdup(astruct->a);
  balph = strdup(astruct->a);
  if (translate) {
    protein_alph = bio_AlphabetStruct("protein/wildcard/stopcodon");
    makeGeneticCode(astruct,protein_alph);
  }

  if (revcomp) forwardstrand=0;
  if (both) forwardstrand = revcomp = 1;

  if (newChar) {
    if (invert) for (k=0; k<astruct->len; ++k) balph[k] = *newChar;
    else for (k=0; k<astruct->len; ++k) ialph[k] = *newChar;
  }
  else if (edit) {
    if (strcmp(edit, "delete")) replace_seq = make_Sequence(edit,NULL,astruct);
    else replace_seq = alloc_Sequence(); // Dummy sequence with length 0 if delete
  }
  if (upper) {
    if (invert) for (k=0; k<astruct->len; ++k) balph[k] = toupper(balph[k]);
    else for (k=0; k<astruct->len; ++k) ialph[k] = toupper(ialph[k]);
  }
  if (lower) {
    if (invert) for (k=0; k<astruct->len; ++k) balph[k] = tolower(balph[k]);
    else for (k=0; k<astruct->len; ++k) ialph[k] = tolower(ialph[k]);
  }

  // Fasta input
  ReadSequenceFileHeader(stdin, 0);
  
  if (!intervalfile) {
    // Operations on complete input sequences
    if (newChar) ERROR("Option -newChar should not be used when no intervals are given.",1);
    if (split) fprintf(stderr,"WARNING: -split has no effect when no intervals are given.\n");
    if (invert) ERROR("Option -invert should not be used when no intervals are given.",1);
    if (skip) fprintf(stderr,"WARNING: -skip has no effect when no intervals are given.\n");
    if (edit) fprintf(stderr,"WARNING: -edit has no effect when no intervals are given.\n");

    while ( (seq = readFasta(stdin,astruct,1024,0,&eof)) ) {
      newseq->s = seq->s;
      newseq->id = seq->id;
      newseq->len = seq->len;
      if ( forwardstrand ) {
	if (translate) {
	  newseq->s = translateDNA(seq,astruct);
	  newseq->id = newid;
	  sprintf(newseq->id,"%s_transl",seq->id);
	  printFasta(stdout,newseq,protein_alph->a,linelen);
	  free(newseq->s);
	}
	else printFasta(stdout,newseq,ialph,linelen);
      }
      if (revcomp) {
	revcompSequence(seq,astruct);
	newseq->id = newid;
	sprintf(newseq->id,"%s_revcomp",seq->id);
	if (translate) {
	  newseq->s = translateDNA(seq,astruct);
	  sprintf(newseq->id,"%s_transl",newseq->id);
	  printFasta(stdout,newseq,protein_alph->a,linelen);
	  free(newseq->s);
	}
	else printFasta(stdout,newseq,ialph,linelen);
      }
      free_Sequence(seq);
      newseq->id = newseq->s = NULL;
      seq=NULL;
    }
  }

  else {
    if (revcomp) ERROR("Revcomp cannot be used with intervals.\nUse intervals first and pipe to editSequence -r",1);
    fp = open_file_read(intervalfile,NULL,NULL);
    newseq->id=newid;
  }

  while ( intervalfile ) {
    f = read_bed_like(fp,normalIntervals,normalIntervals);
    // if (f) fprintf(stderr,"%s %ld %ld\n",f->id, f->start, f->end);

    while ( !seq || !f || strcmp(seq->id,f->id) ) {
      // If there already is a sequence with different id, finish it
      if ( seq ) {
	if ( (!split || invert ) && ( last_end>0 || !skip) ) {
	  // Print intervening of last seq
	  newseq->s = seq->s + last_end;
	  newseq->len = seq->len - last_end;
	  if ( newseq->len>0 && edit && invert ) {
	    if (replace_seq) {
	      newseq->s = replace_seq->s;
	      newseq->len = replace_seq->len;
	    }
	    else newseq->len=0;
	  }
	  if (split) {
	    if (newseq->len>0) {
	      sprintf(newseq->id,"%s_%ld-%ld",seq->id,last_end+1,seq->len);
	      printFasta(stdout,newseq,balph,linelen);
	    }
	  }
	  else {
	    if (newseq->len==0 && last_end<seq->len) printf("\n"); // Happens if "delete" and end of sequence
	    else linepos = continue_printFasta(stdout,newseq,balph,linelen,linepos,1);
	  }
	}
	free_Sequence(seq);
	seq=NULL;
      }
      //if (!f) break;
      // Read sequence
      if (!eof) seq = readFasta(stdin,astruct,1024,0,&eof);
      if (!seq) break;

      // fprintf(stderr,"Read sequence %s of length %ld\n",seq->id, seq->len);
      // Prepare next round
      last_end=0;
      linepos=0;
      if (!split) newseq->id=seq->id;
    }

    if (!f || !seq) break;

    if ( !split || invert ) {
      newseq->s = seq->s + last_end;
      newseq->len = f->start - last_end;
      if ( edit && invert ) {
	if (replace_seq) {
	  newseq->s = replace_seq->s;
	  newseq->len = replace_seq->len;
	}
	else newseq->len = 0;
      }
      if (split) {
	if (newseq->len>0) {
	  sprintf(newseq->id,"%s_%ld-%ld",seq->id,last_end+1,f->start);
	  printFasta(stdout,newseq,balph,linelen);
	}
      }
      else linepos = continue_printFasta(stdout,newseq,balph,linelen,linepos,0);
    }
    if ( !split || !invert ) {
      newseq->s = seq->s + f->start;
      newseq->len = f->len;
      if ( edit && !invert ) {
	newseq->s = replace_seq->s;
	newseq->len = replace_seq->len;
      }
      if (split) {
	sprintf(newseq->id,"%s_%ld-%ld",seq->id,f->start+1,f->end);
	linepos=0;
      }

      if (split || f->end>=seq->len) k=1;
      else k=0;

      if (k && newseq->len==0) printf("\n"); // Happens if "delete" and end of sequence
      else linepos = continue_printFasta(stdout,newseq,ialph,linelen,linepos,k);      
    }

    last_end = f->end;
    free_seqFeature(f);
    f=NULL;
  }

  if (seq) free_Sequence(seq);
  if (replace_seq) free_Sequence(replace_seq);
  newseq->s = newseq->id = NULL;
  free_Sequence(newseq);
  if (f) free_seqFeature(f);
  free_AlphabetStruct(astruct);
  if (protein_alph) free_AlphabetStruct(protein_alph);
  free(ialph);
  free(balph);

}
