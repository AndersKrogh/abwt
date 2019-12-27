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



/*
  p[i]*log2(p[i]/0.25) = p[i]*( log2(p[i]) + log2(4) ) = p[i]*( log2(p[i]) + 2 )
*/
double entropy(double *p){
  int i;
  double E=0.;
  for (i=1; i<5; ++i) { if (p[i]>0.) E += p[i]*( 2.+log2(p[i]) ); }
  return E;
}


static inline void interval_reverse(interval *si, char *compTrans) {
  IndexType l = interval_width(si);
  IndexType i = si->i[0];
  si->i[0]=si->r;
  si->i[1]=si->r+l;
  si->r=i;
  if (compTrans) si->c = compTrans[si->c];
}



// Clean up and keep only the relevant SI
static inline interval *child_interval(interval *si, int refbase, int alen) {
  si->next = si->kids[refbase];
  si->kids[refbase]=NULL;
  // Free unused kids
  interval_free_kids(si, alen);
  return free_interval(si);  // free_interval returns next;
}


/*
  Count number of suffixes for letter 1 to 4 in SI
  Total count is returned (and also in count[0])
  If compTrans in non-null, kids are reverse complemented
 */
static inline int count_nucleotides_in_kids(interval *si, int *count, char *compTrans) {
  int c,r;
  count[0]=0;
  for (r=1; r<=4; ++r) {
    if (compTrans) c = compTrans[r];
    else c=r;
    count[c]=0;
    if (si->kids[r]) count[c]=interval_width(si->kids[r]);
    count[0] += count[c];
    //fprintf(stderr," - %1d/%1d %d\n",r,c,count[c]);
  }
  return count[0];
}


/*
  Caclulate probs from count array
  Regularizes with alpha*p[c]
  alpha sits in p[0]!!
  1 is subtracted from reference count in the array unless refbase=0
 */
static inline void count2probs_nucleotides(int *count, int refbase, double *p) {
  int c;
  double norm;

  if ( refbase && count[refbase]<=0 ) ERROR("count2probs: Reference base has 0 prob",1);
  if ( count[0]>0 ) {
    // Subtract 1 for the reference base
    if (refbase) count[refbase] -= 1;
    // Calc probs
    norm = 1.0/(p[0]+count[0]-1);
    for (c=1; c<=4; ++c) p[c] = (p[0]*p[c]+count[c])*norm;
  }
}



/*
  Markov probabilities are calculated for position k_max of sequence (str)
       P(b| x_0, ...., x_{kmax-1})
  Returned in p (p[1]...p[4]) - assuming alphabet *ACGTN
  Interpolates from kmin to kmax

  Idea is to calculate the suffix interval (SI) for sequence
  x_{kmax-k},...,x_{kmax-1}, while keeping track of the 4 reverse complement
  SIs for sequence x_{kmax-k},...,x_{kmax}

  If subtract1 is set, 1 is subtracted from reference base count.
 */
static int MarkovProbs(BWT *bwt, char *str, int kmin, int kmax, double *p, double alpha, int subtract1) {
  int k, c;
  int alen = bwt->astruct->len;
  interval *si = alloc_interval(0, bwt->f->bwtlen, 0);
  interval *rsi = alloc_interval(0, bwt->f->bwtlen, 0);
  int count[6];
  char ref;

  p[0]=alpha;
  for (c=1; c<=4; ++c) p[c]=0.;

  for (k=1; k<=kmax; ++k) {
    // Backward SIs for sequence x_{kmax-k},...,x_{kmax-1}
    UpdateInterval(bwt->f, si);
    // Because we need the SIs for the revcomp, all 4 SIs are needed.
    // Here the SIs for the reverse complement is calculated for the kids
    interval_update_reverse(si, bwt->astruct->compTrans, alen);
    // Choosing the backward interval only for letter str[kmax-k]
    si = child_interval(si, str[kmax-k],alen);
    if ( si == NULL ) return 0;

    if (k>=kmin) {
      // Start interpolation
      // Make the reverse interval
      rsi->i[0]=si->r;
      rsi->i[1]=si->r + interval_width(si);
      // Update the reverse complement interval
      UpdateInterval(bwt->f, rsi);
      // We now have the SI for the 4 complement bases for k
      count_nucleotides_in_kids(rsi, count, bwt->astruct->compTrans);
      // The reference base is str[kmax]
      ref = str[kmax];
      if (!subtract1) ref=0;
      if (count[0]) count2probs_nucleotides(count, ref, p);

      //for (c=1; c<=4; ++c) fprintf(stderr, "%d %lf ", count[c], p[c]);
      //fprintf(stderr,"tot=%d k=%d\n",count[0],k);

      interval_free_kids(rsi,alen);
    }

  }
  free_interval(si);
  free_interval(rsi);
  return count[0];
}



/* Extend an interval both ways with characters first_char and last_char
   Assumes that si->r has ben updated
 */
static interval *extend_both_ways(BWT *bwt, interval *si, int first_char, int last_char) {

  // Backward children
  UpdateInterval(bwt->f, si);
  interval_update_reverse(si, bwt->astruct->compTrans, bwt->astruct->len);

  // Save only relevant child
  si = child_interval(si, first_char, bwt->astruct->len);
  if (!si) return NULL;

  // Calculate the forward intervals
  interval_reverse(si,bwt->astruct->compTrans);
  UpdateInterval(bwt->f, si);
  interval_update_reverse(si, bwt->astruct->compTrans, bwt->astruct->len);
  // Save only relevant child (for revcomp of last_char)
  si = child_interval(si, bwt->astruct->compTrans[last_char], bwt->astruct->len);
  if (!si) return NULL;

  // Reverse to get back to original strand
  interval_reverse(si,bwt->astruct->compTrans);
  return si;
}


static int CentralProbs(BWT *bwt, char *str, int kmin, int kmax, double *p, double alpha, int subtract1) {
  int k, c;
  int alen = bwt->astruct->len;
  int count[6];
  interval *si[5];
  char ref;

  // Initialize probabilities
  p[0]=alpha;
  for (c=1; c<=4; ++c) p[c]=0.;
  // Get first intervals
  si[0]=alloc_interval(0, bwt->f->bwtlen, 0);
  UpdateInterval(bwt->f, si[0]);
  interval_update_reverse(si[0], bwt->astruct->compTrans, bwt->astruct->len);
  // Initialize one interval per central letter
  for (c=1; c<=4; ++c) {
    si[c]=si[0]->kids[c];
    si[0]->kids[c] = NULL;
  }
  interval_free_kids(si[0], alen);
  free_interval(si[0]);
  si[0]=NULL;

  //  for (c=1; c<=4; ++c) fprintf(stderr,"SI for lett %d: %ld - %ld\n",c,si[c]->i[0],si[c]->i[1]);

  for (k=1; k<=kmax; ++k) {
    // Calculate the 4 new intervals
    for (c=1; c<=4; ++c) if (si[c]) si[c] = extend_both_ways(bwt, si[c], str[kmax-k], str[kmax+k]);
    //for (c=1; c<=4; ++c) fprintf(stderr,"SI for lett %d: %ld - %ld\n",c,si[c]->i[0],si[c]->i[1]);

    if (k>=kmin) {
      // We now have the SI for the 4 bases for k
      count[0]=0;
      for (c=1; c<=4; ++c) {
	count[c]=0;
	if (si[c]) count[0] += ( count[c] = interval_width(si[c]) );
      }
      // The reference base is str[kmax]
      ref = str[kmax];
      if (!subtract1) ref=0;
      if (ref && count[(int)ref]<=0) ERROR("CentralProbs: non-existing sequence",1);
      if (count[0]) count2probs_nucleotides(count, ref, p);

      //for (c=1; c<=4; ++c) fprintf(stderr, "%d %lf ", count[c], p[c]);
      //fprintf(stderr,"tot=%d k=%d\n",count[0],k);
    }
  }
  for (c=1; c<=4; ++c) if (si[c]) free_interval(si[c]);
  return count[0];
}


typedef struct {
  BWT *bwt;
  int kmax;
  int kmin;
  int wsize;
  double alpha;
  int Markov;
  int bothdir;
  int printPos;
  int subtract1;
  IndexType correct;
  IndexType total;
} containerStruct;

static void printVCFheader(FILE *out) {
  fprintf(out,"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER INFO\n");
}

static void printVCFline(FILE *out, char *id, IndexType position, char *seq, int seqlen,
			 int iref, int maxi, char *alphabet, double *p) {
  int i;
  fprintf(out,"%s\t%ld\t",id,position);
  printSeqRaw(out, seq, seqlen+1, alphabet, 0, seqlen);
  fprintf(out,"\t%c\t%c\t%lf\t%s\tE=%lf;P",alphabet[iref],alphabet[maxi],
	 p[iref],
	 (maxi==iref?"PASS":"no"),
	 entropy(p));
  for (i=1; i<=4; ++i) fprintf(out,"%c%lf",(i==1?'=':','),p[i]);
  printf("\n");
}


/*
  Estimate probabilities and print for sequence seq.
  The sequence must have 2*kmax letters (or kmax+1 if Markov and NOT bidirectional)
*/
static void seq_probabilities (char *seq, char *rev, char *id, IndexType pos, containerStruct *c) {
  int N, i, count, maxi, iref;
  int printlen=c->wsize;
  double p[6], p2[6];
  char revcomp[256];

  // Make revcomp
  if (!rev && c->bothdir) {
    rev=revcomp;
    for (i=0; i<c->wsize;++i) rev[c->wsize-1-i] = c->bwt->astruct->compTrans[(int)(seq[i])];
  }

  // Check if seq is appropriate
  for (i=0, N=0; i<c->wsize; ++i) if ( seq[i]<1 || seq[i]>4 ) { N=1; break; }
  if (N==0) {
    //if (c->printPos) printf("%s\t%ld\t",id,pos+c->kmax+1);   // first base has number 1
    if (c->Markov) {
      count = MarkovProbs(c->bwt, seq, c->kmin, c->kmax, p, c->alpha, c->subtract1);
      if (c->bothdir) {
	count += MarkovProbs(c->bwt, rev, c->kmin, c->kmax, p2, c->alpha, c->subtract1);
	for (i=1; i<=4; ++i) p[i] = 0.5*(p[i]+p2[(int)(c->bwt->astruct->compTrans[i])]);
      }
      else printlen=c->kmax+1;
    }
    else {
      count = CentralProbs(c->bwt, seq, c->kmin, c->kmax, p, c->alpha, c->subtract1);
    }
    maxi=1;
    for (i=2; i<=4; ++i) if (p[maxi]<p[i]) maxi=i;
    iref = seq[c->kmax];
    if (iref == maxi) c->correct+=1;
    c->total += 1;

    if (c->printPos) printVCFline(stdout, id, pos+c->kmax+1, seq, printlen,
		 iref, maxi, c->bwt->astruct->a, p);
  }
}




int main (int argc, char **argv) {
  FILE *fp;
  Sequence *seq;
  int i, iseq, dz;
  IndexType l, k, pos;
  const int bufsize=64, bufsize2 = 128;
  char dbuffer[bufsize2], rbuffer[bufsize2], *id;
  uchar cur;

  containerStruct *cs = (containerStruct*)malloc(sizeof(containerStruct));
  
  /* Parsing options and arguments */
#include "predictDNA_vars.inc"
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  // Check args
  if (kmin<1 || kmin>16) { fprintf(stderr, "kmin=%d ", kmin); ERROR("doesn't seem right",1); }
  if (kmax<kmin || kmax>16) { fprintf(stderr, "kmax=%d ", kmax); ERROR("doesn't seem right",1); }

  if (bothdir && !Markov) ERROR("Option -bothdir/-b only works with -Markov/-M",1);
  if (GenomeStats) { cs->printPos=0; genome=1; }
  else cs->printPos=1;

  // Open index file
  fp = open_file_read(indexfile,"fmi",NULL);
  fprintf(stderr,"Reading index from file %s\n",indexfile);
  cs->bwt = readIndexes(fp);
  fclose(fp);
  fprintf(stderr,"BWT of length %ld has been read with %d sequences\n",
	  cs->bwt->len,cs->bwt->s->nseq);

  //  Check that things make sense for DNA (makeabwt -r -a DNA/w)
  fprintf(stderr,"Alphabet in index file: %s\n",cs->bwt->astruct->a);
  if ( strncmp(cs->bwt->astruct->a+1,"ACGT",4) )
    fprintf(stderr,"WARNING: Program assumes ACGT for letter 2 to 5\n");

  if ( 2&(cs->bwt->reverse) ) {   // Bit 1 is reverse, bit 2 is complement
    fprintf(stderr,"Index contains the reverse complement (as it should)\n");
  }
  else {
    fprintf(stderr,"Reverse complement flag (bit 2 must be set): %d\n",cs->bwt->reverse);
    ERROR("Index must contain reverse complement (makeabwt -r)",1);
  }


  cs->kmax=kmax;
  cs->kmin = kmin;
  cs->alpha = alpha;
  cs->Markov = Markov;
  cs->bothdir = bothdir;
  cs->correct = 0;
  cs->total = 0;
  cs->subtract1 = 1;

  //total window width
  if ( !cs->Markov || cs->bothdir ) cs->wsize = 2*cs->kmax+1;
  else cs->wsize=cs->kmax+1;

  /*
  if (sequence) {
    if (strlen(sequence)<kmax+1) ERROR("Sequence too short",1);
    c = sequence[kmax];
    sequence[kmax] = '\0';
    printf("Sequence %s* ref=%c\n", sequence, c);
    sequence[kmax] = c;
    seq = make_Sequence(sequence,strdup("testseq"), cs->bwt->astruct);
    if (subtractRef < 0) cs->subtract1 = 0;
    if (Markov) MarkovProbs(cs->bwt, seq->s, kmin, kmax, p, alpha, cs->subtract1);
    else CentralProbs(cs->bwt, s->s, kmin, kmax, p, alpha, cs->subtract1);
    for (c=1; c<=4; ++c) printf("%lf ",p[c]);
    printf("<- probs for A,C,G and T\n");
  }
  */

  if (sequence) {
    char *s;
    char eof = 0;
    FILE *infile = stdin;
    if ( strncmp(sequence,"stdin",5) ) infile=open_file_read(sequence,NULL,"Problem opening sequence file");
    i = ReadSequenceFileHeader(infile,'>');
    if (subtractRef < 0) cs->subtract1 = 0;
    printVCFheader(stdout);
    while (!eof) {
      seq = readFasta(infile,cs->bwt->astruct,1024,0,&eof);
      if (!seq) break;
      l = seq->len;
      fprintf(stderr,"Analyzing sequence \"%s\" of length %ld\n",seq->id, l);
      s = seq->s;
      for (pos=0; pos<=l-cs->wsize; ++pos) {
	seq_probabilities(s, NULL, seq->id, pos, cs);
	++s;
      }
    }
    fclose(infile);
  }

  if (genome) {

    // For sequence included in index, we assume that 1 is subtracted
    // unless it is set to 0 by user.
    if (subtractRef == 0) cs->subtract1 = 0;

    printVCFheader(stdout);
    for (iseq=0; iseq<cs->bwt->s->nseq; ++iseq) {
      l = cs->bwt->s->seqlengths[iseq];
      if (l>=0) i=1;
      else {
	i=0; 
	if (identifier) i = strcmp(cs->bwt->s->ids[iseq],identifier);
	// The right seq is found if i=0
      }
      if (i==0) {  // We travel backwards in reverse complement sequence
	l=-l;
	id = cs->bwt->s->ids[iseq];
	fprintf(stderr,"Analyzing sequence \"%s\" of length %ld\n",cs->bwt->s->ids[iseq], l);
	// Fill buffer
	k = (IndexType)(cs->bwt->s->seqTermOrder[iseq]);
	for (dz=0; dz<bufsize; ++dz) {
	  k = FMindexCurrent(cs->bwt->f, &cur, k);
	  dbuffer[dz] = cs->bwt->astruct->compTrans[cur];
	  rbuffer[bufsize2-1-dz] = cur;
	}
	for (pos=0; pos<=l-cs->wsize; ++pos) {
	  // pos is the start of the window
	  // dz is start in buffer of direct strand
	  dz = (int)(pos%bufsize);
	  // call function with dbuffer+dz
	  // Last base of reverse: x=bufsize2-1-dz;
	  // First base is then x-2*kmax = bufsize2-dz-wsize
	  seq_probabilities(dbuffer+dz, rbuffer+bufsize2-dz-cs->wsize, id, pos, cs);

	  if (pos+bufsize<l) {
	    k = FMindexCurrent(cs->bwt->f, &cur, k);
	    dbuffer[dz+bufsize] = dbuffer[dz] = cs->bwt->astruct->compTrans[cur];
	    rbuffer[bufsize2-1-dz] = rbuffer[bufsize-1-dz] = cur;
	  }
	}     // for (pos...
      }       // if (l<0)
    }         // for (iseq
    fprintf(stderr,"Total sites %ld, Correct sites %ld, Fraction correct %lf\n",cs->total,
	    cs->correct,(double)cs->correct/(double)cs->total);
  }

}
