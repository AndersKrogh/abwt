int count_infilename=0;
char* infilename = NULL;
int count_outfilename=0;
char* outfilename = NULL;
int count_Alphabet=0;
char* Alphabet = "protein";
int count_nThreads=0;
int nThreads = 2;
int count_length=0;
double length = 0;
int count_checkpoint=0;
int checkpoint = 5;
int count_revComp=0;
int revComp = 0;
int count_reverse=0;
int reverse = 0;
int count_revsort=0;
int revsort = 0;
int count_word_len=0;
int word_len = 0;
int count_help=0;
int help = 0;
#ifndef VERSION
#define VERSION "0"
#endif
int count_OPTversion=0;
int OPTversion = 0;
#ifndef PRGNAME
#define PRGNAME ""
#endif

OPT_STRUCT opt_struct[14] = {
	{OPTTYPE_SWITCH,VARTYPE_int,NULL,NULL,NULL,"---\nmakeabwt takes a fasta file as argument and calculates the BWT and FM index.\nOutput file name given with -o (defaults to the input file name)\n\nExample cmd line\n   makeabwt -a DNA -o outputname infilename.fsa\nor for proteins (default alphabet)\n   makeabwt -o outputname infilename.fsa\nor for some other alphabet\n   makeabwt -a abcdefgHIJK -o outputname infilename.fsa\n\nIt can also take sequences on stdin, in which case you have to give\nthe filesize in millions of letters (rounded up), e.g. -l 3000\ncorresponding to 3 billion letters.\n\nFiles are created with outputname followed by the extension fmi\n\nWhile running, the program creates a couple of temporary files, which\nare deleted at the end.\n\nSee options below\n---\n"},
	{OPTTYPE_ARG,VARTYPE_charS,(void *)&infilename,&count_infilename,"|infilename|","      Name of an input file (stdin if no file is given, in which case you\n      need to give length)"},
	{OPTTYPE_VALUE,VARTYPE_charS,(void *)&outfilename,&count_outfilename,"|outfilename|o|","      Name of output. Several files with different extensions are produced\n      (if not given, input file name is used)."},
	{OPTTYPE_VALUE,VARTYPE_charS,(void *)&Alphabet,&count_Alphabet,"|Alphabet|a|","      Alphabet used. Instead of alphabet you can specify protein, DNA, IUPAC, or\n      RNA, in which case the alphabet is ACDEFGHIKLMNPQRSTVWYX, ACGTN,\n      ACGTNRYSWKMBDHVN, or ACGUN.\n      You can give qualifiers: /wildcard, /casesens, /variants to add\n      wildcard character (N or X), make a case sensitive index and add '|'\n      for variant encoding, repectively. These can be abbreviated.\n      Example: DNA/wild/var will give the alphabet ACGT|N"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&nThreads,&count_nThreads,"|nThreads|t|","      Number of threads"},
	{OPTTYPE_VALUE,VARTYPE_double,(void *)&length,&count_length,"|length|l|","      Length of concatenated sequence in millions (one decimal, round up).\n      Used when reading from stdin. If file name is given, length is estimated\n      from file size and length needs not be specified."},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&checkpoint,&count_checkpoint,"|checkpoint|e|","      Exponent for suffix array checkpoints. There is a checkpoint for every\n      2^e points. Value around 5 is a good compromise between speed and memory."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&revComp,&count_revComp,"|revComp|r|","      Add reverse complement sequences, so both strands are present. Works only\n      for DNA."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&reverse,&count_reverse,"|reverse|R|","      Make the BWT of the reversed sequences."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&revsort,&count_revsort,"|revsort|s|","      The termination symbols sorts as reverse sequences. This will make the\n      BWT more compressible."},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&word_len,&count_word_len,"|word_len|k|","      Word length for buckets in sorting. If k=0, it is set internally to\n      k=4 for alphabets with at most 8 characters, k=3 for alphabets with at most\n      16 characters and k=2 otherwise."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&help,(void *)&count_help,"|help|h|","      Prints summary of options and arguments"},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&OPTversion,(void *)&count_OPTversion,"|version|v|","      Prints version number"},
	{0,0,NULL,NULL,PRGNAME,VERSION}
};

// This is how you call the cmd line parser etc
// int main(int argc, char **argv) {
//     OPT_read_cmdline(opt_struct, argc, argv);
//     if (help) OPT_help(opt_struct);
//     OPT_print_vars(stderr, opt_struct, "# ", 0);
// }
