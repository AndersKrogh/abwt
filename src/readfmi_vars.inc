int count_bwtfile=0;
char* bwtfile = NULL;
int count_Lookup=0;
char* Lookup = NULL;
int count_lookup=0;
char* lookup = NULL;
int count_sarray=2;
int* sarray=NULL;
int count_marray=2;
int* marray=NULL;
int count_format=0;
char* format = "b10i10l10 p5 B s30";
int count_PrintId=0;
int PrintId = 0;
int count_PrintSeq=0;
int PrintSeq = 0;
int count_forward=0;
int forward = 0;
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

OPT_STRUCT opt_struct[13] = {
	{OPTTYPE_SWITCH,VARTYPE_int,NULL,NULL,NULL,"---\nreadfmi reads an fmi file (made by mkbwt and mkfmi) and outputs\nsequences in fasta, seq ids, or a portion of the suffix array.\n\nFirst arg is index name and then options follows.\n\nReconstruct sequences (linelength 70):\nreadfmi <indexname> -p 70\n\nReconstruct suffix array from position 1000 to 1000+200:\nreadfmi <indexname> -s 1000 200\n\nPrint ids and seq info:\nreadfmi <indexname> -i\n\nIt is very slow and mostly for debugging\n\nSee options below\n---\n"},
	{OPTTYPE_ARG,VARTYPE_charS,(void *)&bwtfile,&count_bwtfile,"|bwtfile|","      Name of file containing BWT (output from mkfmi)"},
	{OPTTYPE_VALUE,VARTYPE_charS,(void *)&Lookup,&count_Lookup,"|Lookup|L|","      Look up the suffix interval for a sequence. Outputs the sequence, the start\n      of the SI and the length of the SI. Values of 0 0 means no match.\n      For example \"-l AACC\" will look up the SI for AACC."},
	{OPTTYPE_VALUE,VARTYPE_charS,(void *)&lookup,&count_lookup,"|lookup|l|","      As -L, but also prints the correponding suffix array (as -s). Use -m to add\n      margins."},
	{OPTTYPE_VALUE,VARTYPE_intS,(void *)&sarray,&count_sarray,"|sarray|s|","      Print suffix array interval. Takes two integers: start and length. Print the\n      interval starting at the first value with a length of the second value.\n      For instance \"-s 100 20\" will print it starting at 100 with length 20.\n      Values are ignored if -l or -L are used."},
	{OPTTYPE_VALUE,VARTYPE_intS,(void *)&marray,&count_marray,"|marray|m|","      Extra margin for the array you want to print. Upper and lower. Let us say you\n      want to print and interval from 100 and a length of 20. Then if you have a\n      margin set to 10 and 5 (-m 10 5), you get an interval from 90 with a length\n      of 35. So \"-s 100 20 -m 10 5\" is identical to \"-s 90 35\"."},
	{OPTTYPE_VALUE,VARTYPE_charS,(void *)&format,&count_format,"|format|f|","      Format telling what to print.\n        b is SA index.\n        i is position in sequence.\n        l is seq length (negative if reverse strand).\n        f is fmi value for BWT character at this position.\n        p is prefix (except BWT character)\n        s is suffix\n      All the above can be followed by a field width.\n        B is BWT charater (last char in prefix)\n        n is sequence name (should always be last - there is no field width)\n      Example: \"b9f9  p5 B s50 n\" prints SA index, FMI value (both 9 chars), a\n      prefix of 5 chars (p5), the BWT letter (B), 50 chars of the suffix (s50),\n      and the seq id (n)."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&PrintId,&count_PrintId,"|PrintId|i|","      Print sequence information"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&PrintSeq,&count_PrintSeq,"|PrintSeq|p|","      Print all sequences in fastaformat (value is line length). If value is\n      negative (!), print ID, space, sequence on single line. SLOW!"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&forward,&count_forward,"|forward|a|","      Takes integer suffix array index. Reconstruct sequence by forward recursion\n      from that index."},
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
