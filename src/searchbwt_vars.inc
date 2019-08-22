int count_indexfile=0;
char* indexfile = NULL;
int count_qfile=0;
char* qfile = NULL;
int count_minLen=0;
int minLen = 5;
int count_maxMatch=0;
int maxMatch = 0;
int count_caseSens=0;
int caseSens = 0;
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

OPT_STRUCT opt_struct[9] = {
	{OPTTYPE_SWITCH,VARTYPE_int,NULL,NULL,NULL,"---\nsearchbwt takes index name as argument and searches for longest matches\n\nExample cmd line\n   searchbwt -L 5 database queryfile\n\nwhere -L sets the minimum length of a match\n\nSee options below\n---\n"},
	{OPTTYPE_ARG,VARTYPE_charS,(void *)&indexfile,&count_indexfile,"|indexfile|","      Name of file containing index (output from mkfmi)"},
	{OPTTYPE_ARG,VARTYPE_charS,(void *)&qfile,&count_qfile,"|qfile|","      Name of file with query sequences (stdin if no file is given)"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&minLen,&count_minLen,"|minLen|L|","      Minimum length of a match. Make sure to set it (e.g. 15 for DNA and 5\n      for proteins)"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&maxMatch,&count_maxMatch,"|maxMatch|m|","      Maximum number of matches returned. However, all matches of a certain\n      length are always returned. With a value of 1, ALL maximum length matces\n      are returned."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&caseSens,&count_caseSens,"|caseSens|c|","      The sequences are read case sensitive"},
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
