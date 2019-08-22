int count_indexfile=0;
char* indexfile = NULL;
int count_minlen=0;
int minlen = 25;
int count_maxlen=0;
int maxlen = 45;
int count_min_interval=0;
int min_interval = 10;
int count_maxfrac=0;
double maxfrac = 0.8;
int count_printall=0;
int printall = 0;
int count_prefix=0;
char* prefix = "Repeat";
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

OPT_STRUCT opt_struct[11] = {
	{OPTTYPE_SWITCH,VARTYPE_int,NULL,NULL,NULL,"---\nexactRepeats takes index name as argument and searches for likely\nexact repeats in a length interval specified.\n\nProgram is made with CRISPR repeats in mind.\n\nIndex of reads/contigs/genomes must be built using makeabwt\n\nExample cmd line\n   exactRepeats index_file_name -minlen 20 -maxlen 50\n\nto search for repeats between 20 and 50\n\nSee options below\n---\n"},
	{OPTTYPE_ARG,VARTYPE_charS,(void *)&indexfile,&count_indexfile,"|indexfile|","      Name of file containing BWT/FMI index"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&minlen,&count_minlen,"|minlen|m|","      Integer specifying the minimum length of reapeats"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&maxlen,&count_maxlen,"|maxlen|M|","      Integer specifying the maximum length of reapeats"},
	{OPTTYPE_VALUE,VARTYPE_int,(void *)&min_interval,&count_min_interval,"|min_interval|i|","      Integer specifying the maximum length of reapeats"},
	{OPTTYPE_VALUE,VARTYPE_double,(void *)&maxfrac,&count_maxfrac,"|maxfrac|f|","      Select repeats where the maximum fraction of the neihboring letter is maxfrac."},
	{OPTTYPE_SWITCH,VARTYPE_int,(void *)&printall,&count_printall,"|printall|a|","      Print all repeats found with SUBSTRING keyword in id of those\n      that would not be printed without -a"},
	{OPTTYPE_VALUE,VARTYPE_charS,(void *)&prefix,&count_prefix,"|prefix|p|","      The prefix for the repeats discovered"},
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
