/*
//Options parser for MQM
//Danny Arends (c) 5-Aug-2009
//   -v      Be verbose
//   -d n    Debug level n (default 0)
//   -p s    Phenotype file
//   -g s    Genotype file
//   -m s    Marker file
//   -c s    Cofactor file (optional)
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void printoptionshelp(void){
	printf ("Commandline switches:\n");
	printf ("-h      		This help\n");
	printf ("-v      		Verbose (produce a lot of textoutput)\n");
	printf ("-p(INT) 		DebugLevel -d0,-d1\n");
	printf ("-p(FILE_NAME)	Phenotypes file in plain textformat\n");
	printf ("-g(FILE_NAME)	Genotypes file in plain textformat\n");
	printf ("-m(FILE_NAME)	Marker and Chromosome descriptionfile in plain textformat\n");
	printf ("-c(FILE_NAME)	Cofactors file in plain textformat\n");
 }

//Functions
void exitonerror(const char *msg){
	fprintf(stderr, msg);
	printoptionshelp();
	exit(1);
 }
  
 //Main function
int main (int argc, char **argv){
	//variables
	bool verboseflag = false;
	bool helpflag = false;
	int debuglevel = 0;
	char *phenofile = NULL;
	char *genofile = NULL;
	char *markerfile = NULL;
	char *coffile = NULL;       
	unsigned int index;
	signed int c;

//Parsing of arguments     
while ((c = getopt (argc, argv, "vd:h:p:g:m:c:")) != -1)
switch (c)
{
	case 'v':
		verboseflag = true;
	break;
	case 'h':
		helpflag = true;
	break;	
	case 'd': 
		debuglevel = atoi(optarg);
	break;
	case 'p':
		phenofile = optarg;
	break;
	case 'g':
		genofile = optarg;
	break;             
	case 'm':
		markerfile = optarg;
	break;
	case 'c':
		coffile = optarg;
	break;
	default:
		fprintf (stderr, "Unknown option character '%c'.\n", optopt);
}
if(helpflag){
	printoptionshelp();
	return 0;
}else{
	printf ("Options parsed for MQM:\n");
//Verbose & debug
	printf ("verboseflag = %d, debuglvl = %d\n",verboseflag, debuglevel);
//Needed files
	if(!phenofile) exitonerror("Please supply a phenotypefile.\n");
	printf ("Phenotypes = %s\n",phenofile);
	if(!genofile)  exitonerror("Please supply a genofile.\n");
	printf ("Genotypes = %s\n",genofile);
	if(!markerfile) exitonerror("Please supply a markerfile.\n");
	printf ("Markers = %s\n",markerfile);
//Optional files
	if(coffile != NULL){
		printf ("Cofactors = %s\n",coffile);
	}
//Warn people for non-existing options
	for (index = optind; index < argc; index++){
		printf ("Non-option argument %s\n", argv[index]);
	}
}
return 0;
}