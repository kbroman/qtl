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

int main (int argc, char **argv) {
  bool verboseflag = false;
  int debuglevel = 0;
  char *phenofile = NULL;
  char *genofile = NULL;
  char *markerfile = NULL;
  char *coffile = NULL;
  int index;
  int c;

//Parsing of arguments
  while ((c = getopt (argc, argv, "vd:p:g:m:c:")) != -1)
    switch (c) {
    case 'v':
      verboseflag = true;
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
  printf ("Options parsed for MQM:\n");
//Verbose & debug
  printf ("verboseflag = %d, bflag = %d\n",verboseflag, debuglevel);
//Needed files
  if (phenofile==NULL) {
    fprintf (stderr, "Please supply a phenotypefile.\n");
    return 0;
  }
  printf ("Phenotypes = %s\n",phenofile);
  if (genofile==NULL) {
    fprintf (stderr, "Please supply a genofile.\n");
    return 0;
  }
  printf ("Genotypes = %s\n",genofile);
  if (markerfile==NULL) {
    fprintf (stderr, "Please supply a markerfile.\n");
    return 0;
  }
  printf ("Markers = %s\n",markerfile);
//Optional files
  if (coffile != NULL) {
    printf ("Cofactors = %s\n",coffile);
  }
//Warn people for non-existing options
  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);
  return 0;
}
