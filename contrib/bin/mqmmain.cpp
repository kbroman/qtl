/**********************************************************************
 *
 * mqmmain.cpp - standalone MQM edition
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Apr,2009
 * first written Feb, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 *
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 * Contains: R_scanMQM, scanMQM
 *
 **********************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "mqm.h"
#include <getopt.h>


using namespace std;

bool checkfileexists(const char *filename) {
  ifstream myfile;
  bool exists;
  myfile.open(filename);
  exists = myfile.is_open();
  myfile.close();
  return exists;
}

struct algorithmsettings {
  unsigned int nind;
  unsigned int nmark;
  unsigned int npheno;
  int stepmin;
  int stepmax;
  unsigned int stepsize;
  unsigned int windowsize;
  double alpha;
  unsigned int maxiter;
  char estmap;
  int max_totalaugment;
  int max_indaugment;
  int neglect_unlikely;
  char suggestedcross;
};

struct markersinformation {
  ivector markerchr;
  vector markerdistance;
  ivector markerparent;
};

struct algorithmsettings loadmqmsetting(const char* filename,const algorithmsettings commandline, bool verbose) {
  algorithmsettings runsettings=commandline;
  if (verbose) printf("INFO: Loading settings from file\n");
  ifstream instream(filename, ios::in);
  instream >> runsettings.nind >> runsettings.nmark >> runsettings.npheno;
  //instream >> runsettings.stepmin >> runsettings.stepmax >> runsettings.stepsize;
  //instream >> runsettings.windowsize >> runsettings.alpha;
  //instream >> runsettings.maxiter >> runsettings.estmap;
  //instream >> runsettings.max_totalaugment >> runsettings.max_indaugment >> runsettings.neglect_unlikely;
  instream >> runsettings.suggestedcross;
  if (verbose) {
    Rprintf("number of individuals: %d\n",runsettings.nind);
    Rprintf("number of markers: %d\n",runsettings.nmark);
    Rprintf("number of phenotypes: %d\n",runsettings.npheno);
    //Rprintf("stepmin: %d\n",runsettings.stepmin);
    //Rprintf("stepmax: %d\n",runsettings.stepmax);
    //Rprintf("stepsize: %d\n",runsettings.stepsize);
    //Rprintf("windowsize for dropping qtls: %d\n",runsettings.windowsize);
    //Rprintf("Alpha level considered to be significant: %f\n",runsettings.alpha);
    //Rprintf("Max iterations using EM: %d\n",runsettings.maxiter);
    //Rprintf("Re-estimating map-positions: %c\n",runsettings.estmap);
    Rprintf("Suggested cross: %c\n",runsettings.suggestedcross);
    //Rprintf("Data-augmentation parameters: max:%d maxind:%d neglect:%d\n",runsettings.max_totalaugment,runsettings.max_indaugment,runsettings.neglect_unlikely);
  }
  return runsettings;
}


cmatrix readgenotype(const char* filename,const unsigned int nind,const unsigned int nmar,const bool verbose) {
  unsigned int j = 0;  //current marker
  unsigned int i = 0;  //current individual
  cmatrix genomarkers = newcmatrix(nmar,nind);
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof()) {
    if (j < nmar) {
      myfstream >> genomarkers[j][i];
      j++;
    } else {
      j = 0;
      i++;
    }
  }
  if (verbose) Rprintf("Individuals: %d\n",i);
  myfstream.close();
  return genomarkers;
}

matrix readphenotype(const char* filename,const unsigned int nind,const unsigned int nphe,const bool verbose) {
  unsigned int p = 0;  // current phenotype
  unsigned int i = 0;  //current individual
  matrix phenovalues = newmatrix(nphe,nind);
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof()) {
    if (p < nphe) {
      myfstream >> phenovalues[p][i];
      p++;
    } else {
      p = 0;
      i++;
    }
  }
  if (verbose) Rprintf("Individuals: %d\n",i);
  myfstream.close();
  return phenovalues;
}

struct markersinformation readmarkerfile(const char* filename,const unsigned int nmar,const bool verbose) {
  unsigned int j = 0;  //current marker
  markersinformation info;
  ivector markerchr = newivector(nmar);
  vector markerdistance= newvector(nmar);
  std::string markernames[nmar];
  ivector markerparent = newivector(nmar);		//Parental genotype
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof() && j < nmar) {
    myfstream >> markerchr[j];
    myfstream >> markernames[j];
    myfstream >> markerdistance[j];
    markerparent[j] = 12;
    //if (verbose) Rprintf("Marker %d: %s %d %f\n",j,markernames[j].c_str(),markerchr[j],markerdistance[j]);
    j++;
  }
  if (verbose) Rprintf("Markers: %d\n",j);
  myfstream.close();
  info.markerchr=markerchr;
  info.markerdistance=markerdistance;
  info.markerparent=markerparent;
  return info;
}

unsigned int readcofactorfile(const char* filename,cvector *cofactors,const unsigned int nmar,const bool verbose) {
  //Cofactor is pass by value
  if (checkfileexists(filename)) {
    unsigned int j = 0;     //current marker
    unsigned int num = 0;   //number of co-factors encountered
    ifstream myfstream(filename, ios::in);
    while (!myfstream.eof()) {
      myfstream >> (*cofactors)[j];
      if ((*cofactors)[j]!='0') num++;
      j++;
    }
    myfstream.close();
    if (verbose) Rprintf("Cofactors/Markers: %d/%d\n",num,j);
    return num;
  } else {
    return 0;
  }
}

void printhelp(void) {
  printf ("Commandline switches:\n");
  printf ("-h      		This help.\n");
  printf ("-v      		Verbose (produce a lot of textoutput).\n");
  printf ("-p(INT) 		DebugLevel -d0,-d1.\n");
  printf ("-t(INT) 		Phenotype under analysis.\n");
  printf ("-p(FILE_NAME)	Phenotypes file in plain textformat.\n");
  printf ("-g(FILE_NAME)	Genotypes file in plain textformat.\n");
  printf ("-m(FILE_NAME)	Marker and Chromosome descriptionfile in plain textformat.\n");
  printf ("-s(FILE_NAME)	Settings file in plain textformat.\n");
  printf ("-c(FILE_NAME)	Optional Cofactors file to do backward elimination on in plain textformat.\n");
  printf ("-o(FILE_NAME)	Optional output file to save MQM-QTL mapping results in.\n");
  printf ("--smin(INT)	Start of mapping in Cm.\n");
  printf ("--smax(INT)	End of mapping in Cm.\n");
  printf ("--sstep(INT)	Stepsize of mapping in Cm.\n");
  printf ("--alpha(FLOAT)	Significance level.\n");
  printf ("--window(INT)	Windowsize for dropping QTLs in Cm.\n");
  printf ("--maxiter(INT)	Maximum number of EM iterations.\n");
  printf ("--estmap(CHAR)	Reestimate marker positions y/n?.\n");
  printf ("--maugment(INT)	Maximum size of augmented dataset.\n");
  printf ("--miaugment(INT)	Maximum number of individual replications inside a dataset.\n");
  printf ("--Neglect(INT)	Drop genotypes more unlikely that neglect/1000.\n");
}

//Functions
void exit_on_error(const char *msg) {
  fprintf(stderr, msg);
  printhelp();
  exit(1);
}

static struct option long_options[] = {
  {"smin",  required_argument, 0, 'a'},
  {"smax",  required_argument, 0, 'b'},
  {"sstep",  required_argument, 0, 'n'},
  {"alpha",  required_argument, 0, 'e'},
  {"window",    required_argument, 0, 'f'},
  {"maxiter",    required_argument, 0, 'q'},
  {"estmap",    required_argument, 0, 'i'},
  {"maugment",    required_argument, 0, 'j'},
  {"miaugment",    required_argument, 0, 'k'},
  {"neglect",    required_argument, 0, 'l'},
  {0, 0, 0, 0}
};


int main(int argc,char *argv[]) {
  Rprintf("MQM standalone version\n");
  bool verbose = false;
  bool helpflag = false;
  int debuglevel = 0;
  unsigned int phenotype = 0; //analyse the first phenotype
  char *phenofile = NULL;
  char *genofile = NULL;
  char *markerfile = NULL;
  char *coffile = NULL;
  char *settingsfile = NULL;
  char *outputfile = NULL;
  struct algorithmsettings mqmalgorithmsettings;
  struct markersinformation mqmmarkersinfo;
  unsigned int index;
  signed int c;
  ofstream outstream; //Could be needed when -o is set

  int option_index = 0;
  //Parsing of arguments
  while ((c = getopt_long(argc, argv, "vd:h:p:g:m:c:s:t:o:a:b:e:f:q:i:j:k:l:",long_options, &option_index)) != -1)
    switch (c) {
    case 'v':
      verbose = true;
      break;
    case 'h':
      helpflag = true;
      break;
    case 'd':
      debuglevel = atoi(optarg);
      break;
    case 't':						//1 phenotype at a time
      phenotype = atoi(optarg);
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
    case 's':
      settingsfile = optarg;
      break;
    case 'c':
      coffile = optarg;
      break;
    case 'o':
      outputfile = optarg;
      break;      
    case 'a':
      mqmalgorithmsettings.stepmin = atoi(optarg);
      printf("Option (a) smin: %d\n",mqmalgorithmsettings.stepmin);
    case 'b':
      mqmalgorithmsettings.stepmax = atoi(optarg);
      printf("Option (b) smax: %d\n",mqmalgorithmsettings.stepmax);
    break;
    case 'n':
      mqmalgorithmsettings.stepsize = atoi(optarg);
      printf("Option (n) ssize: %d\n",mqmalgorithmsettings.stepsize);
    break;
    case 'e':
      mqmalgorithmsettings.alpha = atof(optarg);
      printf("Option (e) alpha: %f\n",mqmalgorithmsettings.alpha);
    break;
    case 'f':
      mqmalgorithmsettings.windowsize = atoi(optarg);
      printf("Option (f) window: %d\n",mqmalgorithmsettings.windowsize);
    break;  
    case 'q':
      mqmalgorithmsettings.maxiter = atoi(optarg);
      printf("Option (q) maxiter: %d\n",mqmalgorithmsettings.maxiter);
    break; 
    case 'i':
      mqmalgorithmsettings.estmap = optarg[0];
      printf("Option (i) estmap: %d\n",mqmalgorithmsettings.estmap);
    break; 
    case 'j':
      mqmalgorithmsettings.max_totalaugment = atoi(optarg);
      printf("Option (j) max_totalaugment: %d\n",mqmalgorithmsettings.max_totalaugment);
    break; 
    case 'k':
      mqmalgorithmsettings.max_indaugment = atoi(optarg);
      printf("Option (k) max_indaugment: %d\n",mqmalgorithmsettings.max_indaugment);
    break; 
    case 'l':
      mqmalgorithmsettings.neglect_unlikely = atoi(optarg);
      printf("Option (l) neglect: %d\n",mqmalgorithmsettings.neglect_unlikely);
    break; 
    default:
      fprintf (stderr, "Unknown option character '%c'.\n", optopt);
  }
  if (helpflag) {
    printhelp();
    return 0;
  } else {
    printf ("Options for MQM:\n");
    //Verbose & debug
    printf ("verbose = %d, debuglevel = %d\n",verbose, debuglevel);
    //Needed files
    if (!phenofile) exit_on_error("Please supply a phenotypefile argument.\n");
    if (!checkfileexists(phenofile)) exit_on_error("Phenotypefile not found on your filesystem.\n");
    printf ("Phenotypefile = %s\n",phenofile);
    if (!genofile)  exit_on_error("Please supply a genofile argument.\n");
    if (!checkfileexists(genofile)) exit_on_error("Genotypefile not found on your filesystem.\n");
    printf ("Genotypefile = %s\n",genofile);
    if (!markerfile) exit_on_error("Please supply a markerfile argument.\n");
    if (!checkfileexists(genofile)) exit_on_error("Markerfile not found on your filesystem.\n");
    printf ("Markerfile = %s\n",markerfile);
    if (!settingsfile) exit_on_error("Please supply a settingsfile argument.\n");
    if (!checkfileexists(settingsfile)) exit_on_error("settingsfile not found on your filesystem.\n");
    printf ("settingsfile = %s\n",settingsfile);
    //Optional files
    if (!coffile) {
      if (!checkfileexists(coffile)) {
        printf("Cofactorfile not found on your filesystem.\n");
      } else {
        printf ("Cofactorfile = %s\n",coffile);
      }
    }
    //Check the output file
    if (outputfile) printf("Output file specified: %s\n",outputfile);
    if (checkfileexists(outputfile)) printf("Outputfile exists\n !!! overwriting previous outputfile !!!\n");
    //Warn people for non-existing options
    for (index = optind; index < argc; index++) {
      printf ("Non-option argument %s\n", argv[index]);
    }
    //Read in settingsfile
    mqmalgorithmsettings = loadmqmsetting(settingsfile,mqmalgorithmsettings,verbose);
    //Create large datastructures
    double **QTL;
    ivector f1genotype = newivector(mqmalgorithmsettings.nmark);
    ivector chr = newivector(mqmalgorithmsettings.nmark);
    cvector cofactor = newcvector(mqmalgorithmsettings.nmark);
    vector mapdistance = newvector(mqmalgorithmsettings.nmark);
    vector pos = newvector(mqmalgorithmsettings.nmark);
    matrix pheno_value = newmatrix(mqmalgorithmsettings.npheno,mqmalgorithmsettings.nind);
    cmatrix markers= newcmatrix(mqmalgorithmsettings.nmark,mqmalgorithmsettings.nind);
    ivector INDlist= newivector(mqmalgorithmsettings.nind);
    //Some additional variables
    int set_cofactors=0;			//Markers set as cofactors
    int backwards=0;				//Backward elimination ?
    MQMCrossType crosstype = CUNKNOWN;
    if(mqmalgorithmsettings.suggestedcross=='F'){
      crosstype = CF2;	//Crosstype
    }
    if(mqmalgorithmsettings.suggestedcross=='B'){
      crosstype = CBC;	//Crosstype
    }
    if(mqmalgorithmsettings.suggestedcross=='R'){
      crosstype = CRIL;	//Crosstype
    }
    //Here we know what we need so we can start reading in files with the new loader functions
    markers = readgenotype(genofile,mqmalgorithmsettings.nind,mqmalgorithmsettings.nmark,verbose);

    if (verbose) Rprintf("Genotypefile done\n");

    pheno_value = readphenotype(phenofile,mqmalgorithmsettings.nind,mqmalgorithmsettings.npheno,verbose);

    if (verbose) Rprintf("Phenotypefile done \n");

    mqmmarkersinfo = readmarkerfile(markerfile,mqmalgorithmsettings.nmark,verbose);
    chr = mqmmarkersinfo.markerchr;
    pos = mqmmarkersinfo.markerdistance;
    f1genotype = mqmmarkersinfo.markerparent;

    if (verbose) Rprintf("Markerposition file done\n");

    //Determin how many chromosomes we have
    unsigned int max_chr=0;
    for (int m=0; m < mqmalgorithmsettings.nmark; m++) {
      if (max_chr<chr[m]) {
        max_chr = chr[m];
      }
    }
    if (verbose)  Rprintf("# %d Chromosomes\n",max_chr);
    //Create a QTL object holding all our output location
    int locationsoutput = 3*max_chr*(((mqmalgorithmsettings.stepmax)-(mqmalgorithmsettings.stepmin))/ (mqmalgorithmsettings.stepsize));
    QTL = newmatrix(1,locationsoutput);
    //initialize cofactors to 0 and mapdistances to 999.0 Cm
    for (int i=0; i< mqmalgorithmsettings.nmark; i++) {
      cofactor[i] = '0';
      mapdistance[i]=999.0;
      mapdistance[i]=pos[i];
      //if (verbose) Rprintf("Distance %d, %f\n",i,mapdistance[i]);
    }

    //Danny: Cofactors are now read-in. the output with cofactors.txt set is not equal to MQM_test0.txt
    //MQM_test0.txt says it uses cofactors but it doesn't, because they are not eliminated
    //The message: "INFO: Marker XX is dropped, resulting in logL of reduced model = -8841.452934" is missing
    //Also the result of MQM without cofactors is equal
    set_cofactors = readcofactorfile(coffile,&cofactor,mqmalgorithmsettings.nmark,verbose);
    if (set_cofactors > 0) {
      backwards = 1;
      if (verbose) Rprintf("%d markers with cofactors. Backward elimination enabled\n",set_cofactors);
    }
    
    //Initialize an empty individuals list
    for (int i=0; i< mqmalgorithmsettings.nind; i++) {
      INDlist[i] = i;
    }

    //<dataaugmentation>
    //Variables for the returned augmented markers,phenotype,individualmapping
    cmatrix newmarkerset;
    vector new_y;
    ivector new_ind;
    int nind = mqmalgorithmsettings.nind;
    int augmentednind = mqmalgorithmsettings.nind;
    cvector position = locate_markers(mqmalgorithmsettings.nmark,chr);
    vector r = recombination_frequencies(mqmalgorithmsettings.nmark, position, mapdistance);
    augmentdata(markers, pheno_value[phenotype], &newmarkerset, &new_y, &new_ind, &nind, &augmentednind,  mqmalgorithmsettings.nmark, position, r, mqmalgorithmsettings.max_totalaugment, mqmalgorithmsettings.max_indaugment, mqmalgorithmsettings.neglect_unlikely, crosstype, verbose);
    if (verbose) Rprintf("Settingsnind: %d nind: %d augmentednind: %d\n",mqmalgorithmsettings.nind,nind,augmentednind);
    //Now to set the values we got back into the variables
    pheno_value[phenotype] = new_y;
    INDlist = new_ind;
    markers = newmarkerset;
    //Cleanup dataaugmentation:
    freevector((void *)position);
    freevector((void *)r);
    // </dataaugmentation>
    // Uncomment to inspect the augmented dataset
    //for (int m=0; m < mqmalgorithmsettings.nmark; m++) {
    //  for (int i=0; i < mqmalgorithmsettings.nind; i++) {
    //    if(verbose) Rprintf("%c ",markers[m][i]);
    //  }
    //  if(verbose) Rprintf("\n");
    //}
    
    //Missing values create an augmented set,
    analyseF2(mqmalgorithmsettings.nind, mqmalgorithmsettings.nmark, &cofactor, markers, pheno_value[phenotype], f1genotype, backwards,QTL, &mapdistance,&chr,0,0,mqmalgorithmsettings.windowsize,
              mqmalgorithmsettings.stepsize,mqmalgorithmsettings.stepmin,mqmalgorithmsettings.stepmax,mqmalgorithmsettings.alpha,mqmalgorithmsettings.maxiter,augmentednind,&INDlist,mqmalgorithmsettings.estmap,crosstype,0,verbose);

    // Open outputstream if specified
    if (outputfile){
      outstream.open(outputfile);
    }
    //Write final QTL profile (screen and file)
    for (int q=0; q<locationsoutput; q++) {
      if (outputfile) outstream << q << "\t" << QTL[0][q] << "\n";
      //if (verbose) Rprintf("%5d%10.5f\n",q,QTL[0][q]);
      
    }
    //close the outputstream
    if (outputfile) outstream.close();
    
    //Cleanup
    freevector((void *)f1genotype);
    freevector((void *)cofactor);
    freevector((void *)mapdistance);
    freematrix((void **)markers,mqmalgorithmsettings.nmark);
    freematrix((void **)pheno_value,mqmalgorithmsettings.npheno);
    freevector((void *)chr);
    freevector((void *)INDlist);
    freevector((void *)pos);
    freematrix((void **)QTL,1);
    return 0;
  }
}
