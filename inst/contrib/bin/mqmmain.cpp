/**********************************************************************
 *
 * mqmmain.cpp
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Pjotr Prins and Danny Arends
 * last modified December 2009
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

FILE *redirect_info = stdout;
int debuglevel = 0;

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
  int nmark;
  unsigned int npheno;
  int stepmin;
  int stepmax;
  unsigned int stepsize;
  unsigned int windowsize;
  double alpha;
  unsigned int maxiter;
  char estmap;
  unsigned int max_totalaugment; // always >= 0
  unsigned int max_indaugment;   // always >= 0
  double neglect_unlikely;
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


MQMMarkerMatrix readgenotype(const char* filename,const unsigned int nind,const unsigned int nmar,const bool verbose) {
  unsigned int j = 0;  //current marker
  unsigned int i = 0;  //current individual
  MQMMarkerMatrix genomarkers = newMQMMarkerMatrix(nmar,nind);
  ifstream myfstream(filename, ios::in);
  char c;
  while (!myfstream.eof() && i < nind) {
    if (j < nmar) {
      myfstream >> c;
      genomarkers[j][i] = (MQMMarker)c;
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
  // std::string markernames[nmar];

  ivector markerparent = newivector(nmar);		//Parental genotype
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof() && j < nmar) {
    myfstream >> markerchr[j];
    std::string markername;
    myfstream >> markername;
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
    // No silent failures!!
    Rprintf("File not found %s",filename);
    exit(1);
  }
}

void printhelp(void) {
  printf ("Commandline switches:\n");
  printf ("-h      		This help.\n");
  printf ("-v      		Verbose (produce a lot of textoutput).\n");
  printf ("-d(INT) 		DebugLevel -d0,-d1.\n");
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
  printf ("--minprob(FLOAT)	Drop genotypes more unlikely that minprob.\n");
}

//Functions
void exit_on_error(const char *msg) {
  info("EXIT ERROR: %s",msg);
  printhelp();
  exit(1);
}

void exit_on_error_gracefull(const char *msg) {
  info("EXIT ERROR: %s",msg);
  printhelp();
  exit(0);
}


bool selectivelygenotyped(const MQMMarkerMatrix markers,const ivector chr,const unsigned int nind, const unsigned int nmar){
  int currentchr =0;
  int count = 0;
  int countmissing = 0;
  for(unsigned int i = 0;i < nind; i++){
    for(unsigned int j = 0;j < nmar; j++){
      if(chr[j] > currentchr && currentchr != 0){
        if(count==countmissing){
          return TRUE;
        }else{
          count = 0;
          countmissing = 0;
        }
        currentchr = chr[j];
      }else{
        if(markers[j][i]==9){
          countmissing++;
        }
        count++;
      }
      
    }
    count=0;
    countmissing=0;
    currentchr=0;
  }
  return FALSE;
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
  {"minprob",    required_argument, 0, 'l'},
  {0, 0, 0, 0}
};


int main(int argc,char *argv[]) {
  printf("MQM standalone version\n");
  bool verbose = false;
  bool helpflag = false;
  unsigned int phenotype = 0; //analyse the first phenotype
  char *phenofile = NULL;
  char *genofile = NULL;
  char *markerfile = NULL;
  char *coffile = NULL;
  char *settingsfile = NULL;
  char *outputfile = NULL;
  struct algorithmsettings mqmalgorithmsettings;
  struct markersinformation mqmmarkersinfo;
  int index; // aligned with argc and optind
  signed int c;

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
      debug_trace("Option (a) smin: %d\n",mqmalgorithmsettings.stepmin);
    case 'b':
      mqmalgorithmsettings.stepmax = atoi(optarg);
      debug_trace("Option (b) smax: %d\n",mqmalgorithmsettings.stepmax);
    break;
    case 'n':
      mqmalgorithmsettings.stepsize = atoi(optarg);
      debug_trace("Option (n) ssize: %d\n",mqmalgorithmsettings.stepsize);
    break;
    case 'e':
      mqmalgorithmsettings.alpha = atof(optarg);
      debug_trace("Option (e) alpha: %f\n",mqmalgorithmsettings.alpha);
    break;
    case 'f':
      mqmalgorithmsettings.windowsize = atoi(optarg);
      debug_trace("Option (f) window: %d\n",mqmalgorithmsettings.windowsize);
    break;  
    case 'q':
      mqmalgorithmsettings.maxiter = atoi(optarg);
      debug_trace("Option (q) maxiter: %d\n",mqmalgorithmsettings.maxiter);
    break; 
    case 'i':
      mqmalgorithmsettings.estmap = optarg[0];
      debug_trace("Option (i) estmap: %d\n",mqmalgorithmsettings.estmap);
    break; 
    case 'j':
      mqmalgorithmsettings.max_totalaugment = atoi(optarg);
      debug_trace("Option (j) max_totalaugment: %d\n",mqmalgorithmsettings.max_totalaugment);
    break; 
    case 'k':
      mqmalgorithmsettings.max_indaugment = atoi(optarg);
      debug_trace("Option (k) max_indaugment: %d\n",mqmalgorithmsettings.max_indaugment);
    break; 
    case 'l':
      mqmalgorithmsettings.neglect_unlikely = atof(optarg);
      debug_trace("Option (l) minprob: %f\n",mqmalgorithmsettings.neglect_unlikely);
    break; 
    default:
      fprintf(stderr, "Unknown option character '%c'.\n", optopt);
  }
  if (helpflag) {
    printhelp();
    return 0;
  } else {
    //Check the output file
    if (outputfile) debug_trace("Output file specified: %s\n",outputfile);
    if (checkfileexists(outputfile)) debug_trace("Outputfile exists\n !!! overwriting previous outputfile !!!\n");
    // Open outputstream if specified - using C type for redirection
    FILE *fout = stdout;
    if (outputfile){
      fout = fopen(outputfile,"w");
      redirect_info = fout;
    }
    debug_trace ("Options for MQM:\n");
    //Verbose & debug
    debug_trace ("verbose = %d, debuglevel = %d\n",verbose, debuglevel);
    //Needed files
    if (!phenofile) exit_on_error("Please supply a phenotypefile argument.\n");
    if (!checkfileexists(phenofile)) exit_on_error("Phenotypefile not found on your filesystem.\n");
    debug_trace ("Phenotypefile = %s\n",phenofile);
    if (!genofile)  exit_on_error("Please supply a genofile argument.\n");
    if (!checkfileexists(genofile)) exit_on_error("Genotypefile not found on your filesystem.\n");
    debug_trace ("Genotypefile = %s\n",genofile);
    if (!markerfile) exit_on_error("Please supply a markerfile argument.\n");
    if (!checkfileexists(genofile)) exit_on_error("Markerfile not found on your filesystem.\n");
    debug_trace ("Markerfile = %s\n",markerfile);
    if (!settingsfile) exit_on_error("Please supply a settingsfile argument.\n");
    if (!checkfileexists(settingsfile)) exit_on_error("settingsfile not found on your filesystem.\n");
    debug_trace ("settingsfile = %s\n",settingsfile);
    //Optional files
    if (!coffile) {
      if (!checkfileexists(coffile)) {
        debug_trace("Cofactorfile not found on your filesystem.\n");
      } else {
        debug_trace("Cofactorfile = %s\n",coffile);
      }
    }

    //Warn people for non-existing options
    for (index = optind; index < argc; index++) {
      debug_trace("Non-option argument %s\n", argv[index]);
    }
    
    //Read in settingsfile
    mqmalgorithmsettings = loadmqmsetting(settingsfile,mqmalgorithmsettings,verbose);
    //Create large datastructures
    double **QTL;
    ivector chr = newivector(mqmalgorithmsettings.nmark);
    cvector cofactor = newcvector(mqmalgorithmsettings.nmark);
    vector mapdistance = newvector(mqmalgorithmsettings.nmark);
    vector pos = newvector(mqmalgorithmsettings.nmark);
    matrix pheno_value = newmatrix(mqmalgorithmsettings.npheno,mqmalgorithmsettings.nind);
    MQMMarkerMatrix markers= newMQMMarkerMatrix(mqmalgorithmsettings.nmark,mqmalgorithmsettings.nind);
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

    debug_trace("Genotypefile done\n");

    pheno_value = readphenotype(phenofile,mqmalgorithmsettings.nind,mqmalgorithmsettings.npheno,verbose);

    debug_trace("Phenotypefile done \n");

    mqmmarkersinfo = readmarkerfile(markerfile,mqmalgorithmsettings.nmark,verbose);
    chr = mqmmarkersinfo.markerchr;
    pos = mqmmarkersinfo.markerdistance;

    debug_trace("Markerposition file done\n");

    //Determine how many chromosomes we have
    int max_chr=0;
    for (unsigned int m=0; m < (unsigned int) mqmalgorithmsettings.nmark; m++) {
      if (max_chr<chr[m]) {
        max_chr = chr[m];
      }
    }
    debug_trace("# %d Chromosomes\n",max_chr);
    //Create a QTL object holding all our output location
    int locationsoutput = 3*max_chr*(((mqmalgorithmsettings.stepmax)-(mqmalgorithmsettings.stepmin))/ (mqmalgorithmsettings.stepsize));
    QTL = newmatrix(1,locationsoutput);
    //initialize cofactors to 0 and mapdistances to UNKNOWN Cm
    for (unsigned int i=0; i< (unsigned int) mqmalgorithmsettings.nmark; i++) {
      cofactor[i] = '0';
      mapdistance[i]=POSITIONUNKNOWN;
      mapdistance[i]=pos[i];
    }

    if (coffile)
      set_cofactors = readcofactorfile(coffile,&cofactor,mqmalgorithmsettings.nmark,verbose);
    if (set_cofactors > 0) {
      backwards = 1;
    }
    
    //Initialize an empty individuals list
    for (unsigned int i=0; i< mqmalgorithmsettings.nind; i++) {
      INDlist[i] = i;
    }

    int nind = mqmalgorithmsettings.nind;
    int augmentednind = mqmalgorithmsettings.nind;
    
    //<dataaugmentation>
     if(mqmalgorithmsettings.max_totalaugment <= mqmalgorithmsettings.nind) 
      exit_on_error_gracefull("Augmentation parameter conflict max_augmentation <= individuals");   
    //int testje = calculate_augmentation(mqmalgorithmsettings.nind,mqmalgorithmsettings.nmark,markers,crosstype);
    
    if(selectivelygenotyped(markers,chr,mqmalgorithmsettings.nind,mqmalgorithmsettings.nmark)){
      fprintf(fout,"Warning: Selective genotyped set, only including most likely (neglect_unlikely set to 1)\n");
      mqmalgorithmsettings.neglect_unlikely = 1;
    }

   mqmaugmentfull(&markers,&nind,&augmentednind,&INDlist,mqmalgorithmsettings.neglect_unlikely, mqmalgorithmsettings.max_totalaugment, 
   mqmalgorithmsettings.max_indaugment,&pheno_value,mqmalgorithmsettings.nmark,chr,mapdistance,1,crosstype,verbose);

    // Start scanning for QTLs
    double logL = analyseF2(augmentednind, &mqmalgorithmsettings.nmark, &cofactor, (MQMMarkerMatrix)markers, pheno_value[phenotype], backwards,QTL, &mapdistance,&chr,0,0,mqmalgorithmsettings.windowsize,
              mqmalgorithmsettings.stepsize,mqmalgorithmsettings.stepmin,mqmalgorithmsettings.stepmax,mqmalgorithmsettings.alpha,mqmalgorithmsettings.maxiter,nind,&INDlist,mqmalgorithmsettings.estmap,crosstype,false,verbose);
    // Write final QTL profile (screen and file)
    if (!isinf(logL) && !isnan(logL)) {
      for (int q=0; q<locationsoutput; q++) {
        double qtlvalue = QTL[0][q];
        fprintf(fout,"%5d\t",q);
        // The following prints a 'standardized' value on Windows and Unix for regression tests (for nan and inf)
        if(isnan(qtlvalue)){
          fprintf(fout,"       NAN\n");
        }else{
          if(isinf(qtlvalue)){
            fprintf(fout,"  INFINITE\n");
          }else{
            fprintf(fout,"%.3f\n",ftruncate3(QTL[0][q]));
          }
        }
      }
    }

    freevector((void *)cofactor);
    freevector((void *)mapdistance);
    freematrix((void **)markers,mqmalgorithmsettings.nmark);
    freematrix((void **)pheno_value,mqmalgorithmsettings.npheno);
    freevector((void *)chr);
    freevector((void *)INDlist);
    freevector((void *)pos);
    freematrix((void **)QTL,1);
    if (outputfile) fclose(fout);
    return 0;
  }
}
