/**********************************************************************
 *
 * opt_main.cpp
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Pjotr Prins and Danny Arends
 * last modified August 2009
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
#include <iostream>
#include <fstream>

struct algorithmsettings {
  unsigned int nind;
  unsigned int nmark;
  unsigned int npheno;
  int stepmin;
  int stepmax;
  unsigned int stepsize;
  unsigned int windowsize;
  unsigned int alpha;
  unsigned int maxiter;
};

struct mqmalgorithmsettings loadmqmsetting(const char* filename) {
  algorithmsettings runsettings;
  if (verbose) printf("INFO: Loading settings from file\n");
  ifstream settingsstream(filename, ios::in);
  settingsstream >> runsettings.nind;
  settingsstream >> runsettings.nmark; //NEW we dont want to guess this: Should be Added to testfile
  settingsstream >> runsettings.npheno;
  settingsstream >> runsettings.stepmin;
  settingsstream >> runsettings.stepmax;
  settingsstream >> runsettings.stepsize;
  settingsstream >> runsettings.windowsize;
  settingsstream >> runsettings.alpha;
  settingsstream >> runsettings.maxIter;
  if (verbose) {
    Rprintf("number of individuals: %d\n",runsettings.nind);
    Rprintf("number of markers: %d\n",runsettings.nind);
    Rprintf("number of phenotypes: %d\n",runsettings.npheno);
    Rprintf("stepmin: %d\n",runsettings.stepmin);
    Rprintf("stepmax: %d\n",runsettings.stepmax);
    Rprintf("stepsize: %d\n",runsettings.stepsize);
    Rprintf("windowsize for dropping qtls: %d\n",runsettings.windowsize);
    Rprintf("Alpha level considered to be significant: %f\n",runsettings.alpha);
    Rprintf("Max iterations using EM: %d\n",runsettings.maxiter);
  }
  return runsettings;
}

cmatrix readgenotype(const char* filename,const unsigned int nind,const unsigned int nmar) {
  unsigned int cmarker = 0;
  unsigned int cindividual = 0;
  cmatrix genomarkers = newcmatrix(nphe)(nind);
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof()) {
    if (cmarker < nmar) {
      myfstream >> genomarkers[cmarker][cindividual];
    } else {
      cmarker = 0;
    }
    cindividual++;
  }
  myfstream.close();
  return genomarkers;
}

matrix readphenotype(const char* filename,const unsigned int nind,const unsigned int nphe) {
  unsigned int cphenotype = 0;
  unsigned int cindividual = 0;
  matrix phenovalues = newmatrix(nphe)(nind);
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof()) {
    if (cphenotype < nphe) {
      myfstream >> phenovalues[cphenotype][cindividual];
    } else {
      cphenotype = 0;
    }
    cindividual++;
  }
  myfstream.close();
  //TODO how to pass large structures or multiple arrays that we read in
  return phenovalues;
}

void readmarkerfile(const char* filename,const unsigned int nmar) {
  unsigned int cmarker = 0;
  markerchr = newivector(nmar);			//NEW !!! chr-> should be added to test
  markerdistance= newvector(nmar);		//pos
  markernames = char*[nmar];				//NEW !!!
  markerparent = newivector(nmar);		//f1genotype
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof()) {
    myfstream >> markerchr[x];
    myfstream >> markerdistance[x];
    myfstream >> markernames[x];
    markerparent[x] = 12;
  }
  //TODO get arrays back to main
  myfstream.close();
}

unsigned int readcofactorfile(const char* filename,const unsigned int nmar) {
  unsigned int cmarker = 0;
  unsigned int set_cofactors = 0;
  cofactors = newcvector(nmar);
  ifstream myfstream(filename, ios::in);
  while (!myfstream.eof()) {
    myfstream >> markerchr[x];
    if (markerchr[x]) set_cofactors++;
  }
  myfstream.close();
  //TODO get array back to main
  return set_cofactors;
}



void printoptionshelp(void) {
  printf ("Commandline switches:\n");
  printf ("-h      		This help.\n");
  printf ("-v      		Verbose (produce a lot of textoutput).\n");
  printf ("-p(INT) 		DebugLevel -d0,-d1.\n");
  printf ("-p(FILE_NAME)	Phenotypes file in plain textformat.\n");
  printf ("-g(FILE_NAME)	Genotypes file in plain textformat.\n");
  printf ("-m(FILE_NAME)	Marker and Chromosome descriptionfile in plain textformat.\n");
  printf ("-s(FILE_NAME)	Settings file in plain textformat.\n");
  printf ("-c(FILE_NAME)	Optional Cofactors file  to do backward elimination on in plain textformat.\n");
}

//Functions
void exitonerror(const char *msg) {
  fprintf(stderr, msg);
  printoptionshelp();
  exit(1);
}

bool checkfileexists(const char *filename) {
  ifstream myfile;
  bool exists;
  myfile.open(filename);
  exists = myfile.is_open();
  myfile.close();
  return exists;
}

//Main function
int main (unsigned int argc, char **argv) {
  //variables
  bool verboseflag = false;
  bool helpflag = false;
  int debuglevel = 0;
  char *phenofile = NULL;
  char *genofile = NULL;
  char *markerfile = NULL;
  char *coffile = NULL;
  char *settingsfile = NULL;
  unsigned int index;
  signed int c;

//Parsing of arguments
  while ((c = getopt (argc, argv, "vd:h:p:g:m:c:")) != -1)
    switch (c) {
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
    case 's':
      settings = optarg;
      break;
    case 'c':
      coffile = optarg;
      break;
    default:
      fprintf (stderr, "Unknown option character '%c'.\n", optopt);
    }
  if (helpflag) {
    printoptionshelp();
    return 0;
  } else {
    printf ("Options for MQM:\n");
//Verbose & debug
    printf ("verboseflag = %d, debuglevel = %d\n",verboseflag, debuglevel);
//Needed files
    if (!phenofile) exitonerror("Please supply a phenotypefile argument.\n");
    if (!checkfileexists(phenofile)) exitonerror("Phenotypefile not found on your filesystem.\n");
    printf ("Phenotypefile = %s\n",phenofile);
    if (!genofile)  exitonerror("Please supply a genofile argument.\n");
    if (!checkfileexists(genofile)) exitonerror("Genotypefile not found on your filesystem.\n");
    printf ("Genotypefile = %s\n",genofile);
    if (!markerfile) exitonerror("Please supply a markerfile argument.\n");
    if (!checkfileexists(genofile)) exitonerror("Markerfile not found on your filesystem.\n");
    printf ("Markerfile = %s\n",markerfile);
    if (!settingsfile) exitonerror("Please supply a settingsfile argument.\n");
    if (!checkfileexists(settingsfile)) exitonerror("settingsfile not found on your filesystem.\n");
    printf ("settingsfile = %s\n",settingsfile);
//Optional files
    if (!coffile) {
      if (!checkfileexists(coffile)) {
        printf("Cofactorfile not found on your filesystem.\n");
      } else {
        printf ("Cofactorfile = %s\n",coffile);
      }
    }
//Warn people for non-existing options
    for (index = optind; index < argc; index++) {
      printf ("Non-option argument %s\n", argv[index]);
    }
//Here we know what we need so we can start reading in files with the new loader functions
    mqmalgorithmsettings = loadmqmsetting(settingsfile);
  }
  return 0;
}