/**********************************************************************
 *
 * MQMmain.cpp
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

#include <fstream>
#include <iostream>


#include "mqm.h"

using namespace std;

int count_lines(const char *file) {
  //NUM: number of elements on 1 line
  int cnt=0;
  char line[100];
  ifstream file_stream(file, ios::in);
  while (!file_stream.eof()) {
    file_stream >> line;
    cnt++;
  }
  file_stream.close();
  return cnt;
}

#ifdef STANDALONE

int main(int argc,char *argv[]) {
  Rprintf("MQM standalone version\n");
  int phenotype=0;
  int verbose=0;
  for (int i=1; i<argc; i++) {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') Rprintf("dash needed at argument");
    char c = toupper(argv[i][1]);
    if (c == 'V') {
      verbose =1;
      continue;
    }
    // -argum=value
    if (argv[i][2]!='=') {
      Rprintf("equal symbol needed at argument");
    }

    switch (c) {
    case 'T':
      phenotype = atoi(&argv[i][3]);
      break;
    default:
      Rprintf("Unknown parameter");
    }
  }
  const char *genofile = "geno.dat";
  const char *phenofile = "pheno.dat";
  const char *mposfile = "markerpos.txt";
  const char *chrfile = "chrid.dat";
  const char *setfile = "settings.dat";
  double **QTL;
  ivector f1genotype;
  ivector chr;
  cvector cofactor;
  vector mapdistance;
  vector pos;
  matrix pheno_value;
  cmatrix markers;
  ivector INDlist;
  int stepmin = 0;
  int stepmax = 220;
  int stepsize = 5;

  int cnt=0;
  int cInd=0; //Current individual
  int nInd=0;
  int nMark=0;
  int backwards=0;
  int nPheno=0;
  int windowsize=0;
  cnt = 0;
  char *name;
  int maxIter;
  double alpha;
  nMark=count_lines(chrfile);
  char peek_c;

  if (verbose) {
    Rprintf("INFO: Loading settings from file\n");
  }
  ifstream setstr(setfile, ios::in);
  setstr >> nInd;
  if (verbose) {
    Rprintf("nInd: %d\n",nInd);
  }
  setstr >> nPheno;
  if (verbose) {
    Rprintf("nPheno: %d\n",nPheno);
  }
  setstr >> stepmin;
  if (verbose) {
    Rprintf("SMin: %d\n",stepmin);
  }
  setstr >> stepmax;
  if (verbose) {
    Rprintf("SMax: %d\n",stepmax);
  }
  setstr >> stepsize;
  if (verbose) {
    Rprintf("SSiz: %d\n",stepsize);
  }
  setstr >> windowsize;
  if (verbose) {
    Rprintf("WSiz: %d\n",windowsize);
  }
  setstr >> alpha;
  if (verbose) {
    Rprintf("A: %f\n",alpha);
  }
  setstr >> maxIter;
  if (verbose) {
    Rprintf("Miter: %d\n",maxIter);
  }
  f1genotype = newivector(nMark);
  cofactor= newcvector(nMark);
  mapdistance= newvector(nMark);
  markers= newcmatrix(nMark,nInd);
  pheno_value = newmatrix(nPheno,nInd);
  chr = newivector(nMark);
  INDlist= newivector(nInd);
  pos = newvector(nMark);
  int sum = 0;
  for (int i=0; i< nMark; i++) {
    setstr >> cofactor[i];
    if (cofactor[i] == '1') {
      sum++;
    }
  }

  if (sum > 0) {
    if (verbose) {
      Rprintf("# Starting backward elimination of %d cofactors\n",sum);
    }
    backwards = 1;
  } else {
    backwards = 0;
  }
  setstr.close();
  if (verbose) {
    Rprintf("# of individuals: %d\n",nInd);
  }
  if (verbose) {
    Rprintf("# of markers: %d\n",nMark);
  }
  cnt=0;
  cInd = 0;
  ifstream geno(genofile, ios::in);
  while (!geno.eof()) {
    if (cnt < nMark) {
      geno >> markers[cnt][cInd];
      cnt++;
    } else {
      cnt = 0;
      cInd++;
    }
  }
  geno.close();
  if (verbose) {
    Rprintf("Genotypes done %d %d\n",cInd,cnt);
  }
  cnt = 0;
  cInd = 0;
  ifstream pheno(phenofile, ios::in);
  while (!pheno.eof()) {
    if (cnt < nPheno) {
      pheno >> pheno_value[cnt][cInd];
      //Rprintf("%d,%d\n",cnt,cInd);
      cnt++;
    } else {
      cnt = 0;
      cInd++;
    }
  }
  pheno.close();
  if (verbose) {
    Rprintf("Phenotype done %d %d\n",cInd,cnt);
  }
  cnt = 0;
  ifstream mpos(mposfile, ios::in);
  while (!mpos.eof()) {
    peek_c=mpos.peek();
    if (peek_c=='\t' || peek_c == ' ') {
      mpos >> pos[cnt];
      //  Rprintf("%f\n",pos[cnt]);
      cnt++;
    } else {
      mpos >> peek_c;
    }
  }
  mpos.close();

  if (verbose) {
    Rprintf("Positions done %d\n",cnt);
  }
  cnt = 0;
  ifstream chrstr(chrfile, ios::in);
  int max_chr = 0;
  while (!chrstr.eof()) {
    chrstr >> chr[cnt];
    if (chr[cnt] > max_chr) {
      max_chr = chr[cnt];
    }
    cnt++;
  }
  chrstr.close();
  if (verbose) {
    Rprintf("Chromosomes done %d -> # %d Chromosomes\n",cnt,max_chr);
  }
  int something = 3*max_chr*(((stepmax)-(stepmin))/ (stepsize));
  QTL = newmatrix(1,something);

  for (int i=0; i< nMark; i++) {
    cofactor[i] = '0';
    f1genotype[i] = 12;
    mapdistance[i]=999.0;
    mapdistance[i]=pos[i];
  }
  for (int i=0; i< nInd; i++) {
    INDlist[i] = i;
  }
  char estmap = 'n';
  //reorg_pheno(2*(*chromo) * (((*stepma)-(*stepmi))/ (*steps)),1,qtl,&QTL);

  //Rprintf("INFO: Cofactors %d\n",sum);

  Rprintf("Starting phenotype: %d\n",phenotype);
  analyseF2(nInd, nMark, &cofactor, markers, pheno_value[phenotype], f1genotype, backwards,QTL, &mapdistance,&chr,0,0,windowsize,stepsize,stepmin,stepmax,alpha,maxIter,nInd,&INDlist,estmap,'F',0,verbose);

  // Output marker info
  for (int m=0; m<nMark; m++) {
    Rprintf("%5d%3d%9.3f\n",m,chr[m],mapdistance[m]);
  }
  // Output (augmented) QTL info
  for (int q=0; q<something; q++) {
    Rprintf("%5d%10.5f\n",q,QTL[0][q]);
  }
  freevector((void *)f1genotype);
  freevector((void *)cofactor);
  freevector((void *)mapdistance);
  freematrix((void **)markers,nMark);
  freematrix((void **)pheno_value,nPheno);
  freevector((void *)chr);
  freevector((void *)INDlist);
  freevector((void *)pos);
  freematrix((void **)QTL,1);
  return 0;
}

#else
#error "Is this a STANDALONE version? STANDALONE should be defined in the build system!"
#endif
