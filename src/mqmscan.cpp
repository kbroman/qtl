/**********************************************************************
 *
 * mqmscan.cpp
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

#include "mqm.h"
#include <Rmath.h>
#include <limits>

using namespace std;

inline int mqmmod(int a, int b) {
  return a%b;
}

/*
 * Helper function for truncate
 */
static double ftruncate(double n, double p = 3){
  int sign = 1;           // Assume positive sign
  if(n < 0) sign = -1;    // Test the assumption and change the sign if needed

  double val = fabs((pow(10,p)) * n);
  val = floor(val);
  val /= pow(10,p);
  return (double) sign * val;
}

/*
 * Truncate a floating point to 3 decimal numbers. This is used for output
 * functions, in particular for regression tests - so floating point problems on
 * different platforms are eliminated
 */
double ftruncate3(double n){
  return ftruncate(n,3);
}

double Lnormal(double residual, double variance) {
  //Now using R-library for Lnormal
  double result = dnorm(residual,0,sqrt(variance),0);
  debug_trace("Lnormal result:%f, residual: %f, variance %f\n",result,residual,variance);
  return result;
}

void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno) {
//reorganisation of doubles into a matrix
  int i;
  *Pheno = (double **)R_alloc(n_mar, sizeof(double *));
  (*Pheno)[0] = pheno;
  for (i=1; i< n_mar; i++)
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}


void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno) {
//reorganisation of integers into a matrix
  int i;
  *Pheno = (int **)R_alloc(n_mar, sizeof(int *));
  (*Pheno)[0] = pheno;
  for (i=1; i< n_mar; i++)
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}


/*
 * analyseF2 - analyse one F2/RIL/BC family
 * This is the main controller - called by mqmscan
 *
 * Returns logL
 */

double analyseF2(int Nind, int *nummark, cvector *cofactor, MQMMarkerMatrix marker,
               vector y, int Backwards, double **QTL,vector
               *mapdistance, int **Chromo, int Nrun, int RMLorML, double
               windowsize, double stepsize, double stepmin, double stepmax,
               double alfa, int em, int out_Naug, int **INDlist, char
               reestimate, MQMCrossType crosstype, bool dominance, int verbose) {
  if (verbose) Rprintf("INFO: Starting C-part of the MQM analysis\n");

  int  Naug, Nmark = (*nummark), run = 0;
  bool useREML = true, fitQTL = false;
  bool warned = false;

  ivector chr = newivector(Nmark); // The chr vector contains the chromosome number for every marker
  for(int i = 0; i < Nmark; i++){  // Rprintf("INFO: Receiving the chromosome matrix from R");
    chr[i] = Chromo[0][i];
  }
  if(RMLorML == 1) useREML=false;  // use ML instead

  // Create an array of marker positions - and calculate R[f] based on these locations
  cvector position = relative_marker_position(Nmark,chr);
  vector  r = recombination_frequencies(Nmark, position, (*mapdistance));

  //Rprintf("INFO: Initialize Frun and informationcontent to 0.0");
  const int Nsteps = (int)(chr[Nmark-1]*((stepmax-stepmin)/stepsize+1));
  matrix Frun = newmatrix(Nsteps,Nrun+1);
  vector informationcontent = newvector(Nsteps);
  for (int i = 0; i < (Nrun+1); i++) {
    for (int ii = 0; ii < Nsteps; ii++) {
      if(i==0) informationcontent[ii] = 0.0;
      Frun[ii][i]= 0.0;
    }
  }

  bool dropj = false;
  int jj=0;

  // Rprintf("any triple of non-segregating markers is considered to be the result of:\n");
  // Rprintf("identity-by-descent (IBD) instead of identity-by-state (IBS)\n");
  // Rprintf("no (segregating!) cofactors are fitted in such non-segregating IBD regions\n");
  for (int j=0; j < Nmark; j++) { // WRONG: (Nmark-1) Should fix the out of bound in mapdistance, it does fix, but created problems for the last marker
    dropj = false;
    if(j+1 < Nmark){  // Check if we can look ahead
      if(((*mapdistance)[j+1]-(*mapdistance)[j])==0.0){ dropj=true; }
    }
    if (!dropj) {
      marker[jj]          = marker[j];
      (*cofactor)[jj]     = (*cofactor)[j];
      (*mapdistance)[jj]  = (*mapdistance)[j];
      chr[jj]             = chr[j];
      r[jj]               = r[j];
      position[jj]        = position[j];
      jj++;
    } else{
      if (verbose) Rprintf("INFO: Marker %d at chr %d is dropped\n",j,chr[j]);
      if ((*cofactor)[j]==MCOF) {
        if (verbose) Rprintf("INFO: Cofactor at chr %d is dropped\n",chr[j]);
      }
    }
  }
  //if(verbose) Rprintf("INFO: Number of markers: %d -> %d\n",Nmark,jj);
  Nmark = jj;
  (*nummark) = jj;

  // Update the array of marker positions - and calculate R[f] based on these new locations
  position = relative_marker_position(Nmark,chr);

  r = recombination_frequencies(Nmark, position, (*mapdistance));

  debug_trace("After dropping of uninformative cofactors\n");

  ivector newind; // calculate Traits mean and variance
  vector newy;
  MQMMarkerMatrix newmarker;
  double ymean = 0.0, yvari = 0.0;
  //Rprintf("INFO: Number of individuals: %d Number Aug: %d",Nind,out_Naug);
  int cur = -1;
  for (int i=0; i < Nind; i++){
    if(INDlist[0][i] != cur){
      ymean += y[i];
      cur = INDlist[0][i];
    }
  }
  ymean/= out_Naug;

  for (int i=0; i < Nind; i++){
    if(INDlist[0][i] != cur){
      yvari += pow(y[i]-ymean, 2);
      cur = INDlist[0][i];
    }
  }
  yvari /= (out_Naug-1);

  Naug      = Nind;                             // Fix for not doing dataaugmentation, we just copy the current as the augmented and set Naug to Nind
  Nind      = out_Naug;
  newind    = newivector(Naug);
  newy      = newvector(Naug);
  newmarker = newMQMMarkerMatrix(Nmark,Naug);
  for (int i=0; i<Naug; i++) {
    newy[i]= y[i];
    newind[i]= INDlist[0][i];
    for (int j=0; j<Nmark; j++) {
      newmarker[j][i]= marker[j][i];
    }
  }
  // End fix

  vector newweight = newvector(Naug);

  double max = rmixture(newmarker, newweight, r, position, newind,Nind, Naug, Nmark, mapdistance,reestimate,crosstype,verbose);   //Re-estimation of mapdistances if reestimate=TRUE

  if(max > stepmax){ fatal("ERROR: Re-estimation of the map put markers at: %f Cm, run the algorithm with a step.max larger than %f Cm", max, max); }

  //Check if everything still is correct positions and R[f]
  position = relative_marker_position(Nmark,chr);

  r = recombination_frequencies(Nmark, position, (*mapdistance));

  /* eliminate individuals with missing trait values */
  //We can skip this part iirc because R throws out missing phenotypes beforehand
  int oldNind = Nind;
  for (int i=0; i<oldNind; i++) {
    Nind -= ((y[i]==TRAITUNKNOWN) ? 1 : 0);
  }

  int oldNaug = Naug;
  for (int i=0; i<oldNaug; i++) {
    Naug -= ((newy[i]==TRAITUNKNOWN) ? 1 : 0);
  }

  marker        = newMQMMarkerMatrix(Nmark+1,Naug);
  y             = newvector(Naug);
  ivector ind   = newivector(Naug);
  vector weight = newvector(Naug);
  int newi = 0;
  for (int i=0; i < oldNaug; i++)
    if (newy[i]!=TRAITUNKNOWN) {
      y[newi]= newy[i];
      ind[newi]= newind[i];
      weight[newi]= newweight[i];
      for (int j=0; j<Nmark; j++) marker[j][newi]= newmarker[j][i];
      newi++;
    }
  int diff;
  for (int i=0; i < (Naug-1); i++) {
    diff = ind[i+1]-ind[i];
    if (diff>1) {
      for (int ii=i+1; ii<Naug; ii++){ ind[ii]=ind[ii]-diff+1; }
    }
  }
  //END throwing out missing phenotypes

  double variance=-1.0;
  cvector selcofactor = newcvector(Nmark); /* selected cofactors */
  int dimx   = designmatrixdimensions((*cofactor),Nmark,dominance);
  double F1  = inverseF(1,Nind-dimx,alfa,verbose);
  double F2  = inverseF(2,Nind-dimx,alfa,verbose);
  if (verbose) {
    Rprintf("INFO: dimX: %d, nInd: %d\n",dimx,Nind);
    Rprintf("INFO: F(Threshold, Degrees of freedom 1, Degrees of freedom 2) = Alfa\n");
    Rprintf("INFO: F(%.3f, 1, %d) = %f\n",ftruncate3(F1),(Nind-dimx),alfa);
    Rprintf("INFO: F(%.3f, 2, %d) = %f\n",ftruncate3(F2),(Nind-dimx),alfa);
  }
  F2 = 2.0* F2; // 9-6-1998 using threshold x*F(x,df,alfa)

  weight[0]= -1.0;
  double logL = QTLmixture(marker,(*cofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,useREML,fitQTL,dominance,crosstype, &warned, verbose);
  if(verbose){
    if (!R_finite(logL)) {
      Rprintf("WARNING: Log-likelihood of full model = INFINITE\n");
    }else{
      if (R_IsNaN(logL)) {
        Rprintf("WARNING: Log-likelihood of full model = NOT A NUMBER (NAN)\n");
      }else{
        Rprintf("INFO: Log-likelihood of full model = %.3f\n",ftruncate3(logL));
      }
    }
    Rprintf("INFO: Residual variance = %.3f\n",ftruncate3(variance));
    Rprintf("INFO: Trait mean= %.3f; Trait variation = %.3f\n",ftruncate3(ymean),ftruncate3(yvari));
  }
  if (R_finite(logL) && !R_IsNaN(logL)) {
    if(Backwards==1){    // use only selected cofactors
      logL = backward(Nind, Nmark, (*cofactor), marker, y, weight, ind, Naug, logL,variance, F1, F2, &selcofactor, r,
                      position, &informationcontent, mapdistance,&Frun,run,useREML,fitQTL,dominance, em, windowsize,
                      stepsize, stepmin, stepmax,crosstype,verbose);
    }else{ // use all cofactors
      logL = mapQTL(Nind, Nmark, (*cofactor), (*cofactor), marker, position,(*mapdistance), y, r, ind, Naug, variance,
                    'n', &informationcontent,&Frun,run,useREML,fitQTL,dominance, em, windowsize, stepsize, stepmin,
                    stepmax,crosstype,verbose); // printout=='n'
    }
  }
  // Write output and/or send it back to R
  // Cofactors that made it to the final model
  for (int j=0; j<Nmark; j++) {
    if (selcofactor[j]==MCOF) {
      (*cofactor)[j]=MCOF;
    }else{
      (*cofactor)[j]=MNOCOF;
    }
  }

  if (verbose) Rprintf("INFO: Number of output datapoints: %d\n", Nsteps);  // QTL likelihood for each location
  for (int ii=0; ii<Nsteps; ii++) {
    //Convert LR to LOD before sending back
    QTL[0][ii] = Frun[ii][0] / 4.60517;
    QTL[0][Nsteps+ii] = informationcontent[ii];
  }
  return logL;
}

/**********************************************************************
 *
 * mqmscan
 *
 *
 **********************************************************************/

void mqmscan(int Nind, int Nmark,int Npheno,int **Geno,int **Chromo, double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,
             int Emiter, double Windowsize,double Steps, double Stepmi,double Stepma,int NRUN,int out_Naug,int **INDlist, double **QTL, int re_estimate,
             RqtlCrossType rqtlcrosstype,int domi,int verbose){
  int cof_cnt=0;
  MQMMarkerMatrix markers = newMQMMarkerMatrix(Nmark+1,Nind);
  cvector cofactor        = newcvector(Nmark);
  vector mapdistance      = newvector(Nmark);

  MQMCrossType crosstype = determine_MQMCross(Nmark,Nind,(const int **)Geno,rqtlcrosstype);

  change_coding(&Nmark, &Nind, Geno, markers, crosstype); // Change all the markers from R/qtl format to MQM internal

  for (int i=0; i< Nmark; i++) {
    mapdistance[i] = POSITIONUNKNOWN;  // Mapdistances
    mapdistance[i] = Dist[0][i];
    cofactor[i]    = MNOCOF;           // Cofactors
    if (Cofactors[0][i] == 1) {
      cofactor[i] = MCOF;              // Set cofactor
      cof_cnt++;
    }
    if (Cofactors[0][i] == 2) {
      cofactor[i] = MSEX;
      cof_cnt++;
    }
    if (cof_cnt+10 > Nind){ fatal("Setting %d cofactors would leave less than 10 degrees of freedom.\n", cof_cnt); }
  }

  char reestimate = 'y';
  if(re_estimate == 0) reestimate = 'n';

  if (crosstype != CF2) {  // Determine what kind of cross we have
    if (verbose==1) Rprintf("INFO: Dominance setting ignored (setting dominance to 0)\n"); // Update dominance accordingly
    domi = 0;
  }

  bool dominance=false;
  if(domi != 0){ dominance=true; }

  //WE HAVE EVERYTHING START WITH MAIN SCANNING FUNCTION
  analyseF2(Nind, &Nmark, &cofactor, (MQMMarkerMatrix)markers, Pheno[(Npheno-1)], Backwards, QTL, &mapdistance, Chromo, NRUN, RMLorML, Windowsize,
            Steps, Stepmi, Stepma, Alfa, Emiter, out_Naug, INDlist, reestimate, crosstype, dominance, verbose);

  if (re_estimate) {
    if (verbose==1) Rprintf("INFO: Sending back the re-estimated map used during the MQM analysis\n");
    for (int i=0; i< Nmark; i++) {
      Dist[0][i] = mapdistance[i];
    }
  }
  if (Backwards) {
    if (verbose==1) Rprintf("INFO: Sending back the model\n");
    for (int i=0; i< Nmark; i++) { Cofactors[0][i] = cofactor[i]; }
  }

  if(verbose) Rprintf("INFO: All done in C returning to R\n");
  #ifndef STANDALONE
    R_CheckUserInterrupt(); /* check for ^C */
    R_FlushConsole();
  #endif
  return;
}  /* end of function mqmscan */

/**********************************************************************
 *
 * R_mqmscan
 *
 **********************************************************************/

void R_mqmscan(int *Nind,int *Nmark,int *Npheno,
               int *geno,int *chromo, double *dist, double *pheno,
               int *cofactors, int *backwards, int *RMLorML,double *alfa,int *emiter,
               double *windowsize,double *steps,double *stepmi,double *stepma,
               int *nRun, int *out_Naug, int *indlist, double *qtl, int *reestimate, int *crosstype, int *domi, int *verbose) {
  int **Geno;
  int **Chromo;
  double **Dist;
  double **Pheno;
  double **QTL;
  int **Cofactors;
  int **INDlist;

  // Reorganise the pointers into arrays, singletons are just cast into the function
  reorg_geno(*Nind,*Nmark,geno,&Geno);
  reorg_int(*Nmark,1,chromo,&Chromo);
  reorg_pheno(*Nmark,1,dist,&Dist);
  // Here we have  the assumption that step.min is negative this needs to be split in 2
  reorg_pheno((int)(2*(*chromo) * (((*stepma)-(*stepmi))/ (*steps))),1,qtl,&QTL);
  reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
  reorg_int(*Nmark,1,cofactors,&Cofactors);
  reorg_int(*out_Naug,1,indlist,&INDlist);

  mqmscan(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,*out_Naug,INDlist,QTL, *reestimate,(RqtlCrossType)*crosstype,*domi,*verbose);
}
