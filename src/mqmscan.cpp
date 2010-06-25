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

inline int mqmmod(int a, int b) {
  return a%b;
}

/*
 * Helper function for truncate
 */
static double ftruncate(double n, double p = 3){
  int sign = 0;
  if(n >= 0){
      sign = 1;
  }else{
      sign = -1;
  }
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
               vector y, ivector f1genotype, int Backwards, double **QTL,vector
               *mapdistance, int **Chromo, int Nrun, int RMLorML, double
               windowsize, double stepsize, double stepmin, double stepmax,
               double alfa, int em, int out_Naug, int **INDlist, char
               reestimate, MQMCrossType crosstype, bool dominance, int verbose) {
  if (verbose) {
    Rprintf("INFO: Starting C-part of the MQM analysis\n");
  }
  int Naug;
  int Nmark = (*nummark);
  int run=0;
  cvector position;
  vector informationcontent;
  ivector chr;
  matrix Frun;
  vector r;
  r= newvector(Nmark);
  position= newcvector(Nmark);
  bool useREML=true;
  bool fitQTL=false;

  // The chr vector contains the chromosome number for every marker
  chr= newivector(Nmark);
  //info("Receiving the chromosome matrix from R");
  for (int i=0; i< Nmark; i++) {
    chr[i] = Chromo[0][i];
  }
  if (RMLorML == 1) useREML=false;  // use ML instead

  // Create an array of marker positions - and calculate R[f] based on these locations
  position = relative_marker_position(Nmark,chr);
  r = recombination_frequencies(Nmark, position, (*mapdistance));

  //info("Initialize Frun and informationcontent to 0.0");
  const int Nsteps = (int)(chr[Nmark-1]*((stepmax-stepmin)/stepsize+1));
  Frun= newmatrix(Nsteps,Nrun+1);
  informationcontent= newvector(Nsteps);
  for (int i=0; i<Nrun+1; i++) {
    for (int ii=0; ii<Nsteps; ii++) {
      Frun[ii][i]= 0.0;
    }
  }
  for (int ii=0; ii<Nsteps; ii++) {
    informationcontent[ii]= 0.0;
  }

  bool dropj=true;
  int jj=0;
  //Rprintf("any triple of non-segregating markers is considered to be the result of:\n");
  //Rprintf("identity-by-descent (IBD) instead of identity-by-state (IBS)\n");
  //  Rprintf("no (segregating!) cofactors are fitted in such non-segregating IBD regions\n");
  for (int j=0; j<Nmark; j++) {
    if (mqmmod(f1genotype[j],11)!=0) {
      dropj=false;
      if(((*mapdistance)[j+1]-(*mapdistance)[j])==0.0) dropj=true;
    } else if ((*cofactor)[j]==MNOCOF) {
      dropj=true;
    } else if (position[j]==MLEFT) {
      // (cofactor[j]!=MNOCOF) cofactor at non-segregating marker test whether next segregating marker is nearby (<20cM)
      dropj='y';
      if ((((*mapdistance)[j+1]-(*mapdistance)[j])<windowsize)) dropj=false;
      else if (position[j+1]!=MRIGHT)
        if ((((*mapdistance)[j+2]-(*mapdistance)[j])<windowsize)) dropj=false;
    } else if (position[j]==MMIDDLE) {
      dropj=true;
      if ((((*mapdistance)[j]-(*mapdistance)[j-1])<windowsize)) dropj=false;
      else if ((((*mapdistance)[j+1]-(*mapdistance)[j])<windowsize)) dropj=false;
    } else if (position[j]==MRIGHT) {
      dropj=true;
      if ((((*mapdistance)[j]-(*mapdistance)[j-1])<windowsize)) dropj=false;
      else if (position[j-1]!=MLEFT)
        if ((((*mapdistance)[j]-(*mapdistance)[j-2])<windowsize)) dropj=false;
    }
    if (!dropj) {
      marker[jj]= marker[j];
      (*cofactor)[jj]= (*cofactor)[j];
      (*mapdistance)[jj]= (*mapdistance)[j];
      chr[jj]= chr[j];
      r[jj]= r[j];
      position[jj]= position[j];
      jj++;
    } else{
      if (verbose) info("Marker %d at chr %d is dropped",j,chr[j]);
      if ((*cofactor)[j]==MCOF) {
        if (verbose) info("Cofactor at chr %d is dropped",chr[j]);
      }
    }
  }
  debug_trace("Num markers: %d -> %d\n",Nmark,jj);
  Nmark= jj;
  (*nummark) = jj;
  position = relative_marker_position(Nmark,chr);
  r = recombination_frequencies(Nmark, position, (*mapdistance));

  debug_trace("After dropping of uninformative cofactors\n");
  //calculate Traits mean and variance
  ivector newind;
  vector newy;
  MQMMarkerMatrix newmarker;
  double ymean=0.0, yvari=0.0;
  //info("Number of individuals: %d Number Aug: %d",Nind,out_Naug);
  int cur = -1;
  for (int i=0; i<Nind; i++){
    if(INDlist[0][i] != cur){
      ymean += y[i];
      cur = INDlist[0][i];
    }
  }
  ymean/= out_Naug;
  
  for (int i=0; i<Nind; i++){
    if(INDlist[0][i] != cur){
      yvari += pow(y[i]-ymean,2);
      cur = INDlist[0][i];
    }
  }  
  yvari/= (out_Naug-1);
  //Fix for not doing dataaugmentation, we just copy the current as the augmented and set Naug to Nind
  Naug=Nind;
  Nind=out_Naug;
  newind= newivector(Naug);
  newy= newvector(Naug);
  newmarker= newMQMMarkerMatrix(Nmark,Naug);
  for (int i=0; i<Naug; i++) {
    newy[i]= y[i];
    newind[i]= INDlist[0][i];
    for (int j=0; j<Nmark; j++) {
      newmarker[j][i]= marker[j][i];
    }
  }
  //End fix
  vector newweight;
  newweight= newvector(Naug);
  //Re-estimation of mapdistances if reestimate=TRUE
  double max;
  max = rmixture(newmarker, newweight, r, position, newind,Nind, Naug, Nmark, mapdistance,reestimate,crosstype,verbose);
  if (max > stepmax) {
    info("ERROR: Reestimation of the map put markers at: %f Cm",max);
    info("ERROR: Rerun the algorithm with a step.max larger than %f Cm",max);
    return std::numeric_limits<double>::quiet_NaN();
  } else {
   // if (verbose) info("Reestimation of the map finished. MAX Cm: %f Cm",max);
  }

  //Check if everything still is correct positions and R[f]
  position = relative_marker_position(Nmark,chr);
  r = recombination_frequencies(Nmark, position, (*mapdistance));
  
  /* eliminate individuals with missing trait values */
  //We can skip this part iirc because R throws out missing phenotypes beforehand
  int oldNind=Nind;
  for (int i=0; i<oldNind; i++) {
    Nind-= ((y[i]==TRAITUNKNOWN) ? 1 : 0);
  }

  int oldNaug=Naug;
  for (int i=0; i<oldNaug; i++) {
    Naug-= ((newy[i]==TRAITUNKNOWN) ? 1 : 0);
  }

  vector weight;
  ivector ind;
  marker= newMQMMarkerMatrix(Nmark+1,Naug);
  y= newvector(Naug);
  ind= newivector(Naug);
  weight= newvector(Naug);
  int newi=0;
  for (int i=0; i<oldNaug; i++)
    if (newy[i]!=TRAITUNKNOWN) {
      y[newi]= newy[i];
      ind[newi]= newind[i];
      weight[newi]= newweight[i];
      for (int j=0; j<Nmark; j++) marker[j][newi]= newmarker[j][i];
      newi++;
    }
  int diff;
  for (int i=0; i<Naug-1; i++) {
    diff=ind[i+1]-ind[i];
    if  (diff>1)
      for (int ii=i+1; ii<Naug; ii++) ind[ii]=ind[ii]-diff+1;
  }
  delMQMMarkerMatrix(newmarker,Nmark);
  Free(newy);
  Free(newind);
  Free(newweight);
  //END throwing out missing phenotypes

  double variance=-1.0;
  double logL;
  cvector selcofactor;
  selcofactor= newcvector(Nmark); /* selected cofactors */

  int dimx=designmatrixdimensions((*cofactor),Nmark,dominance);
  double F1, F2;

  F1= inverseF(1,Nind-dimx,alfa,verbose);
  F2= inverseF(2,Nind-dimx,alfa,verbose);
  if (verbose) {
    info("dimX:%d nInd:%d",dimx,Nind);
    info("F(Threshold,Degrees of freedom 1,Degrees of freedom 2)=Alfa");
    info("F(%.3f,1,%d)=%f",ftruncate3(F1),(Nind-dimx),alfa);
    info("F(%.3f,2,%d)=%f",ftruncate3(F2),(Nind-dimx),alfa);
  }
  F2= 2.0* F2; // 9-6-1998 using threshold x*F(x,df,alfa)

  weight[0]= -1.0;
  logL = QTLmixture(marker,(*cofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,useREML,fitQTL,dominance,crosstype,verbose);
  if(verbose)
  {
    if (!R_finite(logL)) {
      info("Log-likelihood of full model= INFINITE");
    } else 
      if (R_IsNaN(logL)) {
        info("Log-likelihood of full model= NOT A NUMBER (NAN)");
      }
      else {
        info("Log-likelihood of full model= %.3f",ftruncate3(logL));
      }
    info("Residual variance= %.3f",ftruncate3(variance));
    info("Trait mean= %.3f; Trait variation= %.3f",ftruncate3(ymean),ftruncate3(yvari));
  }
  if (R_finite(logL) && !R_IsNaN(logL)) {
    if (Backwards==1)    // use only selected cofactors
      logL = backward(Nind, Nmark, (*cofactor), marker, y, weight, ind, Naug, logL,variance, F1, F2, &selcofactor, r, position, &informationcontent, mapdistance,&Frun,run,useREML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype,verbose);
    if (Backwards==0) // use all cofactors
      logL = mapQTL(Nind, Nmark, (*cofactor), (*cofactor), marker, position,(*mapdistance), y, r, ind, Naug, variance, 'n', &informationcontent,&Frun,run,useREML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype,verbose); // printout=='n'
  }
  // Write output and/or send it back to R
  // Cofactors that made it to the final model
  for (int j=0; j<Nmark; j++) {
    if (selcofactor[j]==MCOF) {
      (*cofactor)[j]=MCOF;
    } else {
      (*cofactor)[j]=MNOCOF;
    }
  }
  // QTL likelihood for each location
  if (verbose) info("Number of output datapoints: %d",Nsteps);
  for (int ii=0; ii<Nsteps; ii++) {
    //Convert LR to LOD before sending back
    QTL[0][ii] = Frun[ii][0] / 4.60517;
    QTL[0][Nsteps+ii] = informationcontent[ii];
  }
  //Free used memory
  Free(position);
  Free(weight);
  Free(ind);
  Free(r);
  Free(informationcontent);
  freematrix((void **)Frun,Nsteps);
  delMQMMarkerMatrix(marker,Nmark+1);
  Free(y);
  Free(chr);
  Free(selcofactor);
  //info("Analysis of data finished");
  return logL;
}



/**********************************************************************
 *
 * mqmscan
 *
 *
 **********************************************************************/

void mqmscan(int Nind, int Nmark,int Npheno,int **Geno,int **Chromo,
             double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,int Emiter,
             double Windowsize,double Steps,
             double Stepmi,double Stepma,int NRUN,int out_Naug,int **INDlist, double **QTL, int re_estimate,RqtlCrossType rqtlcrosstype,int domi,int verbose) {

  ivector f1genotype;
  MQMMarkerMatrix markers;
  cvector cofactor;
  vector mapdistance;

  markers= newMQMMarkerMatrix(Nmark+1,Nind);
  f1genotype = newivector(Nmark);
  cofactor= newcvector(Nmark);
  mapdistance= newvector(Nmark);

  int cof_cnt=0;

  MQMCrossType crosstype = determine_MQMCross(Nmark,Nind,(const int **)Geno,rqtlcrosstype);
  //Change all the markers from R/qtl format to MQM internal
  change_coding(&Nmark,&Nind,Geno,markers,crosstype);

  for (int i=0; i< Nmark; i++) {
    f1genotype[i] = 12;               //The parental strain for all markers, this was used to asses marker information
    mapdistance[i]=POSITIONUNKNOWN;   //Mapdistances
    mapdistance[i]=Dist[0][i];  
    cofactor[i] = MNOCOF;             //Cofactors
    if (Cofactors[0][i] == 1) {
      cofactor[i] = MCOF;             //Set cofactor
      cof_cnt++;
    }
    if (Cofactors[0][i] == 2) {
      cofactor[i] = MSEX;
      cof_cnt++;
    }
    if (cof_cnt+10 > Nind) {
      Rprintf("ERROR: Setting %d cofactors would leave less than 10 degrees of freedom.\n",cof_cnt);
      return;
    }
  }

  char reestimate = 'y';
  if (re_estimate == 0) {
    reestimate = 'n';
  }
  //determine what kind of cross we have
  //set dominance accordingly
  if (crosstype != CF2) {
    if (verbose==1) {
      Rprintf("INFO: Dominance setting ignored (dominance=0)\n");
    }
    domi = 0;
  } else {
    domi = domi;
  }

  bool dominance=false;
  if (domi != 0) {
    dominance=true;
  }

  //WE HAVE EVERYTHING START WITH MAIN SCANNING FUNCTION
  analyseF2(Nind, &Nmark, &cofactor, (MQMMarkerMatrix)markers, Pheno[(Npheno-1)], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter,out_Naug,INDlist,reestimate,crosstype,dominance,verbose);

  if (re_estimate) {
    if (verbose==1) {
      Rprintf("INFO: Sending back the reestimated map used during analysis\n");
    }
    for (int i=0; i< Nmark; i++) {
      Dist[0][i]=mapdistance[i];
    }
  }
  if (Backwards) {
    if (verbose==1) {
      Rprintf("INFO: Sending back the model\n");
    }
    for (int i=0; i< Nmark; i++) {
      Cofactors[0][i]=cofactor[i];
    }
  }
  //Rprintf("Starting Cleanup\n");
  //[Danny]: The markers are already DELETED in analyseF2, So this should have always failed... 
  //but in some way it doesn't. R and Memory managment, perhaps its best to use the c++ stdlibs features
  //like vector<int> etc, (would mean again mayor recoding)
  //delMQMMarkerMatrix(markers,Nmark+1);
  Free(f1genotype);
  Free(cofactor);
  Free(mapdistance);
  if (verbose==1) {
    Rprintf("INFO: All done in C returning to R\n");
  }
#ifndef STANDALONE
  R_CheckUserInterrupt(); /* check for ^C */
  // R_ProcessEvents(); /*  Try not to crash windows etc*/
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

