/**********************************************************************
 *
 * mqmscan.cpp
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
 * Contains: R_mqmscan, mqmscan
 *
 **********************************************************************/
#include "mqm.h"
#include <Rmath.h>

inline int mqmmod(int a, int b) {
  return a%b;
}

double Lnormal(double residual, double variance) {
  //double Likelihood,likelyhood;
  //Likelihood=exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
  //if(fabs(Likelihood-likelyhood)>0.05){
  //Rprintf("ERROR: Lnormal error\n");
  //}
  //return likelyhood;
  return dnorm(residual,0,sqrt(variance),0);
}


void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno) {
  int i;

  *Pheno = (double **)R_alloc(n_mar, sizeof(double *));

  (*Pheno)[0] = pheno;
  for (i=1; i< n_mar; i++)
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}


void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno) {
  int i;

  *Pheno = (int **)R_alloc(n_mar, sizeof(int *));

  (*Pheno)[0] = pheno;
  for (i=1; i< n_mar; i++)
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}


/*
 * analyseF2 - analyse one F2/RIL/BC family
 *
 * This is the main controller - called by mqmscan
 *
 */

void analyseF2(int Nind, int Nmark, cvector *cofactor, cmatrix marker,
               vector y, ivector f1genotype, int Backwards, double **QTL,vector
               *mapdistance, int **Chromo, int Nrun, int RMLorML, double
               windowsize, double stepsize, double stepmin, double stepmax,
               double alfa, int em, int out_Naug, int **INDlist, char
               reestimate, MQMCrossType crosstype, char dominance, int verbose) {
  if (verbose) {
    Rprintf("INFO: Starting C-part of the MQM analysis\n");
  }
  int Naug;
  int run=0;
  cvector position;
  vector informationcontent;
  //char dominance='n';
  //char perm_simu=MH;
  ivector chr;
  matrix Frun;
  vector r;
  r= newvector(Nmark);
  position= newcvector(Nmark);
  bool useREML=true;
  char fitQTL='n';

  // The chr vector contains the chromosome number for every marker
  chr= newivector(Nmark);
  info("Receiving the chromosome matrix from R");
  for (int i=0; i< Nmark; i++) {
    chr[i] = Chromo[0][i];
  }
  if (RMLorML == 1) useREML=false;  // use ML instead

  // Create an array of marker positions - each marker is one of LMRU (left,
  // middle, right, unknown/single)
  info("Calculating relative genomepositions of the markers");
  for (int j=0; j<Nmark; j++) {
    if (j==0) {
      // first marker is MLEFT; if single marker MUNLINKED
      if (chr[j]==chr[j+1]) position[j]=MLEFT;
      else position[j]=MUNLINKED;
    } else if (j==(Nmark-1)) {
      // Last marker is MRIGHT; if single marker MUNLINKED
      if (chr[j]==chr[j-1]) position[j]=MRIGHT;
      else position[j]=MUNLINKED;
    } else if (chr[j]==chr[j-1]) { // marker on the left is on the same chromosome
      //  MMIDDLE the marker to the right is on the same chromosome, otherwise MRIGHT
      if (chr[j]==chr[j+1]) position[j]=MMIDDLE;
      else position[j]=MRIGHT;
    } else { // marker on the left is not on the same chromosome
      // MLEFT if marker to the right is on the same chromosome, otherwise MUNLINKED
      if (chr[j]==chr[j+1]) position[j]=MLEFT;
      else position[j]=MUNLINKED;
    }
  }

  info("Estimating recombinant frequencies");
  for (int j=0; j<Nmark; j++) {
    r[j]= 999.0;
    if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
      r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
    }
    //Rprintf("R[j] value: %f\n",r[j]);
  }
  info("Initialize Frun and informationcontent to 0.0");
  const int Nsteps = chr[Nmark-1]*((stepmax-stepmin)/stepsize+1);
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

  char dropj='y';
  int jj=0;
  //Rprintf("any triple of non-segregating markers is considered to be the result of:\n");
  //Rprintf("identity-by-descent (IBD) instead of identity-by-state (IBS)\n");
  //  Rprintf("no (segregating!) cofactors are fitted in such non-segregating IBD regions\n");
  for (int j=0; j<Nmark; j++) {
    if (mqmmod(f1genotype[j],11)!=0) {
      dropj='n';
    } else if ((*cofactor)[j]==MAA) {
      dropj='y';
    } else if (position[j]==MLEFT) {
      // (cofactor[j]!=MAA) cofactor at non-segregating marker
      // test whether next segregating marker is nearby (<20cM)
      dropj='y';
      if ((((*mapdistance)[j+1]-(*mapdistance)[j])<20)&&(mqmmod(f1genotype[j+1],11)!=0)) dropj='n';
      else if (position[j+1]!=MRIGHT)
        if ((((*mapdistance)[j+2]-(*mapdistance)[j])<20)&&(mqmmod(f1genotype[j+2],11)!=0)) dropj='n';
    } else if (position[j]==MMIDDLE) {
      dropj='y';
      if ((((*mapdistance)[j]-(*mapdistance)[j-1])<20)&&(mqmmod(f1genotype[j-1],11)!=0)) dropj='n';
      else if ((((*mapdistance)[j+1]-(*mapdistance)[j])<20)&&(mqmmod(f1genotype[j+1],11)!=0)) dropj='n';
    } else if (position[j]==MRIGHT) {
      dropj='y';
      if ((((*mapdistance)[j]-(*mapdistance)[j-1])<20)&&(mqmmod(f1genotype[j-1],11)!=0)) dropj='n';
      else if (position[j-1]!=MLEFT)
        if ((((*mapdistance)[j]-(*mapdistance)[j-2])<20)&&(mqmmod(f1genotype[j-2],11)!=0)) dropj='n';
    }
    if (dropj=='n') {
      marker[jj]= marker[j];
      (*cofactor)[jj]= (*cofactor)[j];
      (*mapdistance)[jj]= (*mapdistance)[j];
      chr[jj]= chr[j];
      r[jj]= r[j];
      position[jj]= position[j];
      jj++;
    } else if ((*cofactor)[j]==MH) {
      if (verbose) {
        info("Cofactor at chr %d is dropped",chr[j]);
      }
    }
  }
  Nmark= jj;
  if (verbose) {
    Rprintf("Num markers: %d",Nmark);
  }
  // FIXME this is duplication of code above - should be a (unit tested) method
  for (int j=0; j<Nmark; j++) {
    r[j]= 999.0;
    if (j==0) {
      if (chr[j]==chr[j+1]) position[j]=MLEFT;
      else position[j]=MUNLINKED;
    } else if (j==(Nmark-1)) {
      if (chr[j]==chr[j-1]) position[j]=MRIGHT;
      else position[j]=MUNLINKED;
    } else if (chr[j]==chr[j-1]) {
      if (chr[j]==chr[j+1]) position[j]=MMIDDLE;
      else position[j]=MRIGHT;
    } else {
      if (chr[j]==chr[j+1]) position[j]=MLEFT;
      else position[j]=MUNLINKED;
    }
  }
  for (int j=0; j<Nmark; j++) {
    if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
      r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
      if (r[j]<0) {
        Rprintf("ERROR: Position=%d r[j]=%f\n",position[j], r[j]);
        fatal("Recombination frequency is negative");
        return;
      }
    }
  }
  info("After dropping of uninformative cofactors");
  ivector newind;
  vector newy;
  cmatrix newmarker;
  double ymean=0.0, yvari=0.0;
  for (int i=0; i<Nind; i++) ymean += y[i];
  ymean/= Nind;
  for (int i=0; i<Nind; i++) yvari += pow(y[i]-ymean,2);
  yvari/= (Nind-1);
  //Fix for not doing dataaugmentation, we just copy the current as the augmented and set Naug to Nind
  Naug=Nind;
  Nind=out_Naug;
  newind= newivector(Naug);
  newy= newvector(Naug);
  newmarker= newcmatrix(Nmark,Naug);
  for (int i=0; i<Naug; i++) {
    newy[i]= y[i];
    newind[i]= INDlist[0][i];
    for (int j=0; j<Nmark; j++) {
      newmarker[j][i]= marker[j][i];
    }
  }

  vector newweight;
  newweight= newvector(Naug);
  //Re-estimation of recombinant frequencies
  double max;
  max = rmixture(newmarker, newweight, r, position, newind,Nind, Naug, Nmark, mapdistance,reestimate,crosstype,verbose);
  if (max > stepmax) {
    info("ERROR: Reestimation of the map put markers at: %f Cm",max);
    info("ERROR: Rerun the algorithm with a step.max larger than %f Cm",max);
    return;
  } else {
    if (verbose) {
      Rprintf("Reestimation of the map finished. MAX Cm: %f Cm",max);
    }
  }

  //Check if everything still is correct
  for (int j=0; j<Nmark; j++) {
    if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
      r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
      if (r[j]<0) {
        info("ERROR: Recombination frequency is negative");
        info("ERROR: Position=%d r[j]=%f",position[j], r[j]);
        return;
      }
    }
  }
  /* eliminate individuals with missing trait values */
  //We can skip this part iirc because R throws out missing phenotypes beforehand
  int oldNind=Nind;
  for (int i=0; i<oldNind; i++) {
    Nind-= ((y[i]==999.0) ? 1 : 0);
  }

  int oldNaug=Naug;
  for (int i=0; i<oldNaug; i++) {
    Naug-= ((newy[i]==999.0) ? 1 : 0);
  }

  vector weight;
  ivector ind;
  marker= newcmatrix(Nmark+1,Naug);
  y= newvector(Naug);
  ind= newivector(Naug);
  weight= newvector(Naug);
  int newi=0;
  for (int i=0; i<oldNaug; i++)
    if (newy[i]!=999.0) {
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
  delcmatrix(newmarker,Nmark);
  Free(newy);
  Free(newind);
  Free(newweight);

//    vector Fy;
//    Fy= newvector(Naug);
  double variance=-1.0;
  double logLfull;
  cvector selcofactor;
  selcofactor= newcvector(Nmark); /* selected cofactors */

  int dimx=1;
  for (int j=0; j<Nmark; j++) {
    if ((*cofactor)[j]==MH) {
      dimx+= (dominance=='n' ? 1 : 2);  // per QTL only additivity !!
    } else if ((*cofactor)[j]==MBB) {
      dimx+=1;  /* sex of the mouse */
    }
  }
  double F1, F2;

  F1= inverseF(1,Nind-dimx,alfa,verbose);
  F2= inverseF(2,Nind-dimx,alfa,verbose);
  if(verbose)info("dimX:%d nInd:%d",dimx,Nind);
  if(verbose)info("F(Threshold,Degrees of freedom 1,Degrees of freedom 2)=Alfa");
  if(verbose)info("F(%f,1,%d)=%f",F1,(Nind-dimx),alfa);
  if(verbose)info("F(%f,2,%d)=%f",F2,(Nind-dimx),alfa);
  F2= 2.0* F2; // 9-6-1998 using threshold x*F(x,df,alfa)

  weight[0]= -1.0;
  logLfull= QTLmixture(marker,(*cofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,useREML,fitQTL,dominance,crosstype,verbose);
  if(verbose)info("Log-likelihood of full model= %f",logLfull);
  if(verbose)info("Residual variance= %f",variance);
  if(verbose)info("Trait mean= %f; Trait variation= %f",ymean,yvari);
  if (Backwards==1)    // use only selected cofactors
    logLfull= backward(Nind, Nmark, (*cofactor), marker, y, weight, ind, Naug, logLfull,variance, F1, F2, &selcofactor, r, position, &informationcontent, mapdistance,&Frun,run,useREML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype,verbose);
  if (Backwards==0) // use all cofactors
    logLfull= mapQTL(Nind, Nmark, (*cofactor), (*cofactor), marker, position,(*mapdistance), y, r, ind, Naug, variance, 'n', &informationcontent,&Frun,run,useREML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype,verbose); // printout=='n'

  // ---- Write output / send it back to R
  //Cofactors that made it to the final model
  for (int j=0; j<Nmark; j++) {
    if (selcofactor[j]==MH) {
      (*cofactor)[j]=MH;
    } else {
      (*cofactor)[j]=MAA;
    }
  }
  //QTL likelyhood for each location
  if(verbose) info("Number of output datapoints: %d",Nsteps);

  //ofstream fff("MQM.output", ios::out | ios::app);
  for (int ii=0; ii<Nsteps; ii++) {
    //Convert LR to LOD before sending back
    QTL[0][ii] = Frun[ii][0] / 4.60517;
    QTL[0][Nsteps+ii] = informationcontent[ii];
    //char *outline;
    //Rprintf("LOC: %d QTL: %f INFO: %f\n",ii,QTL[0][ii],QTL[0][Nsteps+ii]);
    //fff << outline;
  }
  //fff.close();
  Free(position);
  Free(weight);
  Free(ind);
  Free(r);
  Free(informationcontent);
  freematrix((void **)Frun,Nsteps);
  delcmatrix(marker,Nmark+1);
  Free(y);
  Free(chr);
  Free(selcofactor);
  info("Analysis of data finished");
  return;
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
  cmatrix markers;
  cvector cofactor;
  vector mapdistance;

  markers= newcmatrix(Nmark+1,Nind);
  f1genotype = newivector(Nmark);
  cofactor= newcvector(Nmark);
  mapdistance= newvector(Nmark);

  int cof_cnt=0;

  MQMCrossType crosstype = determine_MQMCross(Nmark,Nind,(const int **)Geno,rqtlcrosstype);
  //Change all the markers from R/qtl format to MQM internal
  change_coding(&Nmark,&Nind,Geno,markers,crosstype);

  for (int i=0; i< Nmark; i++) {
    f1genotype[i] = 12;
    //receiving mapdistances
    mapdistance[i]=999.0;
    mapdistance[i]=Dist[0][i];
    cofactor[i] = MAA;
    if (Cofactors[0][i] == 1) {
      cofactor[i] = MH;
      cof_cnt++;
    }
    if (Cofactors[0][i] == 2) {
      cofactor[i] = MBB;
      cof_cnt++;
    }
    if (cof_cnt+10 > Nind) {
      Rprintf("ERROR: Setting this many cofactors would leave less than 10 degrees of freedom.\n");
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
    domi= domi;
  }

  char dominance='n';
  if (domi != 0) {
    dominance='y';
  }

  //WE HAVE EVERYTHING START WITH MAIN SCANNING FUNCTION
  analyseF2(Nind, Nmark, &cofactor, markers, Pheno[(Npheno-1)], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter,out_Naug,INDlist,reestimate,crosstype,dominance,verbose);

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
  delcmatrix(markers,Nmark+1);
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
               double *windowsize,double *steps,
               double *stepmi,double *stepma, int *nRun,int *out_Naug,int *indlist,  double *qtl,int *reestimate,int *crosstype,int *domi,int *verbose) {
  int **Geno;
  int **Chromo;
  double **Dist;
  double **Pheno;
  double **QTL;
  int **Cofactors;
  int **INDlist;

  //Reorganise the pointers into arrays, singletons are just cast into the function
  reorg_geno(*Nind,*Nmark,geno,&Geno);
  reorg_int(*Nmark,1,chromo,&Chromo);
  reorg_pheno(*Nmark,1,dist,&Dist);
  //Here we have  the assumption that step.min is negative this needs to be split in 2
  reorg_pheno(2*(*chromo) * (((*stepma)-(*stepmi))/ (*steps)),1,qtl,&QTL);
  reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
  reorg_int(*Nmark,1,cofactors,&Cofactors);
  reorg_int(*out_Naug,1,indlist,&INDlist);
  //Done with reorganising lets start executing

  mqmscan(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,*out_Naug,INDlist,QTL, *reestimate,(RqtlCrossType)*crosstype,*domi,*verbose);
} /* end of function R_mqmscan */

