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

double Lnormal(double residual, double variance) {
  //double Likelihood,likelyhood;
  //Likelihood=exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
  //if(fabs(Likelihood-likelyhood)>0.05){
  //Rprintf("ERROR: Lnormal error\n");
  //}
  //return likelyhood;
  return dnorm(residual,0,sqrt(variance),0);
}


/**********************************************************************
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
  int i;

  *Geno = (int **)R_alloc(n_pos, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_pos; i++)
    (*Geno)[i] = (*Geno)[i-1] + n_ind;

}

*/
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
               reestimate, char crosstype, char dominance, int verbose) {
  if (verbose) {
    Rprintf("INFO: Starting C-part of the MQM analysis\n");
  }
  int Naug;
  int run=0;
  cvector position;
  vector informationcontent;
  //char dominance='n';
  //char perm_simu='1';
  ivector chr;
  matrix Frun;
  vector r;
  r= newvector(Nmark);
  position= newcvector(Nmark);
  char REMLorML='0';
  char fitQTL='n';

  chr= newivector(Nmark);
  if (verbose) {
    Rprintf("INFO: Receiving the chromosome matrix from R\n");
  }
  for (int i=0; i< Nmark; i++) {
    chr[i] = Chromo[0][i];
  }

  if (RMLorML == 1) {
    REMLorML='1';
  }

  if (verbose) {
    Rprintf("INFO: Calculating relative genomepositions of the markers\n");
  }
  for (int j=0; j<Nmark; j++) {
    if (j==0) {
      if (chr[j]==chr[j+1]) position[j]='L';
      else position[j]='U';
    } else if (j==(Nmark-1)) {
      if (chr[j]==chr[j-1]) position[j]='R';
      else position[j]='U';
    } else if (chr[j]==chr[j-1]) {
      if (chr[j]==chr[j+1]) position[j]='M';
      else position[j]='R';
    } else {
      if (chr[j]==chr[j+1]) position[j]='L';
      else position[j]='U';
    }
  }

  if (verbose) {
    Rprintf("INFO: Estimating recombinant frequencies\n");
  }
  for (int j=0; j<Nmark; j++) {
    r[j]= 999.0;
    if ((position[j]=='L')||(position[j]=='M')) {
      r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
    }
    //Rprintf("R[j] value: %f\n",r[j]);
  }
  // ---- Initialize Frun and informationcontent to 0.0
  if (verbose) {
    Rprintf("INFO: Initialize Frun and informationcontent to 0.0\n");
  }
  int Nsteps;
  Nsteps= chr[Nmark-1]*((stepmax-stepmin)/stepsize+1);
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
    if (mod(f1genotype[j],11)!=0) {
      dropj='n';
    } else if ((*cofactor)[j]=='0') {
      dropj='y';
    } else if (position[j]=='L') {
      // (cofactor[j]!='0') cofactor at non-segregating marker
      // test whether next segregating marker is nearby (<20cM)
      dropj='y';
      if ((((*mapdistance)[j+1]-(*mapdistance)[j])<20)&&(mod(f1genotype[j+1],11)!=0)) dropj='n';
      else if (position[j+1]!='R')
        if ((((*mapdistance)[j+2]-(*mapdistance)[j])<20)&&(mod(f1genotype[j+2],11)!=0)) dropj='n';
    } else if (position[j]=='M') {
      dropj='y';
      if ((((*mapdistance)[j]-(*mapdistance)[j-1])<20)&&(mod(f1genotype[j-1],11)!=0)) dropj='n';
      else if ((((*mapdistance)[j+1]-(*mapdistance)[j])<20)&&(mod(f1genotype[j+1],11)!=0)) dropj='n';
    } else if (position[j]=='R') {
      dropj='y';
      if ((((*mapdistance)[j]-(*mapdistance)[j-1])<20)&&(mod(f1genotype[j-1],11)!=0)) dropj='n';
      else if (position[j-1]!='L')
        if ((((*mapdistance)[j]-(*mapdistance)[j-2])<20)&&(mod(f1genotype[j-2],11)!=0)) dropj='n';
    }
    if (dropj=='n') {
      marker[jj]= marker[j];
      (*cofactor)[jj]= (*cofactor)[j];
      (*mapdistance)[jj]= (*mapdistance)[j];
      chr[jj]= chr[j];
      r[jj]= r[j];
      position[jj]= position[j];
      jj++;
    } else if ((*cofactor)[j]=='1') {
      if (verbose) {
        Rprintf("INFO: Cofactor at chr %d is dropped\n",chr[j]);
      }
    }
  }
  Nmark= jj;
  if (verbose) {
    Rprintf("INFO: Num markers: %d\n",Nmark);
  }
  for (int j=0; j<Nmark; j++) {
    r[j]= 999.0;
    if (j==0) {
      if (chr[j]==chr[j+1]) position[j]='L';
      else position[j]='U';
    } else if (j==(Nmark-1)) {
      if (chr[j]==chr[j-1]) position[j]='R';
      else position[j]='U';
    } else if (chr[j]==chr[j-1]) {
      if (chr[j]==chr[j+1]) position[j]='M';
      else position[j]='R';
    } else {
      if (chr[j]==chr[j+1]) position[j]='L';
      else position[j]='U';
    }
  }
  for (int j=0; j<Nmark; j++) {
    if ((position[j]=='L')||(position[j]=='M')) {
      r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
      if (r[j]<0) {
        Rprintf("ERROR: Recombination frequency is negative\n");
        Rprintf("ERROR: Position=%d r[j]=%f\n",position[j], r[j]);
        return;
      }
    }
  }
  if (verbose) {
    Rprintf("INFO: After dropping of uninformative cofactors\n");
  }
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
    Rprintf("ERROR: Reestimation of the map put markers at: %f Cm\n",max);
    Rprintf("ERROR: Rerun the algorithm with a step.max larger than %f Cm\n",max);
    return;
  } else {
    if (verbose) {
      Rprintf("INFO: Reestimation of the map finished. MAX Cm: %f Cm\n",max);
    }
  }

  //Check if everything still is correct
  for (int j=0; j<Nmark; j++) {
    if ((position[j]=='L')||(position[j]=='M')) {
      r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
      if (r[j]<0) {
        Rprintf("ERROR: Recombination frequency is negative\n");
        Rprintf("ERROR: Position=%d r[j]=%f\n",position[j], r[j]);
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
  marker= newcmatrix(Nmark,Naug);
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
  delcmatrix(newmarker);
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
    if ((*cofactor)[j]=='1') {
      dimx+= (dominance=='n' ? 1 : 2);  // per QTL only additivity !!
    } else if ((*cofactor)[j]=='2') {
      dimx+=1;  /* sex of the mouse */
    }
  }
  double F1, F2;

  F1= inverseF(1,Nind-dimx,alfa,verbose);
  F2= inverseF(2,Nind-dimx,alfa,verbose);
  if (verbose) {
    Rprintf("INFO: dimX:%d nInd:%d\n",dimx,Nind);
    Rprintf("INFO: F(Threshold,Degrees of freedom 1,Degrees of freedom 2)=Alfa\n");
    Rprintf("INFO: F(%f,1,%d)=%f\n",F1,(Nind-dimx),alfa);
    Rprintf("INFO: F(%f,2,%d)=%f\n",F2,(Nind-dimx),alfa);
  }
  F2= 2.0* F2; // 9-6-1998 using threshold x*F(x,df,alfa)

  weight[0]= -1.0;
  logLfull= QTLmixture(marker,(*cofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance,crosstype,verbose);
  if (verbose) {
    Rprintf("INFO: Log-likelihood of full model= %f\n",logLfull);
    Rprintf("INFO: Residual variance= %f\n",variance);
    Rprintf("INFO: Trait mean= %f \nINFO: Trait variation= %f\n",ymean,yvari);
  }
  if (Backwards==1)    // use only selected cofactors
    logLfull= backward(Nind, Nmark, (*cofactor), marker, y, weight, ind, Naug, logLfull,variance, F1, F2, &selcofactor, r, position, &informationcontent, mapdistance,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype,verbose);
  if (Backwards==0) // use all cofactors
    logLfull= mapQTL(Nind, Nmark, (*cofactor), (*cofactor), marker, position,(*mapdistance), y, r, ind, Naug, variance, 'n', &informationcontent,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype,verbose); // printout=='n'

  // ---- Write output / send it back to R
  //Cofactors that made it to the final model
  for (int j=0; j<Nmark; j++) {
    if (selcofactor[j]=='1') {
      (*cofactor)[j]='1';
    } else {
      (*cofactor)[j]='0';
    }
  }
  //QTL likelyhood for each location
  if (verbose) {
    Rprintf("INFO: Number of output datapoints: %d\n",Nsteps);
  }
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
  Free(Frun);
  delcmatrix(marker);
  Free(y);
  Free(chr);
  Free(selcofactor);
  if (verbose) {
    Rprintf("INFO: Analysis of data finished\n");
  }
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
             double Stepmi,double Stepma,int NRUN,int out_Naug,int **INDlist, double **QTL, int re_estimate,int crosstype,int domi,int verbose) {

  ivector f1genotype;
  cmatrix markers;
  cvector cofactor;
  vector mapdistance;

  markers= newcmatrix(Nmark,Nind);
  f1genotype = newivector(Nmark);
  cofactor= newcvector(Nmark);
  mapdistance= newvector(Nmark);

  int cof_cnt=0;

  //Change all the markers from R/qtl format to MQM internal
  change_coding(&Nmark,&Nind,Geno,markers, crosstype);

  for (int i=0; i< Nmark; i++) {
    f1genotype[i] = 12;
    //receiving mapdistances
    mapdistance[i]=999.0;
    mapdistance[i]=Dist[0][i];
    cofactor[i] = '0';
    if (Cofactors[0][i] == 1) {
      cofactor[i] = '1';
      cof_cnt++;
    }
    if (Cofactors[0][i] == 2) {
      cofactor[i] = '2';
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
  char cross = determin_cross(&Nmark,&Nind,Geno,&crosstype);
  //set dominance accordingly
  if (cross != 'F') {
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
  analyseF2(Nind, Nmark, &cofactor, markers, Pheno[(Npheno-1)], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter,out_Naug,INDlist,reestimate,cross,dominance,verbose);

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
  delcmatrix(markers);
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

  mqmscan(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,*out_Naug,INDlist,QTL, *reestimate,*crosstype,*domi,*verbose);
} /* end of function R_mqmscan */

