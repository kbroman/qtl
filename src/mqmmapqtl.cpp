/**********************************************************************
 *
 * mqmmapqtl.cpp
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

/* 
 * mapQTL moves a QTL along the chromosome and calculated at each map position
 * the QTL likelihood. Uses either all cofactors, or selected cofactors only
 */
double mapQTL(int Nind, int Nmark, cvector cofactor, cvector selcofactor, 
              MQMMarkerMatrix marker, cvector position, vector mapdistance, vector y,
              vector r, ivector ind, int Naug, double variance, char
              printoutput, vector *informationcontent, matrix *Frun, int run,
              char REMLorML, bool fitQTL, bool dominance, int em, double
              windowsize, double stepsize, double stepmin, double stepmax, 
              MQMCrossType crosstype, int verbose) {
  //Rprintf("INFO: mapQTL function called.\n");
  int Nloci, j, jj, jjj=0;
  vector Fy;
  Fy= newvector(Naug);
  cvector QTLcofactor, saveQTLcofactor;
  QTLcofactor= newcvector(Nmark+1);
  saveQTLcofactor= newcvector(Nmark+1);
  double infocontent;
  vector info0, info1, info2, weight;

  info0= newvector(Nind);
  info1= newvector(Nind);
  info2= newvector(Nind);

  weight= newvector(Naug);
  weight[0]= -1.0;

  /* fit QTL on top of markers (but: should also be done with routine QTLmixture() for exact ML) */

  cvector newcofactor;
  cvector direction;
  newcofactor= newcvector(Nmark);
  direction = newcvector(Nmark);
  vector cumdistance;
  double QTLlikelihood=0.0;
  cumdistance= newvector(Nmark+1);
  for (j=0; j<Nmark; j++) {
    if (position[j]==MLEFT)
      cumdistance[j]= -50*log(1-2.0*r[j]);
    else if (position[j]==MMIDDLE)
      cumdistance[j]= cumdistance[j-1]-50*log(1-2.0*r[j]);
  }
  double savelogL=0.0; // log-likelihood of model with all selected cofactors

  /* fit QTL on top of markers (full ML)   fit QTL between markers (full ML) */
  // cout << "please wait (mixture calculus may take quite a lot of time)" << endl;
  /* estimate variance in mixture model with all marker cofactors */
  // cout << "estimate variance in mixture model with all cofactors" << endl;

  variance= -1.0;
  savelogL= 2.0*QTLmixture(marker, cofactor, r, position, y, ind, Nind, Naug, Nmark, &variance, em, &weight, REMLorML, fitQTL, dominance, crosstype, verbose);
  if (verbose==1) {
    Rprintf("INFO: log-likelihood of full model= %f\n", savelogL/2);
  }
  Nloci= Nmark+1;
  // augment data for missing QTL observations (x 3)
  fitQTL=true;
  int newNaug;
  newNaug= 3*Naug;
  Free(weight);
  weight= newvector(newNaug);
  weight[0]= 1.0;
  vector weight0;
  weight0= newvector(newNaug);
  weight0[0]= -1.0;

//       augmentdataforQTL(marker);
  vector QTLr, QTLmapdistance;
  QTLr= newvector(Nloci);
  QTLmapdistance= newvector(Nloci);
  cvector QTLposition;
  QTLposition= newcvector(Nloci);
  MQMMarkerMatrix QTLloci;

  QTLloci = (MQMMarkerMatrix)Calloc(Nloci, MQMMarkerVector);

  //  Rprintf("DEBUG testing_2");
  double moveQTL= stepmin;
  char nextinterval= 'n', firsttime='y';
  double maxF=0.0, savebaseNoQTLModel=0.0;
  int baseNoQTLModel=0, step=0;
  //Rprintf("DEBUG testing_3");
  for (j=0; j<Nmark; j++) {
    /* 	fit a QTL in two steps:
    1. move QTL along marker interval j -> j+1 with steps of stepsize=20 cM, starting from -20 cM up to 220 cM
    2. all marker-cofactors in the neighborhood of the QTL are dropped by using cM='windows' as criterium
    */
    nextinterval= 'n';
#ifndef STANDALONE
    //Rprintf("TEST mqmmapqtl\n");
    R_CheckUserInterrupt(); /* check for ^C */
    R_FlushConsole();
    //R_ProcessEvents(); /*  Try not to crash windows etc*/
#endif
    while (nextinterval=='n') { // step 1:
      // Rprintf("DEBUG testing STEP 1");
      if (position[j]==MLEFT) {
        if (moveQTL<=mapdistance[j]) {
          QTLposition[j]= position[j];
          QTLposition[j+1]= MMIDDLE;
          QTLr[j]= recombination_frequentie((mapdistance[j]-moveQTL));
          QTLr[j+1]= r[j];
          QTLloci[j+1]= marker[j];
          QTLloci[j]= marker[Nloci-1];
          QTLmapdistance[j]= moveQTL;
          QTLmapdistance[j+1]= mapdistance[j];
          if (firsttime=='y') weight[0]= -1.0;
          moveQTL+= stepsize;
        } else if (moveQTL<=mapdistance[j+1]) {
          QTLposition[j]= position[j];
          QTLposition[j+1]= MMIDDLE;
          QTLr[j]= recombination_frequentie((moveQTL-mapdistance[j]));
          QTLr[j+1]= recombination_frequentie((mapdistance[j+1]-moveQTL)); //r[j];
          QTLloci[j]= marker[j];
          QTLloci[j+1]= marker[Nloci-1];
          QTLmapdistance[j]= mapdistance[j];
          QTLmapdistance[j+1]= moveQTL;
          moveQTL+= stepsize;
        } else nextinterval= 'y';
      } else if (position[j]==MMIDDLE) {
        if (moveQTL<=mapdistance[j+1]) {
          QTLposition[j]= position[j];
          QTLposition[j+1]= MMIDDLE;
          QTLr[j]= recombination_frequentie((moveQTL-mapdistance[j])); //0.0;
          QTLr[j+1]= recombination_frequentie((mapdistance[j+1]-moveQTL)); //r[j];
          QTLloci[j]= marker[j];
          QTLloci[j+1]= marker[Nloci-1];
          QTLmapdistance[j]= mapdistance[j];
          QTLmapdistance[j+1]= moveQTL;
          moveQTL+= stepsize;
        } else nextinterval= 'y';
      } else if (position[j]==MRIGHT) {
        if (moveQTL<=stepmax) {
          QTLposition[j]= MMIDDLE;
          QTLposition[j+1]= MRIGHT;
          QTLr[j]= recombination_frequentie((moveQTL-mapdistance[j])); //0.0;
          QTLr[j+1]= r[j]; // note r[j]=999.0
          QTLloci[j]= marker[j];
          QTLloci[j+1]= marker[Nloci-1];
          QTLmapdistance[j]= mapdistance[j];
          QTLmapdistance[j+1]= moveQTL;
          moveQTL+= stepsize;
        } else {
          nextinterval= 'y';
          moveQTL= stepmin;
        }
      } else if (position[j]==MUNLINKED) {
        QTLposition[j]= MLEFT;
        QTLposition[j+1]= MRIGHT; //position[j] ?? MRIGHT ?
        QTLr[j]= 0.0;
        QTLr[j+1]= r[j];
        QTLloci[j+1]= marker[j];
        QTLloci[j]= marker[Nloci-1];
        QTLmapdistance[j]= mapdistance[j];
        QTLmapdistance[j+1]= mapdistance[j];
        if (firsttime=='y') weight[0]= -1.0;
        nextinterval= 'y';
        moveQTL= stepmin;
      }
      if (nextinterval=='n') { // QTLcofactor[j]= MAA;
        // QTLcofactor[j+1]= MAA;
        for (jj=0; jj<j; jj++) {
          QTLposition[jj]= position[jj];
          QTLr[jj]= r[jj];
          QTLloci[jj]= marker[jj];
          QTLmapdistance[jj]= mapdistance[jj];
          QTLcofactor[jj]= selcofactor[jj];
        }
        for (jj=j+1; jj<Nmark; jj++) {
          QTLposition[jj+1]= position[jj];
          QTLr[jj+1]= r[jj];
          QTLloci[jj+1]= marker[jj];
          QTLcofactor[jj+1]= selcofactor[jj];
          QTLmapdistance[jj+1]= mapdistance[jj];
          QTLcofactor[jj+1]= selcofactor[jj];
        }
        // step 2:
        //  Rprintf("DEBUG testing STEP 2");
        if ((position[j]==MLEFT)&&((moveQTL-stepsize)<=mapdistance[j])) {
          QTLcofactor[j]= MNOCOF;
          QTLcofactor[j+1]=
            (((QTLmapdistance[j+1]-QTLmapdistance[j])<windowsize) ? MNOCOF : selcofactor[j]);
        } else {
          QTLcofactor[j+1]= MNOCOF;
          QTLcofactor[j]=
            (((QTLmapdistance[j+1]-QTLmapdistance[j])<windowsize) ? MNOCOF : selcofactor[j]);
        }
        if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
          jjj=j+2;
          while (QTLposition[jjj]==MMIDDLE) {
            if ((position[j]==MLEFT)&&((moveQTL-stepsize)<=mapdistance[j]))
              QTLcofactor[jjj]=
                (((QTLmapdistance[jjj]-QTLmapdistance[j])<windowsize) ? MNOCOF : QTLcofactor[jjj]);
            else
              QTLcofactor[jjj]=
                (((QTLmapdistance[jjj]-QTLmapdistance[j+1])<windowsize) ? MNOCOF : QTLcofactor[jjj]);
            jjj++;
          }
          QTLcofactor[jjj]=
            (((QTLmapdistance[jjj]-QTLmapdistance[j+1])<windowsize) ? MNOCOF : QTLcofactor[jjj]);
        }
        if ((position[j]==MMIDDLE)||(position[j]==MRIGHT)) {
          jjj=j-1;
          while (QTLposition[jjj]==MMIDDLE) {
            QTLcofactor[jjj]= (((QTLmapdistance[j+1]-QTLmapdistance[jjj])<windowsize) ? MNOCOF : QTLcofactor[jjj]);
            jjj--;
          }
          QTLcofactor[jjj]= (((QTLmapdistance[j+1]-QTLmapdistance[jjj])<windowsize) ? MNOCOF : QTLcofactor[jjj]);
        }

        // fit no-QTL model at current map position (cofactors only)
        if (firsttime=='y') {
          for (jj=0; jj<Nloci; jj++) saveQTLcofactor[jj]= QTLcofactor[jj];
          baseNoQTLModel=1;
          firsttime='n';
        } else {
          baseNoQTLModel=0;
          for (jj=0; jj<Nloci; jj++) baseNoQTLModel+= (saveQTLcofactor[jj]==QTLcofactor[jj] ? 0 : 1);
        }
        //  Rprintf("fitting NO-QTL model\n");
        if (baseNoQTLModel!=0) { // new base no-QTL model
          if ((position[j]==MLEFT)&&((moveQTL-stepsize)<=mapdistance[j])) QTLcofactor[j]= MSEX;
          else QTLcofactor[j+1]= MSEX;
          // Rprintf("INFO: Before base model\n", QTLlikelihood/-2);
          QTLlikelihood= -2.0*QTLmixture(QTLloci, QTLcofactor, QTLr, QTLposition, y, ind, Nind, Naug, Nloci, &variance, em, &weight0, REMLorML, fitQTL, dominance, crosstype, verbose);
          // Rprintf("INFO: log-likelihood of NO QTL model= %f\n", QTLlikelihood/-2);
          weight0[0]= -1.0;
          savebaseNoQTLModel= QTLlikelihood;
          if ((position[j]==MLEFT)&&((moveQTL-stepsize)<=mapdistance[j])) QTLcofactor[j]= MNOCOF;
          else QTLcofactor[j+1]= MNOCOF;
          for (jj=0; jj<Nloci; jj++) saveQTLcofactor[jj]= QTLcofactor[jj];
        } else
          QTLlikelihood= savebaseNoQTLModel;

        // fit QTL-model (plus cofactors) at current map position
        // MNOTAA= QTL
        if ((position[j]==MLEFT)&&((moveQTL-stepsize)<=mapdistance[j])) QTLcofactor[j]= MQTL;
        else QTLcofactor[j+1]= MQTL;
        if (REMLorML==MH) weight[0]= -1.0;
        QTLlikelihood+=2.0*QTLmixture(QTLloci, QTLcofactor, QTLr, QTLposition, y, ind, Nind, Naug, Nloci, &variance, em, &weight, REMLorML, fitQTL, dominance, crosstype, verbose);
        //this is the place we error at, because the likelihood is not correct.
        if (QTLlikelihood<-0.05) {
          Rprintf("WARNING: Negative QTLlikelihood=%f versus BASE MODEL: %f\nThis applies to the QTL at %d\n", QTLlikelihood, (savebaseNoQTLModel/-2), j); //return 0;}
        }
        maxF= (maxF<QTLlikelihood ? QTLlikelihood : maxF);
        if (run>0) (*Frun)[step][run]+= QTLlikelihood;
        else (*Frun)[step][0]+= QTLlikelihood;

        /* 	Each individual has condition multilocus probabilities for being 0, 1 or 2 at the QTL.
        Calculate the maximum per individu. Calculate the mean of this maximum, averaging over all individuals
        This is the information content plotted.
        */
        infocontent= 0.0;
        for (int i=0; i<Nind; i++) {
          info0[i]= 0.0; // qq
          info1[i]= 0.0; // Qq
          info2[i]= 0.0; // QQ
        }
        for (int i=0; i<Naug; i++) {
          info0[ind[i]]+= weight[i];
          info1[ind[i]]+= weight[i+Naug];
          info2[ind[i]]+= weight[i+2*Naug];
        }
        for (int i=0; i<Nind; i++)
          if (info0[i]<info1[i]) infocontent+= (info1[i]<info2[i] ? info2[i] : info1[i]);
          else infocontent+= (info0[i]<info2[i] ? info2[i] : info0[i]);
        (*informationcontent)[step]+=infocontent/Nind;

        step++;
      }
    }
  }
  fitQTL=false;

  freevector(direction);
  Free(info0);
  Free(info1);
  Free(info2);
  Free(weight);
  Free(weight0);
  Free(QTLr);
  Free(QTLposition);
  Free(Fy);
  Free(newcofactor);
  Free(QTLcofactor);
  Free(cumdistance);
  Free(QTLmapdistance);
  Free(QTLloci);
  Free(saveQTLcofactor);
  return maxF; //QTLlikelihood;
}

