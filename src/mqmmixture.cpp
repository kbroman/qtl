/**********************************************************************
 *
 * mqmmixture.cpp
 *
 * Copyright (c) 1996-2010 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Pjotr Prins and Danny Arends
 * last modified August 2010
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
 * ML estimation of recombination frequencies via EM; calculation of multilocus
 * genotype probabilities; ignorance of unlikely genotypes. Called by the
 * mqmscan.  maximum-likelihood estimation of recombination frequencies via the
 * EM algorithm, using multilocus information (default: the recombination
 * frequencies are not estimated but taken from mqm.in)
 *
 * When reestimate is 'n' the method is skipped
 */
double rmixture(MQMMarkerMatrix marker, vector weight, vector r,
                cvector position, ivector ind, int Nind, int Naug, int Nmark,
                vector *mapdistance, char reestimate, MQMCrossType crosstype,
                int verbose) {
  int i, j;
  int iem= 0;
  double Nrecom, oldr=0.0, newr, rdelta=1.0;
  double maximum = 0.0;
  float last_step = 0.0;
  vector indweight;
  indweight = newvector(Nind);
  vector distance;
  distance = newvector(Nmark+1);

  if (reestimate=='n') {
    if (verbose==1) {
      Rprintf("INFO: recombination parameters are not re-estimated\n");
    }
    for (j=0; j<Nmark; j++) {
      if (maximum < (*mapdistance)[j]) {
        maximum = (*mapdistance)[j];
      }
    }
  } else {
    if (verbose==1) {
      Rprintf("INFO: recombination parameters are re-estimated\n");
    }
    //Reestimation of map now works
    while ((iem<1000)&&(rdelta>0.0001)) {
      iem+=1;
      rdelta= 0.0;
      /* calculate weights = conditional genotype probabilities */
      for (i=0; i<Naug; i++) weight[i]=1.0;
      for (j=0; j<Nmark; j++) {
        if ((position[j]==MLEFT)||(position[j]==MUNLINKED))
          for (i=0; i<Naug; i++)
            if (marker[j][i]==MH) weight[i]*= 0.5;
            else weight[i]*= 0.25;
        if ((position[j]==MLEFT)||(position[j]==MMIDDLE))
          for (i=0; i<Naug; i++) {
            double calc_i = left_prob(r[j],marker[j][i],marker[j+1][i],crosstype);            //double calc_i = prob(marker, r, i, j, marker[j+1][i], crosstype, 0);
            weight[i]*=calc_i;
          }
      }
      for (i=0; i<Nind; i++) {
        indweight[i]= 0.0;
      }
      for (i=0; i<Naug; i++) {
        indweight[ind[i]]+=weight[i];
      }
      for (i=0; i<Naug; i++) {
        weight[i]/=indweight[ind[i]];
      }
      for (j=0; j<Nmark; j++) {
        if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
          newr= 0.0;
          for (i=0; i<Naug; i++) {
            Nrecom= fabs((double)(marker[j][i]-marker[j+1][i]));
            if ((marker[j][i]==MH)&&(marker[j+1][i]==MH))
              Nrecom= 2.0*r[j]*r[j]/(r[j]*r[j]+(1-r[j])*(1-r[j]));
            newr+= Nrecom*weight[i];
          }
          if (reestimate=='y' && position[j]!=MRIGHT) { //only update if it isn't the last marker of a chromosome ;)
            oldr=r[j];
            r[j]= newr/(2.0*Nind);
            rdelta+=pow(r[j]-oldr, 2.0);
          } else rdelta+=0.0;
        }
      }
    }
    /*   print new estimates of recombination frequencies */
    //Rprintf("INFO: Reestimate? %c\n", reestimate);
    //Rprintf("INFO: looping over all markers %d\n", Nmark);
    for (j=0; j<Nmark; j++) {
      if (position[j+1]==MRIGHT) {
        last_step = (*mapdistance)[j+1]-(*mapdistance)[j];
      }
      if (position[j]!=MLEFT) {
        if (position[j]!=MRIGHT) {
          (*mapdistance)[j]= -50*log(1-2.0*r[j])+(*mapdistance)[j-1];
        } else {
          (*mapdistance)[j]= (*mapdistance)[j-1]+last_step;
        }
      } else {
        (*mapdistance)[j]= -50*log(1-2.0*r[j]);
      }
      if (maximum < (*mapdistance)[j]) {
        maximum = (*mapdistance)[j];
      }
      //Rprintf("r(%d)= %f -> %f\n", j, r[j], (*mapdistance)[j]);
    }
  }
  if (verbose==1) {
    Rprintf("INFO: Re-estimation of the genetic map took %d iterations, to reach a rdelta of %f\n", iem, rdelta);
  }
  Free(indweight);
  freevector(distance);
  return maximum;
}


/* ML estimation of parameters in mixture model via EM; maximum-likelihood
 * estimation of parameters in the mixture model via the EM algorithm, using
 * multilocus information, but assuming known recombination frequencies
*/
double QTLmixture(MQMMarkerMatrix loci, cvector cofactor, vector r, cvector position,
                  vector y, ivector ind, int Nind, int Naug,
                  int Nloci,
                  double *variance, int em, vector *weight, const bool useREML,const bool fitQTL,const bool dominance, MQMCrossType crosstype, int verbose) {
                  
  //debug_trace("QTLmixture called Nloci=%d Nind=%d Naug=%d, REML=%d em=%d fit=%d domi=%d cross=%c\n",Nloci,Nind,Naug,useREML,em,fitQTL,dominance,crosstype);
  //for (int i=0; i<Nloci; i++){
  // debug_trace("loci %d : recombfreq=%f\n",i,r[i]);
  //}
  int iem= 0, newNaug, i, j;
  bool warnZeroDist=false;
  bool varknown;
  bool biasadj=false;
  double oldlogL=-10000, delta=1.0, calc_i, logP=0.0, Pscale=1.75;
  vector indweight, Ploci, Fy;

  indweight= newvector(Nind);
  newNaug= ((!fitQTL) ? Naug : 3*Naug);
  Fy= newvector(newNaug);
  logP= Nloci*log(Pscale);                          // only for computational accuracy
  debug_trace("logP:%f\n",logP);
  varknown= (((*variance)==-1.0) ? false : true );
  Ploci= newvector(newNaug);
  #ifndef STANDALONE
    R_CheckUserInterrupt(); /* check for ^C */
    //R_ProcessEvents();
    R_FlushConsole();
  #endif
  if ((useREML)&&(!varknown)) {
		//info("INFO: REML");
  }
  if (!useREML) {
		//info("INFO: Maximum Likelihood");
    varknown=false;
    biasadj=false;
  }
  for (i=0; i<newNaug; i++) {
    Ploci[i]= 1.0;
  }
  if (!fitQTL) {
    for (j=0; j<Nloci; j++) {
      for (i=0; i<Naug; i++)
        Ploci[i]*= Pscale;
      if ((position[j]==MLEFT)||(position[j]==MUNLINKED)) {
        for (i=0; i<Naug; i++) {
          calc_i = start_prob(crosstype, loci[j][i]);   // calc_i= prob(loci, r, i, j, MH, crosstype, 0, 1);
          Ploci[i]*= calc_i;
          //Als Ploci > 0 en calc_i > 0 then we want to assert Ploci[] != 0
        }
      }
      if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
        for (i=0; i<Naug; i++) {
          
          calc_i =left_prob(r[j],loci[j][i],loci[j+1][i],crosstype); //calc_i = prob(loci, r, i, j, loci[j+1][i], crosstype, 0);
          if(calc_i == 0.0){calc_i=1.0;warnZeroDist=true;}
          Ploci[i]*= calc_i;
        }
      }
    }
  } else {
    for (j=0; j<Nloci; j++) {
      for (i=0; i<Naug; i++) {
        Ploci[i]*= Pscale;           // only for computational accuracy; see use of logP
        Ploci[i+Naug]*= Pscale;      // only for computational accuracy; see use of logP
        Ploci[i+2*Naug]*= Pscale;    // only for computational accuracy; see use of logP
      }
      if ((position[j]==MLEFT)||(position[j]==MUNLINKED)) {
        if (cofactor[j]<=MCOF){
          for (i=0; i<Naug; i++) {
            calc_i = start_prob(crosstype, loci[j][i]);  // calc_i= prob(loci, r, i, j, MH, crosstype, 0, 1);
            Ploci[i] *= calc_i;
            Ploci[i+Naug] *= calc_i;
            Ploci[i+2*Naug] *= calc_i;
          }
        }else{
          for (i=0; i<Naug; i++) {
            Ploci[i]*= start_prob(crosstype, MAA);          //startvalue for MAA for new chromosome
            Ploci[i+Naug]*= start_prob(crosstype, MH);      //startvalue for MH for new chromosome
            Ploci[i+2*Naug] *= start_prob(crosstype, MBB);  //startvalue for MBB for new chromosome
          }
        }
      }
      if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
        if ((cofactor[j]<=MCOF)&&(cofactor[j+1]<=MCOF))
          for (i=0; i<Naug; i++) {
            calc_i =left_prob(r[j],loci[j][i],loci[j+1][i],crosstype);  //calc_i = prob(loci, r, i, j, loci[j+1][i], crosstype, 0);
            if(calc_i == 0.0){calc_i=1.0;warnZeroDist=true;}
            Ploci[i]*= calc_i;
            Ploci[i+Naug]*= calc_i;
            Ploci[i+2*Naug]*= calc_i;
          }
        else if (cofactor[j]<=MCOF) // locus j+1 == QTL
          for (i=0; i<Naug; i++) { // QTL==MAA || MH || MBB means: What is the prob of finding an MAA at J=1    
            calc_i =left_prob(r[j],loci[j][i],MAA,crosstype);     //calc_i = prob(loci, r, i, j, MAA, crosstype, 0);
            Ploci[i]*= calc_i;      
            calc_i = left_prob(r[j],loci[j][i],MH,crosstype);     //calc_i = prob(loci, r, i, j, MH, crosstype, 0);
            Ploci[i+Naug]*= calc_i;
            calc_i = left_prob(r[j],loci[j][i],MBB,crosstype);    //calc_i = prob(loci, r, i, j, MBB, crosstype, 0);
            Ploci[i+2*Naug]*= calc_i;
          }
        else // locus j == QTL
          for (i=0; i<Naug; i++) { // QTL==MQTL
            calc_i = left_prob(r[j],MAA,loci[j+1][i],crosstype);  //calc_i = prob(loci, r, i, j+1, MAA, crosstype, -1);
            Ploci[i]*= calc_i;
            calc_i = left_prob(r[j],MH,loci[j+1][i],crosstype);   //calc_i = prob(loci, r, i, j+1, MH, crosstype, -1);
            Ploci[i+Naug]*= calc_i;
            calc_i = left_prob(r[j],MBB,loci[j+1][i],crosstype);  //calc_i = prob(loci, r, i, j+1, MBB, crosstype, -1);
            Ploci[i+2*Naug]*= calc_i;
          }
      }
    }
  }
  if(warnZeroDist)info("!!! 0.0 from Prob !!! Markers at same Cm but different genotype !!!"); 
//	Rprintf("INFO: Done fitting QTL's\n");
  if ((*weight)[0]== -1.0) {
    for (i=0; i<Nind; i++) indweight[i]= 0.0;
    if (!fitQTL) {
      for (i=0; i<Naug; i++) indweight[ind[i]]+=Ploci[i];
      for (i=0; i<Naug; i++) (*weight)[i]= Ploci[i]/indweight[ind[i]];
    } else {
      for (i=0; i<Naug; i++) indweight[ind[i]]+=Ploci[i]+Ploci[i+Naug]+Ploci[i+2*Naug];
      for (i=0; i<Naug; i++) {
        (*weight)[i]       = Ploci[i]/indweight[ind[i]];
        (*weight)[i+Naug]  = Ploci[i+Naug]/indweight[ind[i]];
        (*weight)[i+2*Naug]= Ploci[i+2*Naug]/indweight[ind[i]];
      }
    }
  }
  debug_trace("Weights done\n");
  debug_trace("Individual->trait,indweight weight Ploci\n");
  //for (int j=0; j<Nind; j++){
  //  debug_trace("%d->%f,%f %f %f\n", j, y[j],indweight[i], (*weight)[j], Ploci[j]);
  //}
  double logL=0;
  vector indL;
  indL= newvector(Nind);
  while ((iem<em)&&(delta>1.0e-5)) {
    iem+=1;
    if (!varknown) *variance=-1.0;
    logL= regression(Nind, Nloci, cofactor, loci, y,
                     weight, ind, Naug, variance, Fy, biasadj, fitQTL, dominance);
    logL=0.0;
    for (i=0; i<Nind; i++) indL[i]= 0.0;
    if (!fitQTL) // no QTL fitted
      for (i=0; i<Naug; i++) {
        (*weight)[i]= Ploci[i]*Fy[i];
        indL[ind[i]]= indL[ind[i]] + (*weight)[i];
      }
    else // QTL moved along the chromosomes
      for (i=0; i<Naug; i++) {
        (*weight)[i]= Ploci[i]*Fy[i];
        (*weight)[i+Naug]  = Ploci[i+Naug]*  Fy[i+Naug];
        (*weight)[i+2*Naug]= Ploci[i+2*Naug]*Fy[i+2*Naug];
        indL[ind[i]]+=(*weight)[i]+(*weight)[i+Naug]+(*weight)[i+2*Naug];
      }
    for (i=0; i<Nind; i++) logL+=log(indL[i])-logP;
    for (i=0; i<Nind; i++) indweight[i]= 0.0;
    if (!fitQTL) {
      for (i=0; i<Naug; i++) indweight[ind[i]]+=(*weight)[i];
      for (i=0; i<Naug; i++) (*weight)[i]/=indweight[ind[i]];
    } else {
      for (i=0; i<Naug; i++)
        indweight[ind[i]]+=(*weight)[i]+(*weight)[i+Naug]+(*weight)[i+2*Naug];
      for (i=0; i<Naug; i++) {
        (*weight)[i]       /=indweight[ind[i]];
        (*weight)[i+Naug]  /=indweight[ind[i]];
        (*weight)[i+2*Naug]/=indweight[ind[i]];
      }
    }
    delta= fabs(logL-oldlogL);
    oldlogL= logL;
  }
  // bias adjustment after finished ML estimation via EM
  if ((useREML)&&(!varknown)) {
    *variance=-1.0;
    biasadj=true;
    logL= regression(Nind, Nloci, cofactor, loci, y,
                     weight, ind, Naug, variance, Fy, biasadj, fitQTL, dominance);
    logL=0.0;
    for (int _i=0; _i<Nind; _i++) indL[_i]= 0.0;
    if (!fitQTL)
      for (i=0; i<Naug; i++) {
        (*weight)[i]= Ploci[i]*Fy[i];
        indL[ind[i]]+=(*weight)[i];
      }
    else
      for (i=0; i<Naug; i++) {
        (*weight)[i]= Ploci[i]*Fy[i];
        (*weight)[i+Naug]= Ploci[i+Naug]*Fy[i+Naug];
        (*weight)[i+2*Naug]= Ploci[i+2*Naug]*Fy[i+2*Naug];
        indL[ind[i]]+=(*weight)[i];
        indL[ind[i]]+=(*weight)[i+Naug];
        indL[ind[i]]+=(*weight)[i+2*Naug];
      }
    for (i=0; i<Nind; i++) logL+=log(indL[i])-logP;
    for (i=0; i<Nind; i++) indweight[i]= 0.0;
    if (!fitQTL) {
      for (i=0; i<Naug; i++) indweight[ind[i]]+=(*weight)[i];
      for (i=0; i<Naug; i++) (*weight)[i]/=indweight[ind[i]];
    } else {
      for (i=0; i<Naug; i++) {
        indweight[ind[i]]+=(*weight)[i];
        indweight[ind[i]]+=(*weight)[i+Naug];
        indweight[ind[i]]+=(*weight)[i+2*Naug];
      }
      for (i=0; i<Naug; i++) {
        (*weight)[i]       /=indweight[ind[i]];
        (*weight)[i+Naug]  /=indweight[ind[i]];
        (*weight)[i+2*Naug]/=indweight[ind[i]];
      }
    }
  }
  //for (i=0; i<Nind; i++){
  //  debug_trace("IND %d Ploci: %f Fy: %f UNLOG:%f LogL:%f LogL-LogP: %f\n", i, Ploci[i], Fy[i], indL[i], log(indL[i]), log(indL[i])-logP);
  //}
  Free(Fy);
  Free(Ploci);
  Free(indweight);
  Free(indL);
  return logL;
}

