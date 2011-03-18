/**********************************************************************
 *
 * mqmaugment.cpp
 *
 * Copyright (c) 1996-2011 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Pjotr Prins and Danny Arends
 * last modified Feb 2011
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
#include "simulate.h"

/*
 * Augment/expand the dataset by adding additional (likely) genotypes.
 * genotypic marker information is incomplete when
 *           - marker scores are missing or dominant, or
 *           - markers are non-segregating but you want to assume a
 *             segregating QTL on top of this marker (~very closely linked),
 *             and all likely possible configurations are generated 
 *
 * Inputs are number of markers Nmark, the marker matrix, position vector,
 * recombinations r.  The minprob parameter drops individuals from the dataset
 * (the value should be between 1..n).
 *
 * The minprob parameter drops genotypes. E.g. for minprob=0.01 eliminate
 * genotypes 100 times less likely than the most likely configuration.
 *
 * A new markerset is created and returned in augmarker, likewise the
 * phenotypes are returned in augy and individuals in augind. Augmentation
 * halts when the number of individuals maxNaug is reached. Individuals are
 * expanded up to imaxNaug.
 *
 * returns 1 on success, 0 on failure. marker, Nind, Naug, augy, augind are
 * changed to reflect the newly augmented dataset. Nind may be smaller because
 * of dropped individuals.
 *
 */
 
int calculate_augmentation(const int Nind, int const Nmark,const MQMMarkerMatrix markers, const MQMCrossType crosstype){
  unsigned int augmentationfactor=2;                  //RIL and/or backcross
  if(crosstype == CF2){
    augmentationfactor=3;                             //F2 population
  }
  for(int i=0; i<Nind; i++) {
    unsigned int augind=1;                            //How many times did we augment this individual
    int missingmarkers=0;                             //How many markers are missing for this individual
    bool outoflimit = false;
    for(int j=0; j<Nmark;j++){
      switch (markers[j][i]) {
        case MMISSING:
          if(!outoflimit) augind=augind*augmentationfactor;
          missingmarkers++;
        break;
        case MNOTAA:
          if(!outoflimit) augind=augind*(augmentationfactor-1);
          missingmarkers++;
        break;
        case MNOTBB:
          if(!outoflimit) augind=augind*(augmentationfactor-1);
          missingmarkers++;
        break;
        default:
        break;
      }
      if(augind >  UINT_MAX/augmentationfactor){
        outoflimit = true;
      }
    }
    if(!outoflimit){
      info("Individual: %d has %d missing markers, leading to %d augmentations",i,missingmarkers,augind);
    }else{
      info("Individual: %d has %d missing markers",i,missingmarkers);
    }
  }
  return 0;
}


MQMMarker randommarker(const MQMCrossType crosstype){
  double randnum;
  switch (crosstype) {
    case CF2:
      randnum = 4*((double)rand()/(double)RAND_MAX);
      if(randnum <= 1){
        return MAA;
      }
      if(randnum <= 3){
        return MH;
      }
      return MBB;
    break;
    case CBC:
      randnum = 2*((double)rand()/(double)RAND_MAX);   
      if(randnum <= 1){
        return MAA;
      }else{
        return MH;
      }
    break;
    case CRIL:
      randnum = 2*((double)rand()/(double)RAND_MAX);    
      if(randnum <= 1){
        return MAA;
      }else{
        return MBB;
      }
    break;
    case CUNKNOWN:
      fatal("Strange: unknown crosstype in mqm augment()");
    break;
  }
  return MMISSING;
}

int mqmaugmentfull(MQMMarkerMatrix* markers,int* nind, int* augmentednind, ivector* INDlist,
                  double neglect_unlikely, int max_totalaugment, int max_indaugment,
                  const matrix* pheno_value, const int nmark, const ivector chr, const vector mapdistance,
                  const int augment_strategy, const MQMCrossType crosstype,const int verbose){
    //Prepare for the first augmentation
    if (verbose) info("Augmentation routine");
    const int nind0 = *nind;
    const vector originalpheno = (*pheno_value)[0];
    MQMMarkerMatrix newmarkerset;
    vector new_y;                   //Because we do a phenotype matrix, we optimize by storing original the R-individual 
    ivector new_ind;                //numbers inside the trait-values, ands use new_ind etc for inside C
    ivector succes_ind;
    cvector position = relative_marker_position(nmark,chr);
    vector r = recombination_frequencies(nmark, position, mapdistance);
    if(verbose) info("Step 1: Augmentation");
    mqmaugment((*markers), (*pheno_value)[0], &newmarkerset, &new_y, &new_ind, &succes_ind, nind, augmentednind,  nmark, position, r, max_totalaugment, max_indaugment, neglect_unlikely, crosstype, verbose);
    //First round of augmentation, check if there are still individuals we need to do
    int ind_still_left=0;
    int ind_done=0;
    for(int i=0; i<nind0; i++){
      debug_trace("Individual:%d Succesfull?:%d",i,succes_ind[i]);
      if(succes_ind[i]==0){
        ind_still_left++;
      }else{
        ind_done++;
      }
    }
    if(ind_still_left && verbose) info("Step 2: Unaugmented individuals");
    if(ind_still_left && augment_strategy != 3){
      //Second round we augment dropped individuals from the first augmentation
      MQMMarkerMatrix left_markerset;
      matrix left_y_input = newmatrix(1,ind_still_left);
      vector left_y;
      ivector left_ind;
      if(verbose) info("Done with: %d/%d individuals still need to do %d",ind_done,nind0,ind_still_left);
      //Create a new markermatrix for the individuals
      MQMMarkerMatrix indleftmarkers= newMQMMarkerMatrix(nmark,ind_still_left);
      int current_leftover_ind=0;
      for(int i=0;i<nind0;i++){
        if(succes_ind[i]==0){
          debug_trace("IND %d -> %d",i,current_leftover_ind);
          left_y_input[0][current_leftover_ind] = originalpheno[i];
          for(int j=0;j<nmark;j++){
            indleftmarkers[j][current_leftover_ind] = (*markers)[j][i];
          }
          current_leftover_ind++;
        }
      }
      mqmaugment(indleftmarkers, left_y_input[0], &left_markerset, &left_y, &left_ind, &succes_ind, &current_leftover_ind, &current_leftover_ind,  nmark, position, r, max_totalaugment, max_indaugment, 1, crosstype, verbose);
      if(verbose) info("Augmentation step 2 returned most likely for %d individuals",current_leftover_ind);
      //Data augmentation done, we need to return both matrices to R
      int numimputations=1;
      if(augment_strategy==2){
        numimputations=max_indaugment;  //If we do imputation, we should generate enough to not increase likelihood for the 'unlikely genotypes'
      }
      MQMMarkerMatrix newmarkerset_all = newMQMMarkerMatrix(nmark,(*augmentednind)+numimputations*current_leftover_ind);
      vector new_y_all = newvector((*augmentednind)+numimputations*current_leftover_ind);
      ivector new_ind_all = newivector((*augmentednind)+numimputations*current_leftover_ind);;
      for(int i=0;i<(*augmentednind)+current_leftover_ind;i++){    
        int currentind;
        double currentpheno;
        if(i < (*augmentednind)){
          // Results from first augmentation step
          currentind = new_ind[i];
          currentpheno = new_y[i];
          for(int j=0;j<nmark;j++){
            newmarkerset_all[j][i] = newmarkerset[j][i];
          }
          new_ind_all[i]= currentind;
          new_y_all[i]= currentpheno;
        }else{
          // Results from second augmentation step
          currentind = ind_done+(i-(*augmentednind));
          currentpheno = left_y[(i-(*augmentednind))];
          debug_trace("Imputation of individual %d %d",currentind,numimputations);
          for(int a=0;a<numimputations;a++){
            int newindex = (*augmentednind)+a+((i-(*augmentednind))*numimputations);
            debug_trace("i=%d,s=%d,i-s=%d index=%d/%d",i,(*augmentednind),(i-(*augmentednind)),newindex,(*augmentednind)+numimputations*current_leftover_ind);
            if(augment_strategy == 2 && a > 0){
              for(int j=0;j<nmark;j++){  
                // Imputed genotype at 1 ... max_indaugment
                if(indleftmarkers[j][(i-(*augmentednind))]==MMISSING){
                  newmarkerset_all[j][newindex] = randommarker(crosstype);
                }else{
                  newmarkerset_all[j][newindex] = left_markerset[j][(i-(*augmentednind))];
                }
              }        
            }else{
              for(int j=0;j<nmark;j++){  
                // Most likely genotype at 0  
                newmarkerset_all[j][newindex] = left_markerset[j][(i-(*augmentednind))];
              }
            }
            new_ind_all[newindex]= currentind;
            new_y_all[newindex]= currentpheno;
            debug_trace("Individual: %d OriginalID:%f Variant:%d",currentind,currentpheno,a);
          }
        }
      }
      //Everything is added together so lets set out return pointers
      (*pheno_value)[0] = new_y_all;
      (*INDlist) = new_ind_all;
      (*markers) = newmarkerset_all;
      (*augmentednind)=(*augmentednind)+(numimputations*current_leftover_ind);
      (*nind)= (*nind)+(current_leftover_ind);
      debug_trace("nind:%d,naugmented:%d",(*nind)+(current_leftover_ind),(*augmentednind)+(current_leftover_ind));
    }else{
      if(ind_still_left && augment_strategy == 3){
        if (verbose) info("Dropping %d augment_strategy individuals from further analysis",ind_still_left);
      }
      //We augmented all individuals in the first go so lets use those
      (*pheno_value)[0] = new_y;
      (*INDlist) = new_ind;
      (*markers) = newmarkerset;
    }
    if(verbose) info("Done with augmentation");
    return 1;
}

int mqmaugment(const MQMMarkerMatrix marker, const vector y, 
               MQMMarkerMatrix* augmarker, vector *augy, 
               ivector* augind, ivector* sucind, int *Nind, int *Naug, const int Nmark, 
               const cvector position, vector r, const int maxNaug, 
               const int imaxNaug, const double minprob, 
               const MQMCrossType crosstype, const int verbose) 
{
  int retvalue = 1;     //[Danny] Assume everything will go right, (it never returned a 1 OK, initialization to 0 and return
  int jj;
  const int nind0 = *Nind;              //Original number of individuals
  (*Naug) = maxNaug;     // sets and returns the maximum size of augmented dataset
  // new variables sized to maxNaug:
  MQMMarkerMatrix newmarker;
  vector newy;
  MQMMarkerVector imarker;
  ivector newind;
  ivector succesind;
  
  double minprobratio = (1.0f/minprob);
  if(minprob!=1){
    minprobratio += 0.00001;
  }
  newmarker = newMQMMarkerMatrix(Nmark+1, maxNaug);  // augmented marker matrix
  newy      = newvector(maxNaug);            // phenotypes
  newind    = newivector(maxNaug);           // individuals index
  succesind = newivector(nind0);              // Tracks if the augmentation is a succes
  imarker   = newMQMMarkerVector(Nmark);             

  int iaug     = 0;     // iaug keeps track of current augmented individual
  double prob0, prob1, prob2, sumprob,
  prob0left, prob1left, prob2left,
  prob0right=0.0, prob1right=0.0, prob2right = 0.0f;
  vector newprob = newvector(maxNaug);
  vector newprobmax = newvector(maxNaug);
  if (verbose) info("Crosstype determined by the algorithm:%c:", crosstype);
  if (verbose) info("Augmentation parameters: Maximum augmentation=%d, Maximum augmentation per individual=%d, Minprob=%f", maxNaug, imaxNaug, minprob);
  // ---- foreach individual create one in the newmarker matrix
 
  int newNind = nind0;                  //Number of unique individuals
  int previaug = 0;                     // previous index in newmarkers
  for (int i=0; i<nind0; i++) {
    //Loop through individuals
    succesind[i] = 1;                   //Assume we succeed in augmentation
    #ifndef STANDALONE
      //R_ProcessEvents(); /*  Try not to crash windows */
      R_FlushConsole();
    #endif
    const int dropped = nind0-newNind;  //How many are dropped
    const int iidx = i - dropped;       //Individuals I's new individual number based on dropped individuals
    newind[iaug]   = iidx;              // iidx corrects for dropped individuals
    newy[iaug]     = y[i];              // cvariance (phenotype)
    newprob[iaug]  = 1.0;               //prop
    double probmax = 1.0;               //current maximum probability

    for (int j=0; j<Nmark; j++){ 
      newmarker[j][iaug]=marker[j][i];    // copy markers into newmarkers for the new indidivudal under investigation
    }
    for (int j=0; j<Nmark; j++) {
      //Loop through markers:
      const int maxiaug = iaug;          // fixate maxiaug
      if ((maxiaug-previaug)<=imaxNaug)  // within bounds for individual?
        for (int ii=previaug; ii<=maxiaug; ii++) {

	  R_CheckUserInterrupt(); /* check for ^C */

          debug_trace("i=%d ii=%d iidx=%d maxiaug=%d previaug=%d,imaxNaug=%d\n",i,ii,iidx,maxiaug,previaug,imaxNaug);
          // ---- walk from previous augmented to current augmented genotype
          //WE HAVE 3 SPECIAL CASES: (1) NOTAA, (2) NOTBB and (3)UNKNOWN, and the std case of a next known marker
          if (newmarker[j][ii]==MNOTAA) {
            //NOTAA augment data to contain AB and BB
            for (jj=0; jj<Nmark; jj++) imarker[jj] = newmarker[jj][ii];

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob1left= start_prob(crosstype, MH);
              prob2left= start_prob(crosstype, MBB);
            } else {
              prob1left= left_prob(r[j-1],newmarker[j-1][ii],MH,crosstype);      //prob1left= prob(newmarker, r, ii, j-1, MH, crosstype, 0);
              prob2left= left_prob(r[j-1],newmarker[j-1][ii],MBB,crosstype);     //prob2left= prob(newmarker, r, ii, j-1, MBB, crosstype, 0);
            }
            switch (crosstype) {
              case CF2:
                prob1right= right_prob_F2(MH, j, imarker, r, position);          //prob1right= probright(MH, j, imarker, r, position, crosstype);
                prob2right= right_prob_F2(MBB, j, imarker, r, position);         //prob2right= probright(MBB, j, imarker, r, position, crosstype);
              break;
              case CBC:
                prob1right= right_prob_BC(MH, j, imarker, r, position);
                prob2right= right_prob_BC(MBB, j, imarker, r, position);                
              break;
              case CRIL:
                prob1right= right_prob_RIL(MH, j, imarker, r, position);
                prob2right= right_prob_RIL(MBB, j, imarker, r, position);                
              break;
              case CUNKNOWN:
                fatal("Strange: unknown crosstype in mqm augment()");
              break;
            }
            prob1= prob1left*prob1right;
            prob2= prob2left*prob2right;

            if (ii==previaug) probmax = (prob2>prob1 ? newprob[ii]*prob2 : newprob[ii]*prob1);
            if (prob1>prob2) {
              if (probmax/(newprob[ii]*prob2)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MBB;
                newprob[iaug]= newprob[ii]*prob2left;
                newprobmax[iaug]= newprob[iaug]*prob2right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MH;
              newprobmax[ii]= newprob[ii]*prob1;
              newprob[ii]= newprob[ii]*prob1left;
            } else {
              if (probmax/(newprob[ii]*prob1)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MH;
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MBB;
              newprobmax[ii]= newprob[ii]*prob2;
              newprob[ii]*= prob2left;
            }
            probmax = (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
          } else if (newmarker[j][ii]==MNOTBB) {
            //NOTBB: augment data can contain MH and MAA 
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob0left= start_prob(crosstype, MAA);
              prob1left= start_prob(crosstype, MH);
            } else {
              prob0left= left_prob(r[j-1],newmarker[j-1][ii],MAA,crosstype);  //prob0left= prob(newmarker, r, ii, j-1, MAA, crosstype, 0);
              prob1left= left_prob(r[j-1],newmarker[j-1][ii],MH,crosstype);   //prob1left= prob(newmarker, r, ii, j-1, MH, crosstype, 0);
            }
            switch (crosstype) {
              case CF2:
                prob0right= right_prob_F2(MAA, j, imarker, r, position);      //prob0right= probright(MAA, j, imarker, r, position, crosstype);
                prob1right= right_prob_F2(MH, j, imarker, r, position);       //prob1right= probright(MH, j, imarker, r, position, crosstype);
              break;
              case CBC:
                prob0right= right_prob_BC(MAA, j, imarker, r, position);
                prob1right= right_prob_BC(MH, j, imarker, r, position);               
              break;
              case CRIL:
                prob0right= right_prob_RIL(MAA, j, imarker, r, position);
                prob1right= right_prob_RIL(MH, j, imarker, r, position);              
              break;
              case CUNKNOWN:
                fatal("Strange: unknown crosstype in mqm augment()");
              break;
            }
            prob0= prob0left*prob0right;
            prob1= prob1left*prob1right;

            if (ii==previaug) probmax= (prob0>prob1 ? newprob[ii]*prob0 : newprob[ii]*prob1);
            if (prob1>prob0) {
              if (probmax/(newprob[ii]*prob0)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MAA;
                newprob[iaug]= newprob[ii]*prob0left;
                newprobmax[iaug]= newprob[iaug]*prob0right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MH;
              newprobmax[ii]= newprob[ii]*prob1;
              newprob[ii]*= prob1left;
            } else {
              if (probmax/(newprob[ii]*prob1)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MH;
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MAA;
              newprobmax[ii]= newprob[ii]*prob0;
              newprob[ii]*= prob0left;
            }
            probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
          } else if (newmarker[j][ii]==MMISSING) {
            //UNKNOWN: augment data to contain AB, AA and BB
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob0left= start_prob(crosstype, MAA);
              prob1left= start_prob(crosstype, MH);
              prob2left= start_prob(crosstype, MBB);
            } else {
              prob0left= left_prob(r[j-1],newmarker[j-1][ii],MAA,crosstype);  //prob0left= prob(newmarker, r, ii, j-1, MAA, crosstype, 0);
              prob1left= left_prob(r[j-1],newmarker[j-1][ii],MH,crosstype);   //prob1left= prob(newmarker, r, ii, j-1, MH, crosstype, 0);
              prob2left= left_prob(r[j-1],newmarker[j-1][ii],MBB,crosstype);  //prob2left= prob(newmarker, r, ii, j-1, MBB, crosstype, 0);
            }
            switch (crosstype) {
              case CF2:
                prob0right= right_prob_F2(MAA, j, imarker, r, position); //prob0right= probright(MAA, j, imarker, r, position, crosstype);
                prob1right= right_prob_F2(MH, j, imarker, r, position);  //prob1right= probright(MH, j, imarker, r, position, crosstype);
                prob2right= right_prob_F2(MBB, j, imarker, r, position); //prob2right= probright(MBB, j, imarker, r, position, crosstype);
              break;
              case CBC:
                prob0right= right_prob_BC(MAA, j, imarker, r, position);
                prob1right= right_prob_BC(MH, j, imarker, r, position);
                prob2right= 0.0;              
              break;
              case CRIL:
                prob0right= right_prob_RIL(MAA, j, imarker, r, position);
                prob1right= 0.0;
                prob2right= right_prob_RIL(MBB, j, imarker, r, position);              
              break;
              case CUNKNOWN:
                fatal("Strange: unknown crosstype in mqm augment()");
              break;
            }            
            prob0= prob0left*prob0right;
            prob1= prob1left*prob1right;
            prob2= prob2left*prob2right;
            if (ii==previaug) {
              if ((prob2>prob1)&&(prob2>prob0)) probmax= newprob[ii]*prob2;
              else if ((prob1>prob0)&&(prob1>prob2)) probmax= newprob[ii]*prob1;
              else probmax= newprob[ii]*prob0;
            }
            if ((prob2>prob1)&&(prob2>prob0)) {
              if (probmax/(newprob[ii]*prob1)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MH;
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              if (probmax/(newprob[ii]*prob0)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MAA;
                newprob[iaug]= newprob[ii]*prob0left;
                newprobmax[iaug]= newprob[iaug]*prob0right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MBB;
              newprobmax[ii]= newprob[ii]*prob2;
              newprob[ii]*= prob2left;

            } else if ((prob1>prob2)&&(prob1>prob0)) {
              if (probmax/(newprob[ii]*prob2)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MBB;
                newprob[iaug]= newprob[ii]*prob2left;
                newprobmax[iaug]= newprob[iaug]*prob2right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              if (probmax/(newprob[ii]*prob0)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MAA;
                newprob[iaug]= newprob[ii]*prob0left;
                newprobmax[iaug]= newprob[iaug]*prob0right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MH;
              newprobmax[ii]= newprob[ii]*prob1;
              newprob[ii]*= prob1left;
            } else {
              if (probmax/(newprob[ii]*prob1)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MH;
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              if (probmax/(newprob[ii]*prob2)<minprobratio) {
                if (++iaug >= maxNaug) goto bailout;
                newmarker[j][iaug]= MBB;
                newprob[iaug]= newprob[ii]*prob2left;
                newprobmax[iaug]= newprob[iaug]*prob2right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=iidx;
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= MAA;
              newprobmax[ii]= newprob[ii]*prob0;
              newprob[ii]*= prob0left;
            }
            probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
          } else {
            //STD case we know what the next marker is nou use probleft to estimate the likelihood of the current location
            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob0left= start_prob(crosstype, newmarker[j][ii]);
            } else {
              prob0left= left_prob(r[j-1],newmarker[j-1][ii],newmarker[j][ii],crosstype); //prob0left= prob(newmarker, r, ii, j-1, newmarker[j][ii], crosstype, 0);
            }
            newprob[ii]*= prob0left;
          }

          if (iaug+3>maxNaug) {
            Rprintf("ERROR: augmentation (this code should not be reached)\n");  
            goto bailout;
          }
        }
      if ((iaug-previaug+1)>imaxNaug) {
        newNind-= 1;
        iaug= previaug-1;
        succesind[i]=0;
        //for(int x=previaug;x<previaug+imaxNaug;x++){
        //  info("Individual: %d, variant: %d, prob: %f",i,x,newprob[x]);
        //}
        if (verbose) info("Individual %d moved to second augmentation round", i);
      }
      sumprob= 0.0;
      for (int ii=previaug; ii<=iaug; ii++) sumprob+= newprob[ii];
      for (int ii=previaug; ii<=iaug; ii++) newprob[ii]/= sumprob;
    }
    if (++iaug >= maxNaug) goto bailout;
    previaug=iaug;
  }
  *Naug = iaug;
  *Nind = newNind;
  *augmarker = newMQMMarkerMatrix(Nmark, *Naug);
  *augy = newvector(*Naug);
  *augind = newivector(*Naug);
  *sucind = newivector(nind0);
  for (int i=0; i<nind0; i++) {
    (*sucind)[i] = succesind[i];
  }
  for (int i=0; i<(*Naug); i++) {
    (*augy)[i]= newy[i];
    (*augind)[i]= newind[i];
    for (int j=0; j<Nmark; j++) (*augmarker)[j][i]= newmarker[j][i];
  }
  goto cleanup;
bailout:
  Rprintf("ERROR: Dataset too large after augmentation\n");
  if (verbose) Rprintf("INFO: Recall procedure with larger value for augmentation parameters or increase the parameter minprob\n");
  retvalue = 0;
cleanup:
  Free(newy);
  Free(newmarker);
  Free(newind);
  Free(newprob);
  Free(newprobmax);
  Free(imarker);
  return retvalue;
}

/*
 * The R interfact to data augmentation
 */

void R_mqmaugment(int *geno, double *dist, double *pheno, int *auggeno, 
               double *augPheno, int *augIND, int *Nind, int *Naug, int *Nmark,
               int *Npheno, int *maxind, int *maxiaug, double *minprob, int
               *chromo, int *rqtlcrosstypep, int *augment_strategy, int *verbosep) {
  int **Geno;
  double **Pheno;
  double **Dist;
  int **NEW;                      //Holds the output for the augmentdata function
  int **Chromo;
  double **NEWPheno;              //New phenotype vector
  int **NEWIND;                   //New list of individuals 
  const int nind0 = *Nind;        //Individuals we start with
  const int verbose = *verbosep;
  const RqtlCrossType rqtlcrosstype = (RqtlCrossType) *rqtlcrosstypep;

  if(verbose) info("Starting C-part of the data augmentation routine");
  ivector new_ind;
  vector mapdistance;
  cvector position;
  MQMMarkerMatrix markers, new_markers;
  ivector chr;

  markers= newMQMMarkerMatrix(*Nmark, nind0);
  new_markers= newMQMMarkerMatrix(*Nmark, *maxind);
  mapdistance = newvector(*Nmark);
  chr= newivector(*Nmark);

  //Reorganise the pointers into arrays, Singletons are just cast into the function
  reorg_geno(nind0, *Nmark, geno, &Geno);
  reorg_int(*Nmark, 1, chromo, &Chromo);
  reorg_pheno(nind0, *Npheno, pheno, &Pheno);
  reorg_pheno(*Nmark, 1, dist, &Dist);

  reorg_int(*maxind, *Nmark, auggeno, &NEW);
  reorg_int((*maxiaug)*nind0, 1, augIND, &NEWIND);
  reorg_pheno((*maxiaug)*nind0, 1, augPheno, &NEWPheno);

  MQMCrossType crosstype = determine_MQMCross(*Nmark, *Nind, (const int **)Geno, rqtlcrosstype);
  //Change all the markers from R/qtl format to MQM internal
  change_coding(Nmark, Nind, Geno, markers, crosstype);

  if(verbose) info("Filling the chromosome matrix");
  for (int i=0; i<(*Nmark); i++) {
    //Set some general information structures per marker
    mapdistance[i]=POSITIONUNKNOWN;
    mapdistance[i]=Dist[0][i];
    chr[i] = Chromo[0][i];
  }
  //Calculate positions of markers and Recombinant frequencies
  position = relative_marker_position(*Nmark,chr);
  vector r = recombination_frequencies(*Nmark, position, mapdistance);
  //ivector succes_ind;
  /*
  if (mqmaugment(markers, Pheno[(*Npheno-1)], &new_markers, &new_y, &new_ind, &succes_ind, Nind, Naug, *Nmark, position, r, *maxind, *maxiaug, *minprob, crosstype, verbose)==1) {
  
  */
  if(mqmaugmentfull(&markers,Nind,Naug,&new_ind,*minprob, *maxind, *maxiaug,&Pheno,*Nmark,chr,mapdistance,*augment_strategy,crosstype,verbose)){
    //Data augmentation finished succesfully
    //Push it back into RQTL format
    for (int i=0; i<(*Nmark); i++) {
      for (int j=0; j<(*Naug); j++) {
        //info("Phenotype after return: %f",NEWPheno[0][j]);
        NEWPheno[0][j] = Pheno[0][j];
        NEWIND[0][j] = new_ind[j];
        NEW[i][j] = 9;
        if (markers[i][j] == MAA) {
          NEW[i][j] = 1;
        }
        if (markers[i][j] == MH) {
          NEW[i][j] = 2;
        }
        if (markers[i][j] == MBB) {  // [karl:] this might need to be changed for RIL
          crosstype==CRIL ? NEW[i][j]=2 : NEW[i][j] = 3;  //[Danny:] This should solve it 
        }
        if (markers[i][j] == MNOTAA) {
          NEW[i][j] = 5;
        }
        if (markers[i][j] == MNOTBB) {
          NEW[i][j] = 4;
        }
      }
    }
    //delMQMMarkerMatrix(new_markers,*Nmark);
    //delMQMMarkerMatrix(markers,*Nmark);
    Free(mapdistance);
    Free(position);
    Free(r);
    Free(chr);
    if (verbose) {
      Rprintf("# Unique individuals before augmentation:%d\n", nind0);
      Rprintf("# Unique selected individuals:%d\n", *Nind);
      Rprintf("# Marker p individual:%d\n", *Nmark);
      Rprintf("# Individuals after augmentation:%d\n", *Naug);
      info("Data augmentation succesfull");
    }
  } else {
    //Unsuccessfull data augmentation exit
    info("This code should not be reached, data corruption could have occured. Please re-run this analysis.");
    *Naug = nind0;
    for (int i=0; i<(*Nmark); i++) {
      for (int j=0; j<(*Naug); j++) {
        NEWPheno[0][j] = Pheno[0][j];
        NEW[i][j] = 9;
        if (markers[i][j] == MAA) {
          NEW[i][j] = 1;
        }
        if (markers[i][j] == MH) {
          NEW[i][j] = 2;
        }
        if (markers[i][j] == MBB) { // [karl:] this might need to be changed for RIL
          crosstype==CRIL ? NEW[i][j]=2 : NEW[i][j] = 3;  //[Danny:] This should solve it 
        }
        if (markers[i][j] == MNOTAA) {
          NEW[i][j] = 5;
        }
        if (markers[i][j] == MNOTBB) {
          NEW[i][j] = 4;
        }
      }
    }
    delMQMMarkerMatrix(new_markers,*Nmark);
    delMQMMarkerMatrix(markers,*Nmark);
    Free(mapdistance);
    Free(position);
    Free(r);
    Free(chr);
    fatal("Data augmentation failed");
  }
  return;
}


