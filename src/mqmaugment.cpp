/**********************************************************************
 *
 * mqmaugment.cpp
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 * last modified Apr, 2009
 * first written Feb, 2009
 *
 * first written <2000 (unknown)
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
 * Data augmentation routines
 *
 **********************************************************************/

#include "mqm.h"

/*
 * Augment/expand the dataset by adding additional (likely) genotypes.
 * genotypic marker information is incomplete when
 *           - marker scores are missing or dominant, or
 *           - markers are non-segregating but you want to assume a
 *             segregating QTL on top of this marker (~very closely linked),
 *             and all likely possible configurations are generated 
 *
 * Inputs are number of markers Nmark, the marker matrix, position vector,
 * recombinations r.  The neglect parameter drops individuals from the dataset
 * (the value should be between 1..n).
 *
 * The neglect parameter drops genotypes. E.g. for neglect=100  eliminate
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
 * FIXME: increasing the buffers for augmentation can automatic
 */

int augmentdata(const cmatrix marker, const vector y, cmatrix* augmarker, vector *augy, 
            ivector* augind, int *Nind, int *Naug, const int Nmark, 
            const cvector position, vector r, const int maxNaug, const int imaxNaug, 
            const double neglect, const MQMCrossType crosstype, const int verbose) {
  int retvalue = 0;
  int jj;
  (*Naug) = maxNaug;     // sets and returns the maximum size of augmented dataset
  // new variables sized to maxNaug:
  cmatrix newmarker;
  vector newy;
  cvector imarker;
  ivector newind;

  newmarker = newcmatrix(Nmark+1, maxNaug);  // augmented marker matrix
  newy      = newvector(maxNaug);            // phenotypes
  newind    = newivector(maxNaug);           // individuals index
  imarker   = newcvector(Nmark);             

  int iaug     = 0;     // iaug keeps track of current augmented individual
  // int maxiaug  = 0;     // highest reached(?)
  // probabilities:
  double prob0, prob1, prob2, sumprob,
  prob0left, prob1left, prob2left,
  prob0right, prob1right, prob2right;
  vector newprob = newvector(maxNaug);
  vector newprobmax = newvector(maxNaug);
  if (verbose) Rprintf(("Crosstype determined by the algorithm:%c:", crosstype);
  if (verbose) Rprintf(("Augmentation parameters: Maximum augmentation=%d, Maximum augmentation per individual=%d, Neglect=%f", maxNaug, imaxNaug, neglect);
  // ---- foreach individual create one in the newmarker matrix
  const int nind0 = *Nind;
  int newNind = nind0;
  int previaug = 0;                    // previous iaug
  for (int i=0; i<nind0; i++) {
    // ---- for every individual:
    const int dropped = nind0-newNind;
    const int iidx = i - dropped;
    newind[iaug]  = iidx;              // iidx corrects for dropped individuals
    newy[iaug]    = y[i];              // cvariance (phenotype)
    newprob[iaug] = 1.0;
    double probmax = 1.0;
    for (int j=0; j<Nmark; j++) 
      newmarker[j][iaug]=marker[j][i]; // align new markers with markers (current iaug)
    for (int j=0; j<Nmark; j++) {
      // ---- for every marker:
      const int maxiaug = iaug;          // fixate maxiaug
      if ((maxiaug-previaug)<=imaxNaug)  // within bounds for individual?
        for (int ii=previaug; ii<=maxiaug; ii++) {
          // ---- walk from previous augmented to current augmented genotype
          if (newmarker[j][ii]==MNOTAA) {
            // augment data to contain AB and BB
            for (jj=0; jj<Nmark; jj++) imarker[jj] = newmarker[jj][ii];

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob1left= start_prob(crosstype, MH);
              prob2left= start_prob(crosstype, MBB);
            } else {
              prob1left= prob(newmarker, r, ii, j-1, MH, crosstype, 0);
              prob2left= prob(newmarker, r, ii, j-1, MBB, crosstype, 0);
            }

            prob1right= probright(MH, j, imarker, r, position, crosstype);
            prob2right= probright(MBB, j, imarker, r, position, crosstype);
            prob1= prob1left*prob1right;
            prob2= prob2left*prob2right;

            if (ii==previaug) probmax = (prob2>prob1 ? newprob[ii]*prob2 : newprob[ii]*prob1);
            if (prob1>prob2) {
              if (probmax/(newprob[ii]*prob2)<neglect) {
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
              if (probmax/(newprob[ii]*prob1)<neglect) {
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
            // augment data to contain AB and AA 
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob0left= start_prob(crosstype, MAA);
              prob1left= start_prob(crosstype, MH);
            } else {
              prob0left= prob(newmarker, r, ii, j-1, MAA, crosstype, 0);
              prob1left= prob(newmarker, r, ii, j-1, MH, crosstype, 0);
            }

            prob0right= probright(MAA, j, imarker, r, position, crosstype);
            prob1right= probright(MH, j, imarker, r, position, crosstype);
            prob0= prob0left*prob0right;
            prob1= prob1left*prob1right;

            if (ii==previaug) probmax= (prob0>prob1 ? newprob[ii]*prob0 : newprob[ii]*prob1);
            if (prob1>prob0) {
              if (probmax/(newprob[ii]*prob0)<neglect) {
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
              if (probmax/(newprob[ii]*prob1)<neglect) {
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
            // augment data to contain AB, AA and BB
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob0left= start_prob(crosstype, MAA);
              prob1left= start_prob(crosstype, MH);
              prob2left= start_prob(crosstype, MBB);
            } else {
              prob0left= prob(newmarker, r, ii, j-1, MAA, crosstype, 0);
              prob1left= prob(newmarker, r, ii, j-1, MH, crosstype, 0);
              prob2left= prob(newmarker, r, ii, j-1, MBB, crosstype, 0);
            }

            prob0right= probright(MAA, j, imarker, r, position, crosstype);
            prob1right= probright(MH, j, imarker, r, position, crosstype);
            prob2right= probright(MBB, j, imarker, r, position, crosstype);
            prob0= prob0left*prob0right;
            prob1= prob1left*prob1right;
            prob2= prob2left*prob2right;
            if (ii==previaug) {
              if ((prob2>prob1)&&(prob2>prob0)) probmax= newprob[ii]*prob2;
              else if ((prob1>prob0)&&(prob1>prob2)) probmax= newprob[ii]*prob1;
              else probmax= newprob[ii]*prob0;
            }
            if ((prob2>prob1)&&(prob2>prob0)) {
              if (probmax/(newprob[ii]*prob1)<neglect) {
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
              if (probmax/(newprob[ii]*prob0)<neglect) {
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
              if (probmax/(newprob[ii]*prob2)<neglect) {
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
              if (probmax/(newprob[ii]*prob0)<neglect) {
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
              if (probmax/(newprob[ii]*prob1)<neglect) {
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
              if (probmax/(newprob[ii]*prob2)<neglect) {
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
          } else { // newmarker[j][ii] is observed

            if ((position[j]==MLEFT||position[j]==MUNLINKED)) {
              prob0left= start_prob(crosstype, newmarker[j][ii]);
            } else {
              prob0left= prob(newmarker, r, ii, j-1, newmarker[j][ii], crosstype, 0);
            }

            newprob[ii]*= prob0left;

          }

          if (iaug+3>maxNaug) {
            Rprintf("ERROR: augmentation (should not be reached)\n");  
            goto bailout;
          }
        }
      if ((iaug-previaug+1)>imaxNaug) {
        newNind-= 1;
        iaug= previaug-1;
        if (verbose) Rprintf("INFO: Individual %d has been dropped\n", i);
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
  *augmarker = newcmatrix(Nmark, *Naug);
  *augy = newvector(*Naug);
  *augind = newivector(*Naug);
  for (int i=0; i<(*Naug); i++) {
    (*augy)[i]= newy[i];
    (*augind)[i]= newind[i];
    for (int j=0; j<Nmark; j++) (*augmarker)[j][i]= newmarker[j][i];
  }
  goto cleanup;
bailout:
  Rprintf("ERROR: Dataset too large after augmentation\n");
  if (verbose) Rprintf("INFO: Recall procedure with larger value for augmentation parameters or lower the parameter neglect\n");
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

void R_augmentdata(int *geno, double *dist, double *pheno, int *auggeno, 
               double *augPheno, int *augIND, int *Nind, int *Naug, int *Nmark,
               int *Npheno, int *maxind, int *maxiaug, double *neglect, int
               *chromo, int *rqtlcrosstypep, int *verbosep) {
  int **Geno;
  double **Pheno;
  double **Dist;
  int **NEW;
  int **Chromo;
  double **NEWPheno;
  int **NEWIND;
  const int nind0 = *Nind;
  int prior = nind0;
  const int verbose = *verbosep;
  const RqtlCrossType rqtlcrosstype = (RqtlCrossType) *rqtlcrosstypep;

  info("Starting C-part of the data augmentation routine");
  ivector new_ind;
  vector new_y, mapdistance;
  cvector position;
  cmatrix markers, new_markers;
  ivector chr;

  markers= newcmatrix(*Nmark, nind0);
  new_markers= newcmatrix(*Nmark, *maxind);
  mapdistance = newvector(*Nmark);
  chr= newivector(*Nmark);

  //Reorganise the pointers into arrays, Singletons are just cast into the function
  reorg_geno(nind0, *Nmark, geno, &Geno);
  reorg_int(*Nmark, 1, chromo, &Chromo);
  reorg_pheno(nind0, *Npheno, pheno, &Pheno);
  reorg_pheno(*Nmark, 1, dist, &Dist);

  reorg_int(*maxind, *Nmark, auggeno, &NEW);
  reorg_int((*maxiaug)*nind0, 1, augIND, &NEWIND);
  reorg_pheno(*maxind, 1, augPheno, &NEWPheno);

  MQMCrossType crosstype = determine_MQMCross(*Nmark, *Nind, (const int **)Geno, rqtlcrosstype);
  //Change all the markers from R/qtl format to MQM internal
  change_coding(Nmark, Nind, Geno, markers, crosstype);

  info("Filling the chromosome matrix");

  for (int i=0; i<(*Nmark); i++) {
    //Set some general information structures per marker
    mapdistance[i]=999.0;
    mapdistance[i]=Dist[0][i];
    chr[i] = Chromo[0][i];
  }

  position = locate_markers(*Nmark,chr);
  vector r = recombination_frequencies(*Nmark, position, mapdistance);

  if (augmentdata(markers, Pheno[(*Npheno-1)], &new_markers, &new_y, &new_ind, Nind, Naug, *Nmark, position, r, *maxind, *maxiaug, *neglect, crosstype, verbose)==1) {
    //Data augmentation finished succesfully
    //Push it back into RQTL format
    for (int i=0; i<(*Nmark); i++) {
      for (int j=0; j<(*Naug); j++) {
        NEWPheno[0][j] = new_y[j];
        NEWIND[0][j] = new_ind[j];
        NEW[i][j] = 9;
        if (new_markers[i][j] == MAA) {
          NEW[i][j] = 1;
        }
        if (new_markers[i][j] == MH) {
          NEW[i][j] = 2;
        }
        if (new_markers[i][j] == MBB) {  // [karl:] this might need to be changed for RIL
          NEW[i][j] = 3;
        }
        if (new_markers[i][j] == MNOTAA) {
          NEW[i][j] = 5;
        }
        if (new_markers[i][j] == MNOTBB) {
          NEW[i][j] = 4;
        }
      }
    }
    delcmatrix(new_markers,*Nmark);
    delcmatrix(markers,*Nmark);
    Free(mapdistance);
    Free(position);
    Free(r);
    Free(chr);
    if (verbose) {
      Rprintf("# Unique individuals before augmentation:%d\n", prior);
      Rprintf("# Unique selected individuals:%d\n", nind0);
      Rprintf("# Marker p individual:%d\n", *Nmark);
      Rprintf("# Individuals after augmentation:%d\n", *Naug);
      info("Data augmentation succesfull");
    }
  } else {
    //Unsuccessfull data augmentation exit
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
          NEW[i][j] = 3;
        }
        if (markers[i][j] == MNOTAA) {
          NEW[i][j] = 5;
        }
        if (markers[i][j] == MNOTBB) {
          NEW[i][j] = 4;
        }
      }
    }
    delcmatrix(new_markers,*Nmark);
    delcmatrix(markers,*Nmark);
    Free(mapdistance);
    Free(position);
    Free(r);
    Free(chr);
    fatal("Data augmentation failed");
  }
  return;
}


