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
 * Augment/expand the dataset by adding additional marker positions. Inputs are
 * number of markers Nmark, the marker matrix, position vector, recombinations
 * r.  The neglect parameter drops individuals from the dataset (the value
 * should be between 1..n).
 *
 * A new markerset is created and returned in augmarker, likewise the
 * phenotypes are returned in augy and individuals in augind. Augmentation
 * halts when the number of individuals maxNaug is reached. Markers are
 * expanded up to imaxNaug.
 *
 */

int augdata(cmatrix marker, vector y, cmatrix* augmarker, vector *augy, 
            ivector* augind, int *Nind, int *Naug, int Nmark, cvector position,
            vector r, int maxNaug, int imaxNaug, double neglect, char
            crosstype, int verbose) {
  int jj;
  int newNind=(*Nind);
  (*Naug)= maxNaug;     // sets and returns the maximum size of augmented dataset
  // new variables sized to maxNaug:
  cmatrix newmarker;
  vector newy;
  cvector imarker;
  ivector newind;

  newmarker= newcmatrix(Nmark+1, *Naug);
  newy= newvector(*Naug);
  newind= newivector(*Naug);
  imarker= newcvector(Nmark);

  int iaug=0;      // iaug keeps track of current augmented individual
  int maxiaug=0;   // highest reached(?)
  int saveiaug=0;  // previous iaug
  double prob0, prob1, prob2, sumprob,
  prob0left, prob1left, prob2left,
  prob0right, prob1right, prob2right;
  double probmax;
  vector newprob, newprobmax;
  newprob= newvector(*Naug);
  newprobmax= newvector(*Naug);
  if (verbose) {
    Rprintf("INFO: Crosstype determined by the algorithm:%c:\n", crosstype);
    Rprintf("INFO: Augmentation parameters: Maximum augmentation=%d, Maximum augmentation per individual=%d, Neglect=%f\n", maxNaug, imaxNaug, neglect);
  }
  // ---- foreach individual create one in the newmarker matrix
  for (int i=0; i<(*Nind); i++) {
    newind[iaug]=i-((*Nind)-newNind);  // index of individuals
    newy[iaug]= y[i];               // cvariance
    newprob[iaug]= 1.0;
    probmax= 1.0;
    for (int j=0; j<Nmark; j++) newmarker[j][iaug]=marker[j][i];
    for (int j=0; j<Nmark; j++) {
      maxiaug=iaug;
      if ((maxiaug-saveiaug)<=imaxNaug)  // within bounds for individual?
        for (int ii=saveiaug; ii<=maxiaug; ii++) {
          if (newmarker[j][ii]=='3') {
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]=='L'||position[j]=='U')) {
              prob1left= start_prob(crosstype, '1');
              prob2left= start_prob(crosstype, '2');
            } else {
              prob1left= prob(newmarker, r, ii, j-1, '1', crosstype, 1, 0, 0);
              prob2left= prob(newmarker, r, ii, j-1, '2', crosstype, 1, 0, 0);
            }

            prob1right= probright('1', j, imarker, r, position, crosstype);
            prob2right= probright('2', j, imarker, r, position, crosstype);
            prob1= prob1left*prob1right;
            prob2= prob2left*prob2right;

            if (ii==saveiaug) probmax= (prob2>prob1 ? newprob[ii]*prob2 : newprob[ii]*prob1);
            if (prob1>prob2) {
              if (probmax/(newprob[ii]*prob2)<neglect) {
                iaug++;
                newmarker[j][iaug]= '2';
                newprob[iaug]= newprob[ii]*prob2left;
                newprobmax[iaug]= newprob[iaug]*prob2right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '1';
              newprobmax[ii]= newprob[ii]*prob1;
              newprob[ii]= newprob[ii]*prob1left;
            } else {
              if (probmax/(newprob[ii]*prob1)<neglect) {
                iaug++;
                newmarker[j][iaug]= '1';
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '2';
              newprobmax[ii]= newprob[ii]*prob2;
              newprob[ii]*= prob2left;
            }
            probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
          } else if (newmarker[j][ii]=='4') {
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]=='L'||position[j]=='U')) {
              prob0left= start_prob(crosstype, '0');
              prob1left= start_prob(crosstype, '1');
            } else {
              prob0left= prob(newmarker, r, ii, j-1, '0', crosstype, 1, 0, 0);
              prob1left= prob(newmarker, r, ii, j-1, '1', crosstype, 1, 0, 0);
            }

            prob0right= probright('0', j, imarker, r, position, crosstype);
            prob1right= probright('1', j, imarker, r, position, crosstype);
            prob0= prob0left*prob0right;
            prob1= prob1left*prob1right;

            if (ii==saveiaug) probmax= (prob0>prob1 ? newprob[ii]*prob0 : newprob[ii]*prob1);
            if (prob1>prob0) {
              if (probmax/(newprob[ii]*prob0)<neglect) {
                iaug++;
                newmarker[j][iaug]= '0';
                newprob[iaug]= newprob[ii]*prob0left;
                newprobmax[iaug]= newprob[iaug]*prob0right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '1';
              newprobmax[ii]= newprob[ii]*prob1;
              newprob[ii]*= prob1left;
            } else {
              if (probmax/(newprob[ii]*prob1)<neglect) {
                iaug++;
                newmarker[j][iaug]= '1';
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '0';
              newprobmax[ii]= newprob[ii]*prob0;
              newprob[ii]*= prob0left;
            }
            probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
          } else if (newmarker[j][ii]=='9') {
            for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];

            if ((position[j]=='L'||position[j]=='U')) {
              prob0left= start_prob(crosstype, '0');
              prob1left= start_prob(crosstype, '1');
              prob2left= start_prob(crosstype, '2');
            } else {
              prob0left= prob(newmarker, r, ii, j-1, '0', crosstype, 1, 0, 0);
              prob1left= prob(newmarker, r, ii, j-1, '1', crosstype, 1, 0, 0);
              prob2left= prob(newmarker, r, ii, j-1, '2', crosstype, 1, 0, 0);
            }

            prob0right= probright('0', j, imarker, r, position, crosstype);
            prob1right= probright('1', j, imarker, r, position, crosstype);
            prob2right= probright('2', j, imarker, r, position, crosstype);
            prob0= prob0left*prob0right;
            prob1= prob1left*prob1right;
            prob2= prob2left*prob2right;
            if (ii==saveiaug) {
              if ((prob2>prob1)&&(prob2>prob0)) probmax= newprob[ii]*prob2;
              else if ((prob1>prob0)&&(prob1>prob2)) probmax= newprob[ii]*prob1;
              else probmax= newprob[ii]*prob0;
            }
            if ((prob2>prob1)&&(prob2>prob0)) {
              if (probmax/(newprob[ii]*prob1)<neglect) {
                iaug++;
                newmarker[j][iaug]= '1';
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              if (probmax/(newprob[ii]*prob0)<neglect) {
                iaug++;
                newmarker[j][iaug]= '0';
                newprob[iaug]= newprob[ii]*prob0left;
                newprobmax[iaug]= newprob[iaug]*prob0right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '2';
              newprobmax[ii]= newprob[ii]*prob2;
              newprob[ii]*= prob2left;

            } else if ((prob1>prob2)&&(prob1>prob0)) {
              if (probmax/(newprob[ii]*prob2)<neglect) {
                iaug++;
                newmarker[j][iaug]= '2';
                newprob[iaug]= newprob[ii]*prob2left;
                newprobmax[iaug]= newprob[iaug]*prob2right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              if (probmax/(newprob[ii]*prob0)<neglect) {
                iaug++;
                newmarker[j][iaug]= '0';
                newprob[iaug]= newprob[ii]*prob0left;
                newprobmax[iaug]= newprob[iaug]*prob0right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '1';
              newprobmax[ii]= newprob[ii]*prob1;
              newprob[ii]*= prob1left;
            } else {
              if (probmax/(newprob[ii]*prob1)<neglect) {
                iaug++;
                newmarker[j][iaug]= '1';
                newprob[iaug]= newprob[ii]*prob1left;
                newprobmax[iaug]= newprob[iaug]*prob1right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              if (probmax/(newprob[ii]*prob2)<neglect) {
                iaug++;
                newmarker[j][iaug]= '2';
                newprob[iaug]= newprob[ii]*prob2left;
                newprobmax[iaug]= newprob[iaug]*prob2right;
                for (jj=0; jj<Nmark; jj++) {
                  if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                }
                newind[iaug]=i-((*Nind)-newNind);
                newy[iaug]=y[i];
              }
              newmarker[j][ii]= '0';
              newprobmax[ii]= newprob[ii]*prob0;
              newprob[ii]*= prob0left;
            }
            probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
          } else { // newmarker[j][ii] is observed

            if ((position[j]=='L'||position[j]=='U')) {
              prob0left= start_prob(crosstype, newmarker[j][ii]);
            } else {
              prob0left= prob(newmarker, r, ii, j-1, newmarker[j][ii], crosstype, 1, 0, 0);
            }

            newprob[ii]*= prob0left;

          }

          if (iaug+3>maxNaug) {
            Rprintf("ERROR: Dataset too large after augmentation - your CROSS data may have been corrupted\n");  // FIXME!
            if (verbose) Rprintf("INFO: Recall procedure with larger value for augmentation parameters or lower the parameter neglect\n");
            // Better not free them, we don't know if the arrays already contain something, perhaps not... then we would segfault in R
            //Free(newy);
            //Free(newmarker);
            //Free(newind);
            //Free(newprob);
            //Free(newprobmax);
            //Free(imarker);
            return 0;
          }
        }
      if ((iaug-saveiaug+1)>imaxNaug) {
        newNind-= 1;
        iaug= saveiaug-1;
        if (verbose) Rprintf("INFO: Individual %d is eliminated\n", i);
      }
      sumprob= 0.0;
      for (int ii=saveiaug; ii<=iaug; ii++) sumprob+= newprob[ii];
      for (int ii=saveiaug; ii<=iaug; ii++) newprob[ii]/= sumprob;
    }
    iaug++;
    saveiaug=iaug;
  }
  *Naug= iaug;
  *Nind= newNind;
  *augmarker= newcmatrix(Nmark, *Naug);
  *augy= newvector(*Naug);
  *augind = newivector(*Naug);
  for (int i=0; i<(*Naug); i++) {
    (*augy)[i]= newy[i];
    (*augind)[i]= newind[i];
    for (int j=0; j<Nmark; j++) (*augmarker)[j][i]= newmarker[j][i];
  }
  Free(newy);
  Free(newmarker);
  Free(newind);
  Free(newprob);
  Free(newprobmax);
  Free(imarker);
  return 1;
}

/*
 * The R interfact to data augmentation
 */

void R_augdata(int *geno, double *dist, double *pheno, int *auggeno, 
               double *augPheno, int *augIND, int *Nind, int *Naug, int *Nmark,
               int *Npheno, int *maxaug, int *maxiaug, double *neglect, int
               *chromo, int *crosstype, int *verbose) {
  int **Geno;
  double **Pheno;
  double **Dist;
  int **NEW;
  int **Chromo;
  double **NEWPheno;
  int **NEWIND;
  int prior = *Nind;

  if (*verbose) Rprintf("INFO: Starting C-part of the data augmentation routine\n");
  ivector new_ind;
  vector new_y, r, mapdistance;
  cvector position;
  cmatrix markers, new_markers;
  ivector chr;

  markers= newcmatrix(*Nmark, *Nind);
  new_markers= newcmatrix(*Nmark, *maxaug);
  r = newvector(*Nmark);
  mapdistance = newvector(*Nmark);
  position= newcvector(*Nmark);
  chr= newivector(*Nmark);

  //Reorganise the pointers into arrays, Singletons are just cast into the function
  reorg_geno(*Nind, *Nmark, geno, &Geno);
  reorg_int(*Nmark, 1, chromo, &Chromo);
  reorg_pheno(*Nind, *Npheno, pheno, &Pheno);
  reorg_pheno(*Nmark, 1, dist, &Dist);

  reorg_int(*maxaug, *Nmark, auggeno, &NEW);
  reorg_int((*maxiaug)*(*Nind), 1, augIND, &NEWIND);
  reorg_pheno(*maxaug, 1, augPheno, &NEWPheno);

  //Change all the markers from R/qtl format to MQM internal
  change_coding(Nmark, Nind, Geno, markers, *crosstype);

  char cross = determin_cross(Nmark, Nind, Geno, crosstype);
  if (*verbose) Rprintf("INFO: Filling the chromosome matrix\n");

  for (int i=0; i<(*Nmark); i++) {
    //Set some general information structures per marker
    mapdistance[i]=999.0;
    mapdistance[i]=Dist[0][i];
    chr[i] = Chromo[0][i];
  }

  if (*verbose) Rprintf("INFO: Calculating relative genomepositions of the markers\n");
  for (int j=0; j<(*Nmark); j++) {
    if (j==0) {
      if (chr[j]==chr[j+1]) position[j]='L';
      else position[j]='U';
    } else if (j==(*Nmark-1)) {
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

  if (*verbose) Rprintf("INFO: Estimating recombinant frequencies\n");
  for (int j=0; j<(*Nmark); j++) {
    r[j]= 999.0;
    if ((position[j]=='L')||(position[j]=='M')) {
      r[j]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-mapdistance[j])));
      if (r[j]<0) {
        Rprintf("ERROR: Recombination frequency is negative\n");
        Rprintf("ERROR: Position=%d r[j]=%f\n", position[j], r[j]);
        return;
      }
    }
    //RRprintf("recomfreq:%d, %f\n", j, r[j]);
  }

  if (augdata(markers, Pheno[(*Npheno-1)], &new_markers, &new_y, &new_ind, Nind, Naug, *Nmark, position, r, *maxaug, *maxiaug, *neglect, cross, *verbose)==1) {
    //Data augmentation finished succesfully
    //Push it back into RQTL format
    for (int i=0; i<(*Nmark); i++) {
      for (int j=0; j<(*Naug); j++) {
        NEWPheno[0][j] = new_y[j];
        NEWIND[0][j] = new_ind[j];
        NEW[i][j] = 9;
        if (new_markers[i][j] == '0') {
          NEW[i][j] = 1;
        }
        if (new_markers[i][j] == '1') {
          NEW[i][j] = 2;
        }
        if (new_markers[i][j] == '2') {  // [karl:] this might need to be changed for RIL
          NEW[i][j] = 3;
        }
        if (new_markers[i][j] == '3') {
          NEW[i][j] = 5;
        }
        if (new_markers[i][j] == '4') {
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
    if (*verbose) {
      Rprintf("INFO: Data augmentation finished succesfull\n");
      Rprintf("# Unique individuals before augmentation:%d\n", prior);
      Rprintf("# Unique selected individuals:%d\n", *Nind);
      Rprintf("# Marker p individual:%d\n", *Nmark);
      Rprintf("# Individuals after augmentation:%d\n", *Naug);
    }
  } else {
    //Unsuccessfull data augmentation exit
    *Naug = *Nind;
    for (int i=0; i<(*Nmark); i++) {
      for (int j=0; j<(*Naug); j++) {
        NEWPheno[0][j] = Pheno[0][j];
        NEW[i][j] = 9;
        if (markers[i][j] == '0') {
          NEW[i][j] = 1;
        }
        if (markers[i][j] == '1') {
          NEW[i][j] = 2;
        }
        if (markers[i][j] == '2') { // [karl:] this might need to be changed for RIL
          NEW[i][j] = 3;
        }
        if (markers[i][j] == '3') {
          NEW[i][j] = 5;
        }
        if (markers[i][j] == '4') {
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
    Rprintf("Data augmentation failed\n");
  }
  return;
}


