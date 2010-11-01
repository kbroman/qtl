/**********************************************************************
 *
 * mqmaugment.h
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Danny Arends and Pjotr Prins
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

#ifdef __cplusplus
  extern "C" {
#endif

/*
Data augmentation routing
*/

void R_mqmaugment(int *geno, double *dist, double *pheno, int *auggeno, 
               double *augPheno, int *augIND, int *Nind, int *Naug, int *Nmark,
               int *Npheno, int *maxind, int *maxiaug, double *minprob, int
               *chromo,int* augment_strategy,int *crosstype, int *verbose);
               
int mqmaugmentfull(MQMMarkerMatrix* markers,int* nind, int* augmentednind, ivector* INDlist,
                  double neglect_unlikely, int max_totalaugment, int max_indaugment,
                  const matrix* pheno_value,const int nmark,const ivector chr,const vector mapdistance,
                  const int augment_strategy, const MQMCrossType crosstype,const int verbose);      

int calculate_augmentation(const int Nind, int const Nmark,const MQMMarkerMatrix markers, const MQMCrossType crosstype);         

int mqmaugment(const MQMMarkerMatrix marker, const vector y, 
               MQMMarkerMatrix* augmarker, vector *augy, 
               ivector* augind, ivector* sucind, int *Nind, int *Naug, const int Nmark, 
               const cvector position, vector r, const int maxNaug, 
               const int imaxNaug, const double minprob, 
               const MQMCrossType crosstype, const int verbose);
               
MQMMarker randommarker(const MQMCrossType crosstype);

#ifdef __cplusplus
  }
#endif

