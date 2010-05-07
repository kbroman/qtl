/**********************************************************************
 *
 * mqmmixture.h
 *
 * Copyright (c) 1996-2010 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Danny Arends and Pjotr Prins
 * last modified May 2010
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


/* ML estimation of recombination frequencies via EM;
    calculation of multilocus genotype probabilities;
    ignorance of unlikely genotypes*/
double rmixture(MQMMarkerMatrix marker, vector weight, vector r,
                cvector position, ivector ind,
                int Nind, int Naug, int Nmark,vector *mapdistance, char reestimate,MQMCrossType crosstype,int verbose);


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(MQMMarkerMatrix loci, cvector cofactor, vector r, cvector position,
                  vector y, ivector ind, int Nind, int Naug,
                  int Nloci,
                  double *variance, int em, vector *weight, const bool useREML,bool fitQTL,bool dominance,MQMCrossType crosstype,int verbose);

       
#ifdef __cplusplus
  }
#endif

/* end of mqmmixture.h */
