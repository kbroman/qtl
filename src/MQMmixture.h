/**********************************************************************
 *
 * MQMmixture.h
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Apr, 2009
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
 * C external functions used by the MQM algorithm
 * Contains:
 *
 **********************************************************************/


#ifdef __cplusplus
  extern "C" {
#endif



/* ML estimation of recombination frequencies via EM;
    calculation of multilocus genotype probabilities;
    ignorance of unlikely genotypes*/
double rmixture(cmatrix marker, vector weight, vector r,
                cvector position, ivector ind,
                int Nind, int Naug, int Nmark,vector *mapdistance, char reestimate,char crosstype,int verbose);


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
                  vector y, ivector ind, int Nind, int Naug,
                  int Nloci,
                  double *variance, int em, vector *weight,char REMLorML,char fitQTL,char dominance,char crosstype,int verbose);

       
#ifdef __cplusplus
  }
#endif

/* end of MQMmixture.c */
