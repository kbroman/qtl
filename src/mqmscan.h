/**********************************************************************
 *
 * mqmscan.h
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
 * C functions for the R/qtl package
 * Contains: R_mqmscan, mqmscan
 *
 **********************************************************************/

#ifdef __cplusplus
  extern "C" {
#endif
     
double Lnormal(double residual, double variance);
int mod(int a, int b);
void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno);
void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno);

/* analyseF2 - analyse one F2 family */

void analyseF2(int Nind, int Nmark, cvector *cofactor, cmatrix marker, 
               vector y, ivector f1genotype, int Backwards, double **QTL,vector
               *mapdistance,int **Chromo,int Nrun,int RMLorML, double
               windowsize,double stepsize, double stepmin,double stepmax,double
               alfa,int em,int out_Naug,int **INDlist,char reestimate,char
               crosstype,char dominance,int verbose);

#ifdef __cplusplus
  }
#endif

/* end of mqmscan.h */
