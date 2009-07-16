/**********************************************************************
 *
 * MQMRegression.h
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
      

/*
 * Used regression (perhaps change it to something faster)
 *  regression of trait on multiple cofactors  y=xb+e with weight w
*							(xtwx)b=(xtw)y
*							b=inv(xtwx)(xtw)y
 */

double regression(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector* weight, ivector ind, int Naug, double *variance, vector Fy, char biasadj,char fitQTL,char dominance);

/*
-----------------------------------------------------------------------
subroutines from book 'Numerical Recipees in C' for calculating F-probabilities and
for generating randomly permuted trait data for other tasks
-----------------------------------------------------------------------*/

void ludcmp(matrix m, int dim, ivector ndx, int *d);
void lusolve(matrix lu, int dim, ivector ndx, vector b);
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
double inverseF(int df1, int df2, double alfa,int verbose);
#ifdef __cplusplus
   }
#endif
     
/* end of MQMRegression.h */
