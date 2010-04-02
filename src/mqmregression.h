/**********************************************************************
 *
 * mqmregression.h
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Danny Arends and Pjotr Prins
 * last modified September 2009
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

int designmatrixdimensions(const cvector cofactor,const unsigned int nmark,const bool dominance);

/*
 *  Used regression (perhaps change it to something faster)
 *  regression of trait on multiple cofactors  y=xb+e with weight w
*	(xtwx)b=(xtw)y
*	b=inv(xtwx)(xtw)y
 */

double regression(int Nind, int Nmark, cvector cofactor, MQMMarkerMatrix marker, vector y, vector* weight, ivector ind, int Naug, double *variance, vector Fy,const bool biasadj,const bool fitQTL,const bool dominance);

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
     
/* end of mqmregression.h */
