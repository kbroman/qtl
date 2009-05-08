/**********************************************************************
 * 
 * MQMsupport.h
 *
 * copyright (c) 2009 Danny Arends
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
 * Contains: analyseF2, backward
 *
 **********************************************************************/

/* analyseF2 - analyse one F2 family */

void analyseF2(int Nind, int Nmark, cvector *cofactor, cmatrix marker, vector y, ivector f1genotype, int Backwards, 
			   double **QTL,vector *mapdistance,int **Chromo,int Nrun,int RMLorML, double windowsize,double stepsize,
			   double stepmin,double stepmax,double alfa,int em,int out_Naug,int **INDlist,char reestimate,char crosstype,char dominance,int verbose);

/* backward elimination in regression of trait on multiple cofactors routine subX haalt uit matrices voor volledige model de submatrices voor submodellen;
   matrices XtWX en Xt van volledig model worden genoemd fullxtwx en fullxt; analoog vector XtWY wordt full xtwy genoemd;
*/
double backward(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector weight, int* ind, int Naug, double logLfull, double variance, double F1, double F2, cvector* newcofactor, vector r, cvector position,vector *informationcontent,vector *mapdistance,matrix *Frun,int run,char REMLorML,char fitQTL,char dominance,int em, double windowsize,double stepsize,
			  double stepmin,double stepmax,char crosstype,int verbose);
/* end of MQMsupport.h */
