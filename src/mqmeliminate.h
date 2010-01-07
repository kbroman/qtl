/**********************************************************************
 *
 * mqmeliminate.h
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
     
/* backward elimination in regression of trait on multiple cofactors routine subX haalt uit matrices voor volledige model de submatrices voor submodellen;
   matrices XtWX en Xt van volledig model worden genoemd fullxtwx en fullxt; analoog vector XtWY wordt full xtwy genoemd;
*/
double backward(int Nind, int Nmark, cvector cofactor, MQMMarkerMatrix marker, vector y, vector weight, int* ind, 
                int Naug, double logLfull, double variance, double F1, double F2, cvector* newcofactor, vector r, 
                cvector position,vector *informationcontent,vector *mapdistance,matrix *Frun,int run,char REMLorML,
                bool fitQTL,bool dominance,int em, double windowsize,double stepsize,
                double stepmin,double stepmax,MQMCrossType crosstype,int verbose);

#ifdef __cplusplus
  }
#endif
     
/* end of MQMsupport.h */
