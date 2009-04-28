/**********************************************************************
 * 
 * MQMmapQTL.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Feb, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License, as
 *     published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version. 
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the
 *     GNU General Public License for more details.
 * 
 *     A copy of the GNU General Public License is available at
 *     http://www.r-project.org/Licenses/
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/

double mapQTL(int Nind, int Nmark, cvector cofactor, cvector selcofactor, cmatrix marker, cvector position, vector mapdistance, vector y, 
			  vector r, ivector ind, int Naug, double variance, char printoutput,vector *informationcontent,matrix *Frun,int run,char REMLorML,char fitQTL,char dominance,int em, double windowsize,double stepsize,
			  double stepmin,double stepmax,char crosstype,Mmatrix MendelM,int verbose);
 
/* end of MQMmapQTL.h */
