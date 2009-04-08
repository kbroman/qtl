/**********************************************************************
 * 
 * MQMmapQTL.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/

double mapQTL(int Nind, int Nmark, cvector cofactor, cvector selcofactor, cmatrix marker, cvector position, vector mapdistance, vector y, 
			  vector r, ivector ind, int Naug, double variance, char printoutput,vector *informationcontent,matrix *Frun,int run,char REMLorML,char fitQTL,char dominance,int em, double windowsize,double stepsize,
			  double stepmin,double stepmax,char crosstype,Mmatrix MendelM,int verbose);
 
/* end of MQMmapQTL.h */
