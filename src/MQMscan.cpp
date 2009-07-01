/**********************************************************************
 * 
 * MQMscan.cpp
 *
 * copyright (c) 2009
 *
 * last modified Apr,2009
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
 * Contains: R_scanMQM, scanMQM
 *
 **********************************************************************/
using namespace std;

#include <fstream>
#include <iostream>

extern "C"
{
#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include <math.h>
#include <Rmath.h>
#include "MQMdata.h"
#include "MQMsupport.h"
#include "MQMinterfaces.h"   /*Testing */
#include "MQMreDefine.h"

#include "util.h"
double absdouble(double x)
{
//{      double z; z= (x<0 ? -x : x); return z;}
	return fabs(x);
}


double Lnormal(double residual, double variance){
	//double Likelihood,likelyhood;
	//Likelihood=exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
	//if(absdouble(Likelihood-likelyhood)>0.05){
	//Rprintf("ERROR: Lnormal error\n");
	//}
	//return likelyhood;
	return dnorm(residual,0,sqrt(variance),0);
}



int mod(int a, int b)
{      
	return a%b;
	//int c;
    //c= a/b;
    //return a-b*c;
}

/**********************************************************************
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
  int i;

  *Geno = (int **)R_alloc(n_pos, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_pos; i++) 
    (*Geno)[i] = (*Geno)[i-1] + n_ind;

}

*/
void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno)
{
  int i;

  *Pheno = (double **)R_alloc(n_mar, sizeof(double *));

  (*Pheno)[0] = pheno;
  for(i=1; i< n_mar; i++) 
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}


void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno)
{
  int i;

  *Pheno = (int **)R_alloc(n_mar, sizeof(int *));

  (*Pheno)[0] = pheno;
  for(i=1; i< n_mar; i++) 
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}



/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(int Nind, int Nmark,int Npheno,int **Geno,int **Chromo, 
			 double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,int Emiter,
			 double Windowsize,double Steps,
			 double Stepmi,double Stepma,int NRUN,int out_Naug,int **INDlist, double **QTL, int re_estimate,int crosstype,int domi,int verbose){
	
	ivector f1genotype;
	cmatrix markers;
	cvector cofactor;
	vector mapdistance;
	
	markers= newcmatrix(Nmark,Nind);
	f1genotype = newivector(Nmark);
	cofactor= newcvector(Nmark);  
	mapdistance= newvector(Nmark);
	
	int cof_cnt=0;
 	
   	//Change all the markers from Karl format to our own
	change_coding(&Nmark,&Nind,Geno,markers, crosstype);

	for(int i=0; i< Nmark; i++){
		f1genotype[i] = 12;
		//receiving mapdistances
		mapdistance[i]=999.0;
		mapdistance[i]=Dist[0][i];
	 	cofactor[i] = '0';
		if(Cofactors[0][i] == 1){
			cofactor[i] = '1';
			cof_cnt++;
		}
		if(Cofactors[0][i] == 2){
			cofactor[i] = '2';
			cof_cnt++;
		}
		if(cof_cnt+10 > Nind){
			Rprintf("ERROR: Setting this many cofactors would leave less than 10 degrees of freedom.\n");
			return;
		}
	}

	char reestimate = 'y';
	if(re_estimate == 0){
		reestimate = 'n';
	}
	//determine what kind of cross we have
	char cross = determin_cross(&Nmark,&Nind,Geno,&crosstype);
	//set dominance accordingly
	if(cross != 'F'){
		if(verbose==1){Rprintf("INFO: Dominance setting ignored (dominance=0)\n");}
		domi = 0;
	}else{
		domi= domi;
	}

	char dominance='n';
	if(domi != 0){
		dominance='y';
	}
	
	//WE HAVE EVERYTHING START WITH MAIN SCANNING FUNCTION
	analyseF2(Nind, Nmark, &cofactor, markers, Pheno[(Npheno-1)], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter,out_Naug,INDlist,reestimate,cross,dominance,verbose);
	
	if(re_estimate){
		if(verbose==1){Rprintf("INFO: Sending back the reestimated map used during analysis\n");}
		for(int i=0; i< Nmark; i++){
			Dist[0][i]=mapdistance[i];
		}
	}
	if(Backwards){
		if(verbose==1){Rprintf("INFO: Sending back the model\n");}
		for(int i=0; i< Nmark; i++){
			Cofactors[0][i]=cofactor[i];
		}
	}	
	//Rprintf("Starting Cleanup\n");
	delcmatrix(markers);
	Free(f1genotype);
	Free(cofactor);
	Free(mapdistance);
	if(verbose==1){Rprintf("INFO: All done in C returning to R\n");}
	 #ifndef STANDALONE
	 R_CheckUserInterrupt(); /* check for ^C */
	// R_ProcessEvents(); /*  Try not to crash windows etc*/
	 R_FlushConsole();
	 #endif
	return;
}  /* end of function scanMQM */



/**********************************************************************
 * 
 * R_scanMQM
 * 
 **********************************************************************/

void R_scanMQM(int *Nind,int *Nmark,int *Npheno,
			   int *geno,int *chromo, double *dist, double *pheno, 
			   int *cofactors, int *backwards, int *RMLorML,double *alfa,int *emiter,
			   double *windowsize,double *steps,
			   double *stepmi,double *stepma, int *nRun,int *out_Naug,int *indlist,  double *qtl,int *reestimate,int *crosstype,int *domi,int *verbose){
   int **Geno;
   int **Chromo;
   double **Dist;  
   double **Pheno;   
   double **QTL;   
   int **Cofactors;
   int **INDlist;
   
   //Reorganise the pointers into arrays, ginletons are just cast into the function
   reorg_geno(*Nind,*Nmark,geno,&Geno);
   reorg_int(*Nmark,1,chromo,&Chromo);   
   reorg_pheno(*Nmark,1,dist,&Dist);
   //Here we have  the assumption that step.min is negative this needs to be split in 2
   reorg_pheno(2*(*chromo) * (((*stepma)-(*stepmi))/ (*steps)),1,qtl,&QTL);
   reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
   reorg_int(*Nmark,1,cofactors,&Cofactors);  
   reorg_int(*out_Naug,1,indlist,&INDlist);  
   //Done with reorganising lets start executing
   
   scanMQM(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,*out_Naug,INDlist,QTL, *reestimate,*crosstype,*domi,*verbose);
} /* end of function R_scanMQM */

}
