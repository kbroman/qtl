/**********************************************************************
 * 
 * MQMscan.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
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
 * C functions for the R/qtl package
 * Contains: R_scanMQM, scanMQM
 *
 **********************************************************************/


/**********************************************************************
 * 
 * R_scanMQM
 * Wrapper for call from R;
 * 
 **********************************************************************/

void R_scanMQM(int *Nind,int *Nmark,int *Npheno,
			   int *geno,int *chromo, double *dist, double *pheno, 
			   int *cofactors, int *backwards, int *RMLorML,double *alfa,int *emiter,
			   double *windowsize,double *steps,
			   double *stepmi,double *stepma,int *NRun, double *qtl,int *re_estimate,int *crosstype,int *domi,int *verbose);


/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(int Nind, int Nmark,int Npheno,int **Geno,int **Chromo, 
			 double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,int Emiter,
			 double Windowsize,double Steps,
			 double Stepmi,double Stepma,int NRUN, double **QTL,int reestimate,int crosstype,int domi,int verbose);

/**********************************************************************
 * 
 * Helper functions
 *
 *
 **********************************************************************/


double Lnormal(double residual, double variance);
double absdouble(double x);
int mod(int a, int b);
/**********************************************************************
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno);
*/
void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno);

void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno);


/* end of scanMQM.h */
