/**********************************************************************
 * 
 * MQMRegression.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/

		   
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

/* end of MQMRegression.h */
