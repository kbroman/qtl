/**********************************************************************
 * 
 * scanone_em_covar.c
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Nov, 2001
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
 * These functions are for performing a genome scan with a 
 * single QTL model by interval mapping (the EM algorithm) in the 
 * presence of covariates.
 *
 * Contains: scanone_em_covar, estep_em_covar, mstep_em_covar
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "scanone_em.h"
#include "scanone_em_covar.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * scanone_em_covar
 *
 * Performs genome scan using interval mapping in the presence of
 * covariates.  (The multipoint genotype probabilities have already 
 * been calculated in calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates indexed as 
 *              Addcov[cov][ind]
 *
 * n_addcov     Number of columns in Addcov
 *
 * Intcov       Matrix of interactive covariates indexed as 
 *              Intcov[cov][ind]
 *
 * n_intcov     Number of columns in Intcov
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Result vector of length n_pos; upon return, contains 
 *              the LOD scores.
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * verbose        If 1, print out log likelihood at each iteration
 *
 **********************************************************************/

void scanone_em_covar(int n_ind, int n_pos, int n_gen, 
		      double ***Genoprob, double **Addcov, int n_addcov,
		      double **Intcov, int n_intcov, double *pheno, 
		      double *weights,
		      double *result, int maxit, double tol, int verbose)
{
  int i, j, k, s, flag=0, n_par;
  double **wts, *param, *oldparam, regtol; 
  double *x, *coef, *resid, *qty, *qraux, *work;
  int *jpvt, ny, ncol0, error_flag;
  double llik, oldllik=0.0, *work1, *work2, temp, sw;

  /* recenter phenotype to have mean 0, for
     possibly increased numerical stability */
  for(i=0, temp=0.0; i<n_ind; i++) temp += pheno[i];
  temp /= (double)n_ind;
  for(i=0; i<n_ind; i++) pheno[i] -= temp;

  n_par = 1 + n_gen + n_addcov + (n_gen-1)*n_intcov;

  /* Allocate space */
  allocate_dmatrix(n_gen, n_ind, &wts);
  param = (double *)R_alloc(n_par, sizeof(double));
  oldparam = (double *)R_alloc(n_par, sizeof(double));
  work1 = (double *)R_alloc((n_par-1)*(n_par-1),sizeof(double));
  work2 = (double *)R_alloc(n_par-1, sizeof(double));

  ncol0 = 1+n_addcov;
  x = (double *)R_alloc(n_ind*ncol0, sizeof(double));
  coef = (double *)R_alloc(ncol0, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(ncol0, sizeof(int));
  qraux = (double *)R_alloc(ncol0, sizeof(double));
  work = (double *)R_alloc(2*ncol0, sizeof(double));
  regtol = TOL;
  ny = 1;

  /* adjust phenotypes and covariates with weights */
  /* Note: weights are actually sqrt(weights) */
  sw = 0.0;
  for(i=0; i<n_ind; i++) {
    pheno[i] *= weights[i];
    for(j=0; j<n_addcov; j++)
      Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++)
      Intcov[j][i] *= weights[i];
    sw += log(weights[i]);   /* sum of log weights */
  }

  /* NULL model is now done in R ********************
     (only do it once!)
  for(i=0; i<n_ind; i++) {
    x[i] = 1.0;
    for(j=0; j<n_addcov; j++)
      x[i+(j+1)*n_ind] = Addcov[j][i];
  }
  j=0;
  F77_CALL(dqrls)(x, &n_ind, &ncol0, pheno, &ny, &regtol, coef,
	          resid, qty, &k, jpvt, qraux, work);
  sigma0 = llik0 = 0.0;
  for(i=0; i<n_ind; i++) sigma0 += (resid[i] * resid[i]);
  sigma0 = sqrt(sigma0/(double)n_ind);
  for(i=0; i<n_ind; i++) 
    llik0 += log10(dnorm(resid[i],0.0,sigma0,0));
  Null model is now done in R ********************/
  
  /* begin genome scan */
  for(i=0; i<n_pos; i++) { /* loop over marker positions */
    /* initial estimates */
    for(j=0; j<n_ind; j++)
      for(k=0; k<n_gen; k++) 
	wts[k][j] = Genoprob[k][i][j];
    mstep_em_covar(n_ind, n_gen, Addcov, n_addcov, Intcov, n_intcov,
		   pheno, weights, wts, oldparam, work1, work2, &error_flag);
    
    if(!error_flag) { /* only proceed if there's no error */

      if(verbose) {
	estep_em_covar(n_ind, n_gen, i, Genoprob, Addcov, n_addcov,
		       Intcov, n_intcov, pheno, weights, wts, oldparam, 0);
	oldllik = 0.0;
	for(k=0; k<n_ind; k++) {
	  temp = 0.0;
	  for(s=0; s < n_gen; s++) temp += wts[s][k];
	  oldllik += log(temp);
	}
	Rprintf("    %3d %12.6lf\n", i+1, oldllik);
      }
      
      /* begin EM iterations */
      for(j=0; j<maxit; j++) {

	R_CheckUserInterrupt(); /* check for ^C */

	estep_em_covar(n_ind, n_gen, i, Genoprob, Addcov, n_addcov,
		       Intcov, n_intcov, pheno, weights, wts, oldparam, 1);
	mstep_em_covar(n_ind, n_gen, Addcov, n_addcov, Intcov, n_intcov,
		       pheno, weights, wts, param, work1, work2, &error_flag);
	
	if(error_flag) { /* error: X'X singular; break out of EM */
	  flag=0; 
	  break; 
	}

	if(verbose) {
	  estep_em_covar(n_ind, n_gen, i, Genoprob, Addcov, n_addcov,
			 Intcov, n_intcov, pheno, weights, wts, param, 0);
	  llik = 0.0;
	  for(k=0; k<n_ind; k++) {
	    temp = 0.0;
	    for(s=0; s < n_gen; s++) temp += wts[s][k];
	    llik += log(temp);
	  }
	  Rprintf("    %3d %4d %12.6lf\n", i+1, j+1, llik-oldllik);
	  oldllik = llik;
	}

	/* check for convergence */
	flag = 0;
	for(k=0; k<n_par; k++) {
	  if(fabs(param[k]-oldparam[k]) > tol*(fabs(oldparam[k])+tol*100.0)) {
	    flag = 1;
	    break; 
	  }
	}
	if(!flag) break;
	for(k=0; k<n_par; k++) oldparam[k] = param[k];
	
      } /* end of EM iterations */

      if(flag) warning("Didn't converge!\n");

      if(!error_flag) { /* skip if there was an error */
	/* calculate log likelihood */
	estep_em_covar(n_ind, n_gen, i, Genoprob, Addcov, n_addcov,
		       Intcov, n_intcov, pheno, weights, wts, param, 0);
	llik = 0.0;
	for(k=0; k<n_ind; k++) {
	  temp = 0.0;
	  for(s=0; s < n_gen; s++) temp += wts[s][k];
	  llik += log(temp);
	}
	result[i] = -(llik+sw)/log(10.0);
      }
      else result[i] = NA_REAL;
      
      if(verbose) {
	if(error_flag) Rprintf("    %3d NA", i+1);
	else {
	  Rprintf("    %3d %12.6lf\n", i+1, result[i]);
	  Rprintf("        ");
	  for(j=0; j<n_par; j++) Rprintf(" %7.4lf", param[j]);
	}
	Rprintf("\n\n");
      }

    } /* no initial error */
  } /* end loop over marker positions */
}

/**********************************************************************
 * 
 * mstep_em_covar: M-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen    Number of possible QTL genotypes
 *
 * Addcov   Additive covariates
 *
 * n_addcov Number of columns in Addcov
 *
 * Intcov   Interactive covariates
 *
 * n_intcov Number of columns in Intcov
 *
 * pheno    Phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * wts      Pr(QTL gen | phenotype, model, multipoint marker data),
 *          indexed as wts[gen][ind]
 *
 * param    On output, the updated parameter estimates
 *
 * work1    Workspace of doubles, of length (n_par-1)*(n_par-1)
 *
 * work2    Workspace of doubles, of length (n_par-1)
 *
 * error_flag  Set to 1 if E(X'X) is singular
 *
 **********************************************************************/

void mstep_em_covar(int n_ind, int n_gen, double **Addcov, int n_addcov, 
		    double **Intcov, int n_intcov, double *pheno, 
		    double *weights,
		    double **wts, double *param, double *work1, 
		    double *work2, int *error_flag)
{
  int i, j, k, s, sk, nparm1, info;
  double rcond;

  *error_flag=0;

  nparm1 = n_gen + n_addcov + (n_gen-1)*n_intcov;

  /* calculate {E(X)}' y */
  for(j=0; j<nparm1; j++) work2[j] = 0.0;
  for(i=0; i<n_ind; i++) {
    for(j=0; j<n_gen; j++) /* QTL effects */
      work2[j] += (wts[j][i]*pheno[i]*weights[i]);
    for(j=0,k=n_gen; j<n_addcov; j++,k++) /* add covar */
      work2[k] += (Addcov[j][i]*pheno[i]);
    for(j=0,s=n_gen+n_addcov; j<n_gen-1; j++) {
      for(k=0; k<n_intcov; k++, s++)  /* int covar */
	work2[s] += (wts[j][i]*Intcov[k][i]*pheno[i]);
    }
  }

  /* calculate E{X'X}; only the upper right triangle is needed */
  for(j=0; j<nparm1*nparm1; j++) work1[j] = 0.0;
  for(i=0; i<n_ind; i++) {
    for(j=0; j<n_gen; j++) /* QTLxQTL */
      work1[j+nparm1*j] += wts[j][i]*weights[i]*weights[i];
    for(j=0, k=n_gen; j<n_addcov; j++, k++) {
      for(s=j, sk=k; s<n_addcov; s++, sk++)  /* add x add */
	work1[k+nparm1*sk] += (Addcov[j][i]*Addcov[s][i]);
      for(s=0; s<n_gen; s++) /* QTL x add */
	work1[s+nparm1*k] += (Addcov[j][i]*wts[s][i]*weights[i]);
    }
    for(j=0, k=n_gen+n_addcov; j<n_gen-1; j++, k += n_intcov) {
      for(s=0; s<n_intcov; s++) {
	for(sk=s; sk<n_intcov; sk++) /* int x int */
	  work1[k+s+(sk+k)*nparm1] += 
	    (Intcov[s][i]*wts[j][i]*Intcov[sk][i]);
	for(sk=0; sk<n_addcov; sk++) /* add x int */
	  work1[sk+n_gen+(k+s)*nparm1] += 
	    (Addcov[sk][i]*wts[j][i]*Intcov[s][i]);
	work1[j+(k+s)*nparm1] += wts[j][i]*Intcov[s][i]*weights[i]; /* qtl x int */
      }
    }
  }

  /* Copy upper triangle into lower triangle.                   */
  /* This isn't be needed, since the Fortran functions below    */
  /* assume a symmetric matrix and use only the upper triangle. */
  /*                                                            */
  /*  for(j=0; j<nparm1-1; j++)                                 */
  /*    for(k=j+1; k<nparm1; k++)                               */
  /*      work1[k+j*nparm1] = work1[j+k*nparm1];                */

  /* solve work1 * beta = work2 for beta */
  F77_CALL(dpoco)(work1, &nparm1, &nparm1, &rcond, param, &info);
  if(fabs(rcond) < TOL || info != 0) { /* error! */
    warning("X'X matrix is singular.");
    *error_flag = 1;
  }
  else {
    for(j=0; j<nparm1; j++) param[j] = work2[j];
    F77_CALL(dposl)(work1, &nparm1, &nparm1, param);

    /* calculate residual SD */
    param[nparm1] = 0.0;
    for(i=0; i<n_ind; i++) param[nparm1] += pheno[i]*pheno[i];
    for(j=0; j<nparm1; j++) param[nparm1] -= (work2[j]*param[j]);
  
    param[nparm1] = sqrt(param[nparm1] / (double)n_ind);
  }
}


/**********************************************************************
 * 
 * estep_em_covar: E-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen    Number of possible QTL genotypes
 *
 * pos      Position of Genoprob[][][] to consider
 *
 * Genoprob Pr(QTL gen | multipoint marker data)
 *
 * Addcov   Additive covariates
 *
 * n_addcov Number of columns in Addcov
 *
 * Intcov   Interactive covariates
 *
 * n_intcov Number of columns in Intcov
 *
 * pheno    Phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * wts      On output, Pr(QTL gen | pheno, model, multipt marker data), 
 *          indexed as wts[gen][ind]
 *
 * param    Current parameter estimates
 *
 * rescale  If 1, rescale weights so that the sum to 1.
 *          This is done so that by taking rescale=0, we can easily
 *          calculate the log likelihood 
 *
 **********************************************************************/

void estep_em_covar(int n_ind, int n_gen, int pos, double ***Genoprob,
		    double **Addcov, int n_addcov, double **Intcov,
		    int n_intcov, double *pheno, double *weights,
		    double **wts, double *param, int rescale)
{
  int i, j, k, s, nparm1;
  double temp;

  /* nparm1 = number of parameters - 1 (ie, location of resid SD) */
  nparm1 = n_gen + n_addcov + (n_gen-1)*n_intcov;

  for(i=0; i<n_ind; i++) {

    /* calculate fitted values */
    /* additive cov effect */
    temp = 0.0;
    for(j=0, k=n_gen; j<n_addcov; j++, k++)
      temp += (Addcov[j][i]*param[k]);

    /* QTL effect + addcov effect*/
    for(j=0; j<n_gen; j++) wts[j][i] = param[j]*weights[i]+temp;

    /* interactive cov effect */
    for(j=0, s=n_gen+n_addcov; j<n_gen-1; j++) {
      /* s gives location in param vector */
      for(k=0; k<n_intcov; k++, s++) 
	wts[j][i] += (Intcov[k][i]*param[s]);
    }
    /* done calculating fitted values */

    /* calculate p(y|fitted,SD) for normal model 
       and multiple by Genoprob */
    temp=0.0;
    for(j=0; j<n_gen; j++) 
      temp += 
	(wts[j][i] = (dnorm(pheno[i],wts[j][i],param[nparm1],0)*
		      Genoprob[j][pos][i]));
    
    /* rescale wts */
    if(rescale)
      for(j=0; j<n_gen; j++) wts[j][i] /= temp;
      
  } /* end loop over individuals */
    
}


/* end of scanone_em_covar.c */

