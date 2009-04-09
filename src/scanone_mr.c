/**********************************************************************
 * 
 * scanone_mr.c
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
 * single QTL model by marker regression (i.e., analysis of variance at 
 * the marker loci)
 *
 * Contains: R_scanone_mr, scanone_mr
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "scanone_mr.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_mr
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls scanone_mr.
 * 
 **********************************************************************/

void R_scanone_mr(int *n_ind, int *n_pos, int *n_gen, int *geno, 
		  double *addcov, int *n_addcov, double *intcov, 
		  int *n_intcov, double *pheno, double *weights, 
		  double *result)
{
  int **Geno;
  double **Addcov, **Intcov;

  reorg_geno(*n_ind, *n_pos, geno, &Geno);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scanone_mr(*n_ind, *n_pos, *n_gen, Geno, Addcov, *n_addcov,
	     Intcov, *n_intcov, pheno, weights, result);
}

/**********************************************************************
 * 
 * scanone_mr
 *
 * Performs genome scan using marker regression.
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Geno         Genotype matrix
 *
 * Addcov       Matrix of additive covariates: Addcov[cov][ind]
 * 
 * n_addcov     Number of columns of Addcov
 *
 * Intcov       Number of interactive covariates: Intcov[cov][ind]
 *
 * n_intcov     Number of columns of Intcov
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Vector of length n_pos, to contain the RSS
 *
 **********************************************************************/

void scanone_mr(int n_ind, int n_pos, int n_gen, int **Geno, 
		double **Addcov, int n_addcov, double **Intcov,
		int n_intcov, double *pheno, double *weights,
		double *result)
{
  int ny, *jpvt, k, i, j, ncol, ncol0, k2, s;
  double *work, *x, *qty, *qraux, *coef, *resid, tol, rss0, *y;
  double rss0_allind=0.0;
  int *which_ind, this_n_ind, done_allind=0; 

  /* tolerance for linear regression */
  tol = TOL;

  ncol = n_gen + (n_gen-1)*n_intcov+n_addcov;
  ncol0 = n_addcov+1;

  /* allocate space and set things up*/
  x = (double *)R_alloc(n_ind*ncol, sizeof(double));
  coef = (double *)R_alloc(ncol, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(ncol, sizeof(int));
  qraux = (double *)R_alloc(ncol, sizeof(double));
  work = (double *)R_alloc(2 * ncol, sizeof(double));
  which_ind = (int *)R_alloc(n_ind, sizeof(int));
  y = (double *)R_alloc(n_ind, sizeof(double));
  ny = 1;

  for(j=0; j<n_ind; j++) 
    pheno[j] *= weights[j];
  /* note: weights are really square-root of weights */

  for(i=0; i<n_pos; i++) {

    R_CheckUserInterrupt(); /* check for ^C */

    /* genotyped individuals at this marker */
    for(j=0, this_n_ind=0; j<n_ind; j++) {
      if(Geno[i][j] > 0) {
	which_ind[this_n_ind] = j;
	y[this_n_ind] = pheno[j];
	this_n_ind++;
      }
    }

    if((this_n_ind < n_ind) || !done_allind) {
      /* the above is to avoid repeatedly doing the null model
	 regression in the case of complete marker data */
      
      /* null model */
      /* fill up X matrix */
      for(j=0; j<this_n_ind; j++) {
        x[j] = weights[which_ind[j]];
        for(k=0; k<n_addcov; k++) 
  	  x[j+(k+1)*this_n_ind] = Addcov[k][which_ind[j]] * weights[which_ind[j]];
      }
      /* linear regression */
      F77_CALL(dqrls)(x, &this_n_ind, &ncol0, y, &ny, &tol, coef, resid,
		      qty, &k, jpvt, qraux, work);

      /* null RSS */
      rss0 = 0.0;
      for(j=0; j<this_n_ind; j++)  rss0 += (resid[j]*resid[j]);

      if(this_n_ind==n_ind) {
	done_allind=1;
	rss0_allind = rss0;
      }
    }
    else /* already ran null model regression with all individuals */
      rss0 = rss0_allind;
      
    /* regression with QTL genotype */
    for(k=0; k<n_gen; k++) jpvt[k] = k;

    /* fill up X matrix */
    for(j=0; j<this_n_ind; j++) {
      for(k=0; k<n_gen; k++) {
	if(Geno[i][which_ind[j]] == k+1) 
	  x[j+k*this_n_ind] = weights[which_ind[j]];
	else x[j+k*this_n_ind] = 0.0;
      }
      for(k=0; k<n_addcov; k++)
	x[j+(k+n_gen)*this_n_ind] = Addcov[k][which_ind[j]] * weights[which_ind[j]];
      for(k=0,s=0; k<n_gen-1; k++) {
	if(Geno[i][which_ind[j]] == k+1) {
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+(n_gen+n_addcov+s)*this_n_ind] = 
	      Intcov[k2][which_ind[j]] * weights[which_ind[j]];
	}
	else {
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+(n_gen+n_addcov+s)*this_n_ind] = 0.0;
	}
      }
    }

    /* linear regression of phenotype on QTL genotype probabilities */
    F77_CALL(dqrls)(x, &this_n_ind, &ncol, y, &ny, &tol, coef, resid,
		    qty, &k, jpvt, qraux, work);

    /* RSS */
    result[i] = 0.0;
    for(j=0; j<this_n_ind; j++) result[i] += (resid[j]*resid[j]);

    /* log10 likelihood */
    result[i] = (double)this_n_ind/2.0* 
      (log10(rss0)-log10(result[i]));

  } /* loop over marker positions */
}

/* end of scanone_mr.c */
