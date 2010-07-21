/**********************************************************************
 * 
 * scanone_hk_binary.c
 *
 * copyright (c) 2010, Karl W Broman
 *
 * last modified Jul, 2010
 * first written Jun, 2010
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
 * single QTL model by Haley-Knott regression
 *
 * Contains: R_scanone_hk_binary, scanone_hk_binary
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
#include "scanone_hk_binary.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_hk_binary
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk_binary(int *n_ind, int *n_pos, int *n_gen,
			 double *genoprob, double *addcov, int *n_addcov, 
			 double *intcov, int *n_intcov, double *pheno,
			 double *result, double *tol, int *maxit, 
			 int *verbose, int *ind_noqtl)
{
  double ***Genoprob, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scanone_hk_binary(*n_ind, *n_pos, *n_gen, Genoprob, Addcov, *n_addcov,
		    Intcov, *n_intcov, pheno, result, *tol, *maxit, 
		    *verbose, ind_noqtl);
}

/**********************************************************************
 * 
 * scanone_hk_binary
 *
 * Performs genome scan using the Haley-Knott regression method
 * for a binary trait
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              Indexed as Genoprob[gen][pos][ind]
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
 * result       vector of length n_ind, to contain the log10 likelihood values
 *
 * tol          tolerance for convergence
 *
 * maxit        maximum number of iterations
 *
 * verbose      if TRUE, give some output
 *
 * ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *              from QTL effects.  
 *
 **********************************************************************/

void scanone_hk_binary(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		       double **Addcov, int n_addcov, double **Intcov, 
		       int n_intcov, double *pheno, 
		       double *result, double tol, int maxit, int verbose,
		       int *ind_noqtl)
{
  int i, j, k, k2, kk, s, ncolx, thepos, flag, ny, *jpvt;
  double *dwork, *x, *x_bk, *coef, *resid, *qty, *qraux, *ests;
  double *z, *nu, *wt, *pi, tol2;
  double curllik, llik=0.0;

  ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov;
  tol2 = TOL;
  ny = 1;

  /* allocate space and set things up*/
  /* lengths: 2*ncolx:     dwork
              n_ind*ncolx: x, x_bk
              ncolx:       coef, qraux, ests
	      n_ind:       resid, qty, z, nu, wt, pi */
             
  dwork = (double *)R_alloc(2*n_ind*ncolx + ncolx*5 + n_ind*6, sizeof(double)); 
  x = dwork + 2*ncolx;
  x_bk = x + n_ind*ncolx;
  coef = x_bk + n_ind*ncolx;
  resid = coef + ncolx;
  qty = resid + n_ind;
  qraux = qty + n_ind;
  z = qraux + ncolx;
  nu = z + n_ind;
  wt = nu + n_ind;
  pi = wt + n_ind;
  ests = pi + n_ind; /* length ncolx */

  jpvt = (int *)R_alloc(ncolx, sizeof(int));

  for(thepos=0; thepos<n_pos; thepos++) { /* loop over positions */

    for(i=0; i<n_ind*ncolx; i++) x[i] = 0.0;

    /* fill up X matrix */
    for(j=0; j<n_ind; j++) {
      if(!ind_noqtl[j]) 
	for(k=0; k<n_gen; k++) 
	  x[j+k*n_ind] = Genoprob[k][thepos][j];
      for(k=0; k<n_addcov; k++)
	x[j+(k+n_gen)*n_ind] = Addcov[k][j];
      if(!ind_noqtl[j]) {
	for(k=0,s=0; k<n_gen-1; k++)
	  for(k2=0; k2<n_intcov; k2++,s++)
	    x[j+(n_gen+n_addcov+s)*n_ind] = Genoprob[k][thepos][j]*Intcov[k2][j];
      }
    }

    /* make a copy of x matrix */
    memcpy(x_bk, x, n_ind*ncolx*sizeof(double));

    /* starting point for IRLS */
    curllik = 0.0;
    for(j=0; j<n_ind; j++) {
      pi[j] = (pheno[j] + 0.5)/2.0;
      wt[j] = sqrt(pi[j] * (1.0-pi[j]));
      nu[j] = log(pi[j]) - log(1.0-pi[j]);
      z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
      curllik += pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]);
    }
    if(verbose>1) Rprintf("\t%-5d %-5d %-10.5lf\n", thepos+1, 0, curllik);


    /* multiply design matrix by current wts */
    for(i=0; i<ncolx; i++) 
      for(j=0; j<n_ind; j++)
	x[i*n_ind + j] *= wt[j];

    flag = 0;
    for(s=0; s<maxit; s++) { /* IRLS iterations */
      R_CheckUserInterrupt(); /* check for ^C */

      /* make jpvt = numbers 0, 1, ..., (ncolx-1) */
      /*      jpvt keeps track of any permutation of X columns */
      for(i=0; i<ncolx; i++) jpvt[i] = i;
  
      /* call dqrls to fit regression model */
      F77_CALL(dqrls)(x, &n_ind, &ncolx, z, &ny, &tol2, coef, resid,
		      qty, &kk, jpvt, qraux, dwork);
  
      /* get ests; need to permute back */
      for(i=0; i<kk; i++) ests[jpvt[i]] = coef[i];
      for(i=kk; i<ncolx; i++) ests[jpvt[i]] = 0.0;

      if(verbose>1) {
	for(i=0; i<ncolx; i++) Rprintf("%10.5lf ", ests[i]);
	Rprintf("\n");
      }
      
      /* re-form design matrix */
      memcpy(x, x_bk, ncolx*n_ind*sizeof(double));
  
      /* calculate fitted values, probs, new wts, new z's */
      llik = 0.0;
      for(j=0; j<n_ind; j++) {
        nu[j] = 0.0;
        for(i=0; i<ncolx; i++) 
	  nu[j] += x[i*n_ind+j] * ests[i];
        pi[j] = exp(nu[j]);
        pi[j] /= (1.0 + pi[j]);
        wt[j] = sqrt(pi[j] * (1.0-pi[j]));
        z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
        llik += (pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]));

	if(verbose>2) 
	  Rprintf("\t\t%-4d %1lf %-7.5lf\n", j, pheno[j], pi[j]);
  
        /* multiply design matrix by new weights */
        for(i=0; i<ncolx; i++) 
	  x[i*n_ind+j] *= wt[j];
      }
  
      if(verbose>1) Rprintf("\t%-5d %-5d %-10.5lf\n", thepos+1, s+1, llik);

      if(fabs(llik - curllik) < tol) { /* converged? */
        flag = 1;
        break;
      }
      curllik = llik;
  
    } /* end of IRLS iterations */
  
    if(!flag)
      warning("Didn't converge.");

    result[thepos] = llik;
    if(verbose) Rprintf("%-5d final %-10.5lf\n", thepos+1, llik);

  } /* end loop over positions */
}

/* end of scanone_hk_binary.c */
