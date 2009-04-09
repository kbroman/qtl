/**********************************************************************
 * 
 * scanone_hk.c
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
 * single QTL model by Haley-Knott regression
 *
 * Contains: R_scanone_hk, scanone_hk
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
#include "lapackutil.h"
#include "scanone_hk.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk(int *n_ind, int *n_pos, int *n_gen,
		  double *genoprob, double *addcov, int *n_addcov, 
                  double *intcov, int *n_intcov, double *pheno, int *nphe,
		  double *weights, double *result)
{
  double ***Genoprob, **Result, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *nphe, result, &Result); 

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scanone_hk(*n_ind, *n_pos, *n_gen, Genoprob, Addcov, *n_addcov,
	     Intcov, *n_intcov, pheno, *nphe, weights, Result);
}

/**********************************************************************
 * 
 * scanone_hk
 *
 * Performs genome scan using the Haley-Knott regression method
 * (regressing phenotypes on conditional genotype probabilities; the
 * multipoint genotype probabilities have already been calculated in
 * calc.genoprob)
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
 * nphe         number of phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [n_pos x (nphe)] containing the
 *              LOD scores for each phenotype
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
                double **Addcov, int n_addcov, double **Intcov, 
		int n_intcov, double *pheno, int nphe, double *weights, 
		double **Result)
{
  int  i, j, k, k2, s, rank, info, nrss, lwork, ncolx, ind_idx,
    multivar=0;
  double *dwork, *x, *x_bk, *singular, *yfit, *rss, *rss_det=0,
    *work, *tmppheno, *coef;
  double alpha=1.0, beta=0.0, tol=TOL, dtmp;

  /* number of rss's, currently multivar is not used so it's always 0 */
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /* allocate memory */
  rss = (double *)R_alloc(nrss, sizeof(double));
  tmppheno = (double *)R_alloc(n_ind*nphe, sizeof(double));

  /* number of columns in design matrix X for full model */
  ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov;
  /*ncol0 = n_addcov+1;*/
  rank = ncolx;

  /* allocate space and set things up*/
  /*  x = (double *)R_alloc(n_ind*ncol, sizeof(double));
  coef = (double *)R_alloc(ncol, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(ncol, sizeof(int));
  qraux = (double *)R_alloc(ncol, sizeof(double));
  work = (double *)R_alloc(2 * ncol, sizeof(double)); */
  lwork = 3*ncolx+ MAX(n_ind, nphe);
  if(multivar == 1)
    dwork = (double *)R_alloc((2*n_ind+1)*ncolx+lwork+(n_ind+nphe+ncolx)*nphe,
      sizeof(double));
  else
    dwork = (double *)R_alloc((2*n_ind+1)*ncolx+lwork+(n_ind+ncolx)*nphe,
      sizeof(double));

  /* split the memory block */
  singular = dwork;
  work = singular + ncolx;
  x = work + lwork;
  x_bk = x + n_ind*ncolx;
  yfit = x_bk + n_ind*ncolx;
  coef = yfit + n_ind*nphe;
  if(multivar == 1) rss_det = coef + ncolx*nphe;

  /* NULL model is now done in R ********************
     (only do it once!)
  for(j=0; j<n_ind; j++) {
    x[j] = 1.0;
    for(k=0; k<n_addcov; k++) 
      x[j+(k+1)*n_ind] = Addcov[k][j];
  }
  F77_CALL(dqrls)(x, &n_ind, &ncol0, pheno, &ny, &tol, coef, resid,
	          qty, &k, jpvt, qraux, work);
  rss0 = 0.0;
  for(j=0; j<n_ind; j++)  rss0 += (resid[j]*resid[j]);
  Null model is now done in R ********************/

  for(j=0; j<n_ind; j++) 
    for(k=0; k<nphe; k++)
      pheno[j+k*n_ind] *= weights[j];
  /* note: weights are really square-root of weights */

  for(i=0; i<n_pos; i++) { /* loop over positions */
    R_CheckUserInterrupt(); /* check for ^C */

    /* fill up X matrix */
    for(j=0; j<n_ind; j++) {
      for(k=0; k<n_gen; k++)
	x[j+k*n_ind] = Genoprob[k][i][j]*weights[j]; 
      for(k=0; k<n_addcov; k++)
	x[j+(k+n_gen)*n_ind] = Addcov[k][j]*weights[j];
      for(k=0,s=0; k<n_gen-1; k++)
	for(k2=0; k2<n_intcov; k2++,s++) 
	  x[j+(n_gen+n_addcov+s)*n_ind] = Genoprob[k][i][j]*Intcov[k2][j]*weights[j];
    }

    /* linear regression of phenotype on QTL genotype probabilities */
    /*    F77_CALL(dqrls)(x, &n_ind, &ncol, pheno, &ny, &tol, coef, resid,
		    qty, &k, jpvt, qraux, work);
    */
    /* make a copy of x matrix, we may need it */
    memcpy(x_bk, x, n_ind*ncolx*sizeof(double));
    /* make a copy of phenotypes. I'm doing this because 
       dgelss will destroy the input rhs array */
    memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
    /* linear regression of phenotype on QTL genotype probabilities */
    mydgelss (&n_ind, &ncolx, &nphe, x, x_bk, pheno, tmppheno, singular,
      &tol, &rank, work, &lwork, &info);

    /* calculate residual sum of squares */
    if(nphe == 1) {
      /* only one phenotype, this is easier */
      /* if the design matrix is full rank */
      if(rank == ncolx) {
        for (k=rank, rss[0]=0.0; k<n_ind; k++)
          rss[0] += tmppheno[k]*tmppheno[k];
      }
      else {
        /* the desigm matrix is not full rank, this is trouble */
        /* calculate the fitted value */
        matmult(yfit, x_bk, n_ind, ncolx, tmppheno, 1);
        /* calculate rss */
        for (k=0, rss[0]=0.0; k<n_ind; k++)
          rss[0] += (pheno[k]-yfit[k]) * (pheno[k]-yfit[k]);
      }
    }
    else { /* multiple phenotypes */
      if(multivar == 1) {
        /* note that the result tmppheno has dimension n_ind x nphe,
	   the first ncolx rows contains the estimates. */
        for (k=0; k<nphe; k++) 
          memcpy(coef+k*ncolx, tmppheno+k*n_ind, ncolx*sizeof(double));
        /* calculate yfit */
        matmult(yfit, x_bk, n_ind, ncolx, coef, nphe);
        /* calculate residual, put the result in tmppheno */
        for (k=0; k<n_ind*nphe; k++)
          tmppheno[k] = pheno[k] - yfit[k];
        mydgemm(&nphe, &n_ind, &alpha, tmppheno, &beta, rss_det);

        /* calculate the determinant of rss */
        /* do Cholesky factorization on rss_det */
        mydpotrf(&nphe, rss_det, &info);
        for(k=0, rss[0]=1.0;k<nphe; k++)
          rss[0] *= rss_det[k*nphe+k]*rss_det[k*nphe+k];
      } 
      else { /* return rss as a vector */
        if(rank == ncolx) {
          for(k=0; k<nrss; k++) {
            ind_idx = k*n_ind;
            for(j=rank, rss[k]=0.0; j<n_ind; j++) {
              dtmp = tmppheno[ind_idx+j];
              rss[k] += dtmp * dtmp;
            }
          }
        }
        else { /* not full rank, this is troubler */
          /* note that the result tmppheno has dimension n_ind x nphe,
          the first ncolx rows contains the estimates. */
          for (k=0; k<nphe; k++) 
            memcpy(coef+k*ncolx, tmppheno+k*n_ind, ncolx*sizeof(double));
          /* calculate yfit */
          matmult(yfit, x_bk, n_ind, ncolx, coef, nphe);
          /* calculate residual, put the result in tmppheno */
          for (k=0; k<n_ind*nphe; k++)
            tmppheno[k] = pheno[k] - yfit[k];
          /* calculate rss */
          for(k=0; k<nrss; k++) {
            ind_idx = k*n_ind;
            for(j=0, rss[k]=0.0; j<n_ind; j++) {
              dtmp = tmppheno[ind_idx+j];
              rss[k] += dtmp * dtmp;
            }
          }
	  
        }
	
      }
    }
    /* make the result */
    /* log10 likelihood */
    for(k=0; k<nrss; k++) 
       Result[k][i] = (double)n_ind/2.0*log10(rss[k]);

  } /* end loop over positions */
}

/* end of scanone_hk.c */
