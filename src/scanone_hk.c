/**********************************************************************
 * 
 * scanone_hk.c
 *
 * copyright (c) 2001-2012, Karl W Broman
 *
 * last modified Apr, 2012
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
		  double *weights, double *result, 
		  int *ind_noqtl)
{
  double ***Genoprob, **Result, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_errlod(*n_pos, *nphe, result, &Result); 

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scanone_hk(*n_ind, *n_pos, *n_gen, Genoprob, Addcov, *n_addcov,
	     Intcov, *n_intcov, pheno, *nphe, weights, Result, 
	     ind_noqtl);
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
 * ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *              from QTL effects.  
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
                double **Addcov, int n_addcov, double **Intcov, 
		int n_intcov, double *pheno, int nphe, double *weights, 
		double **Result, int *ind_noqtl)
{
  int  i, j, k, k2, s, rank, n_dwork, ncolx;
  double *dwork, *x, *x_copy, *rss, *pheno_copy, tol=TOL;
  int *jpvt;

  /* allocate memory */

  /* number of columns in design matrix X for full model */
  ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov;
  /*ncol0 = n_addcov+1;*/
  rank = ncolx;

  /* allocate space and set things up*/
  setup_linreg_rss(n_ind, ncolx, nphe, &n_dwork, &dwork, &jpvt);
  rss = (double *)R_alloc(nphe, sizeof(double));
  pheno_copy = (double *)R_alloc(n_ind*nphe, sizeof(double));
  x = (double *)R_alloc(n_ind*ncolx, sizeof(double));
  x_copy = (double *)R_alloc(n_ind*ncolx, sizeof(double));

  for(j=0; j<n_ind; j++) 
    for(k=0; k<nphe; k++)
      pheno[j+k*n_ind] *= weights[j];
  /* note: weights are really square-root of weights */

  /* make copy of phenotypes */
  memcpy(pheno_copy, pheno, n_ind*nphe*sizeof(double));

  for(i=0; i<n_pos; i++) { /* loop over positions */
    R_CheckUserInterrupt(); /* check for ^C */

    for(k=0; k<n_ind * ncolx; k++) x[k] = 0.0;

    /* fill up X matrix */
    for(j=0; j<n_ind; j++) {
      if(!ind_noqtl[j]) {
	for(k=0; k<n_gen; k++)
	  x[j+k*n_ind] = Genoprob[k][i][j]*weights[j];
      }
      for(k=0; k<n_addcov; k++)
	x[j+(k+n_gen)*n_ind] = Addcov[k][j]*weights[j];
      if(!ind_noqtl[j]) {
	for(k=0,s=0; k<n_gen-1; k++)
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+(n_gen+n_addcov+s)*n_ind] = Genoprob[k][i][j]*Intcov[k2][j]*weights[j];
      }
    }

    /* make a copy of x matrix, we may need it */
    memcpy(x_copy, x, n_ind*ncolx*sizeof(double));

    /* linear regression of phenotype on QTL genotype probabilities */
    linreg_rss(n_ind, ncolx, x, nphe, pheno, rss, n_dwork, dwork, jpvt,
               x_copy, pheno_copy, tol);

    /* make the result */
    /* log10 likelihood */
    for(k=0; k<nphe; k++) 
       Result[k][i] = (double)n_ind/2.0*log10(rss[k]);

  } /* end loop over positions */
}

/* end of scanone_hk.c */
