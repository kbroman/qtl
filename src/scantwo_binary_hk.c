/**********************************************************************
 * 
 * scantwo_binary_hk.c
 *
 * copyright (c) 2010-2014, Karl W Broman
 *
 * last modified Mar, 2014
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
 * These functions are for performing a two-dimensional genome scan with  
 * a two-QTL model by Haley-Knott regression
 *
 * Contains: R_scantwo_1chr_binary_hk, scantwo_1chr_binary_hk, 
 *           R_scantwo_2chr_binary_hk, scantwo_2chr_binary_hk
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
#include "scantwo_binary_hk.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scantwo_1chr_binary_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_binary_hk.
 * 
 **********************************************************************/

void R_scantwo_1chr_binary_hk(int *n_ind, int *n_pos, int *n_gen,
			      double *genoprob, double *pairprob, 
			      double *addcov, int *n_addcov, 
			      double *intcov, int *n_intcov, 
			      double *pheno, 
			      double *result, int *n_col2drop, int *col2drop,
			      double *tol, int *maxit, int *verbose)
{
  double ***Genoprob, **Result, **Addcov=0, **Intcov=0, *****Pairprob;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
  reorg_errlod(*n_pos, *n_pos, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_1chr_binary_hk(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob, 
			 Addcov, *n_addcov, Intcov, *n_intcov, 
			 pheno, Result, *n_col2drop, col2drop, *tol, 
			 *maxit, *verbose);
}

/**********************************************************************
 * 
 * scantwo_1chr_binary_hk
 *
 * Performs a 2-dimensional genome scan using the Haley-Knott 
 * regression method (regressing phenotypes on conditional genotype 
 * probabilities) for a two-QTL model with the two QTL residing on
 * the same chromosome.
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
 * Pairprob     Array of joint genotype probabilities for QTL
 *              pairs; indexed as Pairprob[gen1][gen2][pos1][pos2][ind]
 *              where pos2 > pos1 (for pos2 <= pos1, points to nothing)
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
 * Result       vector to contain LOD scores
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              testing epistasis.
 *              Note: indexed as Result[col][row]
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_1chr_binary_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
			    double *****Pairprob, double **Addcov, int n_addcov, 
			    double **Intcov, int n_intcov, double *pheno, 
			    double **Result, int n_col2drop,
			    int *col2drop, double tol, int maxit, int verbose)
{
  int n_col_a, n_col_f, n_gen_sq, n_col_a_temp, n_col_f_temp;
  int *jpvt, ny=1, flag;
  int i, i2, j, k, k2, k3, kk, s;
  double *dwork, *x, *x_bk, *coef, *resid, *qty, *qraux, *ests;
  double tol2=TOL, *z, *nu, *wt, *pi, curllik, llik=0.0;
  int *allcol2drop;

  n_gen_sq = n_gen*n_gen;
  /* no. param in additive QTL model */
  n_col_a = (n_gen*2-1)+n_addcov+n_intcov*(n_gen-1)*2; 
  /* no. param full model */
  n_col_f = n_gen_sq+n_addcov+n_intcov*(n_gen_sq-1); 

  /* expand col2drop */
  if(n_col2drop) {
    allocate_int(n_col_f, &allcol2drop);
    expand_col2drop(n_gen, n_addcov, n_intcov, 
		    col2drop, allcol2drop);
  }

  /* allocate space and set things up*/
  /* lengths: 2*n_col_f:     dwork
              n_ind*n_col_f: x, x_bk
              n_col_f:       coef, qraux, ests
	      n_ind:       resid, qty, z, nu, wt, pi */
  dwork = (double *)R_alloc(2*n_ind*n_col_f + n_col_f*5 + n_ind*6, sizeof(double)); 
  x = dwork + 2*n_col_f;
  x_bk = x + n_ind*n_col_f;
  coef = x_bk + n_ind*n_col_f;
  resid = coef + n_col_f;
  qty = resid + n_ind;
  qraux = qty + n_ind;
  z = qraux + n_col_f;
  nu = z + n_ind;
  wt = nu + n_ind;
  pi = wt + n_ind;
  ests = pi + n_ind; /* length n_col_f */

  jpvt = (int *)R_alloc(n_col_f, sizeof(int));

  for(i=0; i<n_pos-1; i++) { 
    for(i2=i+1; i2<n_pos; i2++) { /* loop over pairs of positions */

      R_CheckUserInterrupt(); /* check for ^C */

      /* ADDITIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob[k][i][j];  /* s keeps track of column */
	for(k=0; k<n_gen-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob[k][i2][j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k2=0; k2<n_intcov; k2++) {
	  for(k=0; k<n_gen-1; k++, s++) /* interactive x QTL 1 */
	    x[j+s*n_ind] = Genoprob[k][i][j]*Intcov[k2][j];
	  for(k=0; k<n_gen-1; k++, s++) /* interactive x QTL 2 */
	    x[j+s*n_ind] = Genoprob[k][i2][j]*Intcov[k2][j];
	}
      }

      /* drop cols */
      n_col_a_temp = n_col_a;
      if(n_col2drop) 
	dropcol_x(&n_col_a_temp, n_ind, allcol2drop, x);

      /* make a copy of x matrix */
      memcpy(x_bk, x, n_ind*n_col_a_temp*sizeof(double));

      /* starting point for IRLS */
      curllik = 0.0;
      for(j=0; j<n_ind; j++) {
	pi[j] = (pheno[j] + 0.5)/2.0;
	wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	nu[j] = log(pi[j]) - log(1.0-pi[j]);
	z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	curllik += pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]);

      }
      if(verbose>2) Rprintf("\nadd: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, 0, curllik);


      /* multiply design matrix by current wts */
      for(k=0; k<n_col_a_temp; k++) 
	for(j=0; j<n_ind; j++)
	  x[k*n_ind + j] *= wt[j];

      flag = 0;
      for(s=0; s<maxit; s++) { /* IRLS iterations */
	R_CheckUserInterrupt(); /* check for ^C */

        /* make jpvt = numbers 0, 1, ..., (n_col_a_temp-1) */
        /*      jpvt keeps track of any permutation of X columns */
        for(k=0; k<n_col_a_temp; k++) jpvt[k] = k;
  
        /* call dqrls to fit regression model */
        F77_CALL(dqrls)(x, &n_ind, &n_col_a_temp, z, &ny, &tol2, coef, resid,
			qty, &kk, jpvt, qraux, dwork);
  
	/* get ests; need to permute back */
	for(k=0; k<kk; k++) ests[jpvt[k]] = coef[k];
	for(k=kk; k<n_col_a_temp; k++) ests[jpvt[k]] = 0.0;

#ifdef UNDEFINED
	if(kk < n_col_a_temp) {
	  Rprintf("\nkk=%d n_col_a_temp=%d\n", kk, n_col_a_temp);
	  for(k=0; k<n_col_a_temp; k++) 
	    Rprintf("%10.5lf ", ests[k]);
	  Rprintf("\n");
	}	
#endif

	/* re-form design matrix */
	memcpy(x, x_bk, n_col_a_temp*n_ind*sizeof(double));
  
	/* calculate fitted values, probs, new wts, new z's */
	llik = 0.0;
	for(j=0; j<n_ind; j++) {
	  nu[j] = 0.0;
	  for(k=0; k<n_col_a_temp; k++) 
	    nu[j] += x[k*n_ind+j] * ests[k];
	  pi[j] = exp(nu[j]);
	  pi[j] /= (1.0 + pi[j]);
	  wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	  z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	  llik += (pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]));

	  /* multiply design matrix by new weights */
	  for(k=0; k<n_col_a_temp; k++) 
	    x[k*n_ind+j] *= wt[j];
	}
	if(verbose>2) Rprintf("add: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, s+1, llik);
  
	if(fabs(llik - curllik) < tol) { /* converged? */
	  flag = 1;
	  break;
	}
	curllik = llik;
	
      } /* end of IRLS iterations */
  
      if(!flag)
	warning("Didn't converge.");

      Result[i2][i] = -llik;

      /* INTERACTIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob[k][i][j];  /* s keeps track of column */

	for(k=0; k<n_gen-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob[k][i2][j];

	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];

	for(k2=0; k2<n_intcov; k2++) {
	  for(k=0; k<n_gen-1; k++,s++) /* interactive x QTL 1 */
	    x[j+s*n_ind] = Genoprob[k][i][j]*Intcov[k2][j];

	  for(k=0; k<n_gen-1; k++,s++) /* interactive x QTL 2 */
	    x[j+s*n_ind] = Genoprob[k][i2][j]*Intcov[k2][j];
	}

	for(k=0; k<n_gen-1; k++)
	  for(k2=0; k2<n_gen-1; k2++,s++) /* QTL 1 x QTL 2 */
	    x[j+s*n_ind] = Pairprob[k][k2][i][i2][j];

	for(k3=0; k3<n_intcov; k3++)
	  for(k=0; k<n_gen-1; k++) /* interactive x QTL 1 x QTL 2 */
	    for(k2=0; k2<n_gen-1; k2++,s++) 
	      x[j+s*n_ind] = Pairprob[k][k2][i][i2][j]*Intcov[k3][j];
      }

      /* drop x's */
      n_col_f_temp = n_col_f;
      if(n_col2drop) 
	dropcol_x(&n_col_f_temp, n_ind, allcol2drop, x);

      /* linear regression of phenotype on QTL genotype probabilities */
      /* make a copy of x matrix, we may need it */
      memcpy(x_bk, x, n_ind*n_col_f_temp*sizeof(double));

      /* starting point for IRLS */
      curllik = 0.0;
      for(j=0; j<n_ind; j++) {
	pi[j] = (pheno[j] + 0.5)/2.0;
	wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	nu[j] = log(pi[j]) - log(1.0-pi[j]);
	z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	curllik += pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]);
      }
      if(verbose>2) Rprintf("\nint: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, 0, curllik);

      /* multiply design matrix by current wts */
      for(k=0; k<n_col_f_temp; k++) 
	for(j=0; j<n_ind; j++)
	  x[k*n_ind + j] *= wt[j];

      flag = 0;
      for(s=0; s<maxit; s++) { /* IRLS iterations */
	R_CheckUserInterrupt(); /* check for ^C */

        /* make jpvt = numbers 0, 1, ..., (n_col_f_temp-1) */
        /*      jpvt keeps track of any permutation of X columns */
        for(k=0; k<n_col_f_temp; k++) jpvt[k] = k;
  
        /* call dqrls to fit regression model */
        F77_CALL(dqrls)(x, &n_ind, &n_col_f_temp, z, &ny, &tol2, coef, resid,
			qty, &kk, jpvt, qraux, dwork);
  
	/* get ests; need to permute back */
	for(k=0; k<kk; k++) ests[jpvt[k]] = coef[k];
	for(k=kk; k<n_col_f_temp; k++) ests[jpvt[k]] = 0.0;

#ifdef UNDEFINED
	if(kk < n_col_a_temp) {
	  Rprintf("\nkk=%d n_col_a_temp=%d\n", kk, n_col_a_temp);
	  for(k=0; k<n_col_a_temp; k++) 
	    Rprintf("%10.5lf ", ests[k]);
	  Rprintf("\n");
	}	
#endif

	/* re-form design matrix */
	memcpy(x, x_bk, n_col_f_temp*n_ind*sizeof(double));
  
	/* calculate fitted values, probs, new wts, new z's */
	llik = 0.0;
	for(j=0; j<n_ind; j++) {
	  nu[j] = 0.0;
	  for(k=0; k<n_col_f_temp; k++) 
	    nu[j] += x[k*n_ind+j] * ests[k];
	  pi[j] = exp(nu[j]);
	  pi[j] /= (1.0 + pi[j]);
	  wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	  z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	  llik += (pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]));

	  /* multiply design matrix by new weights */
	  for(k=0; k<n_col_f_temp; k++) 
	    x[k*n_ind+j] *= wt[j];
	}
	if(verbose>2) Rprintf("int: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, s+1, llik);
  
	if(fabs(llik - curllik) < tol) { /* converged? */
	  flag = 1;
	  break;
	}
	curllik = llik;
	
      } /* end of IRLS iterations */
  
      if(!flag)
	warning("Didn't converge.");

      Result[i][i2] = -llik;

    } /* end loop over positions */
  } 
}

/**********************************************************************
 * 
 * R_scantwo_2chr_binary_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_binary_hk.
 * 
 **********************************************************************/

void R_scantwo_2chr_binary_hk(int *n_ind, int *n_pos1, int *n_pos2, 
			      int *n_gen1, int *n_gen2,
			      double *genoprob1, double *genoprob2,
			      double *addcov, int *n_addcov, 
			      double *intcov, int *n_intcov, 
			      double *pheno, 
			      double *result_full, double *result_add, 
			      double *tol, int *maxit, int *verbose)
{
  double ***Genoprob1, ***Genoprob2, **Result_full, **Result_add;
  double **Addcov=0, **Intcov=0;

  reorg_genoprob(*n_ind, *n_pos1, *n_gen1, genoprob1, &Genoprob1);
  reorg_genoprob(*n_ind, *n_pos2, *n_gen2, genoprob2, &Genoprob2);
  reorg_errlod(*n_pos1, *n_pos2, result_full, &Result_full);
  reorg_errlod(*n_pos1, *n_pos2, result_add, &Result_add);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_2chr_binary_hk(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
			 Genoprob1, Genoprob2, Addcov, *n_addcov, Intcov, 
			 *n_intcov, pheno, Result_full, Result_add, *tol, 
			 *maxit, *verbose);
}

/**********************************************************************
 * 
 * scantwo_2chr_binary_hk
 *
 * Performs a 2-dimensional genome scan using the Haley-Knott 
 * regression method (regressing phenotypes on conditional genotype 
 * probabilities) for a two-QTL model with the two QTL residing on
 * the different chromosomes.
 * 
 * n_ind        Number of individuals
 *
 * n_pos1       Number of marker positions on first chromosome
 *
 * n_pos2       Number of marker positions on second chromosome
 *
 * n_gen1       Number of different genotypes for first chromosome
 *
 * n_gen2       Number of different genotypes for second chromosome
 *
 * Genoprob1    Array of conditional genotype probs for 1st chr
 *              Indexed as Genoprob[gen][pos][ind]
 *
 * Genoprob2    Array of conditional genotype probs for 2nd chr
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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos1][pos2]
 *
 * Result_add   Result matrix of size [n_pos2 x n_pos1] 
 *              containing the LODs for add've models
 *              also indexed as Result[pos2][pos1]
 *
 **********************************************************************/

void scantwo_2chr_binary_hk(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
			    int n_gen2, double ***Genoprob1, 
			    double ***Genoprob2, 
			    double **Addcov, int n_addcov, 
			    double **Intcov, int n_intcov, double *pheno, 
			    double **Result_full, double **Result_add,
			    double tol, int maxit, int verbose)
{
  int n_col_a, n_col_f, n_gen_sq;
  int *jpvt, ny=1, flag;
  int i, i2, j, k, k2, k3, kk, s;
  double *dwork, *x, *x_bk, *coef, *resid, *qty, *qraux, *ests;
  double tol2=TOL, *z, *nu, *wt, *pi, curllik, llik=0.0;

  n_gen_sq = n_gen1*n_gen2;
  /* no. param in additive QTL model */
  n_col_a = (n_gen1+n_gen2-1)+n_addcov+n_intcov*(n_gen1+n_gen2-2); 
  /* no. param full model */
  n_col_f = n_gen_sq+n_addcov+n_intcov*(n_gen_sq-1); 

  /* allocate space and set things up*/
  /* lengths: 2*n_col_f:     dwork
              n_ind*n_col_f: x, x_bk
              n_col_f:       coef, qraux, ests
	      n_ind:       resid, qty, z, nu, wt, pi */
             
  dwork = (double *)R_alloc(2*n_ind*n_col_f + n_col_f*5 + n_ind*6, sizeof(double)); 
  x = dwork + 2*n_col_f;
  x_bk = x + n_ind*n_col_f;
  coef = x_bk + n_ind*n_col_f;
  resid = coef + n_col_f;
  qty = resid + n_ind;
  qraux = qty + n_ind;
  z = qraux + n_col_f;
  nu = z + n_ind;
  wt = nu + n_ind;
  pi = wt + n_ind;
  ests = pi + n_ind; /* length n_col_f */

  jpvt = (int *)R_alloc(n_col_f, sizeof(int));

  for(i=0; i<n_pos1; i++) { 
    for(i2=0; i2<n_pos2; i2++) { /* loop over pairs of positions */

      R_CheckUserInterrupt(); /* check for ^C */

      /* ADDITIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen1; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob1[k][i][j];  /* s keeps track of column */
	for(k=0; k<n_gen2-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob2[k][i2][j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k=0; k<n_gen1-1; k++) /* interactive x QTL 1 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Intcov[k2][j];
	for(k=0; k<n_gen2-1; k++) /* interactive x QTL 2 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob2[k][i2][j]*Intcov[k2][j];
      }
      /* linear regression of phenotype on QTL genotype probabilities */
      /* make a copy of x matrix, we may need it */
      memcpy(x_bk, x, n_ind*n_col_a*sizeof(double));

      /* starting point for IRLS */
      curllik = 0.0;
      for(j=0; j<n_ind; j++) {
	pi[j] = (pheno[j] + 0.5)/2.0;
	wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	nu[j] = log(pi[j]) - log(1.0-pi[j]);
	z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	curllik += pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]);
      }
      if(verbose>2) Rprintf("\nadd: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, 0, curllik);

      /* multiply design matrix by current wts */
      for(k=0; k<n_col_a; k++) 
	for(j=0; j<n_ind; j++)
	  x[k*n_ind + j] *= wt[j];

      flag = 0;
      for(s=0; s<maxit; s++) { /* IRLS iterations */
	R_CheckUserInterrupt(); /* check for ^C */

        /* make jpvt = numbers 0, 1, ..., (n_col_a-1) */
        /*      jpvt keeps track of any permutation of X columns */
        for(k=0; k<n_col_a; k++) jpvt[k] = k;
  
        /* call dqrls to fit regression model */
        F77_CALL(dqrls)(x, &n_ind, &n_col_a, z, &ny, &tol2, coef, resid,
			qty, &kk, jpvt, qraux, dwork);
  
	/* get ests; need to permute back */
	for(k=0; k<kk; k++) ests[jpvt[k]] = coef[k];
	for(k=kk; k<n_col_a; k++) ests[jpvt[k]] = 0.0;

#ifdef UNDEFINED
	if(kk < n_col_a_temp) {
	  Rprintf("\nkk=%d n_col_a_temp=%d\n", kk, n_col_a_temp);
	  for(k=0; k<n_col_a_temp; k++) 
	    Rprintf("%10.5lf ", ests[k]);
	  Rprintf("\n");
	}	
#endif

	/* re-form design matrix */
	memcpy(x, x_bk, n_col_a*n_ind*sizeof(double));
  
	/* calculate fitted values, probs, new wts, new z's */
	llik = 0.0;
	for(j=0; j<n_ind; j++) {
	  nu[j] = 0.0;
	  for(k=0; k<n_col_a; k++) 
	    nu[j] += x[k*n_ind+j] * ests[k];
	  pi[j] = exp(nu[j]);
	  pi[j] /= (1.0 + pi[j]);
	  wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	  z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	  llik += (pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]));

	  /* multiply design matrix by new weights */
	  for(k=0; k<n_col_a; k++) 
	    x[k*n_ind+j] *= wt[j];
	}
	if(verbose>2) Rprintf("add: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, s+1, llik);
  
	if(fabs(llik - curllik) < tol) { /* converged? */
	  flag = 1;
	  break;
	}
	curllik = llik;
	
      } /* end of IRLS iterations */

      if(!flag)
	warning("Didn't converge.");

      Result_add[i2][i] = -llik;

      /* INTERACTIVE MODEL */
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen1; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob1[k][i][j];  /* s keeps track of column */
	for(k=0; k<n_gen2-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob2[k][i2][j];
	for(k=0; k<n_gen1-1; k++)
	  for(k2=0; k2<n_gen2-1; k2++,s++) /* QTL 1 x QTL 2 */
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Genoprob2[k2][i2][j];
	for(k=0; k<n_addcov; k++, s++) /* additive covariates */
	  x[j+s*n_ind] = Addcov[k][j];
	for(k=0; k<n_gen1-1; k++) /* interactive x QTL 1 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Intcov[k2][j];
	for(k=0; k<n_gen2-1; k++) /* interactive x QTL 2 */
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+s*n_ind] = Genoprob2[k][i2][j]*Intcov[k2][j];
	for(k=0; k<n_gen1-1; k++) /* interactive x QTL 1 x QTL 2 */
	  for(k2=0; k2<n_gen2-1; k2++) 
	    for(k3=0; k3<n_intcov; k3++,s++)
	      x[j+s*n_ind] = Genoprob1[k][i][j]*Genoprob2[k2][i2][j]*
		Intcov[k3][j];
      }

      /* linear regression of phenotype on QTL genotype probabilities */
      /* make a copy of x matrix, we may need it */
      memcpy(x_bk, x, n_ind*n_col_f*sizeof(double));

      /* starting point for IRLS */
      curllik = 0.0;
      for(j=0; j<n_ind; j++) {
	pi[j] = (pheno[j] + 0.5)/2.0;
	wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	nu[j] = log(pi[j]) - log(1.0-pi[j]);
	z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	curllik += pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]);
      }
      if(verbose>2) Rprintf("\nint: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, 0, curllik);

      /* multiply design matrix by current wts */
      for(k=0; k<n_col_f; k++) 
	for(j=0; j<n_ind; j++)
	  x[k*n_ind + j] *= wt[j];

      flag = 0;
      for(s=0; s<maxit; s++) { /* IRLS iterations */
	R_CheckUserInterrupt(); /* check for ^C */

        /* make jpvt = numbers 0, 1, ..., (n_col_f-1) */
        /*      jpvt keeps track of any permutation of X columns */
        for(k=0; k<n_col_f; k++) jpvt[k] = k;
  
        /* call dqrls to fit regression model */
        F77_CALL(dqrls)(x, &n_ind, &n_col_f, z, &ny, &tol2, coef, resid,
			qty, &kk, jpvt, qraux, dwork);
  
	/* get ests; need to permute back */
	for(k=0; k<kk; k++) ests[jpvt[k]] = coef[k];
	for(k=kk; k<n_col_f; k++) ests[jpvt[k]] = 0.0;

#ifdef UNDEFINED
	if(kk < n_col_a_temp) {
	  Rprintf("\nkk=%d n_col_a_temp=%d\n", kk, n_col_a_temp);
	  for(k=0; k<n_col_a_temp; k++) 
	    Rprintf("%10.5lf ", ests[k]);
	  Rprintf("\n");
	}	
#endif

	/* re-form design matrix */
	memcpy(x, x_bk, n_col_f*n_ind*sizeof(double));
  
	/* calculate fitted values, probs, new wts, new z's */
	llik = 0.0;
	for(j=0; j<n_ind; j++) {
	  nu[j] = 0.0;
	  for(k=0; k<n_col_f; k++) 
	    nu[j] += x[k*n_ind+j] * ests[k];
	  pi[j] = exp(nu[j]);
	  pi[j] /= (1.0 + pi[j]);
	  wt[j] = sqrt(pi[j] * (1.0-pi[j]));
	  z[j] = nu[j]*wt[j] + (pheno[j] - pi[j])/wt[j];
	  llik += (pheno[j] * log10(pi[j]) + (1.0-pheno[j]) * log10(1.0 - pi[j]));

  
	  /* multiply design matrix by new weights */
	  for(k=0; k<n_col_f; k++) 
	    x[k*n_ind+j] *= wt[j];
	}
	if(verbose>2) Rprintf("int: %-4d %-4d : %-4d %-10.5lf\n", i+1, i2+1, s+1, llik);
  
	if(fabs(llik - curllik) < tol) { /* converged? */
	  flag = 1;
	  break;
	}
	curllik = llik;
	
      } /* end of IRLS iterations */
  
      if(!flag)
	warning("Didn't converge.");

      Result_full[i2][i] = -llik;

    } /* end loop over positions */
  } 
}

/* end of scantwo_binary_hk.c */
