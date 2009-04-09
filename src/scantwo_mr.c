/**********************************************************************
 * 
 * scantwo_mr.c
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
 * These functions are for performing a two-dimensional genome scan with  
 * a two-QTL model by marker regression.  Individuals missing genotypes 
 * at either of a pair of markers are dropped.  
 *
 * Contains: R_scantwo_1chr_mr, scantwo_1chr_mr, 
 *           R_scantwo_2chr_mr, scantwo_2chr_mr
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
#include "scantwo_mr.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scantwo_1chr_mr
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_mr.
 * 
 **********************************************************************/

void R_scantwo_1chr_mr(int *n_ind, int *n_pos, int *n_gen, int *geno,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights, double *result,
		       int *n_col2drop, int *col2drop)
{
  int **Geno;
  double **Result, **Addcov, **Intcov;

  reorg_geno(*n_ind, *n_pos, geno, &Geno);
  reorg_errlod(*n_pos, *n_pos, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_1chr_mr(*n_ind, *n_pos, *n_gen, Geno, Addcov, *n_addcov, 
		  Intcov, *n_intcov, pheno, weights, Result,
		  *n_col2drop, col2drop);
}

/**********************************************************************
 * 
 * scantwo_1chr_mr
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
 * Geno         Array of marker genotype data, indexed as
 *              Geno[pos][ind]
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
 * Result       Result matrix of size [n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              additve models.
 *              Note: indexed as Result[col][row]
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_1chr_mr(int n_ind, int n_pos, int n_gen, int **Geno,
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights, double **Result,
		     int n_col2drop, int *col2drop)
{
  int ny, *jpvt, i, i2, j, k, s, this_n_ind, done_allind=0;
  int n_col_0, n_col_a, n_col_f, n_col_a_temp, n_col_f_temp, n_gen_sq, *which_ind;
  double *work, *x, *qty, *qraux, *coef, *resid, tol, lrss0, *y;
  double lrss0_allind=0.0;
  int *allcol2drop;

  /* tolerance for linear regression */
  tol = TOL;

  n_gen_sq = n_gen*n_gen;
  /* no. param in null model */
  n_col_0 = n_addcov+1; 
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
  which_ind = (int *)R_alloc(n_ind, sizeof(int));
  y = (double *)R_alloc(n_ind, sizeof(double));
  x = (double *)R_alloc(n_ind*n_col_f, sizeof(double));
  coef = (double *)R_alloc(n_col_f, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(n_col_f, sizeof(int));
  qraux = (double *)R_alloc(n_col_f, sizeof(double));
  work = (double *)R_alloc(2 * n_col_f, sizeof(double));
  ny = 1;

  /* modify pheno, Addcov and Intcov with weights */
  for(j=0; j<n_ind; j++) {
    pheno[j] *= weights[j];
    for(k=0; k<n_addcov; k++) Addcov[k][j] *= weights[j];
    for(k=0; k<n_intcov; k++) Intcov[k][j] *= weights[j];
  }    

  for(i=0; i<n_pos-1; i++) { 
    for(i2=i+1; i2<n_pos; i2++) { /* loop over pairs of positions */

      R_CheckUserInterrupt(); /* check for ^C */

      /* genotyped individuals at this marker */
      for(j=0, this_n_ind=0; j<n_ind; j++) {
	if(Geno[i][j] > 0 && Geno[i2][j] > 0) {
	  which_ind[this_n_ind] = j;
	  y[this_n_ind] = pheno[j];
	  this_n_ind++;
	}
      }

      if(this_n_ind > 0) {

	if((this_n_ind < n_ind) || !done_allind) {
	  /* the above is to avoid repeatedly doing the null model
	     regression in the case of complete marker data */
	  
	  /* NULL MODEL */
	  /* fill up X matrix */
	  for(j=0; j<this_n_ind; j++) {
	    x[j] = weights[which_ind[j]];
	    for(k=0; k<n_addcov; k++) 
	      x[j+(k+1)*this_n_ind] = Addcov[k][which_ind[j]];
	  }
	  /* linear regression */
	  F77_CALL(dqrls)(x, &this_n_ind, &n_col_0, y, &ny, &tol, 
			  coef, resid, qty, &k, jpvt, qraux, work);
	  /* RSS */
	  lrss0 = 0.0;
	  for(j=0; j<this_n_ind; j++)  lrss0 += (resid[j]*resid[j]);
	  lrss0 = log10(lrss0);
	  
	  if(this_n_ind==n_ind) {
	    done_allind=1;
	    lrss0_allind = lrss0;
	  }
	}
	else /* already ran null model with all individuals */
	  lrss0 = lrss0_allind;
	
	/* ADDITIVE MODEL */
	/* zero out matrix */
	for(j=0; j<this_n_ind*n_col_a; j++) x[j] = 0.0;
	
	/* fill up X matrix */
	for(j=0; j<this_n_ind; j++) { 
	  x[j+(Geno[i][which_ind[j]]-1)*this_n_ind] = 
	    weights[which_ind[j]]; /* QTL 1 */
	  s = n_gen;
	  if(Geno[i2][which_ind[j]] < n_gen) /* QTL 2 */
	    x[j+(Geno[i2][which_ind[j]]-1+s)*this_n_ind] = weights[which_ind[j]];
	  s += (n_gen-1);
	  for(k=0; k<n_addcov; k++) /* additive covariates */
	    x[j+(k+s)*this_n_ind] = Addcov[k][which_ind[j]];
	  s += n_addcov;
	  for(k=0; k<n_intcov; k++) {
	    if(Geno[i][which_ind[j]] < n_gen) /* interactive x QTL 1 */
	      x[j+(s+Geno[i][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen-1);
	    if(Geno[i2][which_ind[j]] < n_gen) /* interactive x QTL 2 */
	      x[j+(s+Geno[i2][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen-1);
	  }
	}
	
	n_col_a_temp = n_col_a;
	if(n_col2drop) 
	  dropcol_x(&n_col_a_temp, this_n_ind, allcol2drop, x);

	/* linear regression of phenotype on QTL genotype probabilities */
	F77_CALL(dqrls)(x, &this_n_ind, &n_col_a_temp, y, &ny, &tol, 
			coef, resid, qty, &k, jpvt, qraux, work);

	/* RSS */
	Result[i2][i] = 0.0;
	for(j=0; j<this_n_ind; j++) Result[i2][i] += (resid[j]*resid[j]);
	Result[i2][i] = log10(Result[i2][i]); /* take log base 10*/
	
	R_CheckUserInterrupt(); /* check for ^C */

	/* INTERACTIVE MODEL */
	/* zero out matrix */
	for(j=0; j<this_n_ind*n_col_f; j++) x[j] = 0.0;
	
	/* fill up X matrix */
	for(j=0; j<this_n_ind; j++) { 
	  x[j+(Geno[i][which_ind[j]]-1)*this_n_ind] = 
	    weights[which_ind[j]]; /* QTL 1 */
	  s = n_gen;

	  if(Geno[i2][which_ind[j]] < n_gen) /* QTL 2 */
	    x[j+(Geno[i2][which_ind[j]]-1+s)*this_n_ind] = 
	      weights[which_ind[j]];
	  s += (n_gen-1);

	  for(k=0; k<n_addcov; k++) /* additive covariates */
	    x[j+(k+s)*this_n_ind] = Addcov[k][which_ind[j]];
	  s += n_addcov;

	  for(k=0; k<n_intcov; k++) {
	    if(Geno[i][which_ind[j]] < n_gen) /* interactive x QTL 1 */
	      x[j+(s+Geno[i][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen-1);
	    if(Geno[i2][which_ind[j]] < n_gen) /* interactive x QTL 2 */
	      x[j+(s+Geno[i2][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen-1);
	  }

	  if(Geno[i][which_ind[j]] < n_gen && 
	     Geno[i2][which_ind[j]] < n_gen) /* QTL x QTL */
	    x[j+((Geno[i][which_ind[j]]-1)*(n_gen-1)+s+
		 Geno[i2][which_ind[j]]-1)*this_n_ind] = weights[which_ind[j]];
	  s += (n_gen-1)*(n_gen-1);

	  for(k=0; k<n_intcov; k++) { /* interactive x QTL x QTL */
	    if(Geno[i][which_ind[j]] < n_gen && 
	       Geno[i2][which_ind[j]] < n_gen) /* QTL x QTL */
	      x[j+((Geno[i][which_ind[j]]-1)*(n_gen-1)+s+
		   Geno[i2][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen-1)*(n_gen-1);
	  }
	}

	/* drop x's */
	n_col_f_temp = n_col_f;
	if(n_col2drop) 
	  dropcol_x(&n_col_f_temp, this_n_ind, allcol2drop, x);

	/* linear regression of phenotype on QTL genotype probabilities */
	F77_CALL(dqrls)(x, &this_n_ind, &n_col_f_temp, y, &ny, &tol, 
			coef, resid, qty, &k, jpvt, qraux, work);
	/* RSS */
	Result[i][i2] = 0.0;
	for(j=0; j<this_n_ind; j++) Result[i][i2] += (resid[j]*resid[j]);
	Result[i][i2] = log10(Result[i][i2]); /* take log base 10*/
	
	/* convert to LODs */
	/* interactive LOD */
	Result[i2][i] = (double)this_n_ind/2.0*(lrss0-Result[i2][i]);
	/* joint LOD */
	Result[i][i2] = (double)this_n_ind/2.0*(lrss0-Result[i][i2]);
	
      } /* > 0 individuals with available data */
    } /* end loop over positions */
  } 
}

/**********************************************************************
 * 
 * R_scantwo_2chr_mr
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_mr.
 * 
 **********************************************************************/

void R_scantwo_2chr_mr(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2,
		       int *geno1, int *geno2,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights,
		       double *result_full, double *result_add)
{
  int **Geno1, **Geno2;
  double **Result_full, **Result_add, **Addcov, **Intcov;

  reorg_geno(*n_ind, *n_pos1, geno1, &Geno1);
  reorg_geno(*n_ind, *n_pos2, geno2, &Geno2);
  reorg_errlod(*n_pos1, *n_pos2, result_full, &Result_full);
  reorg_errlod(*n_pos1, *n_pos2, result_add, &Result_add);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_2chr_mr(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
		  Geno1, Geno2, Addcov, *n_addcov, Intcov, 
		  *n_intcov, pheno, weights, Result_full, Result_add);
}

/**********************************************************************
 * 
 * scantwo_2chr_mr
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
 * Geno1        Matrix of marker genotype data for chr 1,
 *              indexed as Geno1[pos][ind]
 *
 * Geno2        Matrix of marker genotype data for chr 2
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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos2][pos1]
 *
 * Result_add   Result matrix of size [n_pos2 x n_pos1] 
 *              containing the LODs for add've models
 *              also indexed as Result[pos2][pos1]
 *
 **********************************************************************/

void scantwo_2chr_mr(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, int **Geno1, int **Geno2,
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights,
		     double **Result_full, double **Result_add)
{
  int ny, *jpvt, i, i2, j, k, s, this_n_ind, done_allind=0;
  int n_col_0, n_col_a, n_col_f, n_gen_sq, *which_ind;
  double *work, *x, *qty, *qraux, *coef, *resid, tol, lrss0, *y;
  double lrss0_allind=0.0;

  /* tolerance for linear regression */
  tol = TOL;

  n_gen_sq = n_gen1*n_gen2;
  /* no. param in null model */
  n_col_0 = n_addcov+1;
  /* no. param in additive QTL model */
  n_col_a = (n_gen1+n_gen2-1)+n_addcov+n_intcov*(n_gen1+n_gen2-2);
  /* no. param full model */
  n_col_f = n_gen_sq+n_addcov+n_intcov*(n_gen_sq-1);

  /* allocate space and set things up*/
  which_ind = (int *)R_alloc(n_ind, sizeof(int));
  y = (double *)R_alloc(n_ind, sizeof(double));
  x = (double *)R_alloc(n_ind*n_col_f, sizeof(double));
  coef = (double *)R_alloc(n_col_f, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(n_col_f, sizeof(int));
  qraux = (double *)R_alloc(n_col_f, sizeof(double));
  work = (double *)R_alloc(2 * n_col_f, sizeof(double));
  ny = 1;

  /* modify pheno, Addcov and Intcov with weights */
  for(j=0; j<n_ind; j++) {
    pheno[j] *= weights[j];
    for(k=0; k<n_addcov; k++) Addcov[k][j] *= weights[j];
    for(k=0; k<n_intcov; k++) Intcov[k][j] *= weights[j];
  }    

  for(i=0; i<n_pos1; i++) { 
    for(i2=0; i2<n_pos2; i2++) { /* loop over pairs of positions */
      R_CheckUserInterrupt(); /* check for ^C */

      /* genotyped individuals at this marker */
      for(j=0, this_n_ind=0; j<n_ind; j++) {
	if(Geno1[i][j] > 0 && Geno2[i2][j] > 0) {
	  which_ind[this_n_ind] = j;
	  y[this_n_ind] = pheno[j];
	  this_n_ind++;
	}
      }

      if(this_n_ind > 0) {
	if((this_n_ind < n_ind) || !done_allind) {
	  /* the above is to avoid repeatedly doing the null model
	     regression in the case of complete marker data */
      
	  /* NULL MODEL */
	  /* fill up X matrix */
	  for(j=0; j<this_n_ind; j++) {
	    x[j] = weights[which_ind[j]];
	    for(k=0; k<n_addcov; k++) 
	      x[j+(k+1)*this_n_ind] = Addcov[k][which_ind[j]];
	  }
	  /* linear regression */
	  F77_CALL(dqrls)(x, &this_n_ind, &n_col_0, y, &ny, &tol, 
			  coef, resid, qty, &k, jpvt, qraux, work);
	  /* RSS */
	  lrss0 = 0.0;
	  for(j=0; j<this_n_ind; j++)  lrss0 += (resid[j]*resid[j]);
	  lrss0 = log10(lrss0);
	  
	  if(this_n_ind==n_ind) {
	    done_allind=1;
	    lrss0_allind = lrss0;
	  }
	}
	else /* already ran null model with all individuals */
	  lrss0 = lrss0_allind;
	
	/* ADDITIVE MODEL */
	/* zero out matrix */
	for(j=0; j<this_n_ind*n_col_a; j++) x[j] = 0.0;
	
	/* fill up X matrix */
	for(j=0; j<this_n_ind; j++) { 
	  x[j+(Geno1[i][which_ind[j]]-1)*this_n_ind] = 
	    weights[which_ind[j]]; /* QTL 1 */
	  s = n_gen1;
	  if(Geno2[i2][which_ind[j]] < n_gen2) /* QTL 2 */
	    x[j+(Geno2[i2][which_ind[j]]-1+s)*this_n_ind] = 
	      weights[which_ind[j]];
	  s += (n_gen2-1);
	  for(k=0; k<n_addcov; k++) /* additive covariates */
	    x[j+(k+s)*this_n_ind] = Addcov[k][which_ind[j]];
	  s += n_addcov;
	  for(k=0; k<n_intcov; k++) {
	    if(Geno1[i][which_ind[j]] < n_gen1) /* interactive x QTL 1 */
	      x[j+(s+Geno1[i][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen1-1);
	    if(Geno2[i2][which_ind[j]] < n_gen2) /* interactive x QTL 2 */
	      x[j+(s+Geno2[i2][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen2-1);
	  }
	}
	/* linear regression of phenotype on QTL genotype probabilities */
	F77_CALL(dqrls)(x, &this_n_ind, &n_col_a, y, &ny, &tol, 
			coef, resid, qty, &k, jpvt, qraux, work);
	/* RSS */
	Result_add[i2][i] = 0.0;
	for(j=0; j<this_n_ind; j++) Result_add[i2][i] += (resid[j]*resid[j]);
	Result_add[i2][i] = log10(Result_add[i2][i]); /* take log base 10*/
	
	R_CheckUserInterrupt(); /* check for ^C */

	/* INTERACTIVE MODEL */
	/* zero out matrix */
	for(j=0; j<this_n_ind*n_col_f; j++) x[j] = 0.0;
	
	/* fill up X matrix */
	for(j=0; j<this_n_ind; j++) { 
	  x[j+(Geno1[i][which_ind[j]]-1)*this_n_ind] = 
	    weights[which_ind[j]]; /* QTL 1 */
	  s = n_gen1;
	  if(Geno2[i2][which_ind[j]] < n_gen2) /* QTL 2 */
	    x[j+(Geno2[i2][which_ind[j]]-1+s)*this_n_ind] = 
	      weights[which_ind[j]];
	  s += (n_gen2-1);
	  if(Geno1[i][which_ind[j]] < n_gen1 && 
	     Geno2[i2][which_ind[j]] < n_gen2) /* QTL x QTL */
	    x[j+((Geno1[i][which_ind[j]]-1)*(n_gen2-1)+s+
		 Geno2[i2][which_ind[j]]-1)*this_n_ind] = 
	      weights[which_ind[j]];
	  s += (n_gen1-1)*(n_gen2-1);
	  for(k=0; k<n_addcov; k++) /* additive covariates */
	    x[j+(k+s)*this_n_ind] = Addcov[k][which_ind[j]];
	  s += n_addcov;
	  for(k=0; k<n_intcov; k++) {
	    if(Geno1[i][which_ind[j]] < n_gen1) /* interactive x QTL 1 */
	      x[j+(s+Geno1[i][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen1-1);
	    if(Geno2[i2][which_ind[j]] < n_gen2) /* interactive x QTL 2 */
	      x[j+(s+Geno2[i2][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen2-1);
	    if(Geno1[i][which_ind[j]] < n_gen1 && 
	       Geno2[i2][which_ind[j]] < n_gen2) /* QTL x QTL */
	      x[j+((Geno1[i][which_ind[j]]-1)*(n_gen2-1)+s+
		   Geno2[i2][which_ind[j]]-1)*this_n_ind] = 
		Intcov[k][which_ind[j]];
	    s += (n_gen1-1)*(n_gen2-1);
	  }
	}
	/* linear regression of phenotype on QTL genotype probabilities */
	F77_CALL(dqrls)(x, &this_n_ind, &n_col_f, y, &ny, &tol, 
			coef, resid, qty, &k, jpvt, qraux, work);
	/* RSS */
	Result_full[i2][i] = 0.0;
	for(j=0; j<this_n_ind; j++) Result_full[i2][i] += (resid[j]*resid[j]);
	Result_full[i2][i] = log10(Result_full[i2][i]); /* take log base 10*/

	/* convert to LODs */
	Result_add[i2][i] = (double)this_n_ind/2.0*
	  (lrss0 - Result_add[i2][i]);
	Result_full[i2][i] = (double)this_n_ind/2.0*(lrss0 - Result_full[i2][i]);

      } /* > 0 individuals with available data */
    } /* end loop over positions */
  } 
}

/* end of scantwo_mr.c */
