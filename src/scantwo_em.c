/**********************************************************************
 * 
 * scantwo_em.c
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
 * These functions are for performing a 2-dimensional genome scan 
 * with a 2-QTL model by interval mapping.(the EM algorithm).
 *
 * Contains: R_scantwo_1chr_em, scantwo_1chr_em, 
 *           R_scantwo_2chr_em, scantwo_2chr_em,
 *           scantwo_em_estep, scantwo_em_mstep, scantwo_em_loglik
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
#include "scantwo_em.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scantwo_1chr_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_em.
 * 
 **********************************************************************/

void R_scantwo_1chr_em(int *n_ind, int *n_pos, int *n_gen,
		       double *pairprob, double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights, double *result,
		       int *maxit, double *tol, int *verbose,
		       int *n_col2drop, int *col2drop)
{
  double **Result, **Addcov, **Intcov, *****Pairprob;

  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
  reorg_errlod(*n_pos, *n_pos, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_1chr_em(*n_ind, *n_pos, *n_gen, Pairprob, 
		  Addcov, *n_addcov, Intcov, *n_intcov, 
		  pheno, weights, Result, *maxit, *tol, *verbose,
		  *n_col2drop, col2drop);
}

/**********************************************************************
 * 
 * scantwo_1chr_em
 *
 * Performs a 2-dimensional genome scan using the EM algorithm
 * for a two-QTL model with the two QTL residing on the same 
 * chromosome.
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
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
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              additive models.
 *              Note: indexed as Result[col][row]
 *
 * maxit        Maximum number of iterations for EM
 *
 * tol          Tolerance for determining convergence of EM
 *
 * verbose        If >0, print any messages when errors occur
 *                 >1, print out log likelihoods at end of EM
 *                     and check that log likelihood doesn't go down
 *                 >2, print out initial and final log likelihoods
 *                 >3, print out log likelihood at each iteration
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_1chr_em(int n_ind, int n_pos, int n_gen, 
		     double *****Pairprob, double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights,
		     double **Result, int maxit, double tol, int verbose,
		     int n_col2drop, int *col2drop)
{
  int error_flag, i, i1, i2, k1, k2, j, m, n_col[2], n_col_rev[2], nit[2], r, flag=0;
  double *param, *oldparam, ***Wts12, **Wts1, **Wts2;
  double *wts, *work1, *work2, temp, ***Probs, oldllik=0.0, llik[2], sw;
  int *allcol2drop;

  n_col[0] = (2*n_gen-1) + n_addcov + 2*(n_gen-1)*n_intcov;
  n_col[1] = n_gen*n_gen + n_addcov + (n_gen*n_gen-1)*n_intcov;

  /* expand col2drop */
  if(n_col2drop) {
    allocate_int(n_col[1], &allcol2drop);
    expand_col2drop(n_gen, n_addcov, n_intcov, 
		    col2drop, allcol2drop);
  }

  /* revised numbers of parameters */
  if(n_col2drop) {
    n_col_rev[0] = 0;
    for(i=0; i<n_col[0]; i++) 
      if(!allcol2drop[i]) n_col_rev[0]++;
    n_col_rev[1] = n_col_rev[0];
    for(i=n_col[0]; i<n_col[1]; i++)
      if(!allcol2drop[i]) n_col_rev[1]++;
  }
  else {
    n_col_rev[0] = n_col[0];
    n_col_rev[1] = n_col[1];
  }

  /* allocate workspaces */
  wts = (double *)R_alloc(2*n_gen*(n_gen+1)*n_ind, sizeof(double));
  reorg_errlod(n_ind, n_gen, wts, &Wts1);
  reorg_errlod(n_ind, n_gen, wts+n_gen*n_ind, &Wts2);
  reorg_genoprob(n_ind, n_gen, n_gen, wts+2*n_gen*n_ind, &Wts12);
  reorg_genoprob(n_ind, n_gen, n_gen, wts+n_gen*(n_gen+2)*n_ind, &Probs);
  work1 = (double *)R_alloc(n_col[1]*n_col[1], sizeof(double));
  work2 = (double *)R_alloc(n_col[1], sizeof(double));
  param = (double *)R_alloc(n_col[1]+1, sizeof(double));
  oldparam = (double *)R_alloc(n_col[1]+1, sizeof(double));
  
  /* recenter phenotype to have mean 0, for 
     possibly increased numerical stability */
  for(j=0, temp=0.0; j<n_ind; j++) temp += pheno[j];
  temp /= (double)n_ind;
  for(j=0; j<n_ind; j++) pheno[j] -= temp;

  /* adjust phenotypes and covariates with weights */
  /* Note: weights are actually sqrt(weights) */
  sw = 0.0;
  for(i=0; i<n_ind; i++) {
    pheno[i] *= weights[i];
    for(j=0; j<n_addcov; j++) Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++) Intcov[j][i] *= weights[i];
    sw += log(weights[i]);   /* sum of log weights */
  }
  sw /= log(10.0); /* make log 10 */

  /* begin loop over pairs of positions */
  for(i1=0; i1<n_pos-1; i1++) {
    for(i2=i1+1; i2<n_pos; i2++) { /* loop over positions */
      nit[0] = nit[1] = 0;
      llik[0] = llik[1] = NA_REAL;

      /* copy the parts from Pairprob into Probs */
      for(j=0; j<n_ind; j++) 
	for(k1=0; k1<n_gen; k1++)
	  for(k2=0; k2<n_gen; k2++)
	    Probs[k1][k2][j] = Pairprob[k1][k2][i1][i2][j];

      for(m=0; m<2; m++) { /* loop over add've model and full model */
	/* initial estimates */
	for(j=0; j<n_ind; j++) {
	  for(k1=0; k1<n_gen; k1++) {
	    Wts1[k1][j] = Wts2[k1][j] = 0.0;
	    for(k2=0; k2<n_gen; k2++) {
	      Wts1[k1][j] += Probs[k1][k2][j];
	      Wts2[k1][j] += Probs[k2][k1][j];
	    }
	  }
	}
	scantwo_em_mstep(n_ind, n_gen, n_gen, Addcov, n_addcov, 
			 Intcov, n_intcov, pheno, weights, Probs, Wts1, Wts2,
			 oldparam, m, work1, work2, &error_flag,
			 n_col2drop, allcol2drop, verbose);

	if(error_flag) {
	  if(verbose>1) 
	    Rprintf("   [%3d %3d] %1d: Initial model had error.\n",
		    i1+1, i2+1, m+1);
	}
	else { /* only proceed if there's no error */
	  oldllik = scantwo_em_loglik(n_ind, n_gen, n_gen, Probs, Wts12, 
				      Wts1, Wts2, Addcov, n_addcov, Intcov, 
				      n_intcov, pheno, weights, oldparam, m,
				      n_col2drop, allcol2drop);

	  if(verbose>2) 
	    Rprintf("   [%3d %3d] %1d %9.3lf\n", 
		    i1+1, i2+1, m+1, oldllik);
	
	  for(r=0; r<maxit; r++) { /* loop over iterations */
	    R_CheckUserInterrupt(); /* check for ^C */

	    scantwo_em_estep(n_ind, n_gen, n_gen, Probs, Wts12, 
			     Wts1, Wts2, Addcov, n_addcov, Intcov, 
			     n_intcov, pheno, weights, oldparam, m, 1,
			     n_col2drop, allcol2drop);

	    scantwo_em_mstep(n_ind, n_gen, n_gen, Addcov, n_addcov, 
			     Intcov, n_intcov, pheno, weights, Wts12, Wts1, 
			     Wts2, param, m, work1, work2, &error_flag,
			     n_col2drop, allcol2drop, verbose);
	    if(error_flag) {
	      flag=0;
	      if(verbose>1)
		Rprintf("   [%3d %3d] %1d %4d: Error in mstep\n",
			i1+1, i2+1, m+1, r+1);
	      break;
	    }

	    llik[m] = scantwo_em_loglik(n_ind, n_gen, n_gen, Probs, Wts12, 
					Wts1, Wts2, Addcov, n_addcov, Intcov, 
					n_intcov, pheno, weights, param, m,
					n_col2drop, allcol2drop);

	    if(verbose>1) { /* print log likelihood */
	      if(verbose>2) 
		Rprintf("   [%3d %3d] %1d %4d %9.6lf\n", 
			i1+1, i2+1, m+1, r+1, (llik[m]-oldllik));
	      if(llik[m] < oldllik-tol) 
		Rprintf("** [%3d %3d] %1d %4d %9.6lf **\n",
			i1+1, i2+1, m+1, r+1, (llik[m]-oldllik));

	      if(verbose>3) { /* print parameters */
		for(j=0; j<n_col_rev[m]; j++) 
		  Rprintf(" %7.3lf", param[j]);
		Rprintf("\n");
	      }
	    }
	    
	    flag = 1;
	    /* use log likelihood only to check for convergence */
	    if(llik[m]-oldllik < tol) { 
	      flag = 0; 
	      break; 
	    }

	    oldllik = llik[m];
	    for(j=0; j<n_col[m]+1; j++) oldparam[j] = param[j];

	  } /* loop over EM iterations */
	  nit[m] = r+1;
	  if(flag) {
	    if(verbose>1)
	      Rprintf("** [%3d %3d] %1d Didn't converge! **\n",
		      i1+1, i2+1, m+1);
	    warning("Didn't converge!\n");
	  }

	} /* no error in getting initial estimates */
      } /* loop over model */ 


      if(verbose>1) { /* print likelihoods */
	Rprintf("   [%3d %3d]   %4d %4d    %9.6lf", 
		i1+1, i2+1, nit[0], nit[1], llik[1]-llik[0]);
	if(llik[1] < llik[0]) Rprintf(" ****");
	Rprintf("\n");
      }

      Result[i2][i1] = -(llik[0]+sw);
      Result[i1][i2] = -(llik[1]+sw); /* sw = sum[log10(weights)] */

    } /* position 2 */
  } /* position 1 */
}

/**********************************************************************
 * 
 * R_scantwo_2chr_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_em.
 * 
 **********************************************************************/

void R_scantwo_2chr_em(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2, double *genoprob1,
		       double *genoprob2, double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights,
		       double *result_full, double *result_add,
		       int *maxit, double *tol, int *verbose)
{
  double **Result_full, **Result_add, **Addcov, **Intcov;
  double ***Genoprob1, ***Genoprob2;

  reorg_genoprob(*n_ind, *n_pos1, *n_gen1, genoprob1, &Genoprob1);
  reorg_genoprob(*n_ind, *n_pos2, *n_gen2, genoprob2, &Genoprob2);
  reorg_errlod(*n_pos1, *n_pos2, result_full, &Result_full);
  reorg_errlod(*n_pos1, *n_pos2, result_add, &Result_add);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_2chr_em(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2,
		  Genoprob1, Genoprob2, Addcov, *n_addcov, 
		  Intcov, *n_intcov, pheno, weights, 
		  Result_full, Result_add, 
		  *maxit, *tol, *verbose);
}

/**********************************************************************
 * 
 * scantwo_2chr_em
 *
 * Performs a 2-dimensional genome scan using the EM algorithm
 * for a two-QTL model with the two QTL residing on the same 
 * chromosome.
 * 
 * n_ind        Number of individuals
 *
 * n_pos1       Number of marker positions on chr 1
 *
 * n_pos2       Number of marker positions on chr 2
 *
 * n_gen1       Number of different genotypes on chr 1
 *
 * n_gen2       Number of different genotypes on chr 2
 *
 * Genoprob1    Array of genotype probabilities for chr 1
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Genoprob2    Array of genotype probabilities for chr 2
 *              indexed as Genoprob[gen][pos][ind]
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
 *              containing the LODs for add've model
 *              also indexed as Result[pos2][pos1]
 *
 * maxit        Maximum number of iterations for EM
 *
 * tol          Tolerance for determining convergence of EM
 *
 * verbose        If >0, print any messages when errors occur
 *                 >1, print out log likelihoods at end of EM
 *                     and check that log likelihood doesn't go down
 *                 >2, print out initial and final log likelihoods
 *                 >3, print out log likelihood at each iteration
 *
 **********************************************************************/

void scantwo_2chr_em(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, double ***Genoprob1, double ***Genoprob2,
		     double **Addcov, int n_addcov, double **Intcov, 
		     int n_intcov, double *pheno, double *weights,
		     double **Result_full, double **Result_add, 
		     int maxit, double tol, int verbose)
{
  int error_flag, i, i1, i2, k1, k2, j, m, n_col[2], nit[2], r, flag=0;
  double *param, *oldparam, ***Wts12, **Wts1, **Wts2;
  double *wts, *work1, *work2, temp, ***Probs, oldllik=0.0, llik[2], sw;
  int n_col2drop=0, *allcol2drop=0;

  n_col[0] = (n_gen1+n_gen2-1) + n_addcov + (n_gen1+n_gen2-2)*n_intcov;
  n_col[1] = n_gen1*n_gen2 + n_addcov + (n_gen1*n_gen2-1)*n_intcov;

  /* allocate workspaces */
  wts = (double *)R_alloc((2*n_gen1*n_gen2+n_gen1+n_gen2)*n_ind, sizeof(double));
  reorg_errlod(n_ind, n_gen1, wts, &Wts1);
  reorg_errlod(n_ind, n_gen2, wts+n_gen1*n_ind, &Wts2);
  reorg_genoprob(n_ind, n_gen2, n_gen1, wts+(n_gen1+n_gen2)*n_ind, &Wts12);
  reorg_genoprob(n_ind, n_gen2, n_gen1, 
		 wts+(n_gen1*n_gen2+n_gen1+n_gen2)*n_ind, &Probs);
  work1 = (double *)R_alloc(n_col[1]*n_col[1], sizeof(double));
  work2 = (double *)R_alloc(n_col[1], sizeof(double));
  param = (double *)R_alloc(n_col[1]+1, sizeof(double));
  oldparam = (double *)R_alloc(n_col[1]+1, sizeof(double));
  
  /* recenter phenotype to have mean 0, for 
     possibly increased numerical stability */
  for(j=0, temp=0.0; j<n_ind; j++) temp += pheno[j];
  temp /= (double)n_ind;
  for(j=0; j<n_ind; j++) pheno[j] -= temp;

  /* adjust phenotypes and covariates with weights */
  /* Note: weights are actually sqrt(weights) */
  sw = 0.0;
  for(i=0; i<n_ind; i++) {
    pheno[i] *= weights[i];
    for(j=0; j<n_addcov; j++) Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++) Intcov[j][i] *= weights[i];
    sw += log(weights[i]);   /* sum of log weights */
  }
  sw /= log(10.0); /* make log 10 */

  /* begin loop over pairs of positions */
  for(i1=0; i1<n_pos1; i1++) {
    for(i2=0; i2<n_pos2; i2++) { /* loop over positions */
      nit[0] = nit[1] = 0;
      llik[0] = llik[1] = NA_REAL;

      /* copy the parts from Pairprob into Probs */
      for(j=0; j<n_ind; j++) 
	for(k1=0; k1<n_gen1; k1++)
	  for(k2=0; k2<n_gen2; k2++)
	    Probs[k1][k2][j] = Genoprob1[k1][i1][j]*Genoprob2[k2][i2][j];

      for(m=0; m<2; m++) { /* loop over add've model and full model */
	/* initial estimates */
	/* marginal probabilities */
	for(j=0; j<n_ind; j++) {
	  for(k1=0; k1<n_gen1; k1++) {
	    Wts1[k1][j] = 0.0;
	    for(k2=0; k2<n_gen2; k2++) 
	      Wts1[k1][j] += Probs[k1][k2][j];
	  }
	  for(k2=0; k2<n_gen2; k2++) {
	    Wts2[k2][j] = 0.0;
	    for(k1=0; k1<n_gen1; k1++)
	      Wts2[k2][j] += Probs[k1][k2][j];
	  }
	}
	scantwo_em_mstep(n_ind, n_gen1, n_gen2, Addcov, n_addcov, 
			 Intcov, n_intcov, pheno, weights, Probs, Wts1, Wts2,
			 oldparam, m, work1, work2, &error_flag,
			 n_col2drop, allcol2drop, verbose);

	if(error_flag) {
	  if(verbose>1)
	    Rprintf("   [%3d %3d] %1d: Initial model had error.\n",
		    i1+1, i2+1, m+1);
	}
	else { /* only proceed if there's no error */
	  oldllik = scantwo_em_loglik(n_ind, n_gen1, n_gen2, Probs, Wts12, 
				      Wts1, Wts2, Addcov, n_addcov, Intcov, 
				      n_intcov, pheno, weights, oldparam, m,
				      n_col2drop, allcol2drop);

	  if(verbose>2)
	    Rprintf("   [%3d %3d] %1d %9.3lf\n", 
		    i1+1, i2+1, m+1, oldllik);
	
	  for(r=0; r<maxit; r++) { /* loop over iterations */
	    R_CheckUserInterrupt(); /* check for ^C */

	    scantwo_em_estep(n_ind, n_gen1, n_gen2, Probs, Wts12, 
			     Wts1, Wts2, Addcov, n_addcov, Intcov, 
			     n_intcov, pheno, weights, oldparam, m, 1,
			     n_col2drop, allcol2drop);

	    scantwo_em_mstep(n_ind, n_gen1, n_gen2, Addcov, n_addcov, 
			     Intcov, n_intcov, pheno, weights, Wts12, Wts1, 
			     Wts2, param, m, work1, work2, &error_flag,
			     n_col2drop, allcol2drop, verbose);
	    if(error_flag) {
	      flag=0;
	      if(verbose>1)
		Rprintf("   [%3d %3d] %1d %4d: Error in mstep\n",
			i1+1, i2+1, m+1, r+1);
	      break;
	    }

	    llik[m] = scantwo_em_loglik(n_ind, n_gen1, n_gen2, Probs, Wts12, 
					Wts1, Wts2, Addcov, n_addcov, Intcov, 
					n_intcov, pheno, weights, param, m,
					n_col2drop, allcol2drop);

	    if(verbose>1) { /* print log likelihood */
	      if(verbose>2)
		Rprintf("   [%3d %3d] %1d %4d %9.6lf\n", 
			i1+1, i2+1, m+1, r+1, (llik[m]-oldllik));
	      if(llik[m] < oldllik-tol) 
		Rprintf("** [%3d %3d] %1d %4d %9.6lf **\n",
			i1+1, i2+1, m+1, r+1, (llik[m]-oldllik));

	      if(verbose>3) { /* print parameters */
		for(j=0; j<n_col[m]; j++) 
		  Rprintf(" %7.3lf", param[j]);
		Rprintf("\n");
	      }
	    }
	    
	    flag = 1;
	    /* use log likelihood only to check for convergence */
	    if(llik[m]-oldllik < tol) { 
	      flag = 0; 
	      break; 
	    }

	    for(j=0; j<n_col[m]+1; j++) oldparam[j] = param[j];
	    oldllik = llik[m];



	  } /* loop over EM iterations */
	  nit[m] = r+1;
	  if(flag) {
	    if(verbose>1)
	      Rprintf("** [%3d %3d] %1d Didn't converge! **\n",
		      i1+1, i2+1, m+1);
	    warning("Didn't converge!\n");
	  }

	} /* no error in getting initial estimates */
      } /* loop over model */ 


      if(verbose>1) { /* print likelihoods */
	Rprintf("   [%3d %3d]    %4d %4d    %9.6lf", 
		i1+1, i2+1, nit[0], nit[1], llik[1]-llik[0]);
	if(llik[1] < llik[0]) Rprintf(" ****");
	Rprintf("\n");
      }

      Result_add[i2][i1] = -(llik[0]+sw);
      Result_full[i2][i1] = -(llik[1]+sw);  /* sw = sum[log10(weights)] */

    } /* position 2 */
  } /* position 1 */
}

/**********************************************************************
 * 
 * scantwo_em_mstep: M-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen1   Number of possible genotypes at QTL 1
 *
 * n_gen2   Number of possible genotypes at QTL 2
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
 * Wts12    Pr(QTL1=v, QTL2=w | phenotype, model, marker data),
 *          indexed as Wts[v][w][ind]
 *
 * Wts1     Marginal weights for QTL 1
 * 
 * Wts2     Marginal weights for QTL 2
 *
 * param    On output, the updated parameter estimates (incl resid SD)
 *
 * full_model   If 1, include QTLxQTL interaction
 * 
 * work1    Workspace of doubles, of length (n_par-1)*(n_par-1)
 *
 * work2    Workspace of doubles, of length (n_par-1)
 *
 * error_flag     Set to 1 if X'X is singular
 *
 **********************************************************************/

void scantwo_em_mstep(int n_ind, int n_gen1, int n_gen2, 
		      double **Addcov, int n_addcov, 
		      double **Intcov, int n_intcov, double *pheno, 
		      double *weights,
		      double ***Wts12, double **Wts1, double **Wts2,
		      double *param, int full_model,
		      double *work1, double *work2, int *error_flag,
		      int n_col2drop, int *allcol2drop, int verbose)
{
  int i, j, j2, k1, k2, s, s2, nparm1, info;
  double rcond, temp;

  *error_flag=0;

  nparm1 = n_gen1 + n_gen2-1 + n_addcov + n_intcov*(n_gen1+n_gen2-2);
  if(full_model) nparm1 += (n_gen1-1)*(n_gen2-1)*(n_intcov+1);

  /* calculate {E(X)}' y */
  for(j=0; j<nparm1; j++) work2[j] = 0.0;

  for(i=0; i<n_ind; i++) {
    for(k1=0; k1<n_gen1; k1++) /* QTL 1 */
      work2[k1] += (Wts1[k1][i]*pheno[i]*weights[i]);
    s = n_gen1;
    for(k2=0; k2<n_gen2-1; k2++) /* QTL 2 */
      work2[k2+s] += (Wts2[k2][i]*pheno[i]*weights[i]);
    s += (n_gen2-1);
    for(j=0; j<n_addcov; j++) /* add covar */
      work2[j+s] += (Addcov[j][i]*pheno[i]);
    s += n_addcov;
    for(j=0; j<n_intcov; j++) {
      for(k1=0; k1<n_gen1-1; k1++)
	work2[s+k1] += (Wts1[k1][i]*Intcov[j][i]*pheno[i]);
      s += (n_gen1-1);
      for(k2=0; k2<n_gen2-1; k2++)
	work2[s+k2] += (Wts2[k2][i]*Intcov[j][i]*pheno[i]);
      s += (n_gen2-1);
    }

    if(full_model) {
      for(k1=0; k1<n_gen1-1; k1++) 
	for(k2=0; k2<n_gen2-1; k2++) 
	  work2[s+k1*(n_gen2-1)+k2] += (Wts12[k1][k2][i]*pheno[i]*weights[i]);
      s += (n_gen1-1)*(n_gen2-1);
      for(j=0; j<n_intcov; j++) {
	for(k1=0; k1<n_gen1-1; k1++) 
	  for(k2=0; k2<n_gen2-1; k2++) 
	    work2[s+k1*(n_gen2-1)+k2] += 
	      (Wts12[k1][k2][i]*Intcov[j][i]*pheno[i]);
	s += (n_gen1-1)*(n_gen2-1);
      }
    }
  } /* end loop over individuals */

  /* calculate E{X'X}; only the upper right triangle is needed */
  for(j=0; j<nparm1*nparm1; j++) work1[j] = 0.0;
  for(i=0; i<n_ind; i++) {
    /* QTL 1 columns */
    for(k1=0; k1<n_gen1; k1++) 
      work1[k1+nparm1*k1] += Wts1[k1][i]*weights[i]*weights[i]; 

    /* QTL 2 columns */
    for(k2=0, s=n_gen1; k2<n_gen2-1; k2++, s++) {
      work1[s+nparm1*s] += Wts2[k2][i]*weights[i]*weights[i]; 
      for(k1=0; k1<n_gen1; k1++)  
	work1[k1+nparm1*s] += (Wts12[k1][k2][i]*weights[i]*weights[i]);
    }
    
    /* add covar columns */
    for(j=0; j<n_addcov; j++, s++) {
      work1[s+nparm1*s] += (Addcov[j][i]*Addcov[j][i]);
      for(k1=0; k1<n_gen1; k1++) 
	work1[k1+nparm1*s] += (Wts1[k1][i]*Addcov[j][i]*weights[i]);
      for(k2=0, s2=n_gen1; k2<n_gen2-1; k2++, s2++) 
	work1[s2+nparm1*s] += (Wts2[k2][i]*Addcov[j][i]*weights[i]);
      for(j2=0; j2<j; j2++, s2++) 
	work1[s2+nparm1*s] += (Addcov[j2][i]*Addcov[j][i]);
    }

    /* QTL x interactive covariates columns */
    for(j=0; j<n_intcov; j++) {
      /* int covar x QTL 1 */
      for(k1=0; k1<n_gen1-1; k1++, s++) { 
	work1[s+nparm1*s] += (Wts1[k1][i]*Intcov[j][i]*Intcov[j][i]);
	work1[k1+nparm1*s] += (Wts1[k1][i]*Intcov[j][i]*weights[i]); /* x QTL 1 */
	for(k2=0, s2=n_gen1; k2<n_gen2-1; k2++, s2++) /* x QTL 2 */
	  work1[s2+nparm1*s] += (Wts12[k1][k2][i]*Intcov[j][i]*weights[i]);
	for(j2=0; j2<n_addcov; j2++, s2++) /* x add covar */
	  work1[s2+nparm1*s] += (Wts1[k1][i]*Intcov[j][i]*Addcov[j2][i]);
	/* prev interactive covar */
	for(j2=0; j2<j; j2++, s2 += (n_gen1+n_gen2-2)) {
	  /* x (QTL 1 x prev inter've covar) */
	  work1[s2+k1+nparm1*s] += (Wts1[k1][i]*Intcov[j][i]*Intcov[j2][i]);
	  /* x (QTL 2 x prev inter've covar) */
	  for(k2=0; k2<n_gen2-1; k2++) 
	    work1[s2+n_gen1-1+k2+nparm1*s] += 
	      (Wts12[k1][k2][i]*Intcov[j][i]*Intcov[j2][i]);
	}
      }

      /* int covar x QTL 2 */
      for(k2=0; k2<n_gen2-1; k2++, s++) { 
	work1[s+nparm1*s] += (Wts2[k2][i]*Intcov[j][i]*Intcov[j][i]);
	/* x QTL 1 */
	for(k1=0; k1<n_gen1; k1++) 
	  work1[k1+nparm1*s] += (Wts12[k1][k2][i]*Intcov[j][i]*weights[i]);
	/* x QTL 2 */
	work1[n_gen1+k2+nparm1*s] += (Wts2[k2][i]*Intcov[j][i]*weights[i]);
	/* x add covar */
	for(j2=0, s2=n_gen1+n_gen2-1; j2<n_addcov; j2++, s2++)
	  work1[s2+nparm1*s] += (Wts2[k2][i]*Intcov[j][i]*Addcov[j2][i]);
	/* x prev interactive covar */
	for(j2=0; j2<j; j2++, s2 += (n_gen1+n_gen2-2)) {
	  /* x (QTL 1 x prev inter've covar) */
	  for(k1=0; k1<n_gen1-1; k1++) 
	    work1[s2+k1+nparm1*s] += 
	      (Wts12[k1][k2][i]*Intcov[j][i]*Intcov[j2][i]);
	  /* x (QTL 2 x prev inter've covar) */
	  work1[s2+n_gen1-1+k2+nparm1*s] += 
	    (Wts2[k2][i]*Intcov[j][i]*Intcov[j2][i]);
	}
	/* x (QTL 1 x inter've covar) */
	for(k1=0; k1<n_gen1-1; k1++, s2++) 
	  work1[s2+nparm1*s] += (Wts12[k1][k2][i]*Intcov[j][i]*Intcov[j][i]);
      }
    } /* end of interactive covariates */

    if(full_model) {
      /* QTL x QTL interactions */
      for(k1=0; k1<n_gen1-1; k1++) {
	for(k2=0; k2<n_gen2-1; k2++, s++) {
	  work1[s+nparm1*s] += Wts12[k1][k2][i]*weights[i]*weights[i];
	  /* x QTL 1 */
	  work1[k1+nparm1*s] += Wts12[k1][k2][i]*weights[i]*weights[i];
	  /* x QTL 2 */
	  work1[n_gen1+k2+nparm1*s] += Wts12[k1][k2][i]*weights[i]*weights[i];
	  /* x add covar */
	  for(j=0, s2 = n_gen1+n_gen2-1; j<n_addcov; j++, s2++) 
	    work1[s2+nparm1*s] += (Wts12[k1][k2][i]*Addcov[j][i]*weights[i]);
	  /* x interactive covariates */
	  for(j=0; j<n_intcov; j++, s2+=n_gen1+n_gen2-2) {
	    work1[s2+k1+nparm1*s] += (Wts12[k1][k2][i]*Intcov[j][i]*weights[i]);
	    work1[s2+n_gen1-1+k2+nparm1*s] +=
	      (Wts12[k1][k2][i]*Intcov[j][i]*weights[i]);
	  }
	}
      } /* end of QTL x QTL interactions */

      /* QTL x QTL x inter've covar */
      for(j=0; j<n_intcov; j++) {
	for(k1=0; k1<n_gen1-1; k1++) {
	  for(k2=0; k2<n_gen2-1; k2++, s++) {
	    temp = Wts12[k1][k2][i]*Intcov[j][i];
	    work1[s+nparm1*s] += temp*Intcov[j][i];
	    work1[k1+nparm1*s] += temp*weights[i]; /* x QTL 1 */
	    work1[n_gen1+k2+nparm1*s] += temp*weights[i]; /* x QTL 2 */
	    /* x add covar */
	    for(j2=0, s2=n_gen1+n_gen2-1; j2<n_addcov; j2++, s2++)
	      work1[s2+nparm1*s] += temp*Addcov[j2][i];
	    /* x int covar */
	    for(j2=0; j2<n_intcov; j2++, s2 += (n_gen1+n_gen2-2)) {
	      work1[s2+k1+nparm1*s] += temp*Intcov[j2][i];
	      work1[s2+n_gen1-1+k2+nparm1*s] += temp*Intcov[j2][i];
	    }
	    /* x (QTL x QTL) */
	    work1[s2+k1*(n_gen2-1)+k2+nparm1*s] += temp*weights[i];
	    s2 += (n_gen1-1)*(n_gen2-1);
	    /* x (QTL x QTL x prev covar) */
	    for(j2=0; j2<j; j2++, s2 += ((n_gen1-1)*(n_gen2-1))) 
	      work1[s2+k1*(n_gen2-1)+k2+nparm1*s] += (temp*Intcov[j2][i]);
	  }
	}
      } /* end QTL x QTL x inter've covar */

    }
  } /* end loop over individuals */
  /* done calculating E(X'X) */

  if(n_col2drop) { /* drop some columns (X chromosome) */
    dropcol_xpy(nparm1, allcol2drop, work2);

    dropcol_xpx(&nparm1, allcol2drop, work1);
  }

  /* solve work1 * beta = work2 for beta */
  F77_CALL(dpoco)(work1, &nparm1, &nparm1, &rcond, param, &info);
  if(fabs(rcond) < TOL || info != 0) { /* error! */
    if(verbose > 1) Rprintf("X'X matrix is singular.\n");
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
 * scantwo_em_estep: E-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen1   Number of possible genotypes at QTL 1
 *
 * n_gen2   Number of possible genotypes at QTL 2
 *
 * Probs    Pr(QTL1=v, QTL2=w | multipoint marker data)
 *          Indexed as Probs[v][w][ind]
 *
 * Wts12    The output:
 *          Pr(QTL1=v, QTL2=w | marker data, phenotype, covar, param)
 *          Indexed as Wts[v][w][ind]
 * 
 * Wts1     Marginal weights for QTL 1
 * 
 * Wts2     Marginal weights for QTL 2
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
 * param    Current parameter estimates (including the resid SD)
 *
 * full_model   If 1, use the full model (with QTLxQTL interaction)
 *
 * rescale  If 1, rescale weights so that the sum to 1.
 *          This is done so that by taking rescale=0, we can easily
 *          calculate the log likelihood 
 *
 **********************************************************************/

void scantwo_em_estep(int n_ind, int n_gen1, int n_gen2, 
		      double ***Probs, double ***Wts12, 
		      double **Wts1, double **Wts2,
		      double **Addcov, int n_addcov, double **Intcov,
		      int n_intcov, double *pheno, double *weights, 
		      double *param, int full_model, int rescale,
		      int n_col2drop, int *allcol2drop)
{
  int i, j, k1, k2, s, ss;
  double temp;

  for(i=0; i<n_ind; i++) {

    /* Get fitted values and put in Wts12 */
    /* additive covar effect */
    if(n_col2drop) {
      for(ss=0, s=0; ss<n_gen1+n_gen2-1; ss++)
	if(!allcol2drop[ss]) s++;
    }
    else s=n_gen1+n_gen2-1;

    temp = 0.0;
    for(j=0; j<n_addcov; j++, s++) 
      temp += (Addcov[j][i]*param[s]);
    
    /* QTL 1 effect */
    for(k1=0, ss=0, s=0; k1<n_gen1; k1++, ss++, s++) { 
      if(!n_col2drop || !allcol2drop[ss]) {
	for(k2=0; k2<n_gen2; k2++) 
	  Wts12[k1][k2][i] = param[s]*weights[i]+temp;
      }
      else s--;
    }

    /* QTL 2 effect */
    for(k2=0; k2<n_gen2-1; k2++, ss++, s++) { 
      if(!n_col2drop || !allcol2drop[ss]) {
	for(k1=0; k1<n_gen1; k1++) 
	  Wts12[k1][k2][i] += param[s]*weights[i];
      }
      else s--;
    }
    s += n_addcov;
    ss += n_addcov;

    /* QTL x interactive covar */
    for(j=0; j<n_intcov; j++) {
      for(k1=0; k1<n_gen1-1; k1++, ss++, s++) { /* QTL1 x intxn */
	if(!n_col2drop || !allcol2drop[ss]) {
	  for(k2=0; k2<n_gen2; k2++)
	    Wts12[k1][k2][i] += param[s]*Intcov[j][i];
	}
	else s--;
      }
      for(k2=0; k2<n_gen2-1; k2++, ss++, s++) { /* QTL2 x intxn */
	if(!n_col2drop || !allcol2drop[ss]) {
	  for(k1=0; k1<n_gen1; k1++)
	    Wts12[k1][k2][i] += param[s]*Intcov[j][i];
	}
	else s--;
      }
    }
      
    if(full_model) {
      /* QTL x QTL interaction */
      for(k1=0; k1<n_gen1-1; k1++) 
	for(k2=0; k2<n_gen2-1; k2++, ss++, s++) {
	  if(!n_col2drop || !allcol2drop[ss]) 
	    Wts12[k1][k2][i] += param[s]*weights[i];
	  else s--;
	}
      
      /* QTL x QTL x interactive covar */
      for(j=0; j<n_intcov; j++) {
	for(k1=0; k1<n_gen1-1; k1++) {
	  for(k2=0; k2<n_gen2-1; k2++, ss++, s++) {
	    if(!n_col2drop || !allcol2drop[ss]) 
	      Wts12[k1][k2][i] += param[s]*Intcov[j][i];
	    else s--;
	  }
	}
      }
    } 
    /* done calculating fitted values */
    /* s should now be at the location of the residual SD */

    /* calculate p(y|fitted,SD) for normal model 
       and multiple by Genoprob */
    temp=0.0;
    for(k1=0; k1<n_gen1; k1++) 
      for(k2=0; k2<n_gen2; k2++) 
	temp += 
	  (Wts12[k1][k2][i] = (dnorm(pheno[i],Wts12[k1][k2][i],param[s],0)*
			Probs[k1][k2][i]));
    
    /* rescale wts */
    if(rescale) {
      for(k1=0; k1<n_gen1; k1++) 
	for(k2=0; k2<n_gen2; k2++) 
	  Wts12[k1][k2][i] /= temp;

      /* marginal wts */
      for(k1=0; k1<n_gen1; k1++) {
	Wts1[k1][i] = 0.0;
	for(k2=0; k2<n_gen2; k2++)
	  Wts1[k1][i] += Wts12[k1][k2][i];
      }
      for(k2=0; k2<n_gen2; k2++) {
	Wts2[k2][i] = 0.0;
	for(k1=0; k1<n_gen1; k1++)
	  Wts2[k2][i] += Wts12[k1][k2][i];
      }
    } /* end rescale */
      
  } /* end loop over individuals */

}

double scantwo_em_loglik(int n_ind, int n_gen1, int n_gen2, double ***Probs,
			 double ***Wts12, double **Wts1, double **Wts2, 
			 double **Addcov, int n_addcov, double **Intcov,
			 int n_intcov, double *pheno, double *weights,
			 double *param, int full_model,
			 int n_col2drop, int *allcol2drop)
{
  double loglik, temp;
  int j, k1, k2;
  
  scantwo_em_estep(n_ind, n_gen1, n_gen2, Probs, Wts12, 
		   Wts1, Wts2, Addcov, n_addcov, Intcov, 
		   n_intcov, pheno, weights, param, full_model, 0,
		   n_col2drop, allcol2drop);


  loglik=0.0;
  for(j=0; j<n_ind; j++) {
    temp=0.0;
    for(k1=0; k1<n_gen1; k1++)
      for(k2=0; k2<n_gen2; k2++)
	temp += Wts12[k1][k2][j];
    loglik += log10(temp);
  }

  return(loglik);
}

/* end of scantwo_em.c */

