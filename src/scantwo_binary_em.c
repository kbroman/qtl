/**********************************************************************
 * 
 * scantwo_binary_em.c
 *
 * copyright (c) 2004-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Dec, 2004
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
 * with a 2-QTL model by interval mapping (the EM algorithm) for  
 * a binary trait.
 *
 * Contains: R_scantwo_1chr_binary_em, scantwo_1chr_binary_em, 
 *           R_scantwo_2chr_binary_em, scantwo_2chr_binary_em,
 *           scantwo_binary_em_estep, scantwo_binary_em_mstep
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
#include "scantwo_binary_em.h"
#define TOL 1e-12


/**********************************************************************
 * 
 * R_scantwo_1chr_binary_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_binary_em.
 * 
 **********************************************************************/

void R_scantwo_1chr_binary_em(int *n_ind, int *n_pos, int *n_gen,
			      double *pairprob, double *addcov, int *n_addcov, 
			      double *intcov, int *n_intcov, int *pheno, 
			      double *start, double *result, int *maxit, 
			      double *tol, int *verbose, int *n_col2drop,
			      int *col2drop)
{
  double **Result, **Addcov, **Intcov, *****Pairprob;

  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
  reorg_errlod(*n_pos, *n_pos, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_1chr_binary_em(*n_ind, *n_pos, *n_gen, Pairprob, 
			 Addcov, *n_addcov, Intcov, *n_intcov, 
			 pheno, start, Result, *maxit, *tol, *verbose,
			 *n_col2drop, col2drop);
}

/**********************************************************************
 * 
 * scantwo_1chr_binary_em
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
 * Result       Result matrix of size [n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              an additive model.
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
 **********************************************************************/

void scantwo_1chr_binary_em(int n_ind, int n_pos, int n_gen, 
			    double *****Pairprob, double **Addcov, int n_addcov, 
			    double **Intcov, int n_intcov, int *pheno, double *start,
			    double **Result, int maxit, double tol, int verbose,
			    int n_col2drop, int *col2drop)
{
  int error_flag, i, i1, i2, k1, k2, j, m, n_col[2], n_col_rev[2], nit[2], r, flag=0;
  double *param, *oldparam, ***Wts12, pardif;
  double *wts, ***Probs, oldllik=0.0, llik[2];
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
  reorg_genoprob(n_ind, n_gen, n_gen, wts+2*n_gen*n_ind, &Wts12);
  reorg_genoprob(n_ind, n_gen, n_gen, wts+n_gen*(n_gen+2)*n_ind, &Probs);
  param = (double *)R_alloc(n_col[1], sizeof(double));
  oldparam = (double *)R_alloc(n_col[1], sizeof(double));
  
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

	for(j=0; j<n_col[m]; j++) oldparam[j] = start[j];

	scantwo_binary_em_mstep(n_ind, n_gen, n_gen, Addcov, n_addcov, 
				Intcov, n_intcov, pheno, Probs, 
				oldparam, m, n_col[m], &error_flag,
				n_col2drop, allcol2drop, verbose);
	if(error_flag) {
	  if(verbose>1) 
	    Rprintf("   [%3d %3d] %1d: Initial model had error.\n",
		    i1+1, i2+1, m+1);
	}
	else { /* only proceed if there's no error */
	  oldllik = scantwo_binary_em_loglik(n_ind, n_gen, n_gen, Probs, Addcov, 
					     n_addcov, Intcov, n_intcov, pheno, 
					     oldparam, m, n_col2drop, allcol2drop);
	  if(verbose>2) 
	    Rprintf("   [%3d %3d] %1d %9.3lf\n", 
		    i1+1, i2+1, m+1, oldllik);
	
	  for(j=0; j<n_col[m]; j++) param[j] = oldparam[j]; 

	  for(r=0; r<maxit; r++) { /* loop over iterations */

	    R_CheckUserInterrupt(); /* check for ^C */

	    scantwo_binary_em_estep(n_ind, n_gen, n_gen, Probs, Wts12, 
				    Addcov, n_addcov, Intcov, 
				    n_intcov, pheno, oldparam, m, 1,
				    n_col2drop, allcol2drop);

	    scantwo_binary_em_mstep(n_ind, n_gen, n_gen, Addcov, n_addcov, 
				    Intcov, n_intcov, pheno, Wts12, 
				    param, m, n_col[m], &error_flag,
				    n_col2drop, allcol2drop, verbose);
	    if(error_flag) {
	      flag=0;
	      if(verbose>1)
		Rprintf("   [%3d %3d] %1d %4d: Error in mstep\n",
			i1+1, i2+1, m+1, r+1);
	      break;
	    }

	    llik[m] = scantwo_binary_em_loglik(n_ind, n_gen, n_gen, Probs, Addcov, 
					       n_addcov, Intcov, n_intcov, pheno, 
					       param, m, n_col2drop, allcol2drop);

	    if(verbose>1) {
	      if(verbose>2) {
		pardif = fabs(param[0]-oldparam[0]);
		for(j=1; j<n_col[m]; j++) 
		  if(pardif <= fabs(param[j]-oldparam[j]))
		    pardif = fabs(param[j]-oldparam[j]);
		
		Rprintf("   [%3d %3d] %1d %4d %9.6lf    %lf\n", 
			i1+1, i2+1, m+1, r+1, (llik[m]-oldllik), pardif);
	      }
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
	    if((llik[m]-oldllik) < tol) { 
	      flag = 0; 
	      break; 
	    }

	    oldllik = llik[m];
	    for(j=0; j<n_col[m]; j++) oldparam[j] = param[j];

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
	Rprintf("   [%3d %3d]   %4d %4d    %9.6lf %9.6lf    %9.6lf", 
		i1+1, i2+1, nit[0], nit[1], llik[0], llik[1], llik[1]-llik[0]);
	if(llik[1] < llik[0]) Rprintf(" ****");
	Rprintf("\n");
      }

      Result[i2][i1] = -llik[0];
      Result[i1][i2] = -llik[1];

    } /* position 2 */
  } /* position 1 */
}

/**********************************************************************
 * 
 * R_scantwo_2chr_binary_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_binary_em.
 * 
 **********************************************************************/

void R_scantwo_2chr_binary_em(int *n_ind, int *n_pos1, int *n_pos2, 
			      int *n_gen1, int *n_gen2, double *genoprob1,
			      double *genoprob2, double *addcov, int *n_addcov, 
			      double *intcov, int *n_intcov, 
			      int *pheno, double *start,
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

  scantwo_2chr_binary_em(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2,
			 Genoprob1, Genoprob2, Addcov, *n_addcov, 
			 Intcov, *n_intcov, pheno, start,
			 Result_full, Result_add, 
			 *maxit, *tol, *verbose);
}

/**********************************************************************
 * 
 * scantwo_2chr_binary_em
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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos2][pos1]
 *
 * Result_add   Result matrix of size [n_pos2 x n_pos1] 
 *              containing the LODs for additive models
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

void scantwo_2chr_binary_em(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
			    int n_gen2, double ***Genoprob1, double ***Genoprob2,
			    double **Addcov, int n_addcov, double **Intcov, 
			    int n_intcov, int *pheno, double *start,
			    double **Result_full, double **Result_add, 
			    int maxit, double tol, int verbose)
{
  int error_flag, i1, i2, k1, k2, j, m, n_col[2], nit[2], r, flag=0;
  double *param, *oldparam, ***Wts12;
  double *wts, ***Probs, oldllik=0.0, llik[2];
  int n_col2drop=0, *allcol2drop=0;

  n_col[0] = (n_gen1+n_gen2-1) + n_addcov + (n_gen1+n_gen2-2)*n_intcov;
  n_col[1] = n_gen1*n_gen2 + n_addcov + (n_gen1*n_gen2-1)*n_intcov;

  /* allocate workspaces */
  wts = (double *)R_alloc((2*n_gen1*n_gen2+n_gen1+n_gen2)*n_ind, sizeof(double));
  reorg_genoprob(n_ind, n_gen2, n_gen1, wts+(n_gen1+n_gen2)*n_ind, &Wts12);
  reorg_genoprob(n_ind, n_gen2, n_gen1, 
		 wts+(n_gen1*n_gen2+n_gen1+n_gen2)*n_ind, &Probs);
  param = (double *)R_alloc(n_col[1], sizeof(double));
  oldparam = (double *)R_alloc(n_col[1], sizeof(double));
  
  /* begin loop over pairs of positions */
  for(i1=0; i1<n_pos1; i1++) {
    for(i2=0; i2<n_pos2; i2++) { /* loop over positions */
      nit[0] = nit[1] = 0;
      llik[0] = llik[1] = NA_REAL;

      /* calculate joint genotype probabilities */
      for(j=0; j<n_ind; j++) 
	for(k1=0; k1<n_gen1; k1++)
	  for(k2=0; k2<n_gen2; k2++)
	    Probs[k1][k2][j] = Genoprob1[k1][i1][j]*Genoprob2[k2][i2][j];

      for(m=0; m<2; m++) { /* loop over add've model and full model */

	for(j=0; j<n_col[m]; j++) oldparam[j] = start[j];

	scantwo_binary_em_mstep(n_ind, n_gen1, n_gen2, Addcov, n_addcov, 
				Intcov, n_intcov, pheno, Probs, 
				oldparam, m, n_col[m], &error_flag,
				n_col2drop, allcol2drop, verbose);

	if(error_flag) {
	  if(verbose>1)
	    Rprintf("   [%3d %3d] %1d: Initial model had error.\n",
		    i1+1, i2+1, m+1);
	}
	else { /* only proceed if there's no error */
	  oldllik = scantwo_binary_em_loglik(n_ind, n_gen1, n_gen2, Probs, 
					     Addcov, n_addcov, Intcov, 
					     n_intcov, pheno, oldparam, m,
					     n_col2drop, allcol2drop);

	  if(verbose>2)
	    Rprintf("   [%3d %3d] %1d %9.3lf\n", 
		    i1+1, i2+1, m+1, oldllik);
	
	  for(j=0; j<n_col[m]; j++) param[j] = oldparam[j]; 

	  for(r=0; r<maxit; r++) { /* loop over iterations */

	    R_CheckUserInterrupt(); /* check for ^C */

	    scantwo_binary_em_estep(n_ind, n_gen1, n_gen2, Probs, Wts12, 
				    Addcov, n_addcov, Intcov, 
				    n_intcov, pheno, oldparam, m, 1,
				    n_col2drop, allcol2drop);

	    scantwo_binary_em_mstep(n_ind, n_gen1, n_gen2, Addcov, n_addcov, 
				    Intcov, n_intcov, pheno, Wts12, 
				    param, m, n_col[m], &error_flag,
				    n_col2drop, allcol2drop, verbose);

	    if(error_flag) {
	      flag=0;
	      if(verbose>1)
		Rprintf("   [%3d %3d] %1d %4d: Error in mstep\n",
			i1+1, i2+1, m+1, r+1);
	      break;
	    }

	    llik[m] = scantwo_binary_em_loglik(n_ind, n_gen1, n_gen2, Probs, 
					       Addcov, n_addcov, Intcov, 
					       n_intcov, pheno, param, m,
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
	    if((llik[m]-oldllik) < tol) { 
	      flag = 0; 
	      break; 
	    }

	    oldllik = llik[m];
	    for(j=0; j<n_col[m]; j++) oldparam[j] = param[j];

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
	Rprintf("   [%3d %3d]   %4d %4d    %9.6lf %9.6lf    %9.6lf", 
		i1+1, i2+1, nit[0], nit[1], llik[0], llik[1], llik[1]-llik[0]);
	if(llik[1] < llik[0]) Rprintf(" ****");
	Rprintf("\n");
      }

      Result_add[i2][i1] = -llik[0];
      Result_full[i2][i1] = -llik[1];  
    } /* position 2 */
  } /* position 1 */
}

/**********************************************************************
 * 
 * scantwo_binary_em_mstep: M-step of the EM algorithm 
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
 * Wts12    Pr(QTL1=v, QTL2=w | phenotype, model, marker data),
 *          indexed as Wts[v][w][ind]
 *
 * param    On output, the updated parameter estimates (incl resid SD)
 *
 * full_model   If 1, include QTLxQTL interaction
 * 
 * n_col
 *
 * error_flag     Set to 1 if X'X is singular
 *

 **********************************************************************/

void scantwo_binary_em_mstep(int n_ind, int n_gen1, int n_gen2, 
			     double **Addcov, int n_addcov, 
			     double **Intcov, int n_intcov, int *pheno, 
			     double ***Wts12, 
			     double *param, int full_model,
			     int n_col, int *error_flag,
			     int n_col2drop, int *allcol2drop,
			     int verbose)
{
  int i, j, j2, k1, k2, s, s2, info, ss;
  double rcond, temp, *junk;
  double *grad, *jac, **Jac;
  double *tf1, ***f1, *tf2, ***f2;
  double *tfitted, **fitted;

  /* allocate space */
  allocate_double(n_col, &junk);
  allocate_double(n_col, &grad);
  allocate_double(n_col*n_col, &jac);
  reorg_errlod(n_col, n_col, jac, &Jac);

  allocate_double(n_gen1*n_gen2*n_ind, &tf1);
  allocate_double(n_gen1*n_gen2*n_ind, &tf2);
  reorg_genoprob(n_ind, n_gen2, n_gen1, tf1, &f1);
  reorg_genoprob(n_ind, n_gen2, n_gen1, tf2, &f2);

  allocate_double(n_gen1*n_gen2, &tfitted);
  reorg_errlod(n_gen2, n_gen1, tfitted, &fitted);

  *error_flag=0;

  for(j=0; j<n_col; j++) {
    grad[j] = 0.0;
    for(k1=0; k1<n_col; k1++) Jac[j][k1] = 0.0;
  }

  /* calculate fitted values */
  for(i=0; i<n_ind; i++) {
    if(n_col2drop) {
      for(ss=0, s=0; ss<n_gen1+n_gen2-1; ss++) 
	if(!allcol2drop[ss]) s++;
    }
    else s = n_gen1+n_gen2-1;

    temp = 0.0;
    for(j=0; j<n_addcov; j++, s++) 
      temp += (Addcov[j][i]*param[s]);
  
    /* QTL 1 effect */
    for(k1=0, ss=0, s=0; k1<n_gen1; k1++, ss++, s++) { 
      if(!n_col2drop || !allcol2drop[ss]) {
	for(k2=0; k2<n_gen2; k2++) 
	  fitted[k1][k2] = param[s]+temp;
      }
      else s--;
    }

    /* QTL 2 effect */
    for(k2=0; k2<n_gen2-1; k2++, ss++, s++) { 
      if(!n_col2drop || !allcol2drop[ss]) {
	for(k1=0; k1<n_gen1; k1++) 
	  fitted[k1][k2] += param[s];
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
	    fitted[k1][k2] += param[s]*Intcov[j][i];
	}
	else s--;
      }

      for(k2=0; k2<n_gen2-1; k2++, ss++, s++) { /* QTL2 x intxn */
	if(!n_col2drop || !allcol2drop[ss]) {
	  for(k1=0; k1<n_gen1; k1++)
	    fitted[k1][k2] += param[s]*Intcov[j][i];
	}
	else s--;
      }
    }
      
    if(full_model) {
      /* QTL x QTL interaction */
      for(k1=0; k1<n_gen1-1; k1++) 
	for(k2=0; k2<n_gen2-1; k2++, ss++, s++) {
	  if(!n_col2drop || !allcol2drop[ss]) 
	    fitted[k1][k2] += param[s];
	  else s--;
	}

      /* QTL x QTL x interactive covar */
      for(j=0; j<n_intcov; j++) {
	for(k1=0; k1<n_gen1-1; k1++) {
	  for(k2=0; k2<n_gen2-1; k2++, s++, ss++) {
	    if(!n_col2drop || !allcol2drop[ss]) 
	      fitted[k1][k2] += param[s]*Intcov[j][i];
	    else s--;
	  }
	}
      }
    }


    for(k1=0; k1<n_gen1; k1++) {
      for(k2=0; k2<n_gen2; k2++) {
	fitted[k1][k2] = exp(fitted[k1][k2]);
	fitted[k1][k2] /= (1.0+fitted[k1][k2]);
	f1[k1][k2][i] = Wts12[k1][k2][i]*((double)pheno[i]-fitted[k1][k2]);
	f2[k1][k2][i] = Wts12[k1][k2][i]*fitted[k1][k2]*(1.0-fitted[k1][k2]);

      }
    }
  }
    
  /* calculate gradient */
  for(i=0; i<n_ind; i++) {
    for(k1=0; k1<n_gen1; k1++)  /* QTL 1 */
      for(k2=0; k2<n_gen2; k2++) 
	grad[k1] += f1[k1][k2][i];
    s = n_gen1;
    for(k2=0; k2<n_gen2-1; k2++)
      for(k1=0; k1<n_gen1; k1++)
	grad[k2+s] += f1[k1][k2][i];
    s += (n_gen2-1);
    for(j=0; j<n_addcov; j++) /* add covar */
      for(k1=0; k1<n_gen1; k1++)
	for(k2=0; k2<n_gen2; k2++)
	  grad[j+s] += Addcov[j][i]*f1[k1][k2][i];
    s += n_addcov;
    for(j=0; j<n_intcov; j++) {
      for(k1=0; k1<n_gen1-1; k1++)
	for(k2=0; k2<n_gen2; k2++)
	  grad[s+k1] += Intcov[j][i]*f1[k1][k2][i];
      s += (n_gen1-1);
      for(k2=0; k2<n_gen2-1; k2++)
	for(k1=0; k1<n_gen1; k1++)
	  grad[s+k2] += Intcov[j][i]*f1[k1][k2][i];
      s += (n_gen2-1);
    }

    if(full_model) {
      for(k1=0; k1<n_gen1-1; k1++) 
	for(k2=0; k2<n_gen2-1; k2++) 
	  grad[s+k1*(n_gen2-1)+k2] += f1[k1][k2][i];
      s += (n_gen1-1)*(n_gen2-1);
      for(j=0; j<n_intcov; j++) {
	for(k1=0; k1<n_gen1-1; k1++) 
	  for(k2=0; k2<n_gen2-1; k2++) 
	    grad[s+k1*(n_gen2-1)+k2] += 
	      Intcov[j][i]*f1[k1][k2][i];
	s += (n_gen1-1)*(n_gen2-1);
      }
    }
  } /* end loop over individuals */

  /* calculate Jacobian; only the upper right triangle is needed */
  for(i=0; i<n_ind; i++) {
    /* QTL 1 columns */
    for(k1=0; k1<n_gen1; k1++) 
      for(k2=0; k2<n_gen2; k2++)
	Jac[k1][k1] += f2[k1][k2][i];

    /* QTL 2 columns */
    for(k2=0, s=n_gen1; k2<n_gen2-1; k2++, s++) {
      for(k1=0; k1<n_gen1; k1++)
	Jac[s][s] += f2[k1][k2][i];
      for(k1=0; k1<n_gen1; k1++)  
	Jac[s][k1] += f2[k1][k2][i];
    }
    
    /* add covar columns */
    for(j=0; j<n_addcov; j++, s++) {
      for(k1=0; k1<n_gen1; k1++) 
	for(k2=0; k2<n_gen2; k2++) 
	  Jac[s][s] += (Addcov[j][i]*Addcov[j][i]*f2[k1][k2][i]);
      for(k1=0; k1<n_gen1; k1++) 
	for(k2=0; k2<n_gen2; k2++)
	  Jac[s][k1] += (Addcov[j][i]*f2[k1][k2][i]);
      for(k2=0, s2=n_gen1; k2<n_gen2-1; k2++, s2++) 
	for(k1=0; k1<n_gen1; k1++) 
	  Jac[s][s2] += Addcov[j][i]*f2[k1][k2][i];
      for(j2=0; j2<j; j2++, s2++) 
	for(k1=0; k1<n_gen1; k1++)
	  for(k2=0; k2<n_gen2; k2++)
	    Jac[s][s2] += (Addcov[j2][i]*Addcov[j][i]*f2[k1][k2][i]);
    }

    /* QTL x interactive covariates columns */
    for(j=0; j<n_intcov; j++) {
      /* int covar x QTL 1 */
      for(k1=0; k1<n_gen1-1; k1++, s++) {
	for(k2=0; k2<n_gen2; k2++) {
	  Jac[s][s] += (Intcov[j][i]*Intcov[j][i]*f2[k1][k2][i]);
	  Jac[s][k1] += (Intcov[j][i]*f2[k1][k2][i]); /* x QTL 1 */
	}
	for(k2=0, s2=n_gen1; k2<n_gen2-1; k2++, s2++) /* x QTL 2 */
	  Jac[s][s2] += (Intcov[j][i]*f2[k1][k2][i]);
	for(j2=0; j2<n_addcov; j2++, s2++) /* x add covar */
	  for(k2=0; k2<n_gen2; k2++)
	  Jac[s][s2] += (Intcov[j][i]*Addcov[j2][i]*f2[k1][k2][i]);
	/* prev interactive covar */
	for(j2=0; j2<j; j2++, s2 += (n_gen1+n_gen2-2)) {
	  /* x (QTL 1 x prev inter've covar) */
	  for(k2=0; k2<n_gen2; k2++)
	    Jac[s][s2+k1] += (Intcov[j][i]*Intcov[j2][i]*f2[k1][k2][i]);
	  /* x (QTL 2 x prev inter've covar) */
	  for(k2=0; k2<n_gen2-1; k2++) 
	    Jac[s][s2+n_gen1-1+k2] += 
	      (Intcov[j][i]*Intcov[j2][i]*f2[k1][k2][i]);
	}
      }

      /* int covar x QTL 2 */
      for(k2=0; k2<n_gen2-1; k2++, s++) { 
	for(k1=0; k1<n_gen1; k1++)
	  Jac[s][s] += (Intcov[j][i]*Intcov[j][i]*f2[k1][k2][i]);
	/* x QTL 1 */
	for(k1=0; k1<n_gen1; k1++) 
	  Jac[s][k1] += (Intcov[j][i]*f2[k1][k2][i]);
	/* x QTL 2 */
	for(k1=0; k1<n_gen1; k1++)
	  Jac[s][n_gen1+k2] += (Intcov[j][i]*f2[k1][k2][i]);
	/* x add covar */
	for(j2=0, s2=n_gen1+n_gen2-1; j2<n_addcov; j2++, s2++)
	  for(k1=0; k1<n_gen1; k1++)
	    Jac[s][s2] += (Intcov[j][i]*Addcov[j2][i]*f2[k1][k2][i]);
	/* x prev interactive covar */
	for(j2=0; j2<j; j2++, s2 += (n_gen1+n_gen2-2)) {
	  /* x (QTL 1 x prev inter've covar) */
	  for(k1=0; k1<n_gen1-1; k1++) 
	    Jac[s][s2+k1] += 
	      (Intcov[j][i]*Intcov[j2][i]*f2[k1][k2][i]);
	  /* x (QTL 2 x prev inter've covar) */
	  for(k1=0; k1<n_gen1; k1++)
	    Jac[s][s2+n_gen1-1+k2] += 
	      (Intcov[j][i]*Intcov[j2][i]*f2[k1][k2][i]);
	}
	/* x (QTL 1 x inter've covar) */
	for(k1=0; k1<n_gen1-1; k1++, s2++) 
	  Jac[s][s2] += (Intcov[j][i]*Intcov[j][i]*f2[k1][k2][i]);
      }
    } /* end of interactive covariates */

    if(full_model) {
      /* QTL x QTL interactions */
      for(k1=0; k1<n_gen1-1; k1++) {
	for(k2=0; k2<n_gen2-1; k2++, s++) {
	  Jac[s][s] += f2[k1][k2][i];
	  /* x QTL 1 */
	  Jac[s][k1] += f2[k1][k2][i];
	  /* x QTL 2 */
	  Jac[s][n_gen1+k2] += f2[k1][k2][i];
	  /* x add covar */
	  for(j=0, s2 = n_gen1+n_gen2-1; j<n_addcov; j++, s2++) 
	    Jac[s][s2] += (f2[k1][k2][i]*Addcov[j][i]);
	  /* x interactive covariates */
	  for(j=0; j<n_intcov; j++, s2+=n_gen1+n_gen2-2) {
	    Jac[s][s2+k1] += (f2[k1][k2][i]*Intcov[j][i]);
	    Jac[s][s2+n_gen1-1+k2] +=
	      (f2[k1][k2][i]*Intcov[j][i]);
	  }
	}
      } /* end of QTL x QTL interactions */

      /* QTL x QTL x inter've covar */
      for(j=0; j<n_intcov; j++) {
	for(k1=0; k1<n_gen1-1; k1++) {
	  for(k2=0; k2<n_gen2-1; k2++, s++) {
	    temp = f2[k1][k2][i]*Intcov[j][i];
	    Jac[s][s] += temp*Intcov[j][i];
	    Jac[s][k1] += temp; /* x QTL 1 */
	    Jac[s][n_gen1+k2] += temp; /* x QTL 2 */
	    /* x add covar */
	    for(j2=0, s2=n_gen1+n_gen2-1; j2<n_addcov; j2++, s2++)
	      Jac[s][s2] += temp*Addcov[j2][i];
	    /* x int covar */
	    for(j2=0; j2<n_intcov; j2++, s2 += (n_gen1+n_gen2-2)) {
	      Jac[s][s2+k1] += temp*Intcov[j2][i];
	      Jac[s][s2+n_gen1-1+k2] += temp*Intcov[j2][i];
	    }
	    /* x (QTL x QTL) */
	    Jac[s][s2+k1*(n_gen2-1)+k2] += temp;
	    s2 += (n_gen1-1)*(n_gen2-1);
	    /* x (QTL x QTL x prev covar) */
	    for(j2=0; j2<j; j2++, s2 += ((n_gen1-1)*(n_gen2-1))) 
	      Jac[s][s2+k1*(n_gen2-1)+k2] += (temp*Intcov[j2][i]);
	  }
	}
      } /* end QTL x QTL x inter've covar */

    }
  } /* end loop over individuals */
  /* done calculating Jacobian */

  if(n_col2drop) { /* drop some columns (X chromosome) */
    dropcol_xpy(n_col, allcol2drop, grad);

    dropcol_xpx(&n_col, allcol2drop, jac);
  }

  /* solve work1 * beta = work2 for beta */
  F77_CALL(dpoco)(jac, &n_col, &n_col, &rcond, junk, &info);
  if(fabs(rcond) < TOL || info != 0) { /* error! */
    if(verbose > 1) Rprintf("X'X matrix is singular.\n");
    *error_flag = 1;
  }
  else {
    F77_CALL(dposl)(jac, &n_col, &n_col, grad);
    for(j=0; j<n_col; j++) param[j] += grad[j];
  }
}

/**********************************************************************
 * 
 * scantwo_binary_em_estep: E-step of the EM algorithm 
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
 * param    Current parameter estimates (including the resid SD)
 *
 * full_model   If 1, use the full model (with QTLxQTL interaction)
 *
 * rescale  If 1, rescale weights so that the sum to 1.
 *          This is done so that by taking rescale=0, we can easily
 *          calculate the log likelihood 
 *
 **********************************************************************/

void scantwo_binary_em_estep(int n_ind, int n_gen1, int n_gen2, 
			     double ***Probs, double ***Wts12, 
			     double **Addcov, int n_addcov, double **Intcov,
			     int n_intcov, int *pheno, 
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
	  Wts12[k1][k2][i] = param[s]+temp;
      }
      else s--;
    }

    /* QTL 2 effect */
    for(k2=0; k2<n_gen2-1; k2++, ss++, s++) { 
      if(!n_col2drop || !allcol2drop[ss]) {
	for(k1=0; k1<n_gen1; k1++) 
	  Wts12[k1][k2][i] += param[s];
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
	    Wts12[k1][k2][i] += param[s];
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

    temp = 0.0;
    for(k1=0; k1<n_gen1; k1++) { 
      for(k2=0; k2<n_gen2; k2++) {
	Wts12[k1][k2][i] = exp(Wts12[k1][k2][i]);

	if(pheno[i]) 
	  temp += (Wts12[k1][k2][i] = Probs[k1][k2][i]*Wts12[k1][k2][i]/
		   (1.0 + Wts12[k1][k2][i]));
	else
	  temp += (Wts12[k1][k2][i] = Probs[k1][k2][i]/(1.0+Wts12[k1][k2][i]));
      }
    }
    
    /* rescale wts */
    if(rescale) 
      for(k1=0; k1<n_gen1; k1++) 
	for(k2=0; k2<n_gen2; k2++) 
	  Wts12[k1][k2][i] /= temp;

  } /* end loop over individuals */

}

double scantwo_binary_em_loglik(int n_ind, int n_gen1, int n_gen2, 
				double ***Probs, double **Addcov, int n_addcov,
				double **Intcov, int n_intcov, int *pheno,
				double *param, int full_model,
				int n_col2drop, int *allcol2drop)
{
  double *wts, ***Wts;
  double loglik, temp;
  int i, k1, k2;

  /* allocate space */
  allocate_double(n_gen1*n_gen2*n_ind, &wts);
  reorg_genoprob(n_ind, n_gen2, n_gen1, wts, &Wts);

  scantwo_binary_em_estep(n_ind, n_gen1, n_gen2, Probs, Wts, 
			  Addcov, n_addcov, Intcov, n_intcov, 
			  pheno, param, full_model, 0,
			  n_col2drop, allcol2drop);
  
  loglik=0.0;
  for(i=0; i<n_ind; i++) {
    temp = 0.0;
    for(k1=0; k1<n_gen1; k1++)
      for(k2=0; k2<n_gen2; k2++)
	temp += Wts[k1][k2][i];
    loglik += log10(temp);
  }
  return(loglik);
}

/* end of scantwo_binary_em.c */

