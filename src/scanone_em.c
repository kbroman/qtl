/**********************************************************************
 * 
 * scanone_em.c
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
 * single QTL model by interval mapping (the EM algorithm).
 *
 * Contains: R_scanone_em, scanone_em
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
#include "scanone_em.h"
#include "scanone_em_covar.h"

/**********************************************************************
 * 
 * R_scanone_em
 *
 * Wrapper for call from R; reorganizes genotype prob 
 * and calls scanone_em.
 * 
 **********************************************************************/

void R_scanone_em(int *n_ind, int *n_pos, int *n_gen, 
		  double *genoprob, double *addcov, int *n_addcov,
		  double *intcov, int *n_intcov, double *pheno,
		  double *weights,
		  double *result, int *std_start, double *start,
		  int *maxit, double *tol, int *verbose)
{
  double ***Genoprob, **work, **Addcov, **Intcov;
  double *means;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  allocate_dmatrix(4,*n_gen, &work);
  allocate_double(*n_gen, &means);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  if(*n_addcov == 0 && *n_intcov == 0) { /* no covariates */
    /* Read R's random seed */
    GetRNGstate();

    scanone_em(*n_ind, *n_pos, *n_gen, Genoprob, pheno, weights, 
	       result, *std_start, start, *maxit, *tol, work, means);

    /* Write R's random seed */
    PutRNGstate();
  }
  else { /* interval mapping with covariates */
    scanone_em_covar(*n_ind, *n_pos, *n_gen, Genoprob, Addcov,
		     *n_addcov, Intcov, *n_intcov, pheno, weights,
		     result, *maxit, *tol, *verbose);
  }
}

/**********************************************************************
 * 
 * scanone_em
 *
 * Performs genome scan using interval mapping.  (The multipoint
 * genotype probabilities have already been calculated in 
 * calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Upon exit, the LOD scores
 *
 * std_start    If 1, use the usual starting points [initial weights as 
 *                    Pr(QTL geno | marker genotypes)]
 *              If -1, use iid Bernoulli(1/2) for initial weights.
 *              If 0, use the specified values for the means and SD as
 *                    the starting point.
 *
 * start        If std_start = 0, use these as initial estimates of the
 *              genotype-specific means and the residual SD.
 * 
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of dimension 4 x n_gen
 *
 * means        Space for the phenotype means at each genotype
 *
 **********************************************************************/

void scanone_em(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		double *pheno, double *weights, 
		double *result, int std_start, double *start,
		int maxit, double tol, double **work, double *means)
{
  int i, j, k, s, flag=0;
  double s1, s2, s3, oldsig, r, sigma=0.0;

  /* turn weights back into usual scale rather than sqrt(weights) */
  for(j=0; j<n_ind; j++) 
    weights[j] *= weights[j];

  for(i=0; i<n_pos; i++) { /* loop over marker positions */

    /* initiate EM */
    s1 = 0.0;

    if(std_start==0) { /* specified starting point */
      for(k=0; k<n_gen; k++) work[1][k] = start[k];
      oldsig = start[n_gen];
    }
    else {
      if(std_start == 1) { /* the usual starting points */
	for(k=0; k<n_gen; k++) {
	  work[1][k] = s2 = s3 = 0.0;
	  for(j=0; j<n_ind; j++) {
	    s2 += Genoprob[k][i][j]*weights[j]; /* count up numbers */
	    work[1][k] += Genoprob[k][i][j]*pheno[j]*weights[j]; /* means */
	    s3 += Genoprob[k][i][j]*pheno[j]*pheno[j]*weights[j]; /* for RSS */
	  }
	  s1 += (s3 - work[1][k]*work[1][k]/s2); /* RSS */
	  work[1][k] /= s2;
	}
	oldsig = sqrt(s1/(double)n_ind);
      }
      else { /* start using random weights */
	for(k=0; k<n_gen; k++) {
	  work[1][k] = s2 = s3 = 0.0;
	  for(j=0; j<n_ind; j++) {
	    r = unif_rand()/(double)(n_gen); 
	    s2 += r*weights[j]; /* count up numbers */
	    work[1][k] += r*pheno[j]*weights[j]; /* means */
	    s3 += r*pheno[j]*pheno[j]*weights[j]; /* RSS */
	  }
	  s1 += (s3 - work[1][k]*work[1][k]/s2);
	  work[1][k] /= s2;
	}
	oldsig = sqrt(s1/(double)n_ind);
      }
    }

    for(s=0; s < maxit; s++) { /* EM iterations */
    
      R_CheckUserInterrupt(); /* check for ^C */

      for(k=0; k<n_gen; k++) 
	means[k] = work[2][k] = work[3][k] = 0.0;
      sigma=0.0;

      for(j=0; j<n_ind; j++) { /* loop over individuals */
	/* E-step */
	s1=0.0;
	for(k=0; k<n_gen; k++) 
	  s1 += (work[0][k] = Genoprob[k][i][j]*
		 dnorm(pheno[j],work[1][k],oldsig/sqrt(weights[j]),0));
	for(k=0; k<n_gen; k++) 
	  work[0][k] /= s1;
	
	/* M-step */
	for(k=0; k<n_gen; k++) {
	  work[2][k] += work[0][k]*weights[j]; /* count up numbers */
	  means[k] += work[0][k] * pheno[j]*weights[j]; /* means */
	  work[3][k] += work[0][k] * pheno[j] * pheno[j]*weights[j]; /* RSS */
	}
      }
      
      /* complete M-step */
      for(k=0; k<n_gen; k++) {
	sigma += (work[3][k] - means[k]*means[k]/work[2][k]);
	means[k] /= work[2][k];
      }
      sigma = sqrt(sigma/(double)n_ind);

      /* check for convergence */
      flag = 0;
      for(k=0; k<n_gen; k++) {
	if(fabs(means[k] - work[1][k]) > tol*(fabs(work[1][k])+tol*100.0)) {
	  flag = 1;
	  break;
	}
      }
      if(fabs(sigma - oldsig) > tol*(oldsig+tol*100.0)) flag = 1;

      if(!flag) break;

      oldsig = sigma;
      for(k=0; k<n_gen; k++) work[1][k] = means[k];

    } /* end of EM iterations */

    if(flag) warning("Didn't converge!\n");

    /* calculate negative log lik */
    result[i] = 0.0;
    for(j=0; j<n_ind; j++) {
      s1 = 0.0;
      for(k=0; k<n_gen; k++) 
	s1 += Genoprob[k][i][j] * dnorm(pheno[j], means[k], 
					sigma/sqrt(weights[j]), 0);
      result[i] -= log10(s1);
    }

  } /* end loop over marker positions */
}

/* end of scanone_em.c */

