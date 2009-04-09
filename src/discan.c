/**********************************************************************
 * 
 * discan.c
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Oct, 2001
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
 * single QTL model
 *
 * Contains: R_discan_mr, discan_mr,
 *           R_discan_im, discan_im
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "discan.h"

/**********************************************************************
 * 
 * R_discan_mr
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls discan_mr.
 * 
 **********************************************************************/

void R_discan_mr(int *n_ind, int *n_pos, int *n_gen,
		    int *geno, int *pheno, double *result)
{
  int **Geno;
  double *means;

  reorg_geno(*n_ind, *n_pos, geno, &Geno);
  allocate_double(*n_gen, &means);

  discan_mr(*n_ind, *n_pos, *n_gen, Geno, pheno, result, means);
}

/**********************************************************************
 * 
 * R_discan_im
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls discan_im.
 * 
 **********************************************************************/

void R_discan_im(int *n_ind, int *n_pos, int *n_gen, 
		 double *genoprob, int *pheno, double *result, 
		 int *maxit, double *tol)
{
  double ***Genoprob, **work;
  double *means;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  allocate_dmatrix(3, *n_gen, &work);
  allocate_double(*n_gen, &means);

  discan_im(*n_ind, *n_pos, *n_gen, Genoprob, pheno, result, 
	    *maxit, *tol, work, means);

}

/**********************************************************************
 * 
 * discan_mr
 *
 * Performs genome scan using marker regression for a dichotomous trait
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Geno         Genotype matrix
 *
 * pheno        Phenotype data, as a vector
 *
 * result       Upon return, to contain the log10 likelihoods
 *
 * means        Space for the phenotype means for each genotype
 *
 **********************************************************************/

void discan_mr(int n_ind, int n_pos, int n_gen, int **Geno, 
		  int *pheno, double *result, double *means)
{
  int i, j, k, n, tp, *ng, *np;

  /* number of individuals in each genotype group */
  allocate_int(n_gen, &ng);

  /* number of individuals with phenotype=1 in each genotype group */
  allocate_int(n_gen, &np);

  for(i=0; i<n_pos; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    result[i] = 0.0; n=tp=0; 
    for(j=0; j< n_gen; j++) {
      ng[j] = np[j] = 0;

      for(k=0; k<n_ind; k++) {
	if(Geno[i][k] == j+1) {
	  if(pheno[k]) {
	    np[j]++;
	    tp++;
	  }
	  n++; ng[j]++;
	}
      }
      
      if(ng[j] > 0) means[j] = (double)np[j]/(double)ng[j];
      else /* no individuals with this genotype */
	means[j] = NA_REAL;
    } /* loop over genotype groups */

    /* calculate LOD score */
    for(j=0; j<n_gen; j++) {
      if(np[j] > 0 && np[j] < ng[j]) 
	result[i] += ((double)np[j]*log10(means[j]) +
			 (double)(ng[j]-np[j])*log10(1.0-means[j]));
    }
    if(tp > 0 && tp < n) 
      result[i] -= ((double)tp*log10((double)tp/(double)n) +
		       (double)(n-tp)*log10((double)(n-tp)/(double)n));

  } /* loop over marker positions */
}

/**********************************************************************
 * 
 * discan_im
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
 * result       Upon return, to contain the log10 likelihoods
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of dimension 3 x n_gen
 *
 * means        Space for the phenotype means for each genotype
 *
 **********************************************************************/

void discan_im(int n_ind, int n_pos, int n_gen, double ***Genoprob,
	       int *pheno, double *result, 
	       int maxit, double tol, double **work, double *means)
{
  int i, j, k, s, flag=0;
  double sw;

  for(i=0; i<n_pos; i++) { /* loop over marker positions */

    /* initiate EM */
    for(k=0; k<n_gen; k++) {
      means[k] = sw=0.0;
      for(j=0; j<n_ind; j++) {
	sw += Genoprob[k][i][j]; /* k=gen, i=pos, j=ind */
	if(pheno[j]) means[k] += Genoprob[k][i][j];
      }
      means[k] /= sw;
    }

    for(s=0; s < maxit; s++) { /* EM iterations */
    
      R_CheckUserInterrupt(); /* check for ^C */

      /* copy over current estimates */
      for(k=0; k<n_gen; k++) {
	work[0][k] = means[k]; 
	means[k] = work[1][k] = 0.0;
      }
      
      for(j=0; j<n_ind; j++) { /* loop over individuals */
	/* E-step */
	sw = 0.0;
	for(k=0; k<n_gen; k++) {
	  work[2][k] = Genoprob[k][i][j];
	  if(pheno[j]) work[2][k] *= work[0][k];
	  else work[2][k] *= (1.0-work[0][k]);
	  sw += work[2][k];
	}
	for(k=0; k<n_gen; k++) work[2][k] /= sw;

	/* M-step */
	for(k=0; k<n_gen; k++) {
	  work[1][k] += work[2][k];
	  if(pheno[j]) means[k] += work[2][k];
	}
      }
      
      /* complete M-step */
      for(k=0; k<n_gen; k++) 
	means[k] /= work[1][k];

      /* check for convergence */
      flag = 0;
      for(k=0; k<n_gen; k++) {
	if(fabs(means[k] - work[0][k]) > tol*(fabs(work[0][k])+tol*100.0)) {
	  flag = 1;
	  break;
	}
      }

      if(!flag) break;
    } /* end of EM iterations */

    if(flag) warning("Didn't converge!\n");

    /* calculate log lik */
    result[i] = 0.0;
    for(j=0; j<n_ind; j++) {
      sw = 0.0;
      if(pheno[j]) 
	for(k=0; k<n_gen; k++) sw += Genoprob[k][i][j] * means[k];
      else
	for(k=0; k<n_gen; k++) sw += Genoprob[k][i][j] * (1.0-means[k]);
      result[i] += log10(sw);
    }

  } /* end loop over marker positions */
}



/* end of discan.c */

