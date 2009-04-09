/**********************************************************************
 * 
 * discan_covar.c
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
 * These functions are for performing a genome scan with a binary 
 * trait and a single QTL model in the presence of covariates.
 *
 * Contains: discan_covar, discan_covar_em, discan_covar_loglik,
 *           R_discan_covar
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include "util.h"
#include "discan_covar.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * discan_covar
 *
 * Performs genome scan using interval mapping in the presence of
 * covariates.  (The multipoint genotype probabilities have already 
 * been calculated in calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates indexed as 
 *              Addcov[cov][ind]
 *
 * n_addcov     Number of columns in Addcov
 *
 * Intcov       Matrix of interactive covariates indexed as 
 *              Intcov[cov][ind]
 *
 * n_intcov     Number of columns in Intcov
 *
 * pheno        Phenotype data (0/1), as a vector
 *
 * start        Starting values; vector of length 
 *              n_gen + n_addcov + n_intcov*(n_gen-1)
 *
 * result       Result vector of length n_pos; upon return, contains 
 *              the LOD scores.
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * verbose      If 1, print out log likelihood at each iteration
 *
 **********************************************************************/

void discan_covar(int n_ind, int n_pos, int n_gen, 
		  double ***Genoprob, double **Addcov, int n_addcov,
		  double **Intcov, int n_intcov, int *pheno, 
		  double *start, double *result, int maxit, double tol, 
		  int verbose)
{
  int i, n_par;

  n_par = n_gen + n_addcov + (n_gen-1)*n_intcov;

  for(i=0; i<n_pos; i++)  
    result[i] = discan_covar_em(n_ind, i, n_gen, n_par, Genoprob,
				Addcov, n_addcov, Intcov, n_intcov,
				pheno, start, maxit, tol, verbose);

}


void R_discan_covar(int *n_ind, int *n_pos, int *n_gen, 
		    double *genoprob, double *addcov, int *n_addcov,
		    double *intcov, int *n_intcov, int *pheno, 
		    double *start, double *result, int *maxit, double *tol, 
		    int *verbose)
{
  double ***Genoprob, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);
  
  discan_covar(*n_ind, *n_pos, *n_gen, Genoprob, 
	       Addcov, *n_addcov, Intcov, *n_intcov, pheno,
	       start, result, *maxit, *tol, *verbose);
}



/**********************************************************************
 * 
 * discan_covar_em
 *
 **********************************************************************/

double discan_covar_em(int n_ind, int pos, int n_gen, int n_par,
		       double ***Genoprob, double **Addcov, int n_addcov,
		       double **Intcov, int n_intcov, int *pheno, 
		       double *start, int maxit, double tol, int verbose)
{
  int i, j, k, kk, s, offset;
  double *jac, **Jac, *grad;
  double *temp, **wts;
  double fit, **temp1, **temp2;
  double *temp1s, *temp2s;
  double *newpar, *curpar;
  double newllik=0.0, curllik, sum;
  int info;
  double rcond, *junk;

  /* allocate space */
  allocate_double(n_par*n_par, &jac);
  reorg_errlod(n_par, n_par, jac, &Jac);
  allocate_double(n_par, &grad);
  allocate_double(n_ind*n_gen, &temp);
  reorg_errlod(n_gen, n_ind, temp, &wts);
  allocate_double(n_ind*n_gen, &temp);
  reorg_errlod(n_gen, n_ind, temp, &temp1);
  allocate_double(n_ind*n_gen, &temp);
  reorg_errlod(n_gen, n_ind, temp, &temp2);
  allocate_double(n_ind, &temp1s);
  allocate_double(n_ind, &temp2s);
  allocate_double(n_par, &newpar);
  allocate_double(n_par, &curpar);
  allocate_double(n_par, &junk);

  /* initial wts */
  for(i=0; i<n_ind; i++) 
    for(j=0; j<n_gen; j++) 
      wts[i][j] = Genoprob[j][pos][i];

  for(i=0; i<n_par; i++) curpar[i] = start[i];

  curllik = discan_covar_loglik(n_ind, pos, n_gen, n_par, curpar, Genoprob,
				Addcov, n_addcov, Intcov, n_intcov, pheno);
  if(verbose) 
    Rprintf("        %10.5f\n", curllik);

  for(s=0; s<maxit; s++) {

    R_CheckUserInterrupt(); /* check for ^C */

    /****** M STEP ******/

    /* 0's in gradient and Jacobian */
    for(j=0; j<n_par; j++) {
      grad[j] = 0.0;
      for(k=0; k<n_par; k++) 
	Jac[j][k] = 0.0;
    }

    /* calculate gradient and Jacobian */
    for(i=0; i<n_ind; i++) {
      temp1s[i] = temp2s[i] = 0.0;
      for(j=0; j<n_gen; j++) {
	fit = curpar[j];

	if(n_addcov > 0) {
	  for(k=0; k<n_addcov; k++) 
	    fit += Addcov[k][i] * curpar[n_gen+k];
	}
      
	if(n_intcov > 0 && j<n_gen-1) {
	  for(k=0; k<n_intcov; k++) 
	    fit += Intcov[k][i] * curpar[n_gen+n_addcov+n_intcov*j+k];
	}
	fit = exp(fit)/(1.0+exp(fit));

	temp1s[i] += (temp1[i][j] = wts[i][j]*((double)pheno[i] - fit));
	temp2s[i] += (temp2[i][j] = wts[i][j]*fit*(1.0-fit));
      }
    }

    for(j=0; j<n_gen; j++) {
      for(i=0; i<n_ind; i++) {
	grad[j] += temp1[i][j];
	Jac[j][j] += temp2[i][j];
      }
    }

    for(k=0; k<n_addcov; k++) {
      for(i=0; i<n_ind; i++) {
	grad[k + n_gen] += Addcov[k][i] * temp1s[i];

	for(kk=k; kk<n_addcov; kk++) 
	  Jac[kk+n_gen][k+n_gen] += Addcov[k][i]*Addcov[kk][i] * temp2s[i];
	
	for(j=0; j<n_gen; j++) 
	  Jac[k+n_gen][j] += Addcov[k][i] * temp2[i][j];
      
      }
    }
    
    for(j=0; j<n_gen-1; j++) {
      offset = n_gen + n_addcov + n_intcov*j;
      for(k=0; k<n_intcov; k++) {

	for(i=0; i<n_ind; i++) {
	  grad[k + offset] += Intcov[k][i] * temp1[i][j];

	  for(kk=k; kk<n_intcov; kk++) 
	    Jac[kk+offset][k+offset] += Intcov[k][i]*Intcov[kk][i]*temp2[i][j];

	  for(kk=0; kk<n_addcov; kk++)
	    Jac[k+offset][kk+n_gen] += Intcov[k][i]*Addcov[kk][i]*temp2[i][j];
	
	  Jac[k+offset][j] += Intcov[k][i]*temp2[i][j];
	}
      }
    }

    if(verbose > 1) {
        Rprintf("grad: ");
    for(j=0; j<n_par; j++) 
      Rprintf("%f ", grad[j]);
    Rprintf("\n");
    Rprintf("Jac:\n");
    for(j=0; j<n_par; j++) {
      for(k=0; k<n_par; k++) 
	Rprintf("%f ", Jac[j][k]);
      Rprintf("\n");
      } 
    Rprintf("\n");
    }
    
    /* dpoco and dposl from Linpack to calculate Jac^-1 %*% grad */
    F77_CALL(dpoco)(jac, &n_par, &n_par, &rcond, junk, &info);
    if(fabs(rcond) < TOL || info != 0) {
      warning("Jacobian matrix is singular\n");
      return(NA_REAL);
    }
    F77_CALL(dposl)(jac, &n_par, &n_par, grad);

    if(verbose > 1) {
      Rprintf(" solution: ");
      for(j=0; j<n_par; j++) 
	Rprintf(" %f", grad[j]);
      Rprintf("\n");
    }

    /* revised estimates */
    for(j=0; j<n_par; j++) 
      newpar[j] = curpar[j] + grad[j];

    if(verbose>1) {
      for(j=0; j<n_par; j++)
	Rprintf("%f ", newpar[j]);
      Rprintf("\n");
    }

    newllik = discan_covar_loglik(n_ind, pos, n_gen, n_par, newpar, Genoprob,
				  Addcov, n_addcov, Intcov, n_intcov, pheno);
    if(verbose) {
      Rprintf("    %3d %10.5f %10.5f", s+1, newllik, newllik - curllik);
      if(newllik < curllik) Rprintf(" ***");
      Rprintf("\n");
    }

    if(newllik - curllik < tol) return(newllik);
    
    for(j=0; j<n_par; j++) 
      curpar[j] = newpar[j];
    curllik = newllik;

    /* e-step */
    for(i=0; i<n_ind; i++) {
      sum=0.0;

      for(j=0; j<n_gen; j++) {
	fit = curpar[j];

	if(n_addcov > 0) {
	  for(k=0; k<n_addcov; k++) 
	    fit += Addcov[k][i] * curpar[n_gen+k];
	}

	if(n_intcov > 0 && j<n_gen-1) {
	  for(k=0; k<n_intcov; k++) 
	    fit += Intcov[k][i] * curpar[n_gen+n_addcov+n_intcov*j+k];
	}
	fit = exp(fit);

	if(pheno[i]) sum += (wts[i][j] = Genoprob[j][pos][i] * fit/(1.0+fit));
	else sum += (wts[i][j] = Genoprob[j][pos][i] / (1.0+fit));
      }
      for(j=0; j<n_gen; j++) 
	wts[i][j] /= sum;
    }

  } /* end of em-step */

  /* didn't converge */
  return(newllik);
}


/**********************************************************************
 * 
 * discan_covar_loglik
 *
 **********************************************************************/

double discan_covar_loglik(int n_ind, int pos, int n_gen, int n_par,
			   double *par, 
			   double ***Genoprob, double **Addcov, int n_addcov,
			   double **Intcov, int n_intcov, int *pheno) 
{
  int i, j, k;
  double loglik=0.0, fit, temp;

  for(i=0; i<n_ind; i++) {
    temp=0.0;

    for(j=0; j<n_gen; j++) {
      fit = par[j];

      if(n_addcov > 0) {
	for(k=0; k<n_addcov; k++) 
	  fit += Addcov[k][i] * par[n_gen+k];
      }

      if(n_intcov > 0 && j<n_gen-1) {
	for(k=0; k<n_intcov; k++) 
	  fit += Intcov[k][i] * par[n_gen+n_addcov+n_intcov*j+k];
      }
      fit = exp(fit);

      if(pheno[i]) temp += Genoprob[j][pos][i] * fit/(1.0+fit);
      else temp += Genoprob[j][pos][i] /(1.0+fit);
    }

    loglik += log10(temp);
  }
  return(loglik);
}

/* end of discan_covar.c */

