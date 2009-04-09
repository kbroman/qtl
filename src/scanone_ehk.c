/**********************************************************************
 * 
 * scanone_ehk.c
 *
 * copyright (c) 2006-7, Karl W Broman
 *
 * last modified Aug, 2007
 * first written Jul, 2006
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
 * single QTL model by the extended Haley-Knott method
 *
 * Contains: R_scanone_ehk, scanone_ehk, calc_mvz
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
#include "scanone_ehk.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_ehk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_ehk.
 * 
 **********************************************************************/

void R_scanone_ehk(int *n_ind, int *n_pos, int *n_gen,
		   double *genoprob, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno,
		   double *weights, double *result, int *maxit,
		   double *tol)
{
  double ***Genoprob, **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scanone_ehk(*n_ind, *n_pos, *n_gen, Genoprob, Addcov, *n_addcov,
	      Intcov, *n_intcov, pheno, weights, result, *maxit, 
	      *tol);
}

/**********************************************************************
 * 
 * scanone_ehk
 *
 * Performs genome scan using the extended Haley-Knott method
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
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Vector of length n_pos, to contain the LOD scores
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 **********************************************************************/

void scanone_ehk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		 double **Addcov, int n_addcov, double **Intcov, 
		 int n_intcov, double *pheno, double *weights, 
		 double *result, int maxit, double tol)
{
  int ny, *jpvt, k, i, j, ncol, k2, k3, k4, s;
  int info, error_flag, flag;
  double *work, *x, *qty, *qraux, *coef, *resid, tol2;
  double sigsq, *coef_cur, sigsq_cur, loglik=0.0, loglik_cur;
  double *m, *v, *z, *y, *wtsq;
  double *rhs, *lhs, **LHS, rcond;

  /* tolerance for linear regression */
  tol2 = TOL;

  ncol = n_gen + (n_gen-1)*n_intcov+n_addcov;

  /* allocate space and set things up*/
  x = (double *)R_alloc(n_ind*ncol, sizeof(double));
  coef = (double *)R_alloc(ncol, sizeof(double));
  coef_cur = (double *)R_alloc(ncol, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(ncol, sizeof(int));
  qraux = (double *)R_alloc(ncol, sizeof(double));
  work = (double *)R_alloc(2 * ncol, sizeof(double));
  m = (double *)R_alloc(n_ind, sizeof(double));
  v = (double *)R_alloc(n_ind, sizeof(double));
  z = (double *)R_alloc(n_ind, sizeof(double));
  y = (double *)R_alloc(n_ind, sizeof(double));
  wtsq = (double *)R_alloc(n_ind, sizeof(double));
  rhs = (double *)R_alloc(ncol, sizeof(double));
  lhs = (double *)R_alloc(ncol*ncol, sizeof(double));

  reorg_errlod(ncol, ncol, lhs, &LHS);
  ny = 1;

  for(j=0; j<n_ind; j++) {
    y[j] = pheno[j]*weights[j];
    wtsq[j] = weights[j]*weights[j];
  }
  /* note: weights are really square-root of weights */

  for(i=0; i<n_pos; i++) { /* loop over positions */
    R_CheckUserInterrupt(); /* check for ^C */

    error_flag = flag = 0;

#ifdef DEBUG
    Rprintf(" ********* position %d *********\n", i+1);
#endif

    /* First do H-K regerssion to get starting values*/
    for(k=0; k<n_gen; k++) jpvt[k] = k;

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
    F77_CALL(dqrls)(x, &n_ind, &ncol, y, &ny, &tol2, coef_cur, resid,
		    qty, &k, jpvt, qraux, work);

    /* estimate sigsq */
    sigsq_cur = 0.0;
    for(j=0; j<n_ind; j++) sigsq_cur += (resid[j]*resid[j]);
    sigsq_cur /= (double)n_ind;

#ifdef DEBUG
    Rprintf("    ");
    for(k=0; k<ncol; k++)
      Rprintf("%.3lf ", coef_cur[k]);
    Rprintf("%.3lf\n", sigsq_cur);
#endif

    calc_mvz(n_ind, i, n_gen, Genoprob, Addcov, n_addcov, 
	     Intcov, n_intcov, pheno, wtsq, coef_cur, 
	     sigsq_cur, m, v, z);

    /* calculate log likelihood */
    loglik_cur=0.0;
    for(j=0; j<n_ind; j++)
      loglik_cur += dnorm(pheno[j], m[j], sqrt(v[j]), 1);

#ifdef DEBUG
    Rprintf("%d %.4lf\n", 0, loglik_cur);
#endif

    for(s=0; s<maxit; s++) {

      R_CheckUserInterrupt(); /* check for ^C */

      /* right-hand side of equations to be solved */
      for(k=0; k<n_gen; k++) {
	rhs[k] = 0.0;
	for(j=0; j<n_ind; j++) 
	  rhs[k] += Genoprob[k][i][j] * pheno[j] / v[j];
      }
      for(k=0; k<n_addcov; k++) {
	rhs[k+n_gen] = 0.0;
	for(j=0; j<n_ind; j++) 
	  rhs[k+n_gen] += Addcov[k][j] * pheno[j] / v[j];
      }
      for(k=0; k<n_gen-1; k++) {
	for(k2=0; k2<n_intcov; k2++) {
	  rhs[n_gen+n_addcov+k*n_intcov+k2] = 0.0;
	  for(j=0; j<n_ind; j++) 
	    rhs[n_gen+n_addcov+k*n_intcov+k2] += 
	      Genoprob[k][i][j] * Intcov[k2][j] * pheno[j] / v[j];
	}
      }

      /* the left-hand sides of the equations to be solve */
      /* Note: ONLY NEED HALF OF THE MATRIX */
      for(k=0; k<ncol; k++)
	for(k2=0; k2<ncol; k2++)
	  LHS[k][k2] = 0.0;




      for(j=0; j<n_ind; j++) {

	for(k=0; k<n_gen; k++) { /* genotype rows */
	  for(k2=k; k2<n_gen; k2++) 
	    LHS[k2][k] += Genoprob[k][i][j] * Genoprob[k2][i][j] * z[j] / v[j];
	  LHS[k][k] -= Genoprob[k][i][j]*(z[j]-1.0)/v[j];

	  for(k2=0; k2<n_addcov; k2++)
	    LHS[k2+n_gen][k] += Genoprob[k][i][j] * Addcov[k2][j] / v[j];
	  
	  for(k3=0; k3<n_intcov; k3++) {
	    for(k2=0; k2<n_gen-1; k2++) 
	      LHS[n_gen+n_addcov+k2*n_intcov+k3][k] += 
		Genoprob[k][i][j] * Genoprob[k2][i][j] * Intcov[k3][j] * z[j]/v[j];
	    if(k < n_gen-1) 
	      LHS[n_gen+n_addcov+k*n_intcov+k3][k] -=
		Genoprob[k][i][j] * Intcov[k3][j] * (z[j]-1.0)/v[j];
	  }
	}

	for(k=0; k<n_addcov; k++) { /* add've covariate rows */
	  for(k2=k; k2<n_addcov; k2++) 
	    LHS[n_gen+k2][n_gen+k] += Addcov[k][j] * Addcov[k2][j] / v[j];
	    
	  for(k2=0; k2<n_gen-1; k2++) 
	    for(k3=0; k3<n_intcov; k3++)
	      LHS[n_gen+n_addcov+k2*n_intcov+k3][n_gen+k] += 
		Addcov[k][j] * Genoprob[k2][i][j] * Intcov[k3][j] / v[j];
	}

	for(k=0; k<n_gen-1; k++) { /* int've covariate rows */
	  for(k2=0; k2<n_intcov; k2++) {
	    
	    for(k4=0; k4<n_intcov; k4++) {
	      for(k3=k; k3<n_gen-1; k3++) {
		LHS[n_gen+n_addcov+k3*n_intcov+k4][n_gen+n_addcov+k*n_intcov+k2] += 
		  Genoprob[k][i][j] * Genoprob[k3][i][j] * Intcov[k2][j] * Intcov[k4][j] * z[j]/v[j];
	      }

	      LHS[n_gen+n_addcov+k*n_intcov+k4][n_gen+n_addcov+k*n_intcov+k2] -= 
		Genoprob[k][i][j] * Intcov[k2][j] * Intcov[k4][j] * (z[j]-1.0)/v[j];
	    }
	  }
	}
      } /* end loop over individuals */
	
#ifdef DEBUG2
      Rprintf("\n");
      for(k=0; k<ncol; k++) {
	for(k2=0; k2<ncol; k2++)
	  Rprintf("%.6lf ", LHS[k][k2]);
	Rprintf("\n");
      }
      Rprintf("\n");
#endif
	

      /* solve rhs * beta = lhs for beta */
      F77_CALL(dpoco)(lhs, &ncol, &ncol, &rcond, coef, &info);
      if(fabs(rcond) < TOL || info != 0) { /* error */
	warning("X'X matrix is singular.");
	error_flag = 1;
	break;
      }
      
      for(k=0; k<ncol; k++) coef[k] = rhs[k];
      F77_CALL(dposl)(lhs, &ncol, &ncol, coef);
      

      calc_mvz(n_ind, i, n_gen, Genoprob, Addcov, n_addcov, 
	       Intcov, n_intcov, pheno, wtsq, coef, 
	       sigsq_cur, m, v, z);

#ifdef DEBUG    
      /* calculate log likelihood */
      loglik=0.0;
      for(j=0; j<n_ind; j++)
	loglik += dnorm(pheno[j], m[j], sqrt(v[j]), 1);

      Rprintf("%d %.4lf\n", s+1, loglik);
#endif

      /* estimate sigsq */
      sigsq = 0.0;
      for(j=0; j<n_ind; j++) 
	sigsq += z[j];
      sigsq *= sigsq_cur/(double)n_ind;
	  
#ifdef DEBUG
    Rprintf("    ");
    for(k=0; k<ncol; k++)
      Rprintf("%.3lf ", coef[k]);
    Rprintf("%.3lf\n", sigsq);
#endif
      calc_mvz(n_ind, i, n_gen, Genoprob, Addcov, n_addcov, 
	       Intcov, n_intcov, pheno, wtsq, coef, 
	       sigsq, m, v, z);

      /* calculate log likelihood */
      loglik=0.0;
      for(j=0; j<n_ind; j++)
	loglik += dnorm(pheno[j], m[j], sqrt(v[j]), 1);

#ifdef DEBUG
      Rprintf("%d %.4lf\n", s+1, loglik);
#endif

      if(fabs(loglik-loglik_cur) < tol) {
	flag = 1;
	break;
      }

      loglik_cur = loglik;
      sigsq_cur = sigsq;
      for(k=0; k<ncol; k++) coef_cur[k] = coef[k];


    } /* end loop over iterations */
    if(error_flag) result[i] = 0.0; 
    else {
      result[i] = loglik;
      if(!flag) warning("Didn't converge at pos %d\n", i+1);
    }
  } /* end loop over positions */

}


/* calc_mvz */
void calc_mvz(int n_ind, int curpos, int n_gen, double ***Genoprob, 
	      double **Addcov, int n_addcov, double **Intcov, int n_intcov, 
	      double *pheno, double *weights, double *coef, double sigsq,
	      double *m, double *v, double *z)
{
  int j, k, k2;
  double temp;
  
  for(j=0; j<n_ind; j++) {

    m[j] = v[j] = 0.0;

    for(k=0; k<n_gen; k++) {
      temp = coef[k];

      if(k < n_gen-1) {
	for(k2=0; k2<n_intcov; k2++) 
	  temp += Intcov[k2][j] * coef[n_gen + n_addcov + n_intcov*k + k2];
      }
	  
      m[j] += Genoprob[k][curpos][j] * temp;
      v[j] += Genoprob[k][curpos][j] * temp*temp;
    }

    v[j] = v[j] - m[j]*m[j] + sigsq/weights[j];

    for(k2=0; k2<n_addcov; k2++) 
      m[j] += Addcov[k2][j] * coef[n_gen + k2];
    
    z[j] = (pheno[j]-m[j])*(pheno[j]-m[j])/v[j];
  }
}


/* end of scanone_ehk.c */
