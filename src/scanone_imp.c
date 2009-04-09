/**********************************************************************
 * 
 * scanone_imp.c
 *
 * copyright (c) 2001-6, Karl W Broman and Hao Wu
 *
 * This file is written by Hao Wu 
 * with slight modifications by Karl Broman.
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
 * single QTL model by imputation.  
 *
 * Contains: R_scanone_imp, scanone_imp, nullRss, altRss
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
#include "scanone_imp.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scanone_imp
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_imp.
 * 
 **********************************************************************/

void R_scanone_imp(int *n_ind, int *n_pos, int *n_gen, int *n_draws, 
		   int *draws, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, 
		   int *nphe, double *weights,
		   double *result)
{
  /* reorganize draws */
  int ***Draws;
  double **Addcov, **Intcov, **Result;
  
  reorg_draws(*n_ind, *n_pos, *n_draws, draws, &Draws);
  reorg_errlod(*n_pos, *nphe, result, &Result); 

  /* reorganize addcov and intcov (if they are not empty) */
  /* currently reorg_errlod function is used to reorganize the data */
  if(*n_addcov != 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov != 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);
      
  scanone_imp(*n_ind, *n_pos, *n_gen, *n_draws, Draws, 
	      Addcov, *n_addcov, Intcov, *n_intcov, pheno, *nphe, weights,
	      Result);
}


/**********************************************************************
 * 
 * scanone_imp
 *
 * Performs genome scan using the pseudomarker algorithm (imputation) 
 * method of Sen and Churchill (2001).
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of impiutations
 *
 * Draws        Array of genotype imputations, indexed as 
 *              Draws[repl][mar][ind]
 *
 * Addcov	Additive covariates matrix, Addcov[mar][ind]
 *
 * n_addcov     Number of additive covariates
 *
 * Intcov	Interacting covariates matrix, Intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector/matrix
 *
 * nphe         Number of phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Matrix of size [n_pos x nphe]; upon return, contains
 *              the "LPD" (log posterior distribution of QTL location).
 * 
 **********************************************************************/

void scanone_imp(int n_ind, int n_pos, int n_gen, int n_draws, 
		 int ***Draws, double **Addcov, int n_addcov, 
		 double **Intcov, int n_intcov, double *pheno, 
		 int nphe, double *weights,
		 double **Result)
{

  /* create local variables */
  int i, j, k, itmp, nrss, sizefull, sizenull, lwork, 
    idx, multivar=0, trim=1;
  double **lrss0, **lrss1, *LOD, *lod_tmp, dtmp,
    *tmppheno, *dwork_null, *dwork_full, tol;


  /* if number of pheno is 1 or do multivariate model, 
  we have only one rss at each position. Otherwise, 
  we have one rss for each phenotype */
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /* number of columns in design matrices for null and full model */
  sizenull = 1 + n_addcov;
  sizefull = n_gen + n_addcov + n_intcov*(n_gen-1);

  /* allocate memory */
  lod_tmp = (double *)R_alloc(n_draws, sizeof(double));
  tmppheno = (double *) R_alloc(n_ind*nphe, sizeof(double));
  /* for null model */
  lwork = 3*sizenull + MAX(n_ind, nphe);
  if(multivar == 1) /* request to do multivariate normal model */
    dwork_null = (double *)R_alloc(sizenull+lwork+2*n_ind*sizenull+n_ind*nphe+nphe*nphe+sizenull*nphe, 
      sizeof(double));
  else /* normal model, don't need to allocate memory for rss_det, which is nphe^2 */
    dwork_null = (double *)R_alloc(sizenull+lwork+2*n_ind*sizenull+n_ind*nphe+sizenull*nphe,
      sizeof(double));

  /* for full model */
  lwork = 3*sizefull + MAX(n_ind, nphe);
  if(multivar == 1) /* request to do multivariate normal model */
    dwork_full = (double *)R_alloc(sizefull+lwork+2*n_ind*sizefull+n_ind*nphe+nphe*nphe+sizefull*nphe,
      sizeof(double));
  else /* normal model, don't need to allocate memory for rss_det, which is nphe^2 */
    dwork_full = (double *)R_alloc(sizefull+lwork+2*n_ind*sizefull+n_ind*nphe+sizefull*nphe,
      sizeof(double));
  /* for rss' and lod scores - we might not need all of this memory */
  lrss0 = (double **)R_alloc(n_draws, sizeof(double*));
  lrss1 = (double **)R_alloc(n_draws, sizeof(double*));
  /*LOD = (double **)R_alloc(n_draws, sizeof(double*));*/
  for(i=0; i<n_draws; i++) {
    lrss0[i] = (double *)R_alloc(nrss, sizeof(double));
    lrss1[i] = (double *)R_alloc(nrss, sizeof(double));
    /*LOD[i] = (double *)R_alloc(nrss, sizeof(double));*/
  }
  /* LOD matrix - allocate LOD matrix as a pointer to double, then I can call wtaverage
  directly using pointer operation without looping. This will save some time if there are 
  lots of phenotypes */
  LOD = (double *)R_alloc(n_draws*nrss, sizeof(double));


  /* tolerance for linear regression */
  tol = TOL;

  /* adjust phenotypes and covariates using weights */
  /* Note: these are actually square-root of weights */
  for(i=0; i<n_ind; i++) {
    for(j=0; j<nphe; j++)
      pheno[i+j*n_ind] *= weights[i];
    for(j=0; j<n_addcov; j++) 
      Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++)
      Intcov[j][i] *= weights[i];
  }

  /* calculate the number of LOD needs to be thrown */
  if(trim) idx = (int) floor( 0.5*log(n_draws)/log(2) );
  else idx=0;

  /* Call nullRss to calculate the RSS for the null model */
  for (i=0; i<n_draws; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

  /* make a copy of phenotypes. I'm doing this because 
    dgelss will destroy the input rhs array */
    memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
    nullRss(tmppheno, pheno, nphe, n_ind, Addcov, n_addcov,
      dwork_null, multivar, lrss0[i], weights);
  }
  
  /* calculate the LOD score for each marker */
  dtmp = (double)n_ind/2.0; /* this will be used in calculating LOD score */
  for(i=0; i<n_pos; i++) { /* loop over positions */

    for(j=0; j<n_draws; j++) { /* loop over imputations */
      R_CheckUserInterrupt(); /* check for ^C */

      /* loop over imputations */
      /* call altRss to calcualte the RSS for alternative model,
      given marker and imputatin number */
      memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
      altRss1(tmppheno, pheno, nphe, n_ind, n_gen, Draws[j][i], Addcov,
        n_addcov, Intcov, n_intcov, dwork_full, multivar, lrss1[j], weights);

      /* calculate the LOD score for this marker in this imputation */
      for(k=0; k<nrss; k++) 
        LOD[j+k*n_draws] = dtmp*(lrss0[j][k]-lrss1[j][k]);

    } /* end loop over imputations */

    /* calculate the weight average on the LOD score vector
    and fill the result matrix. Note that result is a matrix
    by ROW. I figured this is the most efficient way to calculate it.
    On exit, we need to use matrix(..., byrow=T) to get the correct one */
    if(n_draws > 1) {
      for(k=0; k<nrss; k++) 
        Result[k][i] = wtaverage(LOD+k*n_draws, n_draws);

    }
    else { 
      itmp = i*nrss;
      for(k=0;k<nrss; k++)
        Result[k][i] = LOD[k];
    } 
   

  } /* end loop over positions */
}


/* function to calculate the null model RSS for scanone_imp */

/**********************************************************
 *
 * function to calculate the residual sum of squares (RSS)
 * for the null model: pheno ~ u + addcov
 * This function is used by scanone_imp and scantwo_imp 
 *
 * Note that if the input pheno is a matrix (multiple columns),
 * a multivariate normal model will be applied, which means
 * the result rss should be det((y-yfit)'*(y-yfit))
 *
 **********************************************************/
void nullRss(double *tmppheno, double *pheno, int nphe, int n_ind,
             double **Addcov, int n_addcov, double *dwork_null,
             int multivar, double *rss0, double *weights)
{
  /* create local variables */
  int i, j, ncolx0, lwork, info, rank, nrss, ind_idx;
  double alpha=1.0, beta=0.0, tol=TOL, dtmp;
  double *work, *x0, *x02, *s, *yfit, *rss_det=0, *coef;

  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /*rss0 = mxCalloc(nrss, sizeof(double));*/

  /* split the memory block */
  ncolx0 = 1 + n_addcov; /* number of columns in x0 matrix */
  /* init rank to be ncolx, which means X is of full rank.
  If it's not, the value of rank will be changed by dgelss */
  rank = ncolx0;

  /* lwork = 3*ncolx0 + n_ind;*/
  lwork = 3*ncolx0 + MAX(n_ind,nphe);
  /* allocate memory */
  s = dwork_null;
  work = s + ncolx0;
  x0 = work + lwork;
  x02 = x0 + n_ind*ncolx0;
  yfit = x02 + n_ind*ncolx0;
  coef = yfit + n_ind*nphe;
  if(multivar == 1)
    rss_det = coef + ncolx0*nphe;

  /* fill up x0 matrix */ 
  for (i=0; i<n_ind; i++) {
    x0[i] = weights[i]; /* the first row (column in Fortran) are all 1s */
    for(j=0; j<n_addcov; j++)
      x0[(j+1)*n_ind+i] = Addcov[j][i];
  }

  /* make a copy of x0 matrix, we may need it */
  memcpy(x02, x0, n_ind*ncolx0*sizeof(double));

  /* Now we have the design matrix x, the model is pheno = x*b
  call LAPACK routine DGELSS to do the linear regression.
  Note that DGELSS doesn't have the assumption that X is full rank. */
  /* Pass all arguments to Fortran by reference */
  mydgelss (&n_ind, &ncolx0, &nphe, x0, x02, pheno, tmppheno,
    s, &tol, &rank, work, &lwork, &info);

  /* calculate residual sum of squares */
  if(nphe == 1) {
    /* if there are only one phenotype, this is easier */
    /* if the design matrix is full rank */
    if(rank == ncolx0) {
      rss0[0] = 0.0;
      for (i=rank; i<n_ind; i++)
        rss0[0] += tmppheno[i]*tmppheno[i];
    }
    else {
      /* the design matrix is not full rank, this is trouble */
      /* calculate the fitted value using yfit=x02*tmppheno(1:ncolx0) */
      matmult(yfit, x02, n_ind, ncolx0, tmppheno, 1);
      /* calculate rss */
      for (i=0; i<n_ind; i++)
        rss0[0] += (pheno[i]-yfit[i]) * (pheno[i]-yfit[i]);
    }
  }

  else { /* multiple phenotypes, this is troubler */
    if(multivar == 1) { /* multivariate model */
      /* note that the result tmppheno has dimension n_ind x nphe,
      the first ncolx0 rows contains the estimates. */
      for (i=0; i<nphe; i++)
        memcpy(coef+i*ncolx0, tmppheno+i*n_ind, ncolx0*sizeof(double));
      /* calculate yfit */
      matmult(yfit, x02, n_ind, ncolx0, coef, nphe);
      /* calculate residual, put the result in tmppheno */
      for (i=0; i<n_ind*nphe; i++)
        tmppheno[i] = pheno[i] - yfit[i];

      /* calcualte rss_det = tmppheno'*tmppheno. */
      /* Call BLAS routine dgemm.  Note that the result rss_det is a 
      symemetric positive definite matrix */
      /* the dimension of tmppheno is n_ind x nphe */
      mydgemm(&nphe, &n_ind, &alpha, tmppheno, &beta, rss_det);
      /* calculate the determinant of rss */
      /* do Cholesky factorization on rss_det */
      mydpotrf(&nphe, rss_det, &info);
      for(i=0, rss0[0]=1.0;i<nphe; i++)
        rss0[0] *= rss_det[i*nphe+i]*rss_det[i*nphe+i];
    }

    else { /* return on rss for each phenotype */
      if(rank == ncolx0) { /* if the design matrix is of full rank, it's easier */
        for(i=0; i<nrss; i++) {
          rss0[i] = 0.0;
          ind_idx = i*n_ind;
          for (j=rank; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss0[i] += dtmp * dtmp; 
          }
        }
      }
      else { /* not full rank, this is trouble */
        /* note that the result tmppheno has dimension n_ind x nphe,
        the first ncolx0 rows contains the estimates. */
        for (i=0; i<nphe; i++)
          memcpy(coef+i*ncolx0, tmppheno+i*n_ind, ncolx0*sizeof(double));
        /* calculate yfit */
        matmult(yfit, x02, n_ind, ncolx0, coef, nphe);
        /* calculate residual, put the result in tmppheno */
        for (i=0; i<n_ind*nphe; i++)
          tmppheno[i] = pheno[i] - yfit[i];
        /* calculate rss */
        for(i=0; i<nrss; i++) {
          rss0[i] = 0.0;
          ind_idx = i*n_ind;
          for(j=0; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss0[i] += dtmp * dtmp; 
          }
        }
      }
    }

  }

  /* take log10 */ 
  for(i=0; i<nrss; i++)
    rss0[i] = log10(rss0[i]);

} 
    

/**************************************************************
 *
 * function to calculate alternative model RSS for one QTL model. 
 * Model is pheno ~ u + Q + addcov + Q:intcov
 *
 * This function is called by scanone_imp 
 *
 **************************************************************/

void altRss1(double *tmppheno, double *pheno, int nphe, int n_ind, int n_gen,
	     int *Draws, double **Addcov, int n_addcov, double **Intcov,
	     int n_intcov, double *dwork, int multivar, double *rss, 
	     double *weights)
{
  /* create local variables */
  int i, j, s, s2, ncolx, lwork, rank, info, nrss, ind_idx;
  double *x, *x_bk, *singular, *yfit, *work, *coef, *rss_det=0;
  double alpha=1.0, beta=0.0, tol=TOL, dtmp;
  /* for lapack dgelss */

  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe; 

  /* number of columns in design matrix X */
  ncolx = n_gen + n_addcov + n_intcov*(n_gen-1);
  /* init rank to be ncolx, which means X is of full rank.
  If it's not, the value of rank will be changed by dgelss */
  rank = ncolx;

  /* split the memory block */
  lwork = 3*ncolx+ MAX(n_ind, nphe);
  /*lwork = 3*ncolx + n_ind;*/
  singular = dwork;
  work = singular + ncolx;
  x = work + lwork;
  x_bk = x + n_ind*ncolx;
  yfit = x_bk + n_ind*ncolx;
  coef = yfit + n_ind*nphe;
  if(multivar == 1)
    rss_det = coef + ncolx*nphe;

  /* zero out X matrix */
  for(i=0; i<n_ind*ncolx; i++) x[i] = 0.0;

  /* fill up design matrix */
  for(i=0; i<n_ind; i++) {
    /* QTL genotypes */
    for(s=0; s<n_gen; s++) {
      if(Draws[i] == s+1) x[i+s*n_ind] = weights[i];
      else x[i+s*n_ind] = 0.0;
    }

    /* Additive covariates */
    for(s=0, s2=n_gen; s<n_addcov; s++, s2++)
      x[i+s2*n_ind] = Addcov[s][i];

    /* Interactive covariates */
    for(s=0; s<n_intcov; s++) {
      for(j=0; j<n_gen-1; j++, s2++) {
        if(Draws[i] == j+1) x[i+n_ind*s2] = Intcov[s][i];
        else x[i+n_ind*s2] = 0.0;
      }
    }

  } /* end loop over individuals */
  /* Done filling up X matrix */

  /* make a copy of x matrix, we may need it */
  memcpy(x_bk, x, n_ind*ncolx*sizeof(double));
  /* Call LAPACK engine DGELSS to do linear regression.
  Note that DGELSS doesn't have the assumption that X is full rank. */
  /* Pass all arguments to Fortran by reference */
  mydgelss(&n_ind, &ncolx, &nphe, x, x_bk, pheno, tmppheno,
    singular, &tol, &rank, work, &lwork, &info);

  /* calculate residual sum of squares */
  if(nphe == 1) {
    /* only one phenotype, this is easier */
    /* if the design matrix is full rank */
    if(rank == ncolx)
      for (i=rank, rss[0]=0.0; i<n_ind; i++)
        rss[0] += tmppheno[i]*tmppheno[i];
      else {
        /* the desigm matrix is not full rank, this is trouble */
        /* calculate the fitted value */
        matmult(yfit, x_bk, n_ind, ncolx, tmppheno, 1);
        /* calculate rss */

        for (i=0, rss[0]=0.0; i<n_ind; i++)
          rss[0] += (pheno[i]-yfit[i]) * (pheno[i]-yfit[i]);
      }
  }
  else { /* multiple phenotypes */
    if(multivar == 1) { /* multivariate normal model, this is troubler */
      /* multivariate model, rss=det(rss) */
      /* note that the result tmppheno has dimension n_ind x nphe,
      the first ncolx rows contains the estimates. */
      for (i=0; i<nphe; i++)
        memcpy(coef+i*ncolx, tmppheno+i*n_ind, ncolx*sizeof(double));
      /* calculate yfit */
      matmult(yfit, x_bk, n_ind, ncolx, coef, nphe);
      /* calculate residual, put the result in tmppheno */
      for (i=0; i<n_ind*nphe; i++)
        tmppheno[i] = pheno[i] - yfit[i];

      /* calcualte rss_det = tmppheno'*tmppheno. */
      /* clear rss_det */ 
      for (i=0; i<nphe*nphe; i++) rss_det[i] = 0.0;
      /* Call BLAS routine dgemm.  Note that the result rss_det is a 
      symemetric positive definite matrix */
      /* the dimension of tmppheno is n_ind x nphe */
      mydgemm(&nphe, &n_ind, &alpha, tmppheno, &beta, rss_det);
      /* calculate the determinant of rss */
      /* do Cholesky factorization on rss_det */
      mydpotrf(&nphe, rss_det, &info);
      for(i=0, rss[0]=1.0;i<nphe; i++)
        rss[0] *= rss_det[i*nphe+i]*rss_det[i*nphe+i];
    } 

    else { /* return rss as a vector */
      if(rank == ncolx) { /* design matrix is of full rank, this is easier */
        for(i=0; i<nrss; i++) {
          rss[i] = 0.0;
          ind_idx = i*n_ind;
          for (j=rank; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss[i] += dtmp*dtmp;
          }
        }
      }
      else { /* design matrix is singular, this is troubler */
        /* note that the result tmppheno has dimension n_ind x nphe,
        the first ncolx rows contains the estimates. */
        for (i=0; i<nphe; i++)
          memcpy(coef+i*ncolx, tmppheno+i*n_ind, ncolx*sizeof(double));
        /* calculate yfit */
        matmult(yfit, x_bk, n_ind, ncolx, coef, nphe);
        /* calculate residual, put the result in tmppheno */
        for (i=0; i<n_ind*nphe; i++)
          tmppheno[i] = pheno[i] - yfit[i];
        /* calculate rss */
        for(i=0; i<nrss; i++) {
          rss[i] = 0.0;
          ind_idx = i*n_ind;
          for(j=0; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss[i] += dtmp * dtmp;
          }
        }
      }
    }
  }

  /* take log10 */
  for(i=0; i<nrss; i++)
    rss[i] = log10(rss[i]);

} /* end of function */

/* end of scanone_imp.c */
