/**********************************************************************
 *
 * scantwo_imp.c
 *
 * copyright (c) 2001-6, Karl W Broman and Hao Wu
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
 * with a 2-QTL model by imputation.
 *
 * Contains: R_scantwo_imp, scantwo_imp, altRss2
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
#include "scantwo_imp.h"
#include "scanone_imp.h"

#define TOL 1.0e-12

/**********************************************************************
 *
 * R_scantwo_imp
 *
 * Wrapper for call from R; reorganizes genotype prob, additive and 
 * interactive covariates and result matrix. Then calls scantwo_imp.
 *
 **********************************************************************/

void R_scantwo_imp(int *n_ind, int *same_chr, int *n_pos1, int *n_pos2, 
		   int *n_gen1, int *n_gen2, int *n_draws, int *draws1, 
		   int *draws2, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, int *nphe, 
		   double *weights, double *result, int *n_col2drop,
		   int *col2drop)
{
  int ***Draws1, ***Draws2;
  double **Addcov, **Intcov;

  /* reorganize draws */
  reorg_draws(*n_ind, *n_pos1, *n_draws, draws1, &Draws1);
  if(!(*same_chr)) reorg_draws(*n_ind, *n_pos2, *n_draws, draws2, &Draws2);

  /* reorganize addcov and intcov (if they are not empty) */
  /* currently reorg_geno function is used to reorganized the data */
  if(*n_addcov != 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov != 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  /* call the engine function scantwo_imp */
  scantwo_imp(*n_ind, *same_chr, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
	      *n_draws, Draws1, Draws2, Addcov, *n_addcov, 
	      Intcov, *n_intcov, pheno, *nphe, weights, result,
	      *n_col2drop, col2drop);
}

/**********************************************************************
 * 
 * scantwo_imp
 *
 * Performs genotype pair scan using the pseudomarker algorithm 
 * (imputation) method of Sen and Churchill (2001).
 * 
 * n_ind        Number of individuals
 *
 * same_chr     If = 1, work only with Draws1 and do 2-QTL model with
 *              QTLs on the same chromosome.
 *
 * chr2         Chromesome id 2
 *
 * n_pos1       Number of marker positions in chromesome 1
 *
 * n_pos2       Number of marker positions in chromesome 2
 *
 * n_gen1       Number of different genotypes on chr 1
 *
 * n_gen2       Number of different genotypes on chr 2
 *
 * n_draws      Number of impiutations
 *
 * Draws1       Array of genotype imputations in chromesome 1, 
 *              indexed as Draws1[repl][mar][ind]
 * 
 * Draws2       Array of genotype imputations in chromesome 2, 
 *              indexed as Draws2[repl][mar][ind]
 *
 * addcov	Additive covariates matrix, addcov[mar][ind]
 *
 * n_addcov     Number of additive covariates
 *
 * intcov	Interacting covariates matrix, intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector
 *
 * nphe         Number of phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Result vector of length [n_pos1*n_pos2];
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_imp(int n_ind, int same_chr, int n_pos1, int n_pos2, 
		 int n_gen1, int n_gen2, int n_draws, int ***Draws1, 
		 int ***Draws2, double **Addcov, int n_addcov, 
		 double **Intcov, int n_intcov, double *pheno, int nphe,
		 double *weights, double *result, int n_col2drop,
		 int *col2drop)
{

  /* create local variables */
  int i, i1, i2, j, k; /* loop variants */
  double **lrss0, **lrss1, **LODfull, **LODadd,*lod_tmp;
  double *dwork_null, *dwork_add, *dwork_full, *tmppheno, dtmp;
  int nrss, n_col_null, n_col_a, n_col_f, n_gen_sq, idx;
  int lwork, nlod_per_draw, multivar=0;
  int *allcol2drop;

  /* if number of pheno is 1 or do multivariate model, 
  we have only one rss at each position. Otherwise, 
  we have one rss for each phenotype */
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /* constants */
  n_gen_sq = n_gen1*n_gen2;
  /* number of columns of X for null model */
  n_col_null = 1 + n_addcov;
  /* number of columns of X for additive model */
  n_col_a = (n_gen1+n_gen2-1) + n_addcov + n_intcov*(n_gen1+n_gen2-2);
  /* number of columns of X for full model */
  n_col_f = n_gen_sq + n_addcov + n_intcov*(n_gen_sq-1);

  /* expand col2drop */
  if(n_col2drop) {
    allocate_int(n_col_f, &allcol2drop);
    expand_col2drop(n_gen1, n_addcov, n_intcov, 
		    col2drop, allcol2drop);
  }

  /*********************
   * allocate memory
   *********************/
  tmppheno = (double *)R_alloc(n_ind*nphe, sizeof(double));
  /* for rss' and lod scores - we might not need all of this memory */
  lrss0 = (double **)R_alloc(n_draws, sizeof(double*));
  lrss1 = (double **)R_alloc(n_draws, sizeof(double*));
  LODadd = (double **)R_alloc(n_draws, sizeof(double*));
  LODfull = (double **)R_alloc(n_draws, sizeof(double*));
  for(i=0; i<n_draws; i++) {
    lrss0[i] = (double *)R_alloc(nrss, sizeof(double));
    lrss1[i] = (double *)R_alloc(2*nrss, sizeof(double));
    LODadd[i] = (double *)R_alloc(nrss, sizeof(double));
    LODfull[i] = (double *)R_alloc(nrss, sizeof(double));
  }

  /* the working arrays for the calling of dgelss */
  /* allocate memory */
  /* for null model */
  lod_tmp = (double *)R_alloc(n_draws, sizeof(double));
  lwork = 3*n_col_null + MAX(n_ind, nphe);
  if(multivar == 1) /* request to do multivariate normal model */
    dwork_null = (double *)R_alloc(n_col_null+lwork+2*n_ind*n_col_null+n_ind*nphe+
				   nphe*nphe+n_col_null*nphe,
      sizeof(double));
  else
    dwork_null = (double *)R_alloc(n_col_null+lwork+2*n_ind*n_col_null+n_ind*nphe+n_col_null*nphe,
      sizeof(double));
  /* for additive model */
  lwork = 3*n_col_a + MAX(n_ind, nphe);;
  if(multivar == 1) /* request to do multivariate normal model */
    dwork_add = (double *)R_alloc(n_col_a+lwork+2*n_ind*n_col_a+n_ind*nphe+nphe*nphe+n_col_a*nphe,
      sizeof(double));
  else
    dwork_add = (double *)R_alloc(n_col_a+lwork+2*n_ind*n_col_a+n_ind*nphe+n_col_a*nphe,
      sizeof(double));
  /* for full model */
  lwork = 3*n_col_f + MAX(n_ind, nphe);;
  if(multivar == 1) /* request to do multivariate normal model */
    dwork_full = (double *)R_alloc(n_col_f+lwork+2*n_ind*n_col_f+n_ind*nphe+nphe*nphe+n_col_f*nphe,
      sizeof(double));
  else
    dwork_full = (double *)R_alloc(n_col_f+lwork+2*n_ind*n_col_f+n_ind*nphe+n_col_f*nphe,
      sizeof(double));
  /***************************
   * finish memory allocation
   ***************************/

  /* adjust phenotypes and covariates using weights */
  /* Note: these are actually square-root of weights */
  for(i=0; i<n_ind; i++) {
    for(j=0; j<nphe; j++)     pheno[i+j*n_ind] *= weights[i];
    for(j=0; j<n_addcov; j++) Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++) Intcov[j][i] *= weights[i];
  }

  /* Note that lrss0 is the log10(RSS) for model E(Yi) = b0;
     lrss_add is the log10(RSS) for model 
                 E(Yi) = b0 + b1*q1 + b2*q2;
     lrss_full is the log10(RSS) for model 
                 E(Yi) = b0 + b1*q1 + b2*q2 + b3*(q1*q2); 
     Additive and interactive covariates are included (if any) */

  dtmp = (double)n_ind/2.0; /* this will be used in calculating LOD score */

  /* Call nullRss to calculate the RSS for the null model */
  for (i=0; i<n_draws; i++) {
    /* make a copy of phenotypes. I'm doing this because 
       dgelss will destroy the input rhs array */
    memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
    nullRss(tmppheno, pheno, nphe, n_ind, Addcov, n_addcov,
      dwork_null, multivar, lrss0[i], weights);
  }

  /* calculate the LOD score for each pair of markers */
  if(same_chr) { /* if the pair is on the same chromesome */
    /* number of lod scores per draw */
    nlod_per_draw = n_pos1 * n_pos1;
    for(i1=0; i1<n_pos1-1; i1++) {
      for (i2=i1+1; i2<n_pos1; i2++) {
	for(j=0; j<n_draws; j++) { /* loop over imputations */
	  R_CheckUserInterrupt(); /* check for ^C */

          /* calculate rss */
          memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
          altRss2(tmppheno, pheno, nphe, n_ind, n_gen1, n_gen1, Draws1[j][i1],
		  Draws1[j][i2], Addcov, n_addcov, Intcov, n_intcov,
		  lrss1[j], dwork_add, dwork_full, multivar, weights,
		  n_col2drop, allcol2drop);

	  /* calculate 2 different LOD scores */
          for(k=0; k<nrss; k++) {
            LODadd[j][k] = dtmp*(lrss0[j][k]-lrss1[j][k]);
            LODfull[j][k] = dtmp*(lrss0[j][k]-lrss1[j][k+nrss]);
          }
	}
        /* calculate the weight average on the two LOD score vector
        and fill the result matrix */
        if(n_draws > 1) {
          for(k=0; k<nrss; k++) { /* loop for phenotypes */
            /* for full model LOD */                                            
            for(j=0; j<n_draws; j++)
              lod_tmp[j] = LODfull[j][k];                                       
            result[k*nlod_per_draw+i1*n_pos2+i2] = wtaverage(lod_tmp, n_draws); 
            /* for epistasis LOD */
            for(j=0; j<n_draws; j++)  
              lod_tmp[j] = LODadd[j][k];
            result[k*nlod_per_draw+i2*n_pos1+i1] = wtaverage(lod_tmp, n_draws);
          }
        }
        else { /* only one draw */
          for(k=0;k<nrss; k++) {
	    result[k*nlod_per_draw+i1*n_pos2+i2] = LODfull[0][k];
            result[k*nlod_per_draw+i2*n_pos1+i1] = LODadd[0][k];
          }
        }

      } /* end loop over position 1 */
    } /* end loop over position 2 */
  }

  else { /* the pair is for different chromesome */
    nlod_per_draw = n_pos1*n_pos2;
    idx = n_pos1*n_pos2;
    for(i1=0; i1<n_pos1; i1++) { /* loop over markers on chr 1 */
      for(i2=0; i2<n_pos2; i2++) { /* loop over markers on chr 2 */
	for(j=0; j<n_draws; j++) { /* loop over imputations */
	  R_CheckUserInterrupt(); /* check for ^C */

	  /* rss for alternative model */
          altRss2(tmppheno, pheno, nphe, n_ind, n_gen1, n_gen2, Draws1[j][i1],
		  Draws2[j][i2], Addcov, n_addcov, Intcov, n_intcov,
		  lrss1[j], dwork_add, dwork_full,multivar, weights,
		  n_col2drop, allcol2drop);
          /* calculate 2 different LOD scores */
          for(k=0; k<nrss; k++) {
            LODadd[j][k] = dtmp*(lrss0[j][k]-lrss1[j][k]);
            LODfull[j][k] = dtmp*(lrss0[j][k]-lrss1[j][k+nrss]);
          }
	}
        /* calculate the weight average on the two LOD score vector
	   and fill the result matrix */
        if(n_draws > 1) {
          for(k=0; k<nrss; k++) {
            /* for full model LOD */ 
            for(j=0; j<n_draws; j++) 
              lod_tmp[j] = LODfull[j][k];
            result[(k+nrss)*nlod_per_draw+i1*n_pos2+i2] = wtaverage(lod_tmp, n_draws);
            /* for epistasis LOD */
            for(j=0; j<n_draws; j++)  
              lod_tmp[j] = LODadd[j][k];
            result[k*nlod_per_draw+i2*n_pos1+i1] = wtaverage(lod_tmp, n_draws);
          }
        }
        else { /* only one draw */
          for(k=0;k<nrss; k++) {
	    result[(k+nrss)*nlod_per_draw+i1*n_pos2+i2] = LODfull[0][k]; 
            result[k*nlod_per_draw+i2*n_pos1+i1] = LODadd[0][k];
          }
        }
      } /* end loop over chromesome 2 */
    } /* end loop over chromesome 1 */
  } 
  /* end of scantwo_imp() */
}


/**************************************************************
 *  
 * function to calculate alternative model RSS for two QTL model. 
 * Model is:
 * pheno ~ u+Q1+Q2+addcov+Q1:intcov+Q2:intcov+Q1:Q2+Q1:Q2:intcov
 *
 * This function is called by scantwo_imp 
 *
 **************************************************************/
    
 void altRss2(double *tmppheno, double *pheno, int nphe, int n_ind, int n_gen1, int n_gen2,
	      int *Draws1, int *Draws2, double **Addcov, int n_addcov,
	      double **Intcov, int n_intcov, double *lrss,
	      double *dwork_add, double *dwork_full, int multivar, 
	      double *weights, int n_col2drop, int *allcol2drop)
{
  int i, j, k, s, nrss, lwork, rank, info;
  int n_col_a, n_col_f, n_gen_sq, ind_idx;
  double *x, *x_bk, *singular, *yfit, *work, *rss, *rss_det=0, *coef;
  double alpha=1.0, beta=0.0, tol=TOL, dtmp;
  
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;
  
  /* allocate memory */
  rss = (double *)R_alloc(nrss, sizeof(double));
  
  /* constants */
  /* number of columns for Q1*Q2 */
  n_gen_sq = n_gen1*n_gen2;
  /* number of columns of X for additive model */
  n_col_a = (n_gen1+n_gen2-1) + n_addcov + n_intcov*(n_gen1+n_gen2-2);
  /* number of columns of X for full model */
  n_col_f = n_gen_sq + n_addcov + n_intcov*(n_gen_sq-1);
  
  /**************************************
   * Now work on the additive model
   **************************************/
  /* split the memory block */
  lwork = 3*n_col_a + MAX(n_ind, nphe);
  singular = dwork_add;
  work = singular + n_col_a;
  x = work + lwork;
  x_bk = x + n_ind*n_col_a;
  yfit = x_bk + n_ind*n_col_a;
  coef = yfit + n_ind*nphe;
  if(multivar == 1)
    rss_det = yfit + n_col_a*nphe;
  
  /* zero out X matrix */
  for(i=0; i<n_ind*n_col_a; i++) x[i] = 0.0;
  
  rank = n_col_a;
  /* fill up X matrix */ 
  for(i=0; i<n_ind; i++) {
    x[i+(Draws1[i]-1)*n_ind] = weights[i]; /* QTL 1 */
    s = n_gen1;
    if(Draws2[i] < n_gen2) /* QTL 2 */
      x[i+(Draws2[i]-1+s)*n_ind] = weights[i];
    s += (n_gen2-1);
    for(k=0; k<n_addcov; k++) /* add cov */
      x[i+(k+s)*n_ind] = Addcov[k][i];
    s += n_addcov;
    for(k=0; k<n_intcov; k++) {
      if(Draws1[i] < n_gen1) /* QTL1 x int cov */
        x[i+(Draws1[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen1-1);
      if(Draws2[i] < n_gen2) /* QTL 2 x int cov*/
        x[i+(Draws2[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen2-1);
    }
  } /* end loop over individuals */

  /* drop cols */
  if(n_col2drop)
    dropcol_x(&n_col_a, n_ind, allcol2drop, x);
  rank = n_col_a;

  /* make a copy of x matrix, we may need it */
  memcpy(x_bk, x, n_ind*n_col_a*sizeof(double));
  /* copy pheno to tmppheno, because dgelss will destroy 
     the input pheno array */
  memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
  /* Call LAPACK engine DGELSS to do linear regression.
     Note that DGELSS doesn't have the assumption that X is full rank. */
  /* Pass all arguments to Fortran by reference */
  mydgelss (&n_ind, &n_col_a, &nphe, x, x_bk, pheno, tmppheno,
        singular, &tol, &rank, work, &lwork, &info);

  /* calculate residual sum of squares */
  if(nphe == 1) { /*only one phenotype */
    /* if the design matrix is full rank */
    if(rank == n_col_a)
      for (i=rank, rss[0]=0.0; i<n_ind; i++)
        rss[0] += tmppheno[i]*tmppheno[i];
    else {
        /* the desigm matrix is not full rank, this is trouble */
      /* calculate the fitted value */
      matmult(yfit, x_bk, n_ind, n_col_a, tmppheno, 1);
      /* calculate rss */
      for (i=0, rss[0]=0.0; i<n_ind; i++)
        rss[0] += (pheno[i]-yfit[i]) * (pheno[i]-yfit[i]);
    }
  }

  else { /* multiple phenotypes */
    if(multivar == 1) {
      /* note that the result tmppheno has dimension n_ind x nphe,
      the first ncolx rows contains the estimates. */
      for (i=0; i<nphe; i++)
        memcpy(coef+i*n_col_a, tmppheno+i*n_ind, n_col_a*sizeof(double));
      /* calculate yfit */
      matmult(yfit, x_bk, n_ind, n_col_a, coef, nphe);
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
      if(rank == n_col_a) { /* design matrix is of full rank, this is easier */
        for(i=0; i<nrss; i++) {
          ind_idx = i*n_ind;
          for(j=rank, rss[i]=0.0; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss[i] += dtmp * dtmp;
          }
        }
      }
      else { /* design matrix is singular, this is troubler */
        /* note that the result tmppheno has dimension n_ind x nphe,
        the first ncolx rows contains the estimates. */
        for (i=0; i<nphe; i++)
          memcpy(coef+i*n_col_a, tmppheno+i*n_ind, n_col_a*sizeof(double));
        /* calculate yfit */
        matmult(yfit, x_bk, n_ind, n_col_a, coef, nphe);
        /* calculate residual, put the result in tmppheno */
        for (i=0; i<n_ind*nphe; i++)
          tmppheno[i] = pheno[i] - yfit[i];
        for(i=0; i<nrss; i++) {
          ind_idx = i*n_ind;
          for(j=0, rss[i]=0.0; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss[i] += dtmp * dtmp;
          }
        }
      }
    }
  }      

  /* take log10 */
  for(i=0; i<nrss; i++)
    lrss[i] = log10(rss[i]);


  /**************************************
   * Finish additive model
   **************************************/
    
  /*******************
   * INTERACTIVE MODEL
   *******************/
  /* split the memory block */
  lwork = 3*n_col_f + MAX(n_ind, nphe);
  singular = dwork_full;
  work = singular + n_col_f; 
  x = work + lwork;
  x_bk = x + n_ind*n_col_f;
  yfit = x_bk + n_ind*n_col_f; 
  coef =  yfit + n_ind*nphe;
  if(multivar == 1)
    rss_det = coef + n_col_f*nphe;

  /* zero out X matrix */
  for(i=0; i<n_ind*n_col_f; i++) x[i] = 0.0;

  rank = n_col_f;
  /* fill up X matrix */
  for(i=0; i<n_ind; i++) {
    x[i+(Draws1[i]-1)*n_ind] = weights[i]; /* QTL 1 */
    s = n_gen1;
    if(Draws2[i] < n_gen2) /* QTL 2 */
      x[i+(Draws2[i]-1+s)*n_ind] = weights[i]; 
    s += (n_gen2-1);
    for(k=0; k<n_addcov; k++) /* add cov */
      x[i+(k+s)*n_ind] = Addcov[k][i];
    s += n_addcov;
    for(k=0; k<n_intcov; k++) {
      if(Draws1[i] < n_gen1) /* QTL1 x int cov */
        x[i+(Draws1[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen1-1);
      if(Draws2[i] < n_gen2) /* QTL 2 x int cov */
        x[i+(Draws2[i]-1+s)*n_ind] = Intcov[k][i];
      s += (n_gen2-1);
    }
    if(Draws1[i] < n_gen1 && Draws2[i] < n_gen2) /* QTL x QTL */
      x[i+((Draws1[i]-1)*(n_gen2-1)+Draws2[i]-1+s)*n_ind] = weights[i];
    s += ((n_gen1-1)*(n_gen2-1));
    for(k=0; k<n_intcov; k++) {
      /* QTL x QTL x int cov */
      if(Draws1[i] < n_gen1 && Draws2[i] < n_gen2)
	x[i+((Draws1[i]-1)*(n_gen2-1)+Draws2[i]-1+s)*n_ind] = Intcov[k][i];
      s += ((n_gen1-1)*(n_gen2-1));
    }
  } /* end loop over individuals */

  /* drop x's */
  if(n_col2drop) 
    dropcol_x(&n_col_f, n_ind, allcol2drop, x);
  rank = n_col_f;

  /* make a copy of x matrix, we may need it */
  memcpy(x_bk, x, n_ind*n_col_f*sizeof(double));
  /* copy pheno to tmppheno, because dgelss will destroy 
     the input pheno array */
  memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
  /* Call LAPACK engine DGELSS to do linear regression.
     Note that DGELSS doesn't have the assumption that X is full rank. */
  /* Pass all arguments to Fortran by reference */
  mydgelss (&n_ind, &n_col_f, &nphe, x, x_bk, pheno, tmppheno,
        singular, &tol, &rank, work, &lwork, &info);
  /* calculate residual sum of squares */
  if(nphe == 1) { /* one phenotype */
    /* if the design matrix is full rank */
    if(rank == n_col_f)
      for (i=rank, rss[0]=0.0; i<n_ind; i++)
        rss[0] += tmppheno[i]*tmppheno[i];
    else {
        /* the desigm matrix is not full rank, this is trouble */
      /* calculate the fitted value */
      matmult(yfit, x_bk, n_ind, n_col_f, tmppheno, 1);
        /* calculate rss */
        for (i=0, rss[0]=0.0; i<n_ind; i++)
        rss[0] += (pheno[i]-yfit[i]) * (pheno[i]-yfit[i]);
    }
  }
  else { /* mutlple phenotypes */
    if(multivar == 1) {
      /* multivariate model, rss=det(rss) */
      /* note that the result tmppheno has dimension n_ind x nphe,
      the first ncolx rows contains the estimates. */
      for (i=0; i<nphe; i++)
        memcpy(coef+i*n_col_f, tmppheno+i*n_ind, n_col_f*sizeof(double));
      /* calculate yfit */
      matmult(yfit, x_bk, n_ind, n_col_f, coef, nphe);
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
      if(rank == n_col_f) { /*design matrix is of full rank, this is easier */
        for(i=0; i<nrss; i++) {
          ind_idx = i * n_ind;
          for(j=rank, rss[i]=0.0; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss[i] += dtmp * dtmp;
          }
        }
      } 
      else { /* sigular design matrix */
        /* note that the result tmppheno has dimension n_ind x nphe,
        the first ncolx rows contains the estimates. */
        for (i=0; i<nphe; i++)
          memcpy(coef+i*n_col_f, tmppheno+i*n_ind, n_col_f*sizeof(double));
        /* calculate yfit */
        matmult(yfit, x_bk, n_ind, n_col_f, coef, nphe);
        /* calculate residual, put the result in tmppheno */
        for (i=0; i<n_ind*nphe; i++)
          tmppheno[i] = pheno[i] - yfit[i];
        for(i=0; i<nrss; i++) {
          ind_idx = i * n_ind;
          for(j=0, rss[i]=0.0; j<n_ind; j++) {
            dtmp = tmppheno[ind_idx+j];
            rss[i] += dtmp * dtmp;
          }
        }
      }

    }
  }

  /* take log10 */
  for(i=0; i<nrss; i++)
    lrss[i+nrss] = log10(rss[i]);
}
/* end of altRss2 */
  
  
/* end of scantwo_imp.c */
