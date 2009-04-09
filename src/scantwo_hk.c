/**********************************************************************
 * 
 * scantwo_hk.c
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
 * a two-QTL model by Haley-Knott regression
 *
 * Contains: R_scantwo_1chr_hk, scantwo_1chr_hk, 
 *           R_scantwo_2chr_hk, scantwo_2chr_hk
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
#include "scantwo_hk.h"
#define TOL 1e-12

/**********************************************************************
 * 
 * R_scantwo_1chr_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_hk.
 * 
 **********************************************************************/

void R_scantwo_1chr_hk(int *n_ind, int *n_pos, int *n_gen,
		       double *genoprob, double *pairprob, 
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, int* nphe, double *weights, 
		       double *result, int *n_col2drop, int *col2drop)
{
  double ***Genoprob, ***Result, **Addcov, **Intcov, *****Pairprob;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
  reorg_genoprob(*n_pos, *n_pos, *nphe, result, &Result);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_1chr_hk(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob, 
		  Addcov, *n_addcov, Intcov, *n_intcov, 
		  pheno, *nphe, weights, Result, *n_col2drop,
		  col2drop);
}

/**********************************************************************
 * 
 * scantwo_1chr_hk
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
 * nphe         Number of phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [nphe x n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              testing epistasis.
 *              Note: indexed as Result[iphe][col][row]
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_1chr_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		     double *****Pairprob, double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, int nphe,
		     double *weights, double ***Result, int n_col2drop,
		     int *col2drop)
{
  int n_col_a, n_col_f, n_gen_sq, multivar=0, rank=0, n_col_a_temp, n_col_f_temp;
  int itmp, i, i2, j, k, k2, k3, s, nrss, lwork, info, ind_idx;
  /* additive model working arrays */
  /*  double *dwork_add, *x_add, *x_bk_add, *singular_add, *work_add,
   *yfit_add, *coef_add; */
  /* full model working arrays */
    /*  double *dwork_full, *x_full, *x_bk_full, *singular_full, *work_full,
     *yfit_full, *coef_full;*/
  double *dwork, *x, *x_bk, *singular, *work, *yfit, *coef, *tmppheno;
  double tol=TOL, dtmp=0;
  int *allcol2drop;

  /* number of rss */
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /* tolerance for linear regression */
  tol = TOL;
  
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

  /* allocate space and set things up - I will leave multivariate model at this time */
  tmppheno = (double *)R_alloc(n_ind*nphe, sizeof(double));

  /* for full model */
  lwork = 3*n_col_f + MAX(n_ind, nphe);;
  if(multivar == 1) /* request to do multivariate normal model */
    dwork = (double *)R_alloc(n_col_f+lwork+2*n_ind*n_col_f+n_ind*nphe+nphe*nphe+n_col_f*nphe,
				   sizeof(double));
  else
    dwork = (double *)R_alloc(n_col_f + lwork + 2*n_ind*n_col_f + n_ind*nphe + n_col_f*nphe,
				   sizeof(double));
  /* split memory block */
  lwork = 3*n_col_f + MAX(n_ind, nphe);
  singular = dwork;
  work = singular + n_col_f;
  x = work + lwork;
  x_bk = x + n_ind*n_col_f;
  yfit = x_bk + n_ind*n_col_f;
  coef =  yfit + n_ind*nphe;

  /***************************
   * finish memory allocation
   ***************************/

  /* modify pheno, Addcov and Intcov with weights */
  for(i=0; i<n_ind; i++) {
    for(j=0; j<nphe; j++) pheno[i+j*n_ind] *= weights[i];
    for(j=0; j<n_addcov; j++) Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++) Intcov[j][i] *= weights[i];
  }    

  for(i=0; i<n_pos-1; i++) { 
    for(i2=i+1; i2<n_pos; i2++) { /* loop over pairs of positions */

      R_CheckUserInterrupt(); /* check for ^C */

      /* ADDITIVE MODEL */
      rank = n_col_a;
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob[k][i2][j]*weights[j];
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
      rank = n_col_a_temp;

      /* linear regression of phenotype on QTL genotype probabilities */
      /* make a copy of x matrix, we may need it */
      memcpy(x_bk, x, n_ind*n_col_a_temp*sizeof(double));

      /* copy pheno to tmppheno, because dgelss will destroy 
	 the input pheno array */
      memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
      /* Call LAPACK engine DGELSS to do linear regression.
	 Note that DGELSS doesn't have the assumption that X is full rank. */
      /* Pass all arguments to Fortran by reference */
      mydgelss (&n_ind, &n_col_a_temp, &nphe, x, x_bk, pheno, tmppheno,
		singular, &tol, &rank, work, &lwork, &info);
  
      /*
      F77_CALL(dqrls)(x, &n_ind, &n_col_a, pheno, &ny, &tol, coef, resid,
      qty, &k, jpvt, qraux, work);*/
      /* RSS */ 
      /* calculate residual sum of squares */
      if(nphe == 1) { /*only one phenotype */
	/* if the design matrix is full rank */
	if(rank == n_col_a_temp) {
	  for (itmp=rank, dtmp=0.0; itmp<n_ind; itmp++) 
	    dtmp += tmppheno[itmp]*tmppheno[itmp]; 
	}
	else {
	  /* the design matrix is not full rank, this is trouble */
	  /* calculate the fitted value */
	  matmult(yfit, x_bk, n_ind, n_col_a_temp, tmppheno, 1);
	  /* calculate rss */
	  for (itmp=0, dtmp=0.0; itmp<n_ind; itmp++)
	    dtmp += (pheno[itmp]-yfit[itmp]) * (pheno[itmp]-yfit[itmp]);
	}
	Result[0][i2][i] = log10(dtmp);
      }
      else { /* multiple phenotypes */
	if(multivar == 1) { /* multivariate model, I will leave it now */
	}
	else{
	  if(rank == n_col_a_temp) { /* design matrix is of full rank, this is easier */
	    for(itmp=0, ind_idx=0; itmp<nrss; itmp++, ind_idx+=n_ind) { /* loop thru phenotypes */
	      for(j=rank, dtmp=0.0; j<n_ind; j++) 
		dtmp += tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
	      Result[itmp][i2][i] = log10(dtmp); 
	    }
	  }
	  else { /* design matrix is singular, this is troubler */
	    /* note that the result tmppheno has dimension n_ind x nphe,
	       the first ncolx rows contains the estimates. */
	    for (itmp=0; itmp<nphe; itmp++)
	      memcpy(coef+itmp*n_col_a_temp, tmppheno+itmp*n_ind, n_col_a_temp*sizeof(double));
	    /* calculate yfit */
	    matmult(yfit, x_bk, n_ind, n_col_a_temp, coef, nphe);
	    /* calculate residual, put the result in tmppheno */
	    for (itmp=0; itmp<n_ind*nphe; itmp++)
	      tmppheno[itmp] = pheno[itmp] - yfit[itmp];
	    for(itmp=0; itmp<nrss; itmp++) { /* loop thru phenotypes */
	      ind_idx = itmp*n_ind;
	      for(j=0, dtmp=0.0; j<n_ind; j++) 
		dtmp += tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
	      Result[itmp][i2][i] = log10(dtmp); 
	    }
	  }
	}
      }

      /* INTERACTIVE MODEL */
      rank = n_col_f;
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob[k][i][j]*weights[j];  /* s keeps track of column */

	for(k=0; k<n_gen-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob[k][i2][j]*weights[j];

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
	    x[j+s*n_ind] = Pairprob[k][k2][i][i2][j]*weights[j];

	for(k3=0; k3<n_intcov; k3++)
	  for(k=0; k<n_gen-1; k++) /* interactive x QTL 1 x QTL 2 */
	    for(k2=0; k2<n_gen-1; k2++,s++) 
	      x[j+s*n_ind] = Pairprob[k][k2][i][i2][j]*Intcov[k3][j];
      }

      /* drop x's */
      n_col_f_temp = n_col_f;
      if(n_col2drop) 
	dropcol_x(&n_col_f_temp, n_ind, allcol2drop, x);
      rank = n_col_f_temp;

      /* linear regression of phenotype on QTL genotype probabilities */
      /* make a copy of x matrix, we may need it */
      memcpy(x_bk, x, n_ind*n_col_f_temp*sizeof(double));

      /* copy pheno to tmppheno, because dgelss will destroy 
         the input pheno array */
      memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));

      /* Call LAPACK engine DGELSS to do linear regression.
         Note that DGELSS doesn't have the assumption that X is full rank. */
      /* Pass all arguments to Fortran by reference */
      mydgelss (&n_ind, &n_col_f_temp, &nphe, x, x_bk, pheno, tmppheno,
                singular, &tol, &rank, work, &lwork, &info);

      /* calculate residual sum of squares */
      if(nphe == 1) { /*only one phenotype */
        /* if the design matrix is full rank */
        if(rank == n_col_f_temp) { 
          for (itmp=rank, dtmp=0.0; itmp<n_ind; itmp++)
            dtmp += tmppheno[itmp]*tmppheno[itmp];
        }
        else {
          /* the design matrix is not full rank, this is trouble */
          /* calculate the fitted value */

	  /*          matmult(yfit, x_bk, n_ind, n_col_f_temp, tmppheno, 1); */
          matmult(yfit, x_bk, n_ind, n_col_f_temp, tmppheno, 1);

          /* calculate rss */
          for (itmp=0, dtmp=0.0; itmp<n_ind; itmp++)
            dtmp += (pheno[itmp]-yfit[itmp]) * (pheno[itmp]-yfit[itmp]);
        }
        Result[0][i][i2] = log10(dtmp);
      }
      else { /* multiple phenotypes */
        if(multivar == 1) { /* multivariate model */
        }
        else{
          if(rank == n_col_f_temp) { /* design matrix is of full rank, this is easier */
            for(itmp=0; itmp<nrss; itmp++) { /* loop thru phenotypes */
              ind_idx = itmp*n_ind; 
              for(j=rank, dtmp=0.0; j<n_ind; j++) 
                dtmp += tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
	      Result[itmp][i][i2] = log10(dtmp);
            }
          }
          else { /* design matrix is singular, this is troubler */
            /* note that the result tmppheno has dimension n_ind x nphe,
               the first ncolx rows contains the estimates. */
            for (itmp=0; itmp<nphe; itmp++)
              memcpy(coef+itmp*n_col_f_temp, tmppheno+itmp*n_ind, n_col_f_temp*sizeof(double));
            /* calculate yfit */
            matmult(yfit, x_bk, n_ind, n_col_f_temp, coef, nphe);
            /* calculate residual, put the result in tmppheno */
            for (itmp=0; itmp<n_ind*nphe; itmp++)
              tmppheno[itmp] = pheno[itmp] - yfit[itmp];
            for(itmp=0; itmp<nrss; itmp++) { /* loop thru phenotypes */
              ind_idx = itmp*n_ind;
              for(j=0, dtmp=0.0; j<n_ind; j++) 
                dtmp += tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
	      Result[itmp][i][i2] = log10(dtmp);
            }
          }
        }
      }

      /* convert to LODs */
      for(itmp=0; itmp < nphe; itmp++) {
	Result[itmp][i2][i] = (double)n_ind/2.0*Result[itmp][i2][i];
	Result[itmp][i][i2] = (double)n_ind/2.0*
	  Result[itmp][i][i2];
      }

    } /* end loop over positions */
  } 
}

/**********************************************************************
 * 
 * R_scantwo_2chr_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_hk.
 * 
 **********************************************************************/

void R_scantwo_2chr_hk(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2,
		       double *genoprob1, double *genoprob2,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, int *nphe, double *weights,
		       double *result_full, double *result_add)
{
  double ***Genoprob1, ***Genoprob2, ***Result_full, ***Result_add;
  double **Addcov, **Intcov;

  reorg_genoprob(*n_ind, *n_pos1, *n_gen1, genoprob1, &Genoprob1);
  reorg_genoprob(*n_ind, *n_pos2, *n_gen2, genoprob2, &Genoprob2);
  reorg_genoprob(*n_pos2, *n_pos1, *nphe, result_full, &Result_full);
  reorg_genoprob(*n_pos1, *n_pos2, *nphe, result_add, &Result_add);

  /* reorganize addcov and intcov (if they are not empty) */
  if(*n_addcov > 0) reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  if(*n_intcov > 0) reorg_errlod(*n_ind, *n_intcov, intcov, &Intcov);

  scantwo_2chr_hk(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
		  Genoprob1, Genoprob2, Addcov, *n_addcov, Intcov, 
		  *n_intcov, pheno, *nphe, weights, Result_full, Result_add);
}

/**********************************************************************
 * 
 * scantwo_2chr_hk
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
 * nphe         Number of phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result_full  Result matrix of size [nphe x n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[iphe][pos1][pos2]
 *
 * Result_add   Result matrix of size [nphe x n_pos2 x n_pos1] 
 *              containing the LODs for add've models
 *              also indexed as Result[iphe][pos2][pos1]
 *
 **********************************************************************/

void scantwo_2chr_hk(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, double ***Genoprob1, 
		     double ***Genoprob2, 
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     int nphe, double *weights,
		     double ***Result_full, double ***Result_add)
{
  int n_col_a, n_col_f, n_gen_sq, multivar=0, rank=0;
  int itmp, i, i2, j, k, k2, k3, s, nrss, lwork, info, ind_idx;
  double *dwork, *x, *x_bk, *singular, *work, *yfit, *coef, *tmppheno;
  double tol=TOL, dtmp=0;

  /* number of rss */
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /* tolerance for linear regression */
  tol = TOL;

  n_gen_sq = n_gen1*n_gen2;
  /* no. param in additive QTL model */
  n_col_a = (n_gen1+n_gen2-1)+n_addcov+n_intcov*(n_gen1+n_gen2-2); 
  /* no. param full model */
  n_col_f = n_gen_sq+n_addcov+n_intcov*(n_gen_sq-1); 

  /* allocate space and set things up*/
  tmppheno = (double *)R_alloc(n_ind*nphe, sizeof(double));
  lwork = 3*n_col_f + MAX(n_ind, nphe);;
  if(multivar == 1) /* request to do multivariate normal model */
    dwork = (double *)R_alloc(n_col_f+lwork+2*n_ind*n_col_f+
			      n_ind*nphe+nphe*nphe+n_col_f*nphe, sizeof(double));
  else
    dwork = (double *)R_alloc(n_col_f+lwork+2*n_ind*n_col_f+n_ind*nphe+n_col_f*nphe,
			      sizeof(double));
  /* split memory block */
  lwork = 3*n_col_f + MAX(n_ind, nphe);
  singular = dwork;
  work = singular + n_col_f;
  x = work + lwork;
  x_bk = x + n_ind*n_col_f;
  yfit = x_bk + n_ind*n_col_f;
  coef =  yfit + n_ind*nphe;

  /***************************
   * finish memory allocation
   ***************************/

  /* modify pheno, Addcov and Intcov with weights */
  for(i=0; i<n_ind; i++) {
    for(j=0; j<nphe; j++) pheno[i+j*n_ind] *= weights[i];
    for(j=0; j<n_addcov; j++) Addcov[j][i] *= weights[i];
    for(j=0; j<n_intcov; j++) Intcov[j][i] *= weights[i];
  }    

  for(i=0; i<n_pos1; i++) { 
    for(i2=0; i2<n_pos2; i2++) { /* loop over pairs of positions */

      R_CheckUserInterrupt(); /* check for ^C */

      /* ADDITIVE MODEL */
      rank = n_col_a;
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen1; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob1[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen2-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob2[k][i2][j]*weights[j];
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
        if(rank == n_col_a) {
          for (itmp=rank, dtmp=0.0; itmp<n_ind; itmp++)
            dtmp += tmppheno[itmp]*tmppheno[itmp];
        }
        else {
          /* the design matrix is not full rank, this is trouble */
          /* calculate the fitted value */
          matmult(yfit, x_bk, n_ind, n_col_a, tmppheno, 1);
          /* calculate rss */
          for (itmp=0, dtmp=0.0; itmp<n_ind; itmp++)
            dtmp += (pheno[itmp]-yfit[itmp]) * (pheno[itmp]-yfit[itmp]);
        }
        Result_add[0][i2][i] = log10(dtmp);
      }
      else { /* multiple phenotypes */
        if(multivar == 1) { /* multivariate model, I will leave it now */
        }
        else{
          if(rank == n_col_a) { /* design matrix is of full rank, this is easier */
            for(itmp=0, ind_idx=0; itmp<nrss; itmp++, ind_idx+=n_ind) { /* loop thru phenotypes */
              for(j=rank, dtmp=0.0; j<n_ind; j++)
                dtmp += tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
              Result_add[itmp][i2][i] = log10(dtmp);
            }
          }
          else { /* design matrix is singular, this is troubler */
            /* note that the result tmppheno has dimension n_ind x nphe,
               the first ncolx rows contains the estimates. */
            for (itmp=0; itmp<nphe; itmp++)
              memcpy(coef+itmp*n_col_a, tmppheno+itmp*n_ind, n_col_a*sizeof(double));
            /* calculate yfit */
            matmult(yfit, x_bk, n_ind, n_col_a, coef, nphe);
            /* calculate residual, put the result in tmppheno */
            for (itmp=0; itmp<n_ind*nphe; itmp++)
              tmppheno[itmp] = pheno[itmp] - yfit[itmp];
            for(itmp=0; itmp<nrss; itmp++) { /* loop thru phenotypes */
              ind_idx = itmp*n_ind;
              for(j=0, dtmp=0.0; j<n_ind; j++)
                dtmp = tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
              Result_add[itmp][i2][i] = log10(dtmp);
            }
          }
        }
      }

      /* INTERACTIVE MODEL */
      rank = n_col_f;
      /* fill up X matrix */
      for(j=0; j<n_ind; j++) { 
	for(k=0, s=0; k<n_gen1; k++, s++) /* QTL 1 */
	  x[j+s*n_ind] = Genoprob1[k][i][j]*weights[j];  /* s keeps track of column */
	for(k=0; k<n_gen2-1; k++,s++) /* QTL 2 */
	  x[j+s*n_ind] = Genoprob2[k][i2][j]*weights[j];
	for(k=0; k<n_gen1-1; k++)
	  for(k2=0; k2<n_gen2-1; k2++,s++) /* QTL 1 x QTL 2 */
	    x[j+s*n_ind] = Genoprob1[k][i][j]*Genoprob2[k2][i2][j]*weights[j];
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
      /* copy pheno to tmppheno, because dgelss will destroy 
         the input pheno array */
      memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
      /* Call LAPACK engine DGELSS to do linear regression.
         Note that DGELSS doesn't have the assumption that X is full rank. */
      /* Pass all arguments to Fortran by reference */
      mydgelss (&n_ind, &n_col_f, &nphe, x, x_bk, pheno, tmppheno,
                singular, &tol, &rank, work, &lwork, &info);

      /* calculate residual sum of squares */
      if(nphe == 1) { /*only one phenotype */
        /* if the design matrix is full rank */
        if(rank == n_col_f) {
          for (itmp=rank, dtmp=0.0; itmp<n_ind; itmp++)
            dtmp += tmppheno[itmp]*tmppheno[itmp];
        }
        else {
          /* the design matrix is not full rank, this is trouble */
          /* calculate the fitted value */
          matmult(yfit, x_bk, n_ind, n_col_a, tmppheno, 1);
          /* calculate rss */
          for (itmp=0, dtmp=0.0; itmp<n_ind; itmp++)
            dtmp += (pheno[itmp]-yfit[itmp]) * (pheno[itmp]-yfit[itmp]);
        }
        Result_full[0][i][i2] = log10(dtmp);
      }
      else { /* multiple phenotypes */
        if(multivar == 1) { /* multivariate model */
        }
        else{
          if(rank == n_col_f) { /* design matrix is of full rank, this is easier */
            for(itmp=0, ind_idx=0; itmp<nrss; itmp++, ind_idx+=n_ind) { /* loop thru phenotypes */
              for(j=rank, dtmp=0.0; j<n_ind; j++) 
                dtmp += tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
              Result_full[itmp][i][i2] = log10(dtmp);
            }
          }
          else { /* design matrix is singular, this is troubler */
            /* note that the result tmppheno has dimension n_ind x nphe,
               the first ncolx rows contains the estimates. */
            for (itmp=0; itmp<nphe; itmp++)
              memcpy(coef+itmp*n_col_f, tmppheno+itmp*n_ind, n_col_f*sizeof(double));
            /* calculate yfit */
            matmult(yfit, x_bk, n_ind, n_col_f, coef, nphe);
            /* calculate residual, put the result in tmppheno */
            for (itmp=0; itmp<n_ind*nphe; itmp++)
              tmppheno[itmp] = pheno[itmp] - yfit[itmp];
            for(itmp=0; itmp<nrss; itmp++) { /* loop thru phenotypes */
              ind_idx = itmp*n_ind;
              for(j=0, dtmp=0.0; j<n_ind; j++)
                dtmp = tmppheno[ind_idx+j] * tmppheno[ind_idx+j];
              Result_full[itmp][i][i2] = log10(dtmp);
            }
          }
        }
      }

      /* convert to LODs */
      for(itmp=0; itmp < nphe; itmp++) {
	Result_add[itmp][i2][i] = (double)n_ind/2.0*Result_add[itmp][i2][i];
	Result_full[itmp][i][i2] = (double)n_ind/2.0*Result_full[itmp][i][i2];
      }

    } /* end loop over positions */
  } 
}

/* end of scantwo_hk.c */
