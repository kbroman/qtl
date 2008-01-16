/**********************************************************************
 *
 * lapackutil.h
 *
 * copyright (c) 2006, Hao Wu
 *
 * last modified Feb, 2006 
 * first written Jan, 2006 
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These are some wrapper functions for several LAPACK routines.
 *
 * Contains: mydgelss, mydgemm, mydpotrf, mydpotrs
 *
 **********************************************************************/

/* DGELS/DGELSS */
void mydgelss (int *n_ind, int *ncolx0, int *nphe, double *x0, double *x0_bk,
               double *pheno, double *tmppheno, double *s, double *tol, 
               int *rank, double *work, int *lwork, int *info);

/* DGEMM */
void mydgemm(int *nphe, int *n_ind, double *alpha, double *tmppheno, double *beta, 
             double *rss_det) ;

/* DPOTRF */
void mydpotrf(int *nphe1, double *rss_det, int *info);

/* DPOTRS */
void mydpotrs(char *uplo, int *n, int *nrhs, double *A, 
              int *lda, double *B, int *ldb, int *info);

/* end of lapackutil.h */
