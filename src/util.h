/**********************************************************************
 *
 * util.h
 *
 * copyright (c) 2001-2010, Karl W Broman and Hao Wu
 *
 * This file written mostly by Karl Broman with some additions
 * from Hao Wu.
 *
 * last modified Nov, 2010
 * first written Feb, 2001
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
 * These are utility functions, mostly for the HMM engine.
 *
 * Other functions: addlog, subtrlog, reorg_geno, reorg_genoprob,
 *                  reorg_pairprob, allocate_int,
 *                  allocate_alpha, reorg_draws, allocate_double,
 *                  sample_int, allocate_imatrix, allocate_dmatrix
 *                  reorg_errlod, double_permute, int_permute,
 *                  random_int
 *                  wtaverage, comparegeno, R_comparegeno
 *                  R_locate_xo, locate_xo, matmult, expand_col2drop
 *                  dropcol_xpx, dropcol_xpy, dropcol_x,
 *                  reviseMWril, R_reviseMWril, R_calcPermPval,
 *                  calcPermPval
 *
 **********************************************************************/

#ifndef __UTIL_H
  #define __UTIL_H

#ifdef __cplusplus
  extern "C" {
#endif

/* Macro for getting maximum */
#ifndef MAX
  #define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif


/**********************************************************************
 *
 * addlog
 *
 * Calculate addlog(a,b) = log[exp(a) + exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.
 *
 **********************************************************************/
double addlog(double a, double b);

/**********************************************************************
 *
 * subtrlog
 *
 * Calculate subtrlog(a,b) = log[exp(a) - exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.
 *
 **********************************************************************/
double subtrlog(double a, double b);

/**********************************************************************
 *
 * reorg_geno
 *
 * Reorganize the genotype data so that it is a doubly indexed array
 * rather than a single long vector
 *
 * Afterwards, geno indexed like Geno[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno);

/**********************************************************************
 *
 * reorg_genoprob
 *
 * Reorganize the genotype probability data so that it is a triply
 * indexed array rather than a single long vector
 *
 * Afterwards, genoprob indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_genoprob(int n_ind, int n_pos, int n_gen,
                    double *genoprob, double ****Genoprob);

/**********************************************************************
 *
 * reorg_pairprob
 *
 * Reorganize the joint genotype probabilities so that they form a
 * quintuply indexed array rather than a single long vector
 *
 * Afterwards, pairprob indexed like
 *    Pairprob[gen1][gen2][pos1][pos2][ind] with pos2 > pos1
 *
 * You *must* refer to cases with pos2 > pos1, as cases with
 * pos2 <= pos1 point off into the ether.
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_pairprob(int n_ind, int n_pos, int n_gen,
                    double *pairprob, double ******Pairprob);

/**********************************************************************
 *
 * allocate_alpha
 *
 * Allocate space for alpha and beta matrices
 *
 * Afterwards, indexed like alpha[gen][mar]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_alpha(int n_pos, int n_gen, double ***alpha);

/**********************************************************************
 *
 * reorg_draws
 *
 * Reorganize the simulated genotypes so that it is a triply
 * indexed array rather than a single long vector
 *
 * Afterwards, draws indexed like Draws[repl][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_draws(int n_ind, int n_pos, int n_draws,
                 int *draws, int ****Draws);

/**********************************************************************
 *
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector);

/**********************************************************************
 *
 * allocate_int
 *
 * Allocate space for a vector of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_int(int n, int **vector);

/**********************************************************************
 *
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_dmatrix(int n_row, int n_col, double ***matrix);

/**********************************************************************
 *
 * allocate_imatrix
 *
 * Allocate space for a matrix of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_imatrix(int n_row, int n_col, int ***matrix);

/**********************************************************************
 *
 * sample_int
 *
 * Make a single draw from (1, ..., n) with probs (p_0, ..., p_(n-1))
 *
 **********************************************************************/
int sample_int(int n, double *p);

/**********************************************************************
 *
 * reorg_errlod
 *
 * Just like reorg_geno(), only for a matrix of doubles.
 *
 * Afterwards, errlod indexed like Errlod[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod);

/**********************************************************************
 *
 * double_permute
 *
 *   This function randomly permutes a vector of doubles
 *
 * Input:
 *
 *   array = vector of doubles; on output, it contains a random
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/
void double_permute(double *array, int len);

/**********************************************************************
 *
 * int_permute
 *
 *   This function randomly permutes a vector of int
 *
 * Input:
 *
 *   array = vector of int; on output, it contains a random
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/
void int_permute(int *array, int len);

/**********************************************************************
 *
 * random_int
 *
 * Generates a random int integer between "low" and "high", inclusive.
 *
 *  Input:
 *
 *    low
 *
 *    high
 *
 **********************************************************************/
int random_int(int low, int high);

/**********************************************************************
 * wtaverage
 * calculate the weight average of the LOD scores
 *********************************************************************/
double wtaverage(double *LOD, int n_draws);


/**********************************************************************
 * comparegeno
 *
 * Count number of matches in the genotypes for all pairs of
 * individuals.
 *
 * Input:
 *
 **********************************************************************/
void comparegeno(int **Geno, int n_ind, int n_mar,
                 int **N_Match, int **N_Missing);

/**********************************************************************
 * R_comparegeno: wrapper for R
 **********************************************************************/
void R_comparegeno(int *geno, int *n_ind, int *n_mar,
                   int *n_match, int *n_missing);

void R_locate_xo(int *n_ind, int *n_mar, int *type,
		 int *geno, double *map, 
		 double *location, int *nseen,
		 int *ileft, int *iright, double *left, double *right,
		 int *ntyped, int *full_info);

/* Note: type ==0 for backcross and ==1 for intercross */
void locate_xo(int n_ind, int n_mar, int type, int **Geno,
	       double *map, double **Location, 
	       int *nseen, int **iLeft, int **iRight,
	       double **Left, double **Right, int **nTyped, 
	       int full_info);

/* multiply two matrices - I'm using dgemm from lapack here */
void matmult(double *result, double *a, int nrowa,
             int ncola, double *b, int ncolb);
/* multiply two matrices - I'm using dgemm from lapack here */
/* void matmult2(double *result, double *a, int nrowa,
               int ncola, double *b, int ncolb); */


/**********************************************************************
 *
 * expand_col2drop
 *
 * Used in scantwo_1chr_em for the X chromosome, to figure out
 * what columns to drop in the presence of covariates when certain
 * genotype columns must be dropped
 *
 **********************************************************************/

void expand_col2drop(int n_gen, int n_addcov, int n_intcov,
                     int *col2drop, int *allcol2drop);

void dropcol_xpx(int *n_col, int *col2drop, double *xpx);

void dropcol_xpy(int n_col, int *col2drop, double *xpy);

void dropcol_x(int *n_col, int n_row, int *col2drop, double *x);

/**********************************************************************
 *
 * reviseMWril    Revise genotypes for 4- or 8-way RIL
 *                to form encoding the founders' genotypes
 *
 * n_ril     Number of RILs to simulate
 * n_mar     Number of markers
 * n_str     Number of founder strains
 *
 * Parents   SNP data for the founder strains [dim n_mar x n_str]
 *
 * Geno      On entry, the detailed genotype data; on exit, the
 *           SNP data written bitwise. [dim n_ril x n_mar]
 *
 * Crosses   The crosses [n_ril x n_str]
 *
 **********************************************************************/
void reviseMWril(int n_ril, int n_mar, int n_str,
                 int **Parents, int **Geno, int **Crosses,
                 int missingval);

/* wrapper for calling reviseMWril from R */
void R_reviseMWril(int *n_ril, int *n_mar, int *n_str,
                   int *parents, int *geno, int *crosses,
                   int *missingval);

/* wrapper for calcPermPval */
void R_calcPermPval(double *peaks, int *nc_peaks, int *nr_peaks,
		    double *perms, int *n_perms, double *pval);

/* calculate permutation p-values for summary.scanone() */
void calcPermPval(double **Peaks, int nc_peaks, int nr_peaks,
		  double **Perms, int n_perms, double **Pval);
#ifdef __cplusplus
  }
#endif

#endif 
/* end of util.h */
