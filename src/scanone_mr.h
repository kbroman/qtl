/**********************************************************************
 * 
 * scanone_mr.h
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Feb, 2006
 * first written Nov, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for performing a genome scan with a 
 * single QTL model by marker regression (i.e., analysis of variance at 
 * the marker loci)
 *
 * Contains: R_scanone_mr, scanone_mr
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_mr
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls scanone_mr.
 * 
 **********************************************************************/

void R_scanone_mr(int *n_ind, int *n_pos, int *n_gen, int *geno, 
		  double *addcov, int *n_addcov, double *intcov, 
		  int *n_intcov, double *pheno, double *weights,
		  double *result);

/**********************************************************************
 * 
 * scanone_mr
 *
 * Performs genome scan using marker regression.
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Geno         Genotype matrix
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
 * result       Vector of length n_pos, to contain the RSS
 *
 **********************************************************************/

void scanone_mr(int n_ind, int n_pos, int n_gen, int **Geno, 
		double **Addcov, int n_addcov, double **Intcov,
		int n_intcov, double *pheno, double *weights, 
		double *result);

/* end of scanone_mr.h */
