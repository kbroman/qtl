/**********************************************************************
 * 
 * scantwo_mr.h
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Oct, 2006
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
 * a two-QTL model by marker regression.  Individuals missing genotypes 
 * at either of a pair of markers are dropped.  
 *
 * Contains: R_scantwo_1chr_mr, scantwo_1chr_mr, 
 *           R_scantwo_2chr_mr, scantwo_2chr_mr
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scantwo_1chr_mr
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_mr.
 * 
 **********************************************************************/

void R_scantwo_1chr_mr(int *n_ind, int *n_pos, int *n_gen, int *geno,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights, double *result,
		       int *n_col2drop, int *col2drop);

/**********************************************************************
 * 
 * scantwo_1chr_mr
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
 * Geno         Array of marker genotype data, indexed as
 *              Geno[pos][ind]
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
 * Result       Result matrix of size [n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              testing epistasis.
 *              Note: indexed as Result[col][row]
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_1chr_mr(int n_ind, int n_pos, int n_gen, int **Geno,
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights, double **Result,
		     int n_col2drop, int *col2drop);

/**********************************************************************
 * 
 * R_scantwo_2chr_mr
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_mr.
 * 
 **********************************************************************/

void R_scantwo_2chr_mr(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2,
		       int *geno1, int *geno2,
		       double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights,
		       double *result_full, double *result_add);

/**********************************************************************
 * 
 * scantwo_2chr_mr
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
 * Geno1        Matrix of marker genotype data for chr 1,
 *              indexed as Geno1[pos][ind]
 *
 * Geno2        Matrix of marker genotype data for chr 2
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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos2][pos1]
 *
 * Result_add   Result matrix of size [n_pos2 x n_pos1] 
 *              containing the LODs for add've models
 *              also indexed as Result[pos2][pos1]
 *
 **********************************************************************/

void scantwo_2chr_mr(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, int **Geno1, int **Geno2,
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights,
		     double **Result_full, double **Result_add);

/* end of scantwo_mr.h */
