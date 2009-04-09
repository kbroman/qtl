/**********************************************************************
 * 
 * scantwo_hk.h
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
 * a two-QTL model by Haley-Knott regression
 *
 * Contains: R_scantwo_1chr_hk, scantwo_1chr_hk, 
 *           R_scantwo_2chr_hk, scantwo_2chr_hk
 *  
 **********************************************************************/

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
		       double *pheno, int *nphe, double *weights, 
		       double *result, int *n_col2drop, int *col2drop);

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

void scantwo_1chr_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		     double *****Pairprob, double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, int nphe,
		     double *weights, double ***Result, int n_col2drop, 
		     int *col2drop);

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
		       double *result_full, double *result_add);

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
 * Result_full  Result matrix of size [n_pos1 x n_pos2]
 *              containing the joint LODs
 *              Note: indexed as Result[pos2][pos1]
 *
 * Result_add   Result matrix of size [n_pos2 x n_pos1] 
 *              containing the LODs for add've models
 *              also indexed as Result[pos2][pos1]
 *
 **********************************************************************/

void scantwo_2chr_hk(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, double ***Genoprob1, 
		     double ***Genoprob2, 
		     double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     int nphe, double *weights,
		     double ***Result_full, double ***Result_add);

/* end of scantwo_hk.h */
