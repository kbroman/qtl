/**********************************************************************
 * 
 * scanone_em.h
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Feb, 2006
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
 * single QTL model by interval mapping (the EM algorithm).
 *
 * Contains: R_scanone_em, scanone_em
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_em.
 * 
 **********************************************************************/

void R_scanone_em(int *n_ind, int *n_pos, int *n_gen, 
		  double *genoprob, double *addcov, int *n_addcov,
		  double *intcov, int *n_intcov, double *pheno,
		  double *weights,
		  double *result, int *std_start, double *start,
		  int *maxit, double *tol, int *verbose);

/**********************************************************************
 * 
 * scanone_em
 *
 * Performs genome scan using interval mapping.  (The multipoint
 * genotype probabilities have already been calculated in 
 * calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Upon exit, the LOD scores
 *
 * std_start    If 1, use the usual starting points [initial weights as 
 *                    Pr(QTL geno | marker genotypes)]
 *              If -1, use iid Bernoulli(1/2) for initial weights.
 *              If 0, use the specified values for the means and SD as
 *                    the starting point.
 *
 * start        If std_start = 0, use these as initial estimates of the
 *              genotype-specific means and the residual SD.
 * 
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of dimension 4 x n_gen
 *
 * means        Space for the phenotype means at each genotype
 *
 **********************************************************************/

void scanone_em(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		double *pheno, double *weights, double *result, 
		int std_start, double *start,
		int maxit, double tol, double **work, double *means);

/* end of scanone_em.h */

