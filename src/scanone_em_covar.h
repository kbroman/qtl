/**********************************************************************
 * 
 * scanone_em_covar.h
 *
 * copyright (c) 2001-4, Karl W Broman
 *
 * last modified Nov, 2004
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
 * single QTL model by interval mapping (the EM algorithm) in the 
 * presence of covariates.
 *
 * Contains: scanone_em_covar, estep_em_covar, mstep_em_covar
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * scanone_em_covar
 *
 * Performs genome scan using interval mapping in the presence of
 * covariates.  (The multipoint genotype probabilities have already 
 * been calculated in calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates indexed as 
 *              Addcov[cov][ind]
 *
 * n_addcov     Number of columns in Addcov
 *
 * Intcov       Matrix of interactive covariates indexed as 
 *              Intcov[cov][ind]
 *
 * n_intcov     Number of columns in Intcov
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Result vector of length n_pos; upon return, contains 
 *              the LOD scores.
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * verbose        If 1, print out log likelihood at each iteration
 *
 **********************************************************************/

void scanone_em_covar(int n_ind, int n_pos, int n_gen, 
		      double ***Genoprob, double **Addcov, int n_addcov,
		      double **Intcov, int n_intcov, double *pheno, 
		      double *weights,
		      double *result, int maxit, double tol, int verbose);

/**********************************************************************
 * 
 * mstep_em_covar: M-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen    Number of possible QTL genotypes
 *
 * Addcov   Additive covariates
 *
 * n_addcov Number of columns in Addcov
 *
 * Intcov   Interactive covariates
 *
 * n_intcov Number of columns in Intcov
 *
 * pheno    Phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * wts      Pr(QTL gen | phenotype, model, multipoint marker data),
 *          indexed as wts[gen][ind]
 *
 * param    On output, the updated parameter estimates
 *
 * work1    Workspace of doubles, of length (n_par-1)*(n_par-1)
 *
 * work2    Workspace of doubles, of length (n_par-1)
 *
 * error_flag  Will be set to 1 if the X'X matrix is singular 
 *
 **********************************************************************/

void mstep_em_covar(int n_ind, int n_gen, double **Addcov, int n_addcov, 
		    double **Intcov, int n_intcov, double *pheno, 
		    double *weights,
		    double **wts, double *param, double *work1, 
		    double *work2, int *error_flag);

/**********************************************************************
 * 
 * estep_em_covar: E-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen    Number of possible QTL genotypes
 *
 * pos      Position of Genoprob[][][] to consider
 *
 * Genoprob Pr(QTL gen | multipoint marker data)
 *
 * Addcov   Additive covariates
 *
 * n_addcov Number of columns in Addcov
 *
 * Intcov   Interactive covariates
 *
 * n_intcov Number of columns in Intcov
 *
 * pheno    Phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * wts      On output, Pr(QTL gen | pheno, model, multipt marker data), 
 *          indexed as wts[gen][ind]
 *
 * param    Current parameter estimates
 *
 * rescale  If 1, rescale weights so that the sum to 1.
 *          This is done so that by taking rescale=0, we can easily
 *          calculate the log likelihood 
 *
 **********************************************************************/

void estep_em_covar(int n_ind, int n_gen, int pos, double ***Genoprob,
		    double **Addcov, int n_addcov, double **Intcov,
		    int n_intcov, double *pheno, double *weights,
		    double **wts, double *param, int rescale);

/* end of scanone_em_covar.h */

