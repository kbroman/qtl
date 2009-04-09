/**********************************************************************
 * 
 * scantwo_em.h
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
 * These functions are for performing a 2-dimensional genome scan 
 * with a 2-QTL model by interval mapping.(the EM algorithm).
 *
 * Contains: R_scantwo_1chr_em, scantwo_1chr_em, 
 *           R_scantwo_2chr_em, scantwo_2chr_em,
 *           scantwo_em_estep, scantwo_em_mstep, scantwo_em_loglik,
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scantwo_1chr_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_1chr_em.
 * 
 **********************************************************************/

void R_scantwo_1chr_em(int *n_ind, int *n_pos, int *n_gen,
		       double *pairprob, double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights, double *result,
		       int *maxit, double *tol, int *verbose,
		       int *n_col2drop, int *col2drop);

/**********************************************************************
 * 
 * scantwo_1chr_em
 *
 * Performs a 2-dimensional genome scan using the EM algorithm
 * for a two-QTL model with the two QTL residing on the same 
 * chromosome.
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
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
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [n_pos x n_pos]; the lower
 *              triangle (row > col) contains the joint LODs while 
 *              the upper triangle (row < col) contains the LODs for 
 *              testing epistasis.
 *              Note: indexed as Result[col][row]
 *
 * maxit        Maximum number of iterations for EM
 *
 * tol          Tolerance for determining convergence of EM
 *
 * verbose        If >0, print any messages when errors occur
 *                 >1, print out log likelihoods at end of EM
 *                     and check that log likelihood doesn't go down
 *                 >2, print out initial and final log likelihoods
 *                 >3, print out log likelihood at each iteration
 *
 * n_col2drop   For X chromosome, number of columns to drop
 *
 * col2drop     For X chromosome, indicates which columns to drop
 *
 **********************************************************************/

void scantwo_1chr_em(int n_ind, int n_pos, int n_gen, 
		     double *****Pairprob, double **Addcov, int n_addcov, 
		     double **Intcov, int n_intcov, double *pheno, 
		     double *weights,
		     double **Result, int maxit, double tol, int verbose,
		     int n_col2drop, int *col2drop);

/**********************************************************************
 * 
 * R_scantwo_2chr_em
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scantwo_2chr_em.
 * 
 **********************************************************************/

void R_scantwo_2chr_em(int *n_ind, int *n_pos1, int *n_pos2, 
		       int *n_gen1, int *n_gen2, double *genoprob1,
		       double *genoprob2, double *addcov, int *n_addcov, 
		       double *intcov, int *n_intcov, 
		       double *pheno, double *weights,
		       double *result_full, double *result_add,
		       int *maxit, double *tol, int *verbose);

/**********************************************************************
 * 
 * scantwo_2chr_em
 *
 * Performs a 2-dimensional genome scan using the EM algorithm
 * for a two-QTL model with the two QTL residing on the same 
 * chromosome.
 * 
 * n_ind        Number of individuals
 *
 * n_pos1       Number of marker positions on chr 1
 *
 * n_pos2       Number of marker positions on chr 2
 *
 * n_gen1       Number of different genotypes on chr 1
 *
 * n_gen2       Number of different genotypes on chr 2
 *
 * Genoprob1    Array of genotype probabilities for chr 1
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Genoprob2    Array of genotype probabilities for chr 2
 *              indexed as Genoprob[gen][pos][ind]
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
 * maxit        Maximum number of iterations for EM
 *
 * tol          Tolerance for determining convergence of EM
 *
 * verbose        If >0, print any messages when errors occur
 *                 >1, print out log likelihoods at end of EM
 *                     and check that log likelihood doesn't go down
 *                 >2, print out initial and final log likelihoods
 *                 >3, print out log likelihood at each iteration
 *
 **********************************************************************/

void scantwo_2chr_em(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
		     int n_gen2, double ***Genoprob1, double ***Genoprob2,
		     double **Addcov, int n_addcov, double **Intcov, 
		     int n_intcov, double *pheno, double *weights,
		     double **Result_full, double **Result_add, 
		     int maxit, double tol, int verbose);

/**********************************************************************
 * 
 * scantwo_em_mstep: M-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen1   Number of possible genotypes at QTL 1
 *
 * n_gen2   Number of possible genotypes at QTL 2
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
 * Wts12    Pr(QTL1=v, QTL2=w | phenotype, model, marker data),
 *          indexed as Wts[v][w][ind]
 *
 * Wts1     Marginal weights for QTL 1
 * 
 * Wts2     Marginal weights for QTL 2
 *
 * param    On output, the updated parameter estimates (incl resid SD)
 *
 * full_model   If 1, include QTLxQTL interaction
 * 
 * work1    Workspace of doubles, of length (n_par-1)*(n_par-1)
 *
 * work2    Workspace of doubles, of length (n_par-1)
 *
 * error_flag     Set to 1 if X'X is singular
 *
 **********************************************************************/

void scantwo_em_mstep(int n_ind, int n_gen1, int n_gen2, 
		      double **Addcov, int n_addcov, 
		      double **Intcov, int n_intcov, double *pheno, 
		      double *weights,
		      double ***Wts12, double **Wts1, double **Wts2,
		      double *param, int full_model,
		      double *work1, double *work2, int *error_flag,
		      int n_col2drop, int *allcol2drop, int verbose);

/**********************************************************************
 * 
 * scantwo_em_estep: E-step of the EM algorithm 
 *
 * n_ind    Number of individuals
 *
 * n_gen1   Number of possible genotypes at QTL 1
 *
 * n_gen2   Number of possible genotypes at QTL 2
 *
 * Probs    Pr(QTL1=v, QTL2=w | multipoint marker data)
 *          Indexed as Probs[v][w][ind]
 *
 * Wts12    The output:
 *          Pr(QTL1=v, QTL2=w | marker data, phenotype, covar, param)
 *          Indexed as Wts[v][w][ind]
 * 
 * Wts1     Marginal weights for QTL 1
 * 
 * Wts2     Marginal weights for QTL 2
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
 * param    Current parameter estimates (including the resid SD)
 *
 * full_model   If 1, use the full model (with QTLxQTL interaction)
 *
 * rescale  If 1, rescale weights so that the sum to 1.
 *          This is done so that by taking rescale=0, we can easily
 *          calculate the log likelihood 
 *
 **********************************************************************/

void scantwo_em_estep(int n_ind, int n_gen1, int n_gen2, 
		      double ***Probs, double ***Wts12, 
		      double **Wts1, double **Wts2,
		      double **Addcov, int n_addcov, double **Intcov,
		      int n_intcov, double *pheno, double *weights, 
		      double *param, int full_model, int rescale,
		      int n_col2drop, int *allcol2drop);

double scantwo_em_loglik(int n_ind, int n_gen1, int n_gen2, double ***Probs,
			 double ***Wts12, double **Wts1, double **Wts2, 
			 double **Addcov, int n_addcov, double **Intcov,
			 int n_intcov, double *pheno, double *weights,
			 double *param, int full_model,
			 int n_col2drop, int *allcol2drop);


/* end of scantwo_em.h */

