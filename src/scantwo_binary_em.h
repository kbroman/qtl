/**********************************************************************
 * 
 * scantwo_binary_em.h
 *
 * copyright (c) 2004-6, Karl W Broman
 *
 * last modified Oct, 2006
 * first written Dec, 2004
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
 * with a 2-QTL model by interval mapping (the EM algorithm) for  
 * a binary trait.
 *
 * Contains: R_scantwo_1chr_binary_em, scantwo_1chr_binary_em, 
 *           R_scantwo_2chr_binary_em, scantwo_2chr_binary_em,
 *           scantwo_binary_em_estep, scantwo_binary_em_mstep
 *  
 **********************************************************************/

void R_scantwo_1chr_binary_em(int *n_ind, int *n_pos, int *n_gen,
			      double *pairprob, double *addcov, int *n_addcov, 
			      double *intcov, int *n_intcov, int *pheno, 
			      double *start, double *result, int *maxit, 
			      double *tol, int *verbose, int *n_col2drop,
			      int *col2drop);

void scantwo_1chr_binary_em(int n_ind, int n_pos, int n_gen, 
			    double *****Pairprob, double **Addcov, int n_addcov, 
			    double **Intcov, int n_intcov, int *pheno, double *start,
			    double **Result, int maxit, double tol, int verbose,
			    int n_col2drop, int *col2drop);

void R_scantwo_2chr_binary_em(int *n_ind, int *n_pos1, int *n_pos2, 
			      int *n_gen1, int *n_gen2, double *genoprob1,
			      double *genoprob2, double *addcov, int *n_addcov, 
			      double *intcov, int *n_intcov, 
			      int *pheno, double *start,
			      double *result_full, double *result_add,
			      int *maxit, double *tol, int *verbose);

void scantwo_2chr_binary_em(int n_ind, int n_pos1, int n_pos2, int n_gen1, 
			    int n_gen2, double ***Genoprob1, double ***Genoprob2,
			    double **Addcov, int n_addcov, double **Intcov, 
			    int n_intcov, int *pheno, double *start,
			    double **Result_full, double **Result_add, 
			    int maxit, double tol, int verbose);

void scantwo_binary_em_mstep(int n_ind, int n_gen1, int n_gen2, 
			     double **Addcov, int n_addcov, 
			     double **Intcov, int n_intcov, int *pheno, 
			     double ***Wts12, 
			     double *param, int full_model,
			     int n_col, int *error_flag, int n_col2drop,
			     int *allcol2drop, int verbose);

void scantwo_binary_em_estep(int n_ind, int n_gen1, int n_gen2, 
			     double ***Probs, double ***Wts12, 
			     double **Addcov, int n_addcov, double **Intcov,
			     int n_intcov, int *pheno, 
			     double *param, int full_model, int rescale,
			     int n_col2drop, int *allcol2drop);

double scantwo_binary_em_loglik(int n_ind, int n_gen1, int n_gen2, 
				double ***Probs, double **Addcov, int n_addcov,
				double **Intcov, int n_intcov, int *pheno,
				double *param, int full_model,
				int n_col2drop, int *allcol2drop);

/* end of scantwo_binary_em.h */

