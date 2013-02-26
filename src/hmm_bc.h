/**********************************************************************
 * 
 * hmm_bc.h
 * 
 * copyright (c) 2001-7, Karl W Broman
 *
 * last modified Oct, 2007
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
 * Contains: init_bc, emit_bc, step_bc, nrec_bc, calc_genoprob_bc,
 *           calc_genoprob_special_bc, 
 *           sim_geno_bc, est_map_bc, argmax_geno_bc, errorlod_bc,
 *           calc_errorlod_bc, est_rf_bc, calc_pairprob_bc, 
 *           marker_loglik_bc
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the backcross.
 *
 * Genotype codes:  0=AA; 1=AB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
 *
 **********************************************************************/

double init_bc(int true_gen, int *cross_scheme);

double emit_bc(int obs_gen, int true_gen, double error_prob, int *cross_scheme);

double step_bc(int gen1, int gen2, double rf, double junk, int *cross_scheme);

double nrec_bc(int gen1, int gen2, double rf, int *cross_scheme);

void calc_genoprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob);

void calc_genoprob_special_bc(int *n_ind, int *n_mar, int *geno, 
			      double *rf, double *error_prob, double *genoprob);

void sim_geno_bc(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws);

void est_map_bc(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose);

void argmax_geno_bc(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax);

double errorlod_bc(int obs, double *prob, double error_prob);

void calc_errorlod_bc(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod);

void est_rf_bc(int *n_ind, int *n_mar, int *geno, double *rf);

void calc_pairprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob);

void marker_loglik_bc(int *n_ind, int *geno,
		      double *error_prob, double *loglik);

/* end of hmm_bc.h */
