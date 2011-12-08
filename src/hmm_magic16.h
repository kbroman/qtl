/**********************************************************************
 * 
 * hmm_magic16.h
 * 
 * copyright (c) 2011, Karl W Broman
 *
 * last modified Dec, 2011
 * first written Dec, 2011
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
 * Contains: init_magic16, emit_magic16, step_magic16, 
 *           calc_genoprob_magic16, calc_genoprob_special_magic16,
 *           argmax_geno_magic16, sim_geno_magic16,
 *           est_map_magic16, 
 *           marker_loglik_magic16, calc_pairprob_magic16, 
 *           errorlod_magic16, calc_errorlod_magic16
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for 8-way RIL by selfing.
 *
 * Genotype codes:    1-16
 * "Phenotype" codes: 0=missing; otherwise binary 1-65536, with bit i
 *                    indicating SNP compatible with parent i
 *
 **********************************************************************/

double init_magic16(int true_gen);
double emit_magic16(int obs_gen, int true_gen, double error_prob);
double step_magic16(int gen1, int gen2, double rf, double junk);

void calc_genoprob_magic16(int *n_ind, int *n_mar, int *geno, 
			   double *rf, double *error_prob, double *genoprob);

void calc_genoprob_special_magic16(int *n_ind, int *n_mar, int *geno, 
				   double *rf, double *error_prob, double *genoprob);

void argmax_geno_magic16(int *n_ind, int *n_pos, int *geno,
			 double *rf, double *error_prob, int *argmax);

void sim_geno_magic16(int *n_ind, int *n_pos, int *n_draws, int *geno, 
		      double *rf, double *error_prob, int *draws);

void est_map_magic16(int *n_ind, int *n_mar, int *geno, double *rf, 
		     double *error_prob, double *loglik, int *maxit, 
		     double *tol, int *verbose);

void marker_loglik_magic16(int *n_ind, int *geno,
			   double *error_prob, double *loglik);

void calc_pairprob_magic16(int *n_ind, int *n_mar, int *geno, 
			   double *rf, double *error_prob, 
			   double *genoprob, double *pairprob);

double errorlod_magic16(int obs, double *prob, double error_prob);

void calc_errorlod_magic16(int *n_ind, int *n_mar, int *geno, 
			   double *error_prob, double *genoprob, 
			   double *errlod);

/* end of hmm_magic16.h */
