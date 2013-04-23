/**********************************************************************
 * 
 * hmm_bcsft.h
 * 
 * copyright (c) 2001-7, Karl W Broman 2011 modified by Brian S Yandell and Laura M Shannon
 *
 * last modified Mar, 2011 BSY, LMS
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
 * Contains: init_bcsft,step_bcsft, init_bcsftb, step_bcsftb,
 *           calc_genoprob_bcsft, calc_genoprob_special_bcsft, sim_geno_bcsft, est_map_bcsft, 
 *           argmax_geno_bcsft, errorlod_bcsft, nrec2_bcsft,
 *           logprec_bcsft, est_rf_bcsft, calc_pairprob_bcsft, marker_loglik_bcsft
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the BCSFT intercross.
 *
 * Genotype codes:  0=AA; 1=AB; 2=BB
 * Phenotype codes: 0=missing; 1=AA; 2=AB; 3=BB; 4=not BB; 5=not AA
 *
 **********************************************************************/

double init_bcsft(int true_gen, int *cross_scheme);

double emit_bcsft(int obs_gen, int true_gen, double error_prob, int *cross_scheme);

double step_bcsft(int gen1, int gen2, double rf, double junk, int *cross_scheme);

double init_bcsftb(int true_gen, int *cross_scheme);

double emit_bcsftb(int obs_gen, int true_gen, double error_prob, int *cross_scheme);

double step_bcsftb(int gen1, int gen2, double rf, double junk, int *cross_scheme);

double nrec_bcsftb(int gen1, int gen2, double rf, int *cross_scheme);

void calc_genoprob_bcsft(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob);
  
void calc_genoprob_special_bcsft(int *n_ind, int *n_mar, int *geno, 
			      double *rf, double *error_prob, double *genoprob);

void sim_geno_bcsft(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws);

void est_map_bcsft(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose);

void argmax_geno_bcsft(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax);

double errorlod_bcsft(int obs, double *prob, double error_prob);

double nrec2_bcsft(int obs1, int obs2, double rf, int *cross_scheme);

double logprec_bcsft(int obs1, int obs2, double rf, int *cross_scheme);

void est_rf_bcsft(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol);

void calc_pairprob_bcsft(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob);

void marker_loglik_bcsft(int *n_ind, int *geno,
		      double *error_prob, double *loglik);

/* end of hmm_bcsft.h */
