/**********************************************************************
 * 
 * hmm_bci.h
 * 
 * copyright (c) 2006-7, Karl W Broman
 *         (Some code adapted from code from Nicola Armstrong)
 *
 * last modified Mar, 2007
 * first written Aug, 2006
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
 * Contains: est_map_bci, 
 *           R_est_map_bci, 
 *           emit_bci, nrec_bci, step_bci,
 *           tm_bci, fms_bci, distinct_tm_bci
 *
 * These are functions for the HMM under the Stahl model
 * (with chiasmata coming from two mechanisms: one following a 
 * chi-square model and one following a no interference model).
 * m = interference parameter in the chi-square model (m=0 == NI)
 * p = proportion of chiasmata from the NI model (p=1 == NI)
 *
 * Code for is for a backcross.
 *
 * BACKCROSS:
 * Genotype codes:  0, ..., 2(m+1) - 1, with the first (m+1) 
 *                  corresponding to AA and the others to AB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
 *
 **********************************************************************/

/**********************************************************************
 * 
 * est_map_bci
 *
 * This function re-estimates the genetic map for a chromosome
 * with the Stahl model, taking m and p known, for a backcross
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * d            inter-marker distances in cM
 *              (on exit, contains the new estimates)
 *
 * m            Interference parameter (non-negative integer)
 *
 * p            Proportion of chiasmata from the NI mechanism
 *
 * error_prob   Genotyping error probability
 *
 * loglik       Loglik at final estimates of recombination fractions
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 **********************************************************************/

void est_map_bci(int n_ind, int n_mar, int *geno, double *d, 
		 int m, double p, double error_prob, 
		 double *loglik, int maxit, double tol, int verbose);

/**********************************************************************
 * emit_bci: log Pr(obs_gen | true_gen)
 **********************************************************************/
double emit_bci(int obs_gen, int true_gen, double error_prob,
		int m);

/**********************************************************************
 * nrec_bci: proportion of recombinantion events
 **********************************************************************/
double nrec_bci(int gen1, int gen2, int m);

/* R wrapper for est_map_stahl for backcross */
void R_est_map_bci(int *n_ind, int *n_mar, int *geno, double *d, 
		   int *m, double *p, double *error_prob, 
		   double *loglik, int *maxit, double *tol, int *verbose);

/**********************************************************************
 * step_bci
 * 
 * Calculate transition probabilities (for all intervals) for
 * the Stahl model
 **********************************************************************/
void step_bci(int n_mar, int n_states, double ***tm, double *d, 
	      int m, double p, int maxit, double tol);

/*****************************************************************************
 * tm_bci: this function calculates the required transition probability for the
 * backcross case
 ****************************************************************************/
double tm_bci(int i, int j, double *the_distinct_tm, int m);

/*****************************************************************************
 * fms_bci: this function calculates the sum to infinity part of the
 * transition probabilities for a given lambda_t
 *
 * f should have length 2m+1
 ****************************************************************************/
void fms_bci(double lambda, double *f, int m, double tol, int maxit);

/*****************************************************************************
 * distinct_tm_bci: this function calculates the 3m+2 distinct transition 
 * probabilities for a given lambda_t
 ****************************************************************************/
void distinct_tm_bci(double lambda, double *the_distinct_tm, int m, 
		     double *fms_bci_result);

/* end of hmm_bci.h */
