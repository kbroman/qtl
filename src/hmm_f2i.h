/**********************************************************************
 * 
 * hmm_f2i.h
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
 * Contains: est_map_f2i, R_est_map_f2i, 
 *           emit_f2i, nrec_f2i, step_f2i,
 *
 * These are functions for the HMM under the Stahl model
 * (with chiasmata coming from two mechanisms: one following a 
 * chi-square model and one following a no interference model).
 * m = interference parameter in the chi-square model (m=0 == NI)
 * p = proportion of chiasmata from the NI model (p=1 == NI)
 *
 * Code for is for an intercross.
 *
 * INTERCROSS::
 * Genotype codes:  [0, ..., 2(m+1) - 1] x [1, ..., 2*(m+1)], 
 *                  with the first (m+1) corresponding to A and the 
 *                  others to B, and then for the two chromosomes crossed.
 * Phenotype codes: 0=missing; 1=AA; 2=AB, 3=BB, 4=not BB, 5=not AA
 *
 **********************************************************************/

/**********************************************************************
 * 
 * est_map_f2i
 *
 * This function re-estimates the genetic map for a chromosome
 * with the Stahl model, taking m and p known, for an intercross
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

void est_map_f2i(int n_ind, int n_mar, int *geno, double *d, 
		  int m, double p, double error_prob, 
		  double *loglik, int maxit, double tol, int verbose);

/**********************************************************************
 * emit_f2i: log Pr(obs_gen | true_gen)
 **********************************************************************/
double emit_f2i(int obs_gen, int true_gen, double error_prob,
		int m, int n_bcstates);

/**********************************************************************
 * nrec_f2i: proportion of recombinantion events
 **********************************************************************/
double nrec_f2i(int gen1, int gen2, int m, int n_bcstates);

/* R wrapper for est_map_stahl for intercross */
void R_est_map_f2i(int *n_ind, int *n_mar, int *geno, double *d, 
		   int *m, double *p, double *error_prob, 
		   double *loglik, int *maxit, double *tol, int *verbose);

/**********************************************************************
 * step_f2i
 * 
 * Calculate transition probabilities for Stahl model in an intercross,
 * on the basis of the results for a BC.
 **********************************************************************/
double step_f2i(int g1, int g2, int pos, double ***tm, int n_bcstates);

/* end of hmm_f2i.h */
