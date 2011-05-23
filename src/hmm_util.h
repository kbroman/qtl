/**********************************************************************
 * 
 * hmm_util.c
 * 
 * copyright (c) 2001-9, Karl W Broman
 * modified from hmm_main.c by Brian S Yandell and Laura M Shannon (c) 2011
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
 * Contains: init_stepf, stepfc, forward, backward, golden
 *
 * These are used in hmm_bcsft to simplify calculations.
 * They could be used in hmm_main.c
 *
 **********************************************************************/

void init_stepf(double *rf, double *rf2, int n_gen, int n_mar, int *cross_scheme, 
		double stepf(int, int, double, double, int *),
		double **probmat);

double stepfc(int obs1, int obs2, int mar, double **probmat);

void forward_prob(int i, int n_mar, int n_gen, int curpos, int *cross_scheme, double error_prob,
		  int **Geno, double **probmat, double **alpha,
		  double initf(int, int *), 
		  double emitf(int, int, double, int *));

void backward_prob(int i, int n_mar, int n_gen, int curpos, int *cross_scheme, double error_prob,
		   int **Geno, double **probmat, double **beta,
		   double initf(int, int *), 
		   double emitf(int, int, double, int *));

void calc_probfb(int i, int n_mar, int n_gen, int curpos, double **alpha, double **beta,
		 double ***Genoprob);

double golden_search(double *countmat, int n_gen, int maxit, double tol, int *cross_scheme,
		     double comploglik(double, int, double *, int *));

/* end of hmm_util.h */
