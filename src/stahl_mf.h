/**********************************************************************
 * 
 * stahl_mf.h
 * 
 * copyright (c) 2006-7, Karl W Broman
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
 * Contains: mf_stahl, R_mf_stahl, imf_stahl, R_imf_stahl,
 *           imf_stahl_sub
 *
 * These are functions for the calculating the map function for
 *  the Stahl model for crossover interference (with chiasmata
 * coming from two mechanisms: one following a chi-square model
 * and one following a no interference model).
 * m = interference parameter in the chi-square model (m=0 == NI)
 * p = proportion of chiasmata from the NI model (p=1 == NI)
 *
 **********************************************************************/

/* structure used by imf_stahl and imf_stahl_sub */
struct imf_stahl_data {
  double r;
  int m;
  double p;
};

/***********************************************************************
 * R_mf_stahl: wrapper for R
 * 
 * n_d = length of vector d
 * d   = genetic distances (in Morgans)
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * result = vector of length n_d to contain the results
 **********************************************************************/
void R_mf_stahl(int *n_d, double *d, int *m, double *p, double *result);
  
/**********************************************************************
 * mf_stahl: map function for Stahl model
 * 
 * d   = genetic distances (in cM)
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 **********************************************************************/
double mf_stahl(double d, int m, double p);

/**********************************************************************
 * R_imf_stahl: wrapper for R
 * 
 * n_r = length of vector r
 * r   = recombination fractions
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * result = vector of length n_r to contain the results
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
void R_imf_stahl(int *n_r, double *r, int *m, double *p,
		 double *result, double *tol, int *maxit);

/**********************************************************************
 * imf_stahl: inverse map function for chi-square model
 * 
 * r   = recombination fraction
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
double imf_stahl(double r, int m, double p, double tol, int maxit);

/**********************************************************************
 * imf_stahl_sub: utility function for imf_stahl
 * 
 * r   = recombination fraction
 * m   = interference parameter (non-negative integer)
 * p   = proportion of chiasmata from the NI mechanism (double in [0,1])
 * tol = tolerance for convergence
 * maxit = number of interations
 **********************************************************************/
double imf_stahl_sub(double d, struct imf_stahl_data *info);

/* end of stahl_mf.h */
