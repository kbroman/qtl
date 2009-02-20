/**********************************************************************
 * 
 * vbscan.h
 *
 * copyright (c) 2001, Karl W Broman
 *
 * last modified July, 2001
 * first written May, 2001
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License, as
 *     published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version. 
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the
 *     GNU General Public License for more details.
 * 
 *     A copy of the GNU General Public License is available at
 *     http://www.r-project.org/Licenses/
 *
 * C functions for the R/qtl package
 *
 * 
 * C function for performing QTL mapping with the model
 *     p_g = Pr(pheno unobserved | genotype g)
 *     y | pheno observed, g ~ N(mu_g, sigma) 
 *
 * Contains: vbscan, R_vbscan
 *
 **********************************************************************/

void vbscan(int n_pos, int n_ind, int n_gen, double *genoprob, 
	    double *pheno, int *survived, double *lod, int maxit, 
	    double tol);

void R_vbscan(int *n_pos, int *n_ind, int *n_gen, double *genoprob, 
              double *pheno, int *survived, double *lod, int *maxit, 
              double *tol);

/* end of vbscan.h */
