/**********************************************************************
 * 
 * vbscan.h
 *
 * copyright (c) 2001, Karl W Broman
 *
 * last modified July, 2001
 * first written May, 2001
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
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
