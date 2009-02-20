/**********************************************************************
 * 
 * scanone_np.h
 *
 * copyright (c) 2001, Karl W Broman
 *
 * last modified Nov, 2001
 * first written Nov, 2001
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
 * These functions are for performing a non-parametric genome scan 
 * (with a single QTL model), an extension of the Kruskal-Wallis test.
 *
 * Contains: R_scanone_np, scanone_np
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_np
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_np.
 * 
 **********************************************************************/

void R_scanone_np(int *n_ind, int *n_pos, int *n_gen, 
		  double *genoprob, double *pheno,
		  double *result);

/**********************************************************************
 * 
 * scanone_np
 *
 * Performs genome scan using a non-parametric version of 
 * interval mapping.  (The multipoint genotype probabilities have
 * already been calculated in calc.genoprob).
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              (indexed as Genoprob[gen][pos][ind]
 *
 * pheno        Phenotype data, as a vector of ranks
 *
 * result       Result vector of length n_pos (the lod score)
 *
 **********************************************************************/

void scanone_np(int n_ind, int n_pos, int n_gen, 
		double ***Genoprob, double *pheno,
		double *result);

/* end of scanone_np.h */
