/**********************************************************************
 * 
 * info.h
 *
 * copyright (c) 2001, Karl W Broman
 *
 * last modified Oct, 2001
 * first written Oct, 2001
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
 * This function is for calculating the information contained in the 
 * genotype data on a particular chromosome
 *
 * Contains: R_info 
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_info: calculates information contained in the genotype data for
 *         a particular chromosome.
 *
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * genoprob     Conditional genotype probabilities
 *
 * info1        Vector of length n_pos, to contain the output 
 *              (using the entropy version of the information)
 *
 * info2        Same as info1 (for the prop'n variance version 
 *              of the information)
 *
 * which        0 = only entropy version
 *              1 = only variance version
 *              2 = both
 *
 **********************************************************************/

void R_info(int *n_ind, int *n_pos, int *n_gen, double *genoprob, 
	    double *info1, double *info2, int *which);

/* end of info.h */

