/**********************************************************************
 * 
 * countXO.h
 *
 * copyright (c) 2008-9, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Feb, 2008
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
 * These functions are for comparing marker orders by counts of 
 * obligate crossovers
 *
 * Contains: R_countXO_bc, R_countXO_f2, R_countXO_4way, countXO,
 *           R_countXO_ril48
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * countXO
 *
 * This function counts the number of obligate crossovers for each 
 * individual on a chromosome
 *
 * Input:
 *
 * n_ind    = no. individuals
 *
 * n_mar    = no. markers
 *
 * n_gen    = no. possible genotypes
 *
 * geno     = genotype data [n_ind x n_mar]
 *
 * nxo      = the output; the number of obligate crossovers for each
 *            individual
 *
 * countxo  = function to count the number of obligate crossovers in
 *            an interval and to update the current inferred genotype
 *            (specific for backcross, intercross, and four-way cross)
 *
 **********************************************************************/

void countXO(int n_ind, int n_mar, int n_gen, int *geno,
	     int *nxo, int countxo(int *curgen, int nextgen));


/**********************************************************************
 * 
 * R_countXO_bc
 *
 * Wrapper for call from R for a backcross
 * 
 **********************************************************************/

void R_countXO_bc(int *n_ind, int *n_mar, int *geno, 
		  int *nxo);


/**********************************************************************
 * 
 * R_countXO_f2
 *
 * Wrapper for call from R for an intercross
 * 
 **********************************************************************/

void R_countXO_f2(int *n_ind, int *n_mar, int *geno, 
		  int *nxo);


/**********************************************************************
 * 
 * R_countXO_4way
 *
 * Wrapper for call from R for a four-way cross
 * 
 **********************************************************************/

void R_countXO_4way(int *n_ind, int *n_mar, int *geno, 
		    int *nxo);

/**********************************************************************
 * 
 * R_countXO_ril48
 *
 * Wrapper for call from R for 4- or 8-way RIL by selfing or sibmating
 * 
 **********************************************************************/

void R_countXO_ril48(int *n_ind, int *n_mar, int *geno, 
		     int *nxo);

/* end of countXO.h */

