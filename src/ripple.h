/**********************************************************************
 * 
 * ripple.h
 *
 * copyright (c) 2002-9, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Mar, 2002
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
 * Contains: R_ripple_bc, R_ripple_f2, R_ripple_4way, ripple, 
 *           countxo_bc, countxo_f2, countxo_4way
 *           R_ripple_ril48, countxo_ril48
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * ripple
 *
 * This function inspects each of a set of marker orders and counts 
 * the number of obligate crossovers for each order.
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
 * n_orders = no. different orders
 *
 * orders   = matrix of marker orders [n_orders x n_mar]
 *            (Note: each row contains {0, 1, ..., n_mar-1}
 *
 * nxo      = the output; the number of obligate crossovers for each
 *            order (a vector of length n_orders)
 *
 * print_by = How often to print out information?
 *
 * countxo  = function to count the number of obligate crossovers in
 *            an interval and to update the current inferred genotype
 *            (specific for backcross, intercross, and four-way cross)
 *
 **********************************************************************/

void ripple(int n_ind, int n_mar, int *geno,
	    int n_orders, int *orders, int *nxo, 
	    int print_by, int countxo(int *curgen, int nextgen));

/**********************************************************************
 * 
 * R_ripple_bc
 *
 * Wrapper for call from R for a backcross
 * 
 **********************************************************************/

void R_ripple_bc(int *n_ind, int *n_mar, int *geno, 
		 int *n_orders, int *orders,
		 int *nxo, int *print_by);

/**********************************************************************
 * 
 * countxo_bc
 * 
 * count no. obligate crossovers in a backcross
 *
 **********************************************************************/

int countxo_bc(int *curgen, int nextgen);

/**********************************************************************
 * 
 * R_ripple_f2
 *
 * Wrapper for call from R for an intercross
 * 
 **********************************************************************/

void R_ripple_f2(int *n_ind, int *n_mar, int *geno, 
		 int *n_orders, int *orders,
		 int *nxo, int *print_by);

/**********************************************************************
 * 
 * countxo_f2
 * 
 * count no. obligate crossovers in an intercross
 *
 **********************************************************************/

int countxo_f2(int *curgen, int nextgen);

/**********************************************************************
 * 
 * R_ripple_4way
 *
 * Wrapper for call from R for a four-way cross
 * 
 **********************************************************************/

void R_ripple_4way(int *n_ind, int *n_mar, int *geno, 
		   int *n_orders, int *orders,
		   int *nxo, int *print_by);

/**********************************************************************
 * 
 * countxo_4way
 * 
 * count no. obligate crossovers in a four-way cross
 *
 **********************************************************************/

int countxo_4way(int *curgen, int nextgen);

/**********************************************************************
 * 
 * R_ripple_ril48
 *
 * Wrapper for call from R for 4- or 8-way RIL
 * 
 **********************************************************************/

void R_ripple_ril48(int *n_ind, int *n_mar, int *geno, 
		    int *n_orders, int *orders,
		    int *nxo, int *print_by);

/**********************************************************************
 * 
 * countxo_ril48
 * 
 * count no. obligate crossovers in 4- or 8-way RIL
 *
 **********************************************************************/

int countxo_ril48(int *curgen, int nextgen);

/* end of ripple.h */

