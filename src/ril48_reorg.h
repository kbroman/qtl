/**********************************************************************
 * 
 * ril48_reorg.h
 *
 * copyright (c) 2009, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Apr, 2009
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
 * These functions are for reorganizing results of argmax.geno, calc.genoprob 
 * and sim.geno, for 4- and 8-way RIL
 *
 * Contains: R_reorgRIgenoprob, reorgRIgenoprob, 
 *           R_reorgRIdraws, reorgRIdraws,
 *           R_reorgRIpairprob, reorgRIpairprob
 *  
 **********************************************************************/

/**********************************************************************
 * reorgRIgenoprob
 * 
 * For 4- and 8-way RIL, reorganize the QTL genotype probabilities
 * using the information on the order of the founder strains in each
 * cross.
 **********************************************************************/
void reorgRIgenoprob(int n_ind, int n_mar, int n_str,
		     double ***Prob, int **Crosses);

/* wrapper for R */
void R_reorgRIgenoprob(int *n_ind, int *n_mar, int *n_str,
		       double *prob, int *crosses);

/**********************************************************************
 * reorgRIdraws
 * 
 * For 4- and 8-way RIL, reorganize the imputed QTL genotypes
 * using the information on the order of the founder strains in each
 * cross.
 **********************************************************************/
void reorgRIdraws(int n_ind, int n_mar, int n_str, int n_draws,
		  int ***Draws, int **Crosses);

/* wrapper for R */
void R_reorgRIdraws(int *n_ind, int *n_mar, int *n_str, int *n_draws,
		    int *draws, int *crosses);

/**********************************************************************
 * reorgRIpairprob
 * 
 * For 4- and 8-way RIL, reorganize the QTL the results of calc.pairprob
 * using the information on the order of the founder strains in each
 * cross.
 **********************************************************************/
void reorgRIpairprob(int n_ind, int n_mar, int n_str,
		     double *****PairProb, int **Crosses);

/* wrapper for R */
void R_reorgRIpairprob(int *n_ind, int *n_mar, int *n_str,
		       double *pairprob, int *crosses);

/* end of ril48_reorg.h */
