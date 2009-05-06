/**********************************************************************
 * 
 * ril48_reorg.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "ril48_reorg.h"
#include "util.h"

/**********************************************************************
 * reorgRIgenoprob
 * 
 * For 4- and 8-way RIL, reorganize the QTL genotype probabilities
 * using the information on the order of the founder strains in each
 * cross.
 **********************************************************************/
void reorgRIgenoprob(int n_ind, int n_mar, int n_str,
		     double ***Prob, int **Crosses)
{
  int i, j, k;
  double *temp;

  temp = (double *)R_alloc(n_str, sizeof(double));

  for(i=0; i<n_ind; i++) {
    for(j=0; j<n_mar; j++) {
      for(k=0; k<n_str; k++) 
	temp[k] = Prob[k][j][i];
      for(k=0; k<n_str; k++) 
	Prob[Crosses[k][i]-1][j][i] = temp[k];
    }
  }
}


/* wrapper for R */
void R_reorgRIgenoprob(int *n_ind, int *n_mar, int *n_str,
		       double *prob, int *crosses)
{
  double ***Prob;
  int **Crosses;

  reorg_genoprob(*n_ind, *n_mar, *n_str, prob, &Prob);
  reorg_geno(*n_ind, *n_str, crosses, &Crosses);
  
  reorgRIgenoprob(*n_ind, *n_mar, *n_str, Prob, Crosses);
}


/**********************************************************************
 * reorgRIdraws
 * 
 * For 4- and 8-way RIL, reorganize the imputed QTL genotypes
 * using the information on the order of the founder strains in each
 * cross.
 **********************************************************************/
void reorgRIdraws(int n_ind, int n_mar, int n_str, int n_draws,
		  int ***Draws, int **Crosses)
{
  int i, j, k;

  for(i=0; i<n_ind; i++) 
    for(j=0; j<n_mar; j++)
      for(k=0; k<n_draws; k++) 
	Draws[k][j][i] = Crosses[Draws[k][j][i]-1][i];
}

/* wrapper for R */
void R_reorgRIdraws(int *n_ind, int *n_mar, int *n_str, int *n_draws,
		    int *draws, int *crosses)
{
  int **Crosses, ***Draws;

  reorg_draws(*n_ind, *n_mar, *n_draws, draws, &Draws);
  reorg_geno(*n_ind, *n_str, crosses, &Crosses);
  
  reorgRIdraws(*n_ind, *n_mar, *n_str, *n_draws, Draws, Crosses);
}

/**********************************************************************
 * reorgRIpairprob
 * 
 * For 4- and 8-way RIL, reorganize the QTL the results of calc.pairprob
 * using the information on the order of the founder strains in each
 * cross.
 **********************************************************************/
void reorgRIpairprob(int n_ind, int n_mar, int n_str,
		     double *****PairProb, int **Crosses)
{
  int i, j1, j2, k1, k2;
  double **temp;

  allocate_dmatrix(n_str, n_str, &temp);

  for(i=0; i<n_ind; i++) {
    for(j1=0; j1<n_mar-1; j1++) {
      for(j2=(j1+1); j2<n_mar; j2++) {

	for(k1=0; k1<n_str; k1++) 
	  for(k2=0; k2<n_str; k2++) 
	    temp[k1][k2] = PairProb[k1][k2][j1][j2][i];

	for(k1=0; k1<n_str; k1++) 
	  for(k2=0; k2<n_str; k2++) 
	    PairProb[Crosses[k1][i]-1][Crosses[k2][i]-1][j1][j2][i] = temp[k1][k2];
      }
    }
  }
}


/* wrapper for R */
void R_reorgRIpairprob(int *n_ind, int *n_mar, int *n_str,
		       double *pairprob, int *crosses)
{
  double *****PairProb;
  int **Crosses;

  reorg_pairprob(*n_ind, *n_mar, *n_str, pairprob, &PairProb);
  reorg_geno(*n_ind, *n_str, crosses, &Crosses);
  
  reorgRIpairprob(*n_ind, *n_mar, *n_str, PairProb, Crosses);
}



/* end of ril48_reorg.c */
