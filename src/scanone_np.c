/**********************************************************************
 * 
 * scanone_np.c
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Nov, 2001
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
 * These functions are for performing a non-parametric genome scan 
 * (with a single QTL model), an extension of the Kruskal-Wallis test.
 *
 * Contains: R_scanone_np, scanone_np
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "scanone_np.h"

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
		  double *result)
{
  double ***Genoprob;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);

  scanone_np(*n_ind, *n_pos, *n_gen, Genoprob, pheno, result); 
}

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
		double *result)
{
  int i, j, k;
  double sp, ssp, sr, temp;

  for(i=0; i<n_pos; i++) {

    R_CheckUserInterrupt(); /* check for ^C */

    result[i] = 0.0;
    for(k=0; k<n_gen; k++) {
      sp = ssp = sr = 0.0;

      for(j=0; j<n_ind; j++) {
	sp += Genoprob[k][i][j];
	ssp += (Genoprob[k][i][j]*Genoprob[k][i][j]);
	sr += (Genoprob[k][i][j]*pheno[j]);
      }

      temp = (sr/sp-(double)(n_ind+1)/2.0);
      result[i] += ((double)n_ind-sp)*sp*sp*6.0*temp*temp/
	((double)n_ind*ssp-sp*sp);
    }
    result[i] /= ((double)(n_ind*(n_ind+1))*log(10.0));
  }
}

/* end of scanone_np.c */
