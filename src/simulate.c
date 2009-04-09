/**********************************************************************
 * 
 * simulate.c
 *
 * copyright (c) 2006, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Jul, 2006
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
 * These functions are for simulating backcross genotype data
 *
 * Contains: sim_bc_ni, sim_bc, R_sim_bc, R_sim_bc_ni
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include "util.h"
#include "simulate.h"

/**********************************************************************
 * 
 * R_sim_bc_ni   Wrapper for sim_bc_ni
 *
 * geno is empty, of size n_mar * n_ind
 *
 **********************************************************************/

void R_sim_bc_ni(int *n_mar, int *n_ind, double *rf, int *geno)
{
  int **Geno;

  reorg_geno(*n_ind, *n_mar, geno, &Geno);

  GetRNGstate();
  sim_bc_ni(*n_mar, *n_ind, rf, Geno);
  PutRNGstate();
}

/**********************************************************************
 * 
 * sim_bc_ni    Simulate backcross under no interference
 *
 * n_mar    Number of markers
 * n_ind    Number of individuals
 * rf       recombination fractions (length n_mar-1)
 * Geno     Matrix of size n_ind x n_mar to contain genotype data
 *
 **********************************************************************/

void sim_bc_ni(int n_mar, int n_ind, double *rf, int **Geno)
{
  int i, j;

  for(i=0; i<n_ind; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    if(unif_rand() < 0.5) Geno[0][i] = 1;
    else Geno[0][i] = 2;

    for(j=1; j<n_mar; j++) {
      if(unif_rand() < rf[j-1])
	Geno[j][i] = 3 - Geno[j-1][i];
      else
	Geno[j][i] = Geno[j-1][i];
    }
  }
}


/**********************************************************************
 * 
 * R_sim_bc   Wrapper for sim_bc
 *
 * geno is empty, of size n_mar * n_ind
 *
 **********************************************************************/

void R_sim_bc(int *n_mar, int *n_ind, double *pos,
	      int *m, double *p, int *geno)
{
  int **Geno;

  reorg_geno(*n_ind, *n_mar, geno, &Geno);

  GetRNGstate();
  sim_bc(*n_mar, *n_ind, pos, *m, *p, Geno);
  PutRNGstate();
}

/**********************************************************************
 * 
 * sim_bc    Simulate backcross under Stahl's interference model
 *
 * n_mar    Number of markers
 * n_ind    Number of individuals
 * pos      Positions of markers (in cM)
 * m        Interference parameter (integer > 0)
 * p        Probability chiasma comes from no interference mechanism
 * Geno     Matrix of size n_ind x n_mar to contain genotype data
 *
 **********************************************************************/

void sim_bc(int n_mar, int n_ind, double *pos, int m, double p, int **Geno)
{
  int i, j, k, first;
  int n_chi, max_chi, n_ni_xo;
  double *chi, L;
  
  L = pos[n_mar-1]; /* length of chromosome in cM */

  /* space to place the crossover locations */
  max_chi = qpois(1e-10, L/50.0*(double)(m+2), 0, 0);
  chi = (double *)Calloc(max_chi, double);

  for(i=0; i<n_ind; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    /* genotype at first marker */
    if(unif_rand() < 0.5) Geno[0][i] = 1;
    else Geno[0][i] = 2;

    /* simulate number of chiasmata and intermediate points */
    n_chi = rpois(L/50.0*(double)(m+1)*(1.0-p));

    /* simulate number of crossovers from ni model */
    if(p > 0)
      n_ni_xo = rpois(L/100.0*p);
    else n_ni_xo = 0;

    if(n_chi + n_ni_xo > max_chi) { /* need more space */
      max_chi = n_chi + n_ni_xo;
      chi = (double *)Realloc(chi, max_chi, double);
    }
      
    /* simulate locations */
    for(j=0; j<n_chi; j++) 
      chi[j] = L*unif_rand();
    R_rsort(chi, n_chi);
    
    /* pull out the locations of chiasmata */
    first = random_int(0, m);
    if(first >= n_chi) n_chi = 0;
    else {
      for(j=first, k=0; j<n_chi; j += (m+1), k++)
        chi[k] = chi[j];
      n_chi = k;
    }

    if(n_chi > 0) {
      /* thin with probability 1/2 */
      for(j=0, k=0; j<n_chi; j++) {
        if(unif_rand() < 0.5) {
  	  chi[k] = chi[j];
	  k++;
	}
      }
      n_chi = k;
    }

    /* add additional crossovers */
    for(j=0; j<n_ni_xo; j++) 
      chi[n_chi + j] = L*unif_rand();
    n_chi += n_ni_xo;

    /* re-sort */
    R_rsort(chi, n_chi);

    /* finally, fill in the genotype data */
    first = 0;
    for(j=1; j<n_mar; j++) {
      while(first < n_chi && chi[first] < pos[j-1]) first++;
      k = 0; /* count crossover in interval */
      while(first < n_chi && chi[first] < pos[j]) {
	k++; first++; 
      }
      first--;
      if(first < 0) { first = 0; }
      
      if(k % 2) /* recombination */
	Geno[j][i] = 3 - Geno[j-1][i];
      else
	Geno[j][i] = Geno[j-1][i];
    }
  }
  Free(chi);
}



/* end of simulate.c */

