/**********************************************************************
 * 
 * hmm_f2i.c
 * 
 * copyright (c) 2006-2012, Karl W Broman
 *         (Some code adapted from code from Nicola Armstrong)
 *
 * last modified Aug, 2012
 * first written Aug, 2006
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
 * Contains: est_map_f2i, R_est_map_f2i, 
 *           emit_f2i, nrec_f2i, step_f2i,
 *
 * These are functions for the HMM under the Stahl model
 * (with chiasmata coming from two mechanisms: one following a 
 * chi-square model and one following a no interference model).
 * m = interference parameter in the chi-square model (m=0 == NI)
 * p = proportion of chiasmata from the NI model (p=1 == NI)
 *
 * Code for is for an intercross.
 *
 * INTERCROSS::
 * Genotype codes:  [0, ..., 2(m+1) - 1] x [1, ..., 2*(m+1)], 
 *                  with the first (m+1) corresponding to A and the 
 *                  others to B, and then for the two chromosomes crossed.
 * Phenotype codes: 0=missing; 1=AA; 2=AB, 3=BB, 4=not BB, 5=not AA
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include "util.h"
#include "hmm_main.h"
#include "hmm_bci.h"
#include "hmm_f2i.h"
#include "stahl_mf.h"

/**********************************************************************
 * 
 * est_map_f2i
 *
 * This function re-estimates the genetic map for a chromosome
 * with the Stahl model, taking m and p known, for an intercross
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * d            inter-marker distances in cM
 *              (on exit, contains the new estimates)
 *
 * m            Interference parameter (non-negative integer)
 *
 * p            Proportion of chiasmata from the NI mechanism
 *
 * error_prob   Genotyping error probability
 *
 * loglik       Loglik at final estimates of recombination fractions
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 **********************************************************************/

void est_map_f2i(int n_ind, int n_mar, int *geno, double *d, 
		  int m, double p, double error_prob, 
		  double *loglik, int maxit, double tol, int verbose)
{
  int i, j, j2, v, v2, it, flag=0, **Geno, n_states, n_bcstates;
  double s, **alpha, **beta, **gamma, *cur_d, *rf;
  double ***tm, *temp;
  double curloglik;
  double initprob;
  double maxdif, tempdif;
  int ndigits;
  char pattern[100], text[200];
  
  n_bcstates = 2*(m+1);
  n_states = n_bcstates*n_bcstates;  
  initprob = -log((double)n_states);

  /* allocate space for beta and reorganize geno */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  allocate_alpha(n_mar, n_states, &alpha);
  allocate_alpha(n_mar, n_states, &beta);
  allocate_dmatrix(n_states, n_states, &gamma);
  allocate_double(n_mar-1, &cur_d);
  allocate_double(n_mar-1, &rf);

  /* allocate space for the [backcross] transition matrices */
  /* size n_states x n_states x (n_mar-1) */
  /* tm[state1][state2][interval] */
  tm = (double ***)R_alloc(n_bcstates, sizeof(double **));
  tm[0] = (double **)R_alloc(n_bcstates * n_bcstates, sizeof(double *));
  for(i=1; i<n_bcstates; i++) tm[i] = tm[i-1] + n_bcstates;
  tm[0][0] = (double *)R_alloc(n_bcstates * n_bcstates * (n_mar - 1), 
			       sizeof(double));
  temp = tm[0][0];
  for(i=0; i < n_bcstates; i++) {
    for(j=0; j < n_bcstates; j++) {
      tm[i][j] = temp;
      temp += n_mar-1;
    }
  }

  /* digits in verbose output */
  if(verbose) {
    ndigits = (int)ceil(-log10(tol));
    if(ndigits > 16) ndigits=16;
    sprintf(pattern, "%s%d.%df", "%", ndigits+3, ndigits+1);
  }

  for(j=0; j<n_mar-1; j++) d[j] /= 100.0; /* convert to Morgans */

  /* begin EM algorithm */
  for(it=0; it<maxit; it++) {

    for(j=0; j<n_mar-1; j++) {
      cur_d[j] = d[j];
      rf[j] = 0.0;
    }

    /* calculate the transition matrices [for BC] */
    step_bci(n_mar, n_bcstates, tm, cur_d, m, p, maxit, tol);

    for(i=0; i<n_ind; i++) { /* i = individual */

      R_CheckUserInterrupt(); /* check for ^C */

      /* initialize alpha and beta */
      for(v=0; v<n_states; v++) {
	alpha[v][0] = initprob + emit_f2i(Geno[0][i], v, error_prob, m, n_bcstates);
	beta[v][n_mar-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	
	for(v=0; v<n_states; v++) {
	  alpha[v][j] = alpha[0][j-1] + step_f2i(0, v, j-1, tm, n_bcstates);
	  
	  beta[v][j2] = beta[0][j2+1] + step_f2i(v, 0, j2, tm, n_bcstates) +
	    emit_f2i(Geno[j2+1][i], 0, error_prob, m, n_bcstates);
	  
	  for(v2=1; v2<n_states; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 step_f2i(v2, v, j-1, tm, n_bcstates));
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				 step_f2i(v, v2, j2, tm, n_bcstates) +
				 emit_f2i(Geno[j2+1][i], v2, error_prob, m, n_bcstates));
	  }
	  
	  alpha[v][j] += emit_f2i(Geno[j][i], v, error_prob, m, n_bcstates);
		 
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_states; v++) {
	  for(v2=0; v2<n_states; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + 
	      emit_f2i(Geno[j+1][i], v2, error_prob, m, n_bcstates) +
	      step_f2i(v, v2, j, tm, n_bcstates);

	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}

	for(v=0; v<n_states; v++) {
	  for(v2=0; v2<n_states; v2++) {
	    rf[j] += nrec_f2i(v, v2, m, n_bcstates) * exp(gamma[v][v2] - s);
	  }
	}
      }

    } /* loop over individuals */

    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf[j] /= (double)n_ind;

      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;

    }

    /* use map function to convert back to distances */
    for(j=0; j<n_mar-1; j++)
      d[j] = imf_stahl(rf[j], m, p, 1e-10, 1000);

    if(verbose>1) {
      /* print estimates as we go along*/
      Rprintf("   %4d ", it+1);
      maxdif=0.0;
      for(j=0; j<n_mar-1; j++) {
	tempdif = fabs(d[j] - cur_d[j])/(cur_d[j]+tol*100.0);
	if(maxdif < tempdif) maxdif = tempdif;
      }
      sprintf(text, "%s%s\n", "  max rel've change = ", pattern);
      Rprintf(text, maxdif);
    }

    /* check convergence */
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(d[j] - cur_d[j]) > tol*(cur_d[j]+tol*100.0)) {
	flag = 1; 
	break;
      }
    }

    if(!flag) break;

  } /* end EM algorithm */
  
  if(flag) warning("Didn't converge!\n");

  /* re-calculate transition matrices */
  step_bci(n_mar, n_bcstates, tm, d, m, p, maxit, tol);

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */
    /* initialize alpha */
    for(v=0; v<n_states; v++) 
      alpha[v][0] = initprob + emit_f2i(Geno[0][i], v, error_prob, m, n_bcstates);

    /* forward equations */
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_states; v++) {
	alpha[v][j] = alpha[0][j-1] + step_f2i(0, v, j-1, tm, n_bcstates);
	for(v2=1; v2<n_states; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       step_f2i(v2, v, j-1, tm, n_bcstates));
	alpha[v][j] += emit_f2i(Geno[j][i], v, error_prob, m, n_bcstates);
      }
    }

    curloglik = alpha[0][n_mar-1];
    for(v=1; v<n_states; v++) 
      curloglik = addlog(curloglik, alpha[v][n_mar-1]);
    *loglik += curloglik;
  }

  if(verbose) {
    if(verbose < 2) {
      /* print final estimates */
      Rprintf("  no. iterations = %d\n", it+1);
      maxdif=0.0;
      for(j=0; j<n_mar-1; j++) {
	tempdif = fabs(d[j] - cur_d[j])/(cur_d[j]+tol*100.0);
	if(maxdif < tempdif) maxdif = tempdif;
      }
      sprintf(text, "%s%s\n", "  max rel've change at last step = ", pattern);
      Rprintf(text, maxdif);
    }
    
    Rprintf("  loglik: %10.4lf\n\n", *loglik);
  }

  /* convert distances back to cM */
  for(j=0; j<n_mar-1; j++) d[j] *= 100.0;
}


/**********************************************************************
 * emit_f2i: log Pr(obs_gen | true_gen)
 **********************************************************************/
double emit_f2i(int obs_gen, int true_gen, double error_prob,
		int m, int n_bcstates)
{
  if(obs_gen==0) return(0.0);

  /* genotype as 1,2,3 */
  true_gen = ((true_gen / n_bcstates) / (m+1)) + ((true_gen % n_bcstates) / (m+1)) + 1;

  switch(obs_gen) {
  case 1: case 2: case 3:
    
    if(true_gen == obs_gen) return(log(1.0-error_prob));
    else return(log(error_prob)-M_LN2);

  case 4:
    if(true_gen != 3) return(log(1.0 - error_prob / 2.0));
    else return(log(error_prob)-M_LN2);

  case 5:
    if(true_gen != 1) return(log(1.0 - error_prob / 2.0));
    else return(log(error_prob)-M_LN2);
  }

  return(0.0); /* shouldn't get here */
}



/**********************************************************************
 * nrec_f2i: proportion of recombinantion events
 **********************************************************************/
double nrec_f2i(int gen1, int gen2, int m, int n_bcstates)
{
  int mom1, dad1, mom2, dad2;

  mom1 = (gen1 / n_bcstates) / (m+1);
  mom2 = (gen2 / n_bcstates) / (m+1);
  
  dad1 = (gen1 % n_bcstates) / (m+1);
  dad2 = (gen2 % n_bcstates) / (m+1);

  return((double)((mom1 != mom2) + (dad1 != dad2)) / 2.0);
}


/* R wrapper for est_map_stahl for intercross */
void R_est_map_f2i(int *n_ind, int *n_mar, int *geno, double *d, 
		   int *m, double *p, double *error_prob, 
		   double *loglik, int *maxit, double *tol, int *verbose)
{

  est_map_f2i(*n_ind, *n_mar, geno, d, *m, *p,
	      *error_prob, loglik, *maxit, *tol, *verbose);
}


/**********************************************************************
 * step_f2i
 * 
 * Calculate transition probabilities for Stahl model in an intercross,
 * on the basis of the results for a BC.
 **********************************************************************/
double step_f2i(int g1, int g2, int pos, double ***tm, int n_bcstates)
{
  return(tm[g1 % n_bcstates][g2 % n_bcstates][pos] +
	 tm[g1 / n_bcstates][g2 / n_bcstates][pos]);
}


/* end of hmm_f2i.c */
