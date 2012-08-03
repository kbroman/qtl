/**********************************************************************
 * 
 * hmm_bci.c
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
 * Contains: est_map_bci, 
 *           R_est_map_bci, 
 *           emit_bci, nrec_bci, step_bci,
 *           tm_bci, fms_bci, distinct_tm_bci
 *
 * These are functions for the HMM under the Stahl model
 * (with chiasmata coming from two mechanisms: one following a 
 * chi-square model and one following a no interference model).
 * m = interference parameter in the chi-square model (m=0 == NI)
 * p = proportion of chiasmata from the NI model (p=1 == NI)
 *
 * Code for is for a backcross.
 *
 * BACKCROSS:
 * Genotype codes:  0, ..., 2(m+1) - 1, with the first (m+1) 
 *                  corresponding to AA and the others to AB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
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
#include "hmm_bci.h"
#include "stahl_mf.h"

/**********************************************************************
 * 
 * est_map_bci
 *
 * This function re-estimates the genetic map for a chromosome
 * with the Stahl model, taking m and p known, for a backcross
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

void est_map_bci(int n_ind, int n_mar, int *geno, double *d, 
		 int m, double p, double error_prob, 
		 double *loglik, int maxit, double tol, int verbose)
{
  int i, j, j2, v, v2, it, flag=0, **Geno, n_states;
  double s, **alpha, **beta, **gamma, *cur_d, *rf;
  double ***tm, *temp;
  double curloglik;
  double initprob;
  int ndigits;
  double maxdif, tempdif;
  char pattern[100], text[200];
  
  n_states = 2*(m+1);  
  initprob = -log((double)n_states);

  /* allocate space for beta and reorganize geno */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  allocate_alpha(n_mar, n_states, &alpha);
  allocate_alpha(n_mar, n_states, &beta);
  allocate_dmatrix(n_states, n_states, &gamma);
  allocate_double(n_mar-1, &cur_d);
  allocate_double(n_mar-1, &rf);

  /* allocate space for the transition matrices */
  /* size n_states x n_states x (n_mar-1) */
  /* tm[state1][state2][interval] */
  tm = (double ***)R_alloc(n_states, sizeof(double **));
  tm[0] = (double **)R_alloc(n_states * n_states, sizeof(double *));
  for(i=1; i<n_states; i++) tm[i] = tm[i-1] + n_states;
  tm[0][0] = (double *)R_alloc(n_states * n_states * (n_mar - 1), 
			       sizeof(double));
  temp = tm[0][0];
  for(i=0; i < n_states; i++) {
    for(j=0; j < n_states; j++) {
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

    /* calculate the transition matrices */
    step_bci(n_mar, n_states, tm, cur_d, m, p, maxit, tol);

    for(i=0; i<n_ind; i++) { /* i = individual */

      R_CheckUserInterrupt(); /* check for ^C */

      /* initialize alpha and beta */
      for(v=0; v<n_states; v++) {
	alpha[v][0] = initprob + emit_bci(Geno[0][i], v, error_prob, m);
	beta[v][n_mar-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	
	for(v=0; v<n_states; v++) {
	  alpha[v][j] = alpha[0][j-1] + tm[0][v][j-1];
	  
	  beta[v][j2] = beta[0][j2+1] + tm[v][0][j2] +
	    emit_bci(Geno[j2+1][i], 0, error_prob, m);
	  
	  for(v2=1; v2<n_states; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 tm[v2][v][j-1]);
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				 tm[v][v2][j2] +
				 emit_bci(Geno[j2+1][i], v2, error_prob, m));
	  }
	  
	  alpha[v][j] += emit_bci(Geno[j][i], v, error_prob, m);
		 
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_states; v++) {
	  for(v2=0; v2<n_states; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + 
	      emit_bci(Geno[j+1][i], v2, error_prob, m) +
	      tm[v][v2][j];

	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}

	for(v=0; v<n_states; v++) {
	  for(v2=0; v2<n_states; v2++) {
	    rf[j] += nrec_bci(v, v2, m) * exp(gamma[v][v2] - s);
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
  step_bci(n_mar, n_states, tm, d, m, p, maxit, tol);

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */

    /* initialize alpha */
    for(v=0; v<n_states; v++) 
      alpha[v][0] = initprob + emit_bci(Geno[0][i], v, error_prob, m);

    /* forward equations */
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_states; v++) {
	alpha[v][j] = alpha[0][j-1] + tm[0][v][j-1];
	for(v2=1; v2<n_states; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       tm[v2][v][j-1]);
	alpha[v][j] += emit_bci(Geno[j][i], v, error_prob, m);
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
 * emit_bci: log Pr(obs_gen | true_gen)
 **********************************************************************/
double emit_bci(int obs_gen, int true_gen, double error_prob,
		int m)
{
  true_gen = true_gen / (m+1) + 1;

  switch(obs_gen) {
  case 0: return(0.0);
  case 1: case 2:
    if(obs_gen == true_gen) return(log(1.0-error_prob));
    else return(log(error_prob));
  }
  return(0.0); /* shouldn't get here */
}

/**********************************************************************
 * nrec_bci: proportion of recombinantion events
 **********************************************************************/
double nrec_bci(int gen1, int gen2, int m)
{
  gen1 /= (m+1); /* =0 if < m+1;   =1 otherwise */
  gen2 /= (m+1);

  if(gen1==gen2) return(0.0);
  else return(1.0);
}


/* R wrapper for est_map_stahl for backcross */
void R_est_map_bci(int *n_ind, int *n_mar, int *geno, double *d, 
		   int *m, double *p, double *error_prob, 
		   double *loglik, int *maxit, double *tol, int *verbose)
{
  est_map_bci(*n_ind, *n_mar, geno, d, *m, *p,
	      *error_prob, loglik, *maxit, *tol, *verbose);
}


/**********************************************************************
 * step_bci
 * 
 * Calculate transition probabilities (for all intervals) for
 * the Stahl model
 **********************************************************************/
void step_bci(int n_mar, int n_states, double ***tm, double *d, 
	      int m, double p, int maxit, double tol)
{
  int i, v1, v2;
  double *the_distinct_tm;
  double *fms_bci_result;
  double lambda1, lambda2, rfp;

  allocate_double(2*m+1, &fms_bci_result);
  allocate_double(3*m+2, &the_distinct_tm);

  for(i=0; i<n_mar-1; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    lambda1 = d[i]*(1-p)*(double)(m+1)*2.0;
    lambda2 = d[i]*p*2.0;
    rfp = 0.5*(1.0 - exp(-lambda2));

    fms_bci(lambda1, fms_bci_result, m, tol, maxit);
    distinct_tm_bci(lambda1, the_distinct_tm, m, fms_bci_result);

    for(v1=0; v1<n_states; v1++) {
      for(v2=0; v2<n_states; v2++) {
	tm[v1][v2][i] = tm_bci(v1, v2, the_distinct_tm, m);
	if(p > 0) 
	  tm[v1][v2][i] = (1.0-rfp)*tm[v1][v2][i] + 
	    rfp*tm_bci(v1, (v2+m+1) % (2*m+2), the_distinct_tm, m);
	tm[v1][v2][i] = log(tm[v1][v2][i]);
      }
    }
  }
} 


/*****************************************************************************
 * tm_bci: this function calculates the required transition probability for the
 * backcross case
 ****************************************************************************/

double tm_bci(int i, int j, double *the_distinct_tm, int m)
{
  int s, tempi, tempj;
  
  if ((i<=m && j<=m) || (i>m && j>m)) {
    s=j-i;
    if (s>=0) {
      return(the_distinct_tm[s]);
    }
    else {
      return(the_distinct_tm[abs(s)+2*m+1]);
    }
  }
  else if (i<=m && j>m) {
    if (j>(i+m)) {
      s=j-i;
      return(the_distinct_tm[s]);
    }
    else /* j <=i+m */ {
      s=j-i-(m+1);
      return(the_distinct_tm[abs(s)+2*m+1]);
    }
  }
  else /* i>m && j<=m */ {
    tempi=i-(m+1);
    tempj=j+(m+1);
    if (tempj>(tempi+m)) {
      s=tempj-tempi;
      return(the_distinct_tm[s]);
    }
    else /* tempj <=tempi+m */ {
      s=tempj-tempi-(m+1);
      return(the_distinct_tm[abs(s)+2*m+1]);
    }
  }
}

/*****************************************************************************
 * fms_bci: this function calculates the sum to infinity part of the
 * transition probabilities for a given lambda_t
 *
 * f should have length 2m+1
 ****************************************************************************/
void fms_bci(double lambda, double *f, int m, double tol, int maxit)
{
  int i,k;
  double diff;

  for (i=0; i<2*m+1; i++) {
    k=1;
    f[i]=0;
    if (i <= m) {
      f[i] = dpois((double)(k*(m+1)+i), lambda, 0);

      for(k=2; k<maxit; k++) {
        diff = dpois((double)(k*(m+1)+i), lambda, 0);
	f[i] += diff;

	if(diff < tol) break;
      }
    }
    else /* i<m */ {
      f[i] += dpois((double)(k*(m+1)+(m-i)), lambda, 0);

      for(k=2; k<maxit; k++) {
	diff = dpois((double)(k*(m+1)+(m-i)), lambda, 0);
	f[i] += diff;

	if(diff < tol) break;
      }  
    }  
    f[i] *= 0.5;
  }  
}


/*****************************************************************************
 * distinct_tm_bci: this function calculates the 3m+2 distinct transition 
 * probabilities for a given lambda_t
 ****************************************************************************/

void distinct_tm_bci(double lambda, double *the_distinct_tm, int m, 
		     double *fms_bci_result)
{
  int i;
  
  for (i=0;i<3*m+2;i++) {
    if (i<=m) 
      the_distinct_tm[i] = fms_bci_result[i]+ dpois((double)i, lambda, 0);

    else 
      the_distinct_tm[i] = fms_bci_result[i-(m+1)];

  }

}

/* end of hmm_bci.c */
