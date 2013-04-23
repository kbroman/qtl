/**********************************************************************
 * 
 * hmm_main.c
 *
 * copyright (c) 2001-2010, Karl W Broman
 *
 * last modified Aug, 2010
 * first written Feb, 2001
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
 * These functions are for the main HMM engine
 *
 * Contains: calc_genoprob, calc_genoprob_special, sim_geno, est_map, argmax_geno
 *           calc_errorlod, est_rf, calc_pairprob, calc_pairprob_condindep,
 *           R_calc_pairprob_condindep, marker_loglik
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "hmm_main.h"
#include "util.h"

/**********************************************************************
 * 
 * calc_genoprob
 *
 * This function uses the hidden Markov model technology to calculate 
 * the genotype probabilities at each of marker and (optionally) at 
 * points between markers, conditional on all marker data for a 
 * chromosome.  This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * genoprob     Genotype probabilities (the output); a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              genotype
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_genoprob(int n_ind, int n_pos, int n_gen, int *geno, 
		   double *rf, double *rf2, 
		   double error_prob, double *genoprob, 
		   double initf(int, int *), 
		   double emitf(int, int, double, int *),
		   double stepf(int, int, double, double, int *)) 
{
  int i, j, j2, v, v2;
  double s, **alpha, **beta;
  int **Geno;
  double ***Genoprob;
  int cross_scheme[2];

  /* cross scheme hidden in genoprob argument; used by hmm_bcsft */
  cross_scheme[0] = genoprob[0];
  cross_scheme[1] = genoprob[1];
  genoprob[0] = 0.0;
  genoprob[1] = 0.0;

  /* allocate space for alpha and beta and 
     reorganize geno and genoprob */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_genoprob(n_ind, n_pos, n_gen, genoprob, &Genoprob);
  allocate_alpha(n_pos, n_gen, &alpha);
  allocate_alpha(n_pos, n_gen, &beta);

  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */

    /* initialize alpha and beta */
    for(v=0; v<n_gen; v++) {
      alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
      beta[v][n_pos-1] = 0.0;
    }

    /* forward-backward equations */
    for(j=1,j2=n_pos-2; j<n_pos; j++, j2--) {
      
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, rf[j-1], rf2[j-1], cross_scheme);
	beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,rf[j2], rf2[j2], cross_scheme) + 
	  emitf(Geno[j2+1][i],1,error_prob, cross_scheme);

	for(v2=1; v2<n_gen; v2++) {
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,rf[j-1],rf2[j-1], cross_scheme));
	  beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
			       stepf(v+1,v2+1,rf[j2],rf2[j2], cross_scheme) +
			       emitf(Geno[j2+1][i],v2+1,error_prob, cross_scheme));
	}

	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob, cross_scheme);
      }
    }

    /* calculate genotype probabilities */
    for(j=0; j<n_pos; j++) {
      s = Genoprob[0][j][i] = alpha[0][j] + beta[0][j];
      for(v=1; v<n_gen; v++) {
	Genoprob[v][j][i] = alpha[v][j] + beta[v][j];
	s = addlog(s, Genoprob[v][j][i]);
      }
      for(v=0; v<n_gen; v++) 
	Genoprob[v][j][i] = exp(Genoprob[v][j][i] - s);
    } 

    /* the following is the old version */
    /*    for(j=0; j<n_pos; j++) {
      s = 0.0;
      for(v=0; v<n_gen; v++) 
	s += (Genoprob[v][j][i] = exp(alpha[v][j] + beta[v][j]));

      for(v=0; v<n_gen; v++) 
	Genoprob[v][j][i] /= s;
	} */


  } /* loop over individuals */
  

}



/**********************************************************************
 * 
 * calc_genoprob_special
 *
 * This is a special version of calc_genoprob, rerun specially for
 * each individual at each marker, allowing that genotype to 
 * be in error but assuming all others are without error
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_genoprob_special(int n_ind, int n_pos, int n_gen, int *geno, 
			   double *rf, double *rf2, 
			   double error_prob, double *genoprob, 
			   double initf(int, int *), 
			   double emitf(int, int, double, int *),
			   double stepf(int, int, double, double, int *)) 
{
  int i, j, j2, v, v2, curpos;
  double s, **alpha, **beta;
  int **Geno;
  double ***Genoprob;
  int cross_scheme[2];

  /* cross scheme hidden in genoprob argument; used by hmm_bcsft */
  cross_scheme[0] = genoprob[0];
  cross_scheme[1] = genoprob[1];
  genoprob[0] = 0.0;
  genoprob[1] = 0.0;
  
  /* allocate space for alpha and beta and 
     reorganize geno and genoprob */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_genoprob(n_ind, n_pos, n_gen, genoprob, &Genoprob);
  allocate_alpha(n_pos, n_gen, &alpha);
  allocate_alpha(n_pos, n_gen, &beta);

  for(i=0; i<n_ind; i++) { /* i = individual */

    for(curpos=0; curpos < n_pos; curpos++) {

      if(!Geno[curpos][i]) continue;

      R_CheckUserInterrupt(); /* check for ^C */

      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	if(curpos==0) 
	  alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
	else
	  alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, TOL, cross_scheme);
	beta[v][n_pos-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_pos-2; j<n_pos; j++, j2--) {
      
	for(v=0; v<n_gen; v++) {
	  alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, rf[j-1], rf2[j-1], cross_scheme);
	
	  if(curpos==j2+1)
	    beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,rf[j2], rf2[j2], cross_scheme) + 
	      emitf(Geno[j2+1][i],1,error_prob, cross_scheme);
	  else 
	    beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,rf[j2], rf2[j2], cross_scheme) + 
	      emitf(Geno[j2+1][i],1,TOL, cross_scheme);

	  for(v2=1; v2<n_gen; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 stepf(v2+1,v+1,rf[j-1],rf2[j-1], cross_scheme));
	    if(curpos==j2+1)
	      beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				   stepf(v+1,v2+1,rf[j2],rf2[j2], cross_scheme) +
				   emitf(Geno[j2+1][i],v2+1,error_prob, cross_scheme));
	    else
	      beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				   stepf(v+1,v2+1,rf[j2],rf2[j2], cross_scheme) +
				   emitf(Geno[j2+1][i],v2+1,TOL, cross_scheme));

	  }

	  if(curpos==j)
	    alpha[v][j] += emitf(Geno[j][i],v+1,error_prob, cross_scheme);
	  else
	    alpha[v][j] += emitf(Geno[j][i],v+1,TOL, cross_scheme);
	}
      }

      /* calculate genotype probabilities */
      s = Genoprob[0][curpos][i] = alpha[0][curpos] + beta[0][curpos];
      for(v=1; v<n_gen; v++) {
	Genoprob[v][curpos][i] = alpha[v][curpos] + beta[v][curpos];
	s = addlog(s, Genoprob[v][curpos][i]);
      }
      for(v=0; v<n_gen; v++) 
	Genoprob[v][curpos][i] = exp(Genoprob[v][curpos][i] - s);

    } /* end loop over current position */

  } /* loop over individuals */
}



/**********************************************************************
 * 
 * sim_geno
 *
 * This function simulates from the joint distribution Pr(g | O)
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              simulate genotypes)
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of simulation replicates
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps
 *
 * error_prob   Genotyping error probability
 *
 * draws        Simulated genotypes (the output), a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              draws
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void sim_geno(int n_ind, int n_pos, int n_gen, int n_draws,
	      int *geno, double *rf, double *rf2, 
	      double error_prob, int *draws,
	      double initf(int, int *), 
	      double emitf(int, int, double, int *),
	      double stepf(int, int, double, double, int *)) 
{
  int i, k, j, v, v2;
  double s, **beta, *probs;
  int **Geno, ***Draws, curstate;
  int cross_scheme[2];

  /* cross scheme hidden in draws argument; used by hmm_bcsft */
  cross_scheme[0] = draws[0];
  cross_scheme[1] = draws[1];
  draws[0] = 0;
  draws[1] = 0;
  
  /* allocate space for beta and 
     reorganize geno and draws */
  /* Geno indexed as Geno[pos][ind] */
  /* Draws indexed as Draws[rep][pos][ind] */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_draws(n_ind, n_pos, n_draws, draws, &Draws);
  allocate_alpha(n_pos, n_gen, &beta);
  allocate_double(n_gen, &probs);

  /* Read R's random seed */
  GetRNGstate();

  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */

    /* do backward equations */
    /* initialize beta */
    for(v=0; v<n_gen; v++) beta[v][n_pos-1] = 0.0;

    /* backward equations */
    for(j=n_pos-2; j>=0; j--) {
      
      for(v=0; v<n_gen; v++) {
	beta[v][j] = beta[0][j+1] + stepf(v+1,1,rf[j], rf2[j], cross_scheme) + 
	  emitf(Geno[j+1][i],1,error_prob, cross_scheme);

	for(v2=1; v2<n_gen; v2++) 
	  beta[v][j] = addlog(beta[v][j], beta[v2][j+1] + 
			      stepf(v+1,v2+1,rf[j],rf2[j], cross_scheme) +
			      emitf(Geno[j+1][i],v2+1,error_prob, cross_scheme));
      }
    }

    for(k=0; k<n_draws; k++) { /* k = simulation replicate */

      /* first draw */
      /* calculate probs */
      s = (probs[0] = initf(1, cross_scheme)+emitf(Geno[0][i],1,error_prob, cross_scheme)+beta[0][0]);
      for(v=1; v<n_gen; v++) {
	probs[v] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme) +
	  beta[v][0];
	s = addlog(s, probs[v]);
      }
      for(v=0; v<n_gen; v++) probs[v] = exp(probs[v] - s);

      /* make draw: returns a value from {1, 2, ..., n_gen} */
      curstate = Draws[k][0][i] = sample_int(n_gen, probs);
      
      /* move along chromosome */
      for(j=1; j<n_pos; j++) {
	/* calculate probs */
	for(v=0; v<n_gen; v++) 
	  probs[v] = exp(stepf(curstate,v+1,rf[j-1],rf2[j-1], cross_scheme) +
			 emitf(Geno[j][i],v+1,error_prob, cross_scheme) +
			 beta[v][j] - beta[curstate-1][j-1]);
	/* make draw */
	curstate = Draws[k][j][i] = sample_int(n_gen, probs);
      }

    } /* loop over replicates */

  } /* loop over individuals */
  
  /* write R's random seed */
  PutRNGstate();

}



/**********************************************************************
 * 
 * est_map
 *
 * This function re-estimates the genetic map for a chromosome
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * n_gen        Number of different genotypes
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          Second set of recombination fractions (may not be needed)
 *
 * error_prob   Genotyping error probability
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 * nrecf1       Function returning number of recombinations associated
 *              with (g_1, g_2)
 *
 * nrecf2       Another such function, used only in the case of a sex-
 *              specific map
 *
 * loglik       Loglik at final estimates of recombination fractions
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 * sexsp        Indicates whether sex-specific maps should be estimated
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void est_map(int n_ind, int n_mar, int n_gen, int *geno, double *rf, 
	     double *rf2, double error_prob, double initf(int, int *), 
	     double emitf(int, int, double, int *),
	     double stepf(int, int, double, double, int *), 
	     double nrecf1(int, int, double, int*), double nrecf2(int, int, double, int*), 
	     double *loglik, int maxit, double tol, int sexsp, 
	     int verbose)
{
  int i, j, j2, v, v2, it, flag=0, **Geno, ndigits;
  double s, **alpha, **beta, **gamma, *cur_rf, *cur_rf2;
  double curloglik, maxdif, temp;
  char pattern[100], text[200];
  int cross_scheme[2];

  /* cross scheme hidden in loglik argument; used by hmm_bcsft */
  cross_scheme[0] = (int) ftrunc(*loglik / 1000.0);
  cross_scheme[1] = ((int) *loglik) - 1000 * cross_scheme[0];
  *loglik = 0.0;
  
  /* allocate space for beta and reorganize geno */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  allocate_alpha(n_mar, n_gen, &alpha);
  allocate_alpha(n_mar, n_gen, &beta);
  allocate_dmatrix(n_gen, n_gen, &gamma);
  allocate_double(n_mar-1, &cur_rf);
  allocate_double(n_mar-1, &cur_rf2);

  /* digits in verbose output */
  if(verbose) {
    ndigits = (int)ceil(-log10(tol));
    if(ndigits > 16) ndigits=16;
    sprintf(pattern, "%s%d.%df", "%", ndigits+3, ndigits+1);
  }

  /* begin EM algorithm */
  for(it=0; it<maxit; it++) {

    for(j=0; j<n_mar-1; j++) {
      cur_rf[j] = cur_rf2[j] = rf[j];
      rf[j] = 0.0;
      if(sexsp) {
	cur_rf2[j] = rf2[j];
	rf2[j] = 0.0;
      }
    }

    for(i=0; i<n_ind; i++) { /* i = individual */

      R_CheckUserInterrupt(); /* check for ^C */

      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
	beta[v][n_mar-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	
	for(v=0; v<n_gen; v++) {
	  alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, cur_rf[j-1], cur_rf2[j-1], cross_scheme);
	  beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,cur_rf[j2], cur_rf2[j2], cross_scheme) + 
	    emitf(Geno[j2+1][i],1,error_prob, cross_scheme);
	  
	  for(v2=1; v2<n_gen; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 stepf(v2+1,v+1,cur_rf[j-1],cur_rf2[j-1], cross_scheme));
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				 stepf(v+1,v2+1,cur_rf[j2],cur_rf2[j2], cross_scheme) +
				 emitf(Geno[j2+1][i],v2+1,error_prob, cross_scheme));
	  }
	  
	  alpha[v][j] += emitf(Geno[j][i],v+1,error_prob, cross_scheme);
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + 
	      emitf(Geno[j+1][i], v2+1, error_prob, cross_scheme) +
	      stepf(v+1, v2+1, cur_rf[j], cur_rf2[j], cross_scheme);

	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}

	for(v=0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    rf[j] += nrecf1(v+1,v2+1, cur_rf[j], cross_scheme) * exp(gamma[v][v2] - s);
	    if(sexsp) rf2[j] += nrecf2(v+1,v2+1, cur_rf[j], cross_scheme) * exp(gamma[v][v2] - s);
	  }
	}
      }

    } /* loop over individuals */

    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/1000.0) rf[j] = tol/1000.0;
      else if(rf[j] > 0.5-tol/1000.0) rf[j] = 0.5-tol/1000.0;
      
      if(sexsp) {
	rf2[j] /= (double)n_ind;
	if(rf2[j] < tol/1000.0) rf2[j] = tol/1000.0;
	else if(rf2[j] > 0.5-tol/1000.0) rf2[j] = 0.5-tol/1000.0;
      }
      else rf2[j] = rf[j];
    }

    if(verbose>1) {
      /* print estimates as we go along*/
      Rprintf("   %4d ", it+1);
      maxdif=0.0;
      for(j=0; j<n_mar-1; j++) {
	temp = fabs(rf[j] - cur_rf[j])/(cur_rf[j]+tol*100.0);
	if(maxdif < temp) maxdif = temp;
	if(sexsp) {
	  temp = fabs(rf2[j] - cur_rf2[j])/(cur_rf2[j]+tol*100.0);
	  if(maxdif < temp) maxdif = temp;
	}
	/* bsy add */
	if(verbose > 2)
	  Rprintf("%d %f %f\n", j+1, cur_rf[j], rf[j]); 
	/* bsy add */
      }
      sprintf(text, "%s%s\n", "  max rel've change = ", pattern);
      Rprintf(text, maxdif);
    }

    /* check convergence */
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(rf[j] - cur_rf[j]) > tol*(cur_rf[j]+tol*100.0) || 
	 (sexsp && fabs(rf2[j] - cur_rf2[j]) > tol*(cur_rf2[j]+tol*100.0))) {
	flag = 1; 
	break;
      }
    }

    if(!flag) break;

  } /* end EM algorithm */
  
  if(flag) warning("Didn't converge!\n");

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */
    /* initialize alpha */
    for(v=0; v<n_gen; v++) {
      alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
    }
    /* forward equations */
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + 
	  stepf(1, v+1, rf[j-1], rf2[j-1], cross_scheme);
	
	for(v2=1; v2<n_gen; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,rf[j-1],rf2[j-1], cross_scheme));

	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob, cross_scheme);
      }
    }

    curloglik = alpha[0][n_mar-1];
    for(v=1; v<n_gen; v++)
      curloglik = addlog(curloglik, alpha[v][n_mar-1]);

    *loglik += curloglik;
  }

  if(verbose) {
    if(verbose < 2) {
      /* print final estimates */
      Rprintf("  no. iterations = %d\n", it+1);
      maxdif=0.0;
      for(j=0; j<n_mar-1; j++) {
	temp = fabs(rf[j] - cur_rf[j])/(cur_rf[j]+tol*100.0);
	if(maxdif < temp) maxdif = temp;
	if(sexsp) {
	  temp = fabs(rf2[j] - cur_rf2[j])/(cur_rf2[j]+tol*100.0);
	  if(maxdif < temp) maxdif = temp;
	}
      }
      sprintf(text, "%s%s\n", "  max rel've change at last step = ", pattern);
      Rprintf(text, maxdif);
    }
    
    Rprintf("  loglik: %10.4lf\n\n", *loglik);
  }

}



/**********************************************************************
 * 
 * argmax_geno
 *
 * This function uses the Viterbi algorithm to calculate the most 
 * likely sequence of underlying genotypes, given the observed marker
 * data for a chromosome.
 * This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              find most likely genotypes)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * argmax       Matrix of most likely genotypes (the output); a single 
 *              vector stored by columns (ind moves fastest, then pos)
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void argmax_geno(int n_ind, int n_pos, int n_gen, int *geno, 
		 double *rf, double *rf2, 
		 double error_prob, int *argmax, 
		 double initf(int, int *), 
		 double emitf(int, int, double, int *),
		 double stepf(int, int, double, double, int *)) 
{
  int i, j, v, v2;
  double s, t, *gamma, *tempgamma, *tempgamma2;
  int **Geno, **Argmax, **traceback;
  int cross_scheme[2];

  /* cross scheme hidden in argmax argument; used by hmm_bcsft */
  cross_scheme[0] = argmax[0];
  cross_scheme[1] = argmax[1];
  argmax[0] = geno[0];
  argmax[1] = geno[1];

  /* Read R's random seed */
  /* in the case of multiple "most likely" genotype sequences, 
     we pick from them at random */
  GetRNGstate();

  /* allocate space and 
     reorganize geno and argmax */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_geno(n_ind, n_pos, argmax, &Argmax);
  allocate_imatrix(n_pos, n_gen, &traceback);
  allocate_double(n_gen, &gamma);
  allocate_double(n_gen, &tempgamma);
  allocate_double(n_gen, &tempgamma2);

  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */

    /* begin viterbi algorithm */
    if(n_pos > 1) { /* multiple markers */
      for(v=0; v<n_gen; v++) 
	gamma[v] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
    
      for(j=0; j<n_pos-1; j++) {
	for(v=0; v<n_gen; v++) {
	  tempgamma[v] = s = gamma[0] + stepf(1, v+1, rf[j], rf2[j], cross_scheme);
	  traceback[j][v] = 0;
	  
	  for(v2=1; v2<n_gen; v2++) {
	    t = gamma[v2] + stepf(v2+1, v+1, rf[j], rf2[j], cross_scheme);
	    if(t > s || (fabs(t-s) < TOL && unif_rand() < 0.5)) {
	      tempgamma[v] = s = t;
	      traceback[j][v] = v2;
	    }
	  }
	  tempgamma2[v] = tempgamma[v] + emitf(Geno[j+1][i], v+1, error_prob, cross_scheme);
	}
	for(v=0; v<n_gen; v++) gamma[v] = tempgamma2[v];
      }
    
      /* finish off viterbi and then traceback to get most 
	 likely sequence of genotypes */
      Argmax[n_pos-1][i] = 0;
      s = gamma[0];
      for(v=1; v<n_gen; v++) {
	if(gamma[v] > s || (fabs(gamma[v]-s) < TOL && 
			    unif_rand() < 0.5)) {
	  s = gamma[v];
	  Argmax[n_pos-1][i] = v;
	}
      }
      for(j=n_pos-2; j >= 0; j--) 
	Argmax[j][i] = traceback[j][Argmax[j+1][i]];
    }
    else {  /* for exactly one marker */
      s = initf(1, cross_scheme) + emitf(Geno[0][i], 1, error_prob, cross_scheme);
      Argmax[0][i] = 0;
      for(v=1; v<n_gen; v++) {
	t = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
	if(t > s || (fabs(t-s) < TOL && unif_rand() < 0.5)) {
	  s = t;
	  Argmax[0][i] = v;
	}
      }
    }
    
    /* code genotypes as 1, 2, ... */
    for(j=0; j<n_pos; j++) Argmax[j][i]++;
    
  } /* loop over individuals */
  
  
  /* write R's random seed */
  PutRNGstate();
}



/**********************************************************************
 * 
 * calc_errorlod
 *
 * Uses the results of calc_genoprob to calculate a LOD score for 
 * each genotype, indicating whether it is likely to be in error.
 *
 * n_ind, n_mar, n_gen, geno        These are all as in the above funcs
 * error_prob, genoprob 
 *
 * errlod          The output, as a single vector stored by columns, 
 *                 of size n_ind x n_mar
 * 
 * errorlod        Function taking observed genotype, genotype probs,
 *                 and error probability, and returning the error LOD
 *
 **********************************************************************/

void calc_errorlod(int n_ind, int n_mar, int n_gen, int *geno, 
		   double error_prob, double *genoprob, double *errlod, 
		   double errorlod(int, double *, double))
{
  int i, j, k, **Geno;
  double *p, ***Genoprob, **Errlod;

  /* reorganize geno, genoprob and errlod */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  reorg_genoprob(n_ind, n_mar, n_gen, genoprob, &Genoprob);
  reorg_errlod(n_ind, n_mar, errlod, &Errlod);
  allocate_double(n_gen, &p);

  for(i=0; i<n_ind; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    for(j=0; j<n_mar; j++) {
      for(k=0; k<n_gen; k++) p[k] = Genoprob[k][j][i];
      Errlod[j][i] = errorlod(Geno[j][i], p, error_prob);
    }
  }

}



/**********************************************************************
 * 
 * est_rf
 *
 * Estimate sex-averaged recombination fractions for all pairs of loci
 *
 * This is for f2 and 4way crosses; backcrosses don't need the EM 
 * algorithm, since there is no partially missing data.
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers
 *
 * geno         Matrix of genotype data (n_ind x n_mar), stored as a 
 *              single vector (by columns)
 * 
 * rf           The output: matrix of doubles (n_mar x n_mar), stored
 *              as a single vector (by columns).  The diagonal will 
 *              contain the number of meioses, the lower triangle will
 *              contain the est'd rec fracs, and the upper triangle
 *              will contain the LOD scores (testing rf=0.5)
 *
 * erec         Function returning the expected number of recombination
 *              events given observed marker genotypes
 * 
 * logprec      Function returning the log probability of a pair of 
 *              observed genotypes, given the recombination fraction
 *              (for calculating the LOD score)
 * 
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence of the EM 
 *
 * meioses_per  No. meioses per individual
 *
 **********************************************************************/

void est_rf(int n_ind, int n_mar, int *geno, double *rf, 
	    double erec(int, int, double, int *), 
	    double logprec(int, int, double, int *), 
	    int maxit, double tol, int meioses_per)
{
  int i, j1, j2, s, **Geno, n_mei=0, flag=0;
  double **Rf, next_rf=0.0, cur_rf=0.0;
  int cross_scheme[2];

  /* cross scheme hidden in rf argument; used by hmm_bcsft */
  cross_scheme[0] = rf[0];
  cross_scheme[1] = rf[1];
  rf[0] = 0.0;
  rf[1] = 0.0;

  /* reorganize geno and rf */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  reorg_errlod(n_mar, n_mar, rf, &Rf);

  for(j1=0; j1<n_mar; j1++) {

    /* count number of meioses */
    for(i=0, n_mei=0; i<n_ind; i++) 
      if(Geno[j1][i] != 0) n_mei += meioses_per;
    Rf[j1][j1] = (double) n_mei;
    
    R_CheckUserInterrupt(); /* check for ^C */

    for(j2=j1+1; j2<n_mar; j2++) {
      
      /* count meioses */
      n_mei = flag = 0;
      for(i=0; i<n_ind; i++) {
	if(Geno[j1][i] != 0 && Geno[j2][i] != 0) {
	  n_mei += meioses_per;
	  /* check if informatve */
	  if(fabs(logprec(Geno[j1][i], Geno[j2][i], 0.5, cross_scheme) -
		  logprec(Geno[j1][i], Geno[j2][i], TOL, cross_scheme)) > TOL) flag = 1;
	}
      }
      if(n_mei != 0 && flag == 1) {
	flag = 0;
	/* begin EM algorithm; start with cur_rf = 0.01 */
	for(s=0, cur_rf=0.01; s < maxit; s++) {
	  next_rf = 0.0; 
	  for(i=0; i<n_ind; i++) {
	    if(Geno[j1][i] != 0 && Geno[j2][i] != 0)
	      next_rf += erec(Geno[j1][i], Geno[j2][i], cur_rf, cross_scheme);
	  }

	  next_rf /= (double) n_mei;
	  if(fabs(next_rf - cur_rf) < tol*(cur_rf+tol*100.0)) { 
	    flag = 1;
	    break;
	  }
	  cur_rf = next_rf;
	}
	if(!flag) warning("Markers (%d,%d) didn't converge\n", j1+1, j2+1);

	/* calculate LOD score */
	Rf[j1][j2] = next_rf;
	Rf[j2][j1] = 0.0;
	for(i=0; i<n_ind; i++) {
	  if(Geno[j1][i] != 0 && Geno[j2][i] != 0) {
	    Rf[j2][j1] += logprec(Geno[j1][i],Geno[j2][i], next_rf, cross_scheme);
	    Rf[j2][j1] -= logprec(Geno[j1][i],Geno[j2][i], 0.5, cross_scheme);
	  }
	}
	Rf[j2][j1] /= log(10.0);
	
      }
      else { /* no informative meioses */
	Rf[j1][j2] = NA_REAL;
	Rf[j2][j1] = 0.0;
      }
	  
    } /* end loops over markers */
  }
}

/**********************************************************************
 * 
 * calc_pairprob
 *
 * This function uses the hidden Markov model technology to calculate 
 * the joint genotype probabilities for all pairs of putative QTLs.
 * This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * genoprob     Genotype probabilities (the output); a single vector, 
 *              of length n_ind x n_pos x n_gen, stored by columns 
 *              (ind moves fastest, then mar, then genotype
 *
 * pairprob     Joint genotype probabilities for pairs of positions.
 *              A single vector of length n_ind x n_pos x (n_pos-1)/2 x
 *              n_gen^2.  We only calculate probabilities for 
 *              pairs (i,j) with i < j.
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_pairprob(int n_ind, int n_pos, int n_gen, int *geno, 
		   double *rf, double *rf2, 
		   double error_prob, double *genoprob, 
		   double *pairprob, 
		   double initf(int, int *), 
		   double emitf(int, int, double, int *),
		   double stepf(int, int, double, double, int *)) 
{
  int i, j, j2, v, v2, v3;
  double s=0.0, **alpha, **beta;
  int **Geno;
  double ***Genoprob, *****Pairprob;
  int cross_scheme[2];

  /* cross scheme hidden in genoprob argument; used by hmm_bcsft */
  cross_scheme[0] = genoprob[0];
  cross_scheme[1] = genoprob[1];
  genoprob[0] = 0.0;
  genoprob[1] = 0.0;
  
  /* n_pos must be at least 2, or there are no pairs! */
  if(n_pos < 2) error("n_pos must be > 1 in calc_pairprob");

  /* allocate space for alpha and beta and 
     reorganize geno, genoprob, and pairprob */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_genoprob(n_ind, n_pos, n_gen, genoprob, &Genoprob);
  reorg_pairprob(n_ind, n_pos, n_gen, pairprob, &Pairprob);
  allocate_alpha(n_pos, n_gen, &alpha);
  allocate_alpha(n_pos, n_gen, &beta);

  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */

    /* initialize alpha and beta */
    for(v=0; v<n_gen; v++) {
      alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
      beta[v][n_pos-1] = 0.0;
    }

    /* forward-backward equations */
    for(j=1,j2=n_pos-2; j<n_pos; j++, j2--) {
      
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, rf[j-1], rf2[j-1], cross_scheme);
	
	beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,rf[j2], rf2[j2], cross_scheme) + 
	  emitf(Geno[j2+1][i],1,error_prob, cross_scheme);

	for(v2=1; v2<n_gen; v2++) {
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,rf[j-1],rf2[j-1], cross_scheme));
	  beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
			       stepf(v+1,v2+1,rf[j2],rf2[j2], cross_scheme) +
			       emitf(Geno[j2+1][i],v2+1,error_prob, cross_scheme));
	}

	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob, cross_scheme);
      }
    }

    /* calculate genotype probabilities */
    for(j=0; j<n_pos; j++) {
      s = Genoprob[0][j][i] = alpha[0][j] + beta[0][j];
      for(v=1; v<n_gen; v++) {
	Genoprob[v][j][i] = alpha[v][j] + beta[v][j];
	s = addlog(s, Genoprob[v][j][i]);
      }
      for(v=0; v<n_gen; v++) 
	Genoprob[v][j][i] = exp(Genoprob[v][j][i] - s);
    }

    /* calculate Pr(G[j], G[j+1] | marker data) for i = 1...n_pos-1 */
    for(j=0; j<n_pos-1; j++) {
      for(v=0; v<n_gen; v++) {
	for(v2=0; v2<n_gen; v2++) {
	  Pairprob[v][v2][j][j+1][i] = alpha[v][j] + beta[v2][j+1] +
	    stepf(v+1,v2+1,rf[j],rf2[j], cross_scheme) + 
	    emitf(Geno[j+1][i],v2+1,error_prob, cross_scheme);
	  if(v==0 && v2==0) s=Pairprob[v][v2][j][j+1][i];
	  else s = addlog(s,Pairprob[v][v2][j][j+1][i]);
	}
      }
      /* scale to sum to 1 */
      for(v=0; v<n_gen; v++) 
	for(v2=0; v2<n_gen; v2++) 
	  Pairprob[v][v2][j][j+1][i] = 
	    exp(Pairprob[v][v2][j][j+1][i] - s);
    } 

    /* now calculate Pr(G[i], G[j] | marker data) for j > i+1 */
    for(j=0; j<n_pos-2; j++) {
      for(j2=j+2; j2<n_pos; j2++) {

	for(v=0; v<n_gen; v++) { /* genotype at pos'n j */
	  for(v2=0; v2<n_gen; v2++) { /* genotype at pos'n j2 */

	    Pairprob[v][v2][j][j2][i] = 0.0;

	    for(v3=0; v3<n_gen; v3++) { /* genotype at pos'n j2-1 */
	      s = Genoprob[v3][j2-1][i];
	      if(fabs(s) > TOL) /* avoid 0/0 */
		Pairprob[v][v2][j][j2][i] += Pairprob[v][v3][j][j2-1][i]*
		  Pairprob[v3][v2][j2-1][j2][i]/s;
	    }
		      
	  }
	} /* end loops over genotypes */

      }
    } /* end loops over pairs of positions */
	    
  } /* end loop over individuals */
}


/**********************************************************************
 * 
 * calc_pairprob_condindep
 *
 * This function calculates the joint genotype probabilities assuming
 * conditional independence of QTL genotypes given the marker data
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * genoprob     QTL genotype probabilities given the marker data
 *
 * pairprob     Joint genotype probabilities for pairs of positions.
 *              A single vector of length n_ind x n_pos x (n_pos-1)/2 x
 *              n_gen^2.  We only calculate probabilities for 
 *              pairs (i,j) with i < j.
 *
 **********************************************************************/

void calc_pairprob_condindep(int n_ind, int n_pos, int n_gen, 
			     double ***Genoprob, double *****Pairprob) 
{
  int i, j, j2, v, v2;
  
  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */

    /* calculate Pr(G[j], G[j+1] | marker data) for i = 1...n_pos-1 */
    for(j=0; j<n_pos-1; j++) 
      for(j2=j+1; j2<n_pos; j2++)
	for(v=0; v<n_gen; v++) 
	  for(v2=0; v2<n_gen; v2++) 
	    Pairprob[v][v2][j][j2][i] = Genoprob[v][j][i] * Genoprob[v2][j2][i];
  }

}

/* wrapper for calc_pairprob_condindep */
void R_calc_pairprob_condindep(int *n_ind, int *n_pos, int *n_gen, 
			       double *genoprob, double *pairprob)
{
  double ***Genoprob, *****Pairprob;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);

  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);

  calc_pairprob_condindep(*n_ind, *n_pos, *n_gen,
			  Genoprob, Pairprob);
}



/**********************************************************************
 * 
 * marker_loglik
 *
 * This function calculates the log likelihood for a fixed marker
 *
 * n_ind        Number of individuals
 *
 * n_gen        Number of different genotypes
 *
 * geno         Genotype data, as a single vector
 *
 * error_prob   Genotyping error probability
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * loglik       Loglik at return
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2 */

void marker_loglik(int n_ind, int n_gen, int *geno, 
		   double error_prob, double initf(int, int *), 
		   double emitf(int, int, double, int *),
		   double *loglik)
{
  int i, v;
  double temp;
  int cross_scheme[2];

  /* cross scheme hidden in loglik argument; used by hmm_bcsft */
  cross_scheme[0] = (int) ftrunc(*loglik / 1000.0);
  cross_scheme[1] = ((int) *loglik) - 1000 * cross_scheme[0];
  
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */

    R_CheckUserInterrupt(); /* check for ^C */
    
    temp = initf(1, cross_scheme) + emitf(geno[i], 1, error_prob, cross_scheme);
    for(v=1; v<n_gen; v++) 
      temp = addlog(temp, initf(v+1, cross_scheme) + emitf(geno[i], v+1, error_prob, cross_scheme));

    (*loglik) += temp;
  }
}


/* end of hmm_main.c */
