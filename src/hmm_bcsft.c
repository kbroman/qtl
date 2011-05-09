/**********************************************************************
 * 
 * hmm_bcsft.c
 * 
 * copyright (c) 2001-9, Karl W Broman
 * modified from hmm_f2.c by Brian S Yandell and Laura M Shannon (c) 2010
 *
 * modified Jun, 2010
 * last modified Apr, 2009
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
 * Contains: init_bcsft, emit_bcsft, step_bcsft, init_bcsftb, emit_bcsftb, step_bcsftb,
 *           calc_genoprob_bcsft, calc_genoprob_special_bcsft, sim_geno_bcsft, est_map_bcsft, 
 *           argmax_geno_bcsft, errorlod_bcsft, calc_errorlod_bcsft, nrec2_bcsft,
 *           logprec_bcsft, est_rf_bcsft, calc_pairprob_bcsft, marker_loglik_bcsft
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the F2 intercross.
 *
 * Genotype codes:  0=AA; 1=AB; 2=BB
 * Phenotype codes: 0=missing; 1=AA; 2=AB; 3=BB; 4=not BB; 5=not AA
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "hmm_bcsft.h"
#include "hmm_f2.h"
#include "hmm_bc.h"
#include "util.h"

/* ref: Jiang and Zeng (1997 Genetics) */

void prob_bcsft(double rf, int s, int t, double *transpr);
void count_bcsft(double rf, int s, int t, double *transct);
void expect_bcsft(double rf, int s, int t, double *transexp);

double assign_bcsft(int gen1, int gen2, double *transpr) 
{
  /* joint probability with known genotype, phase unknown */

  switch(gen1) {
  case 1: case 3:
    { /* AA and aa for gen1 */
      if(gen2 == gen1) {
	if(gen1 == 1) return(transpr[0]);            /* 1,1: A1 */
	return(transpr[5]);                          /* 3,3: A0 */
      }
      if(gen2 + gen1 == 4) return(transpr[2]);       /* 1,3: C */
      break;
    }
  case 2:
    { /* Aa and aA for gen1 */
      if(gen2 == gen1) return(transpr[3]);           /* 2,2: D or E */
    }
  }
  if((gen1 == 1) || (gen2 == 1)) return(transpr[1]); /* 1,2: B1 */
  return(transpr[6]);                                /* 2,3: B0 */
}
double assign_bcsftb(int gen1, int gen2, double *transpr)
{
  /* joint probability with known genotype and unknown */

  switch(gen1) {
  case 1: case 4: /* AA and aa for gen1 */
    {
      if(gen2 == gen1) {
	if(gen1 == 1) return(transpr[0]);            /* 1,1: A1 */
	return(transpr[5]);                          /* 4,4: A0 */
      }
      if(gen2 + gen1 == 5) return(transpr[2]);       /* 1,4: C */
      break;
    }
  case 2: case 3: /* Aa and aA for gen1 */
    {
      if(gen2 == gen1) return(transpr[3]);           /* 2,2: D */
      if(gen2 + gen1 == 5) return(transpr[4]);       /* 2,3: E */
    }
  }
  if((gen1 == 1) || (gen2 == 1)) return(transpr[1]); /* 1,2|3: B1 */
  return(transpr[6]);                                /* 2|3,4: B0 */
}
double assign_bcsftc(int obs1, int obs2, double *transval)
{
  /* joint probability of obs2 and obs1, allowing for partially informative genos */

  if((obs1 == 0) || (obs2 == 0)) return(0.0); /* shouldn't get here */

  int temp;

  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    temp = obs2;
    obs2 = obs1;
    obs1 = temp;
  }

  switch(obs1) {
  case 1: case 3: { /* AA and aa for obs1 */
    if(obs2 == obs1) {
      if(obs1 == 1) return(transval[0]);                                    /* 1,1: A1 */
      return(transval[5]);                                                  /* 3,3: A0 */
    }
    if(obs2 + obs1 == 4) return(transval[2]);                               /* 1,3: C */
    if(obs1 == 1) { /* B1 */
      if(obs1 + obs2 == 3) return(transval[1]);                             /* 1,2: B1 */
      if(obs1 + obs2 == 5) return(transval[0] + transval[1]);               /* 1,4: A1 or B1 */
      return(transval[2] + transval[1]);                                    /* 1,5: B1 or C */
    }
    { /* B0 */
      if(obs1 + obs2 == 7) return(transval[2] + transval[6]);               /* 3,4: A1 or B0 */
      return(transval[5] + transval[6]);                                    /* 3,5: A0 or B0 */
    }
  }
  case 2: /* Aa and aA for obs1 */
    {
      if(obs2 == obs1) return(transval[3]);                                 /* 2,2: D or E */
      if(obs1 + obs2 == 5) return(transval[6]);                             /* 2,3: B0 */
      if(obs1 + obs2 == 6) return(transval[1] + transval[3]);               /* 2,4: A1 or B0 */
      return(transval[6] + transval[3]);                                    /* 2,5: A0 or B0 */
    }
  case 4: /* AA or Aa for obs1 */
    {
      if(obs1 == obs2) return(transval[0] + 2 * transval[1] + transval[3]); /* 4,4: 1 or 2 */
      break;
    }
  case 5: /* Aa or aa for obs1 */
    {
      if(obs1 == obs2) return(transval[3] + 2 * transval[6] + transval[5]); /* 5,5: 2 or 3 */
    }
  }
  return(transval[1] + transval[2] + transval[3] + transval[6]);            /* 4,5: 1 or 2 x 2 or 3 */
}
double assign_bcsftd(int n_gen, int obs1, int obs2, double *transval)
{
  if(n_gen == 5) return(assign_bcsftc(obs1, obs2, transval));
  return(assign_bcsftb(obs1, obs2, transval));
}

double init_bcsft(int true_gen, int *cross_scheme)
{
  static double init1 = 0;
  static double init2 = 0;
  static double init3 = 0;
  static int s = -1;
  static int t = -1;
  
  if(s != cross_scheme[0] || t != cross_scheme[1] || init1 == 0) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    /* static variables used frequently */
    if(s == 0) {  /* Ft */
      init2 = (1 - t) * M_LN2;                  /* Aa: log(2 ^ (1-t)) */
      init1 = log1p(-exp(init2)) - M_LN2;       /* AA: log((1 - 2^(1-t)) / 2) */
      init3 = init1;                            /* aa: */ 
    }
    if(s > 0) {
      if(t == 0) { /* BCs */
	init2 = -s * M_LN2;                     /* Aa: log(2 ^ -s) */
	init1 = log1p(-exp(init2));             /* AA: log(1 - 2^-s) */
      }
      if(t > 0) {  /* BCsFt */
	double sm2,tm2;
	sm2 = -s * M_LN2;
	tm2 = -t * M_LN2;
	init2 = sm2 + tm2;                      /* Aa: log(2 ^ -(s+t)) */
	init3 = sm2 + log1p(-exp(tm2)) - M_LN2; /* aa: log(2^-s * (1 - 2^-t) / 2) */
	init1 = log1p(exp(init3) - exp(sm2));   /* AA: log((1 - 2^-s) + 2^-s * (1 - 2^-t)) */
      }
    }
  }
  switch(true_gen) {
  case 1: return(init1);
  case 2: return(init2);
  case 3: return(init3);
  }
  return(0.0); /* should not get here */
}

double emit_bcsft(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  if(cross_scheme[1] > 0) return(emit_f2(obs_gen, true_gen, error_prob,cross_scheme));
  return(emit_bc(obs_gen, true_gen, error_prob,cross_scheme));
}

double step_bcsft(int gen1, int gen2, double rf, double junk, int *cross_scheme) 
{
  static double transpr[10];
  static double oldrf = -1.0;
  static int s = -1;
  static int t = -1;
  
  if(s != cross_scheme[0] || t != cross_scheme[1] || fabs(rf - oldrf) > TOL) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    oldrf = rf;
    if(rf < TOL) rf = TOL;
    
    prob_bcsft(rf, s, t, transpr);

    /* collapse when phase is unknown */
    if(t > 0) { /* only if Ft in play */
      transpr[3] += transpr[4]; /* D or E */
    }

    /* put probabilities on log scale */
    int k;
    for(k=0; k<7; k++) {
      /*      if(transpr[k] > 0.0) */
      transpr[k] = log(transpr[k]);
    }
  }

  double out;
  /* Find joint probability pr(gen1,gen2). */
  out = assign_bcsft(gen1, gen2, transpr);

  /* Divide by marginal prob to get pr(gen2|gen1). */
  out -= transpr[6+gen1];

  return(out);
}

/****************************************************************************/

/* The following are the init, emit and step functions 
   when considering phase-known F2 genotypes 
   (i.e. the 4-state chain: AA, AB, BA, BB             */

double init_bcsftb(int true_gen, int *cross_scheme)
{
  static double init1 = 0;
  static double init2 = 0;
  static double init3 = 0;
  static double init4 = 0;
  static int s = -1;
  static int t = -1;
  
  /* static variables used frequently */
  if(s != cross_scheme[0] || t != cross_scheme[1] || init1 == 0) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    if(s == 0) {  /* Ft */
      init2 = - t * M_LN2;                           /* Aa: log(2 ^ -t) */
      init1 = log1p(-exp(init2 + M_LN2)) - M_LN2;    /* AA: log((1 - 2^(1-t)) / 2) */
      init3 = init2;                                 /* aA: */
      init4 = init1;                                 /* aa: */
    }
    if(s > 0) {
      if(t == 0) { /* BCs */
	init2 = -s * M_LN2;                          /* Aa: log(2 ^ -s) */
	init1 = log1p(-exp(init2));                  /* AA: log(1 - 2^-s) */
	init3 = 0;
	init4 = 0;
      }
      if(t > 0) {  /* BCsFt */
	double sm2,t1m2;
	sm2 = -s * M_LN2;                            /* -s * log(2) = log(2 ^ -s) */
	t1m2 = -(1 + t) * M_LN2;                     /* -2t * log(2) = log(2 ^ -(t+1)) */

	init2 = sm2 + t1m2;                          /* Aa: log(2^-(s+t+1)) */
	init3 = init2;                               /* aA: log(2^-(s+t+1)) */
	init4 = subtrlog(sm2 - M_LN2, init2);        /* aa: log(2^-(s+1) - 2^-(s+t+1)) */
	init1 = addlog(log1p(-exp(sm2)), init4);     /* AA: log((1-2^-s) + (2^-(s+1) - 2^-(s+t+1))) */
      }
    }
  }

  switch(true_gen) {
  case 1: return(init1);
  case 2: return(init2);
  case 3: return(init3);
  case 4: return(init4);
  }
  return(0.0); /* should not get here */
}

double emit_bcsftb(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  if(cross_scheme[1] > 0) return(emit_f2b(obs_gen, true_gen, error_prob,cross_scheme));
  return(emit_bc(obs_gen, true_gen, error_prob,cross_scheme));
}


double step_bcsftb(int gen1, int gen2, double rf, double junk, int *cross_scheme) 
{
  static double oldrf = -1.0;
  static double transpr[10];
  static int s = -1;
  static int t = -1;
  
  if(s != cross_scheme[0] || t != cross_scheme[1] || fabs(rf - oldrf) > TOL) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    oldrf = rf;
    if(rf < TOL) rf = TOL;

    prob_bcsft(rf, s, t, transpr);

    /* expand when phase is known */
    if(t > 0) { /* only if Ft in play */
      transpr[1] /= 2.0; /* B1 split */
      transpr[6] /= 2.0; /* B0 split */
      transpr[3] /= 2.0; /* D split */
      transpr[4] /= 2.0; /* E split */
      transpr[8] -= M_LN2; /* log(pr(gen1=2)) = log(pr(gen2=3)) */
    }

    /* put probabilities on log scale */
    int k;
    for(k=0; k<7; k++) {
      /*      if(transpr[k] > 0.0)  */
      transpr[k] = log(transpr[k]);
    }
  }

  double out;
  /* Find joint probability pr(gen1,gen2). */
  out = assign_bcsftb(gen1, gen2, transpr);

  /* Divide by marginal prob to get pr(gen2|gen1). */
  if(gen1 > 2) gen1--;
  out -= transpr[6+gen1];

  return(out);
}

double nrec_bcsftb(int gen1, int gen2, double rf, int *cross_scheme)
{
  static double oldrf = -1.0;
  static double transexp[10];
  static int s = -1;
  static int t = -1;
  
  if(s != cross_scheme[0] || t != cross_scheme[1] || fabs(rf - oldrf) > TOL) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    oldrf = rf;
    if(rf < TOL) rf = TOL;

    expect_bcsft(rf, s, t, transexp);
    
    /* reduce by half if t>0 *** NOT SURE IF THIS IS RIGHT THING TO DO? */
    if(t > 0) {
      int k;
      for(k=0; k<7; k++)
	transexp[k] /= 2;
    }
  }

  /* Return expected count. */
  return(assign_bcsftb(gen1, gen2, transexp));
}

/****************************************************************************/

void calc_genoprob_bcsft(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob,
                      double *genoprob) 
{
  /* cross_scheme is hidden in genoprob */
  int n_gen;
  n_gen = 2;
  if(genoprob[1] > 0) n_gen = 3;

  calc_genoprob(*n_ind, *n_mar, n_gen, geno, rf, rf, *error_prob,
                genoprob, init_bcsft, emit_bcsft, step_bcsft);
}
   
void calc_genoprob_special_bcsft(int *n_ind, int *n_mar, int *geno, 
			      double *rf, double *error_prob, double *genoprob) 
{
  /* cross_scheme is hidden in genoprob */
  int n_gen;
  n_gen = 2;
  if(genoprob[1] > 0) n_gen = 3;

  calc_genoprob_special(*n_ind, *n_mar, n_gen, geno, rf, rf, *error_prob, genoprob,
			init_bcsft, emit_bcsft, step_bcsft);
}
 
void sim_geno_bcsft(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws)
{
  /* cross_scheme is hidden in draws */
  int n_gen;
  n_gen = 2;
  if(draws[1] > 0) n_gen = 3;

  sim_geno(*n_ind, *n_pos, n_gen, *n_draws, geno, rf, rf, *error_prob,
	   draws, init_bcsft, emit_bcsft, step_bcsft);
}

double comploglik(double rf, int n_gen, double *countmat, int *cross_scheme)
{
  static double transpr[10];
  static double probmat[15];
  static double oldrf = -1.0;
  static int s = -1;
  static int t = -1;
  int obs1,obs2,tmp1;

  if(s != cross_scheme[0] || t != cross_scheme[1] || fabs(rf - oldrf) > TOL) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    oldrf = rf;
    if(rf < TOL) rf = TOL;

    /* compute probabilities */
    prob_bcsft(rf, s, t, transpr);
    transpr[3] += transpr[4];

    for(obs2=1; obs2<=n_gen; obs2++) {
      tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
      for(obs1=1; obs1<=obs2; obs1++)
	probmat[obs1 + tmp1] = assign_bcsftd(n_gen, obs1, obs2, transpr);
    }
  }

  double lod,temp;
 
  /* compute log likelihood */
  lod = 0.0;
  for(obs2=1; obs2<=n_gen; obs2++) {
    tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
    for(obs1=1; obs1<=obs2; obs1++) {
      temp = countmat[obs1 + tmp1];
      if(temp > 0.0)
	lod += temp * log(probmat[obs1 + tmp1]);
    }
  }
  return(lod);
}

double golden(double *countmat, int n_gen, int maxit, double tol, int *cross_scheme)
{
  /* Golden section search. */
  /* en.wikipedia.org/wiki/Golden_section_search */

  static double resphi = 0.0;
  double x[4],y[4];
  int iter;

  if(resphi == 0.0)
    resphi = 1.5 - sqrt(5.0) / 2.0;

  x[0] = 0.0;
  x[2] = 0.5;
  y[0] = comploglik(0.0, n_gen, countmat, cross_scheme);
  y[2] = comploglik(0.5, n_gen, countmat, cross_scheme);

  if(y[2] < y[0]) {
    x[1] = x[0];
    x[0] = x[2];
    x[2] = x[1];
    y[1] = y[0];
    y[0] = y[2];
    y[2] = y[1];
  }
  
  x[1] = x[0] + resphi * (x[2] - x[0]);
  y[1] = comploglik(x[1], n_gen, countmat, cross_scheme);

  /* x0 and x2 are the current bounds; the minimum is between them.
   * x1 is the center point, which is closer to x0 than to x2. */

  for(iter=0; iter<maxit; iter++) {
    /* Create a new possible center in the area between x1 and x2, closer to x1. */
    x[3] = x[1] + resphi * (x[2] - x[1]);

    /* Evaluate termination criterion */
    if(fabs(x[2] - x[0]) < tol)
      break;
 
    y[3] = comploglik(x[3], n_gen, countmat, cross_scheme);

    if(y[3] >= y[1]) {
      x[0] = x[1];
      x[1] = x[3];
      y[0] = y[1];
      y[1] = y[3];
    }
    else {
      x[2] = x[0];
      x[0] = x[3];
      y[2] = y[0];
      y[0] = y[3];
    }
  }
  /* handle boundary situations cleanly */
  if((x[0] == 0.0 && y[0] >= y[1]) || (x[2] == 0.0 && y[2] >= y[1])) return(0.0);
  if((x[0] == 0.5 && y[0] >= y[1]) || (x[2] == 0.5 && y[2] >= y[1])) return(0.5);

  x[1] = (x[2] + x[0]) / 2.0;
  /* make negative if does not converge */
  if(iter >= maxit)
    x[1] = - x[1];
  return(x[1]);
}

/* est_map_bcsft maps correctly for BC and F2 */
/* but not for higher order crosses BCsFt */
/* need instead to use golden section search for rf */

void est_map0_bcsft(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose)
{
  /* cross scheme hidden in loglik argument; used by hmm_bcsft */
  int n_gen,s,t;
  n_gen = 2;
  s = (int) ftrunc(*loglik / 1000.0);
  t = ((int) *loglik) - 1000 * s;
  if(t > 0) n_gen = 4;

  est_map(*n_ind, *n_mar, n_gen, geno, rf, rf, *error_prob, 
	  init_bcsftb, emit_bcsftb, step_bcsftb, nrec_bcsftb, nrec_bcsftb,
	  loglik, *maxit, *tol, 0, *verbose);
}

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void init_step_map(double *cur_rf, double *cur_rf2, int n_gen, int n_mar, int *cross_scheme, double **probmat)
{
  int j,obs1,obs2,tmp1;
  
  for(j=0; j<n_mar; j++) {
    for(obs2=1; obs2<=n_gen; obs2++) {
      tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
      for(obs1=1; obs1<=obs2; obs1++)
	probmat[j][obs1 + tmp1] = step_bcsftb(obs1, obs2, cur_rf[j], cur_rf2[j], cross_scheme);
    }
  }
}
double step_bcsftc(int obs1, int obs2, int mar, double **probmat)
{
  int tmp1;
  
  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    tmp1 = obs2;
    obs2 = obs1;
    obs1 = tmp1;
  }
  tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
  return(probmat[mar][obs1 + tmp1]);
}
void est_map_bcsft(int *n_ind, int *n_mar, int *geno, double *rf, 
		   double *error_prob, double *loglik, int *maxit, 
		   double *tol, int *verbose)
{
  int i, j, j2, v, v2, it, flag=0, **Geno, ndigits,tmp1,tmp2,sexsp;
  double s, **alpha, **beta, **gamma, *cur_rf, *cur_rf2, *rf2;
  double curloglik, maxdif, temp;
  char pattern[100], text[200];
  int cross_scheme[2];
  double **countmat,**countmat2,**probmat;

  /* Force sex specific to be 0 for now. */
  sexsp = 0;
  rf2 = rf;

  /* cross scheme hidden in loglik argument; used by hmm_bcsft */
  cross_scheme[0] = (int) ftrunc(*loglik / 1000.0);
  cross_scheme[1] = ((int) *loglik) - 1000 * cross_scheme[0];
  *loglik = 0.0;

  /* n_gen inferred from cross scheme */
  int n_gen;
  n_gen = 2;
  if(cross_scheme[1] > 0) n_gen = 4;

  /* allocate space for beta and reorganize geno */
  reorg_geno(*n_ind, *n_mar, geno, &Geno);
  allocate_alpha(*n_mar, n_gen, &alpha);
  allocate_alpha(*n_mar, n_gen, &beta);
  allocate_dmatrix(n_gen, n_gen, &gamma);
  allocate_double(*n_mar-1, &cur_rf);
  allocate_double(*n_mar-1, &cur_rf2);
  allocate_dmatrix(*n_mar, 10, &countmat);
  if(sexsp)
    allocate_dmatrix(*n_mar, 10, &countmat2);
  allocate_dmatrix(*n_mar, 10, &probmat);

  /* digits in verbose output */
  if(*verbose) {
    ndigits = (int)ceil(-log10(*tol));
    if(ndigits > 16) ndigits=16;
    sprintf(pattern, "%s%d.%df", "%", ndigits+3, ndigits+1);
  }

  /* begin EM algorithm */
  for(it=0; it<*maxit; it++) {

    for(j=0; j<*n_mar-1; j++) {
      cur_rf[j] = cur_rf2[j] = rf[j];
      rf[j] = 0.0;
      if(sexsp) {
	cur_rf2[j] = rf2[j];
	rf2[j] = 0.0;
      }
    }

    /* initialize stepf calculations */
    init_step_map(cur_rf, cur_rf2, n_gen, *n_mar, cross_scheme, probmat);
    
    /*** reset counts for countmat ***/
    for(j=0; j<*n_mar-1; j++) {
      for(v2=0; v2<n_gen; v2++) {
	tmp1 = (v2 * (v2 + 1)) / 2;
	for(v=0; v<=v2; v++) {
	  countmat[j][v + tmp1] = 0.0;
	  if(sexsp)
	    countmat2[j][v + tmp1] = 0.0;
	}
      }
    }
    
    for(i=0; i<*n_ind; i++) { /* i = individual */
      
      R_CheckUserInterrupt(); /* check for ^C */
      
      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	alpha[v][0] = init_bcsftb(v+1, cross_scheme) +
	  emit_bcsftb(Geno[0][i], v+1, *error_prob, cross_scheme);
	beta[v][*n_mar-1] = 0.0;
      }
      
      /* forward-backward equations */
      
      /* reorder evaluation to reduce switching of rf in calculations */
      
      /* alpha calcs with cur_rf[j-1] */
      for(j=1; j<*n_mar; j++) {
	for(v=0; v<n_gen; v++) {
	  alpha[v][j] = alpha[0][j-1] + step_bcsftc(1, v+1, j-1, probmat);
	  for(v2=1; v2<n_gen; v2++)
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + step_bcsftc(v2+1, v+1, j-1, probmat));
	  alpha[v][j] += emit_bcsftb(Geno[j][i],v+1,*error_prob, cross_scheme);
	}
      }
      
      /* beta calcs with cur_rf[j2] */
      for(j2=*n_mar-2; j2>=0; j2--) {
	for(v=0; v<n_gen; v++) {
	  beta[v][j2] = beta[0][j2+1] + step_bcsftc(v+1, 1, j2, probmat) + 
	    emit_bcsftb(Geno[j2+1][i],1,*error_prob, cross_scheme);
	  for(v2=1; v2<n_gen; v2++)
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + step_bcsftc(v+1, v2+1, j2, probmat) +
				 emit_bcsftb(Geno[j2+1][i],v2+1,*error_prob, cross_scheme));
	}
      }
      
      for(j=0; j<*n_mar-1; j++) {
	
	/* calculate gamma = log Pr(v, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + step_bcsftc(v+1, v2+1, j, probmat) +
	      emit_bcsftb(Geno[j+1][i], v2+1, *error_prob, cross_scheme);
	    
	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}
	
	for(v=0; v<n_gen; v++) {
	  tmp2 = (v * (v + 1)) / 2;
	  for(v2=0; v2<n_gen; v2++) {
	    /* rf[j] += nrec_bcsftb(v+1,v2+1, cur_rf[j], cross_scheme) * exp(gamma[v][v2] - s); */
	    temp = exp(gamma[v][v2] - s);
	    if(v2 <= v)
	      tmp1 = v2 + tmp2;
	    else
	      tmp1 = v + (v2 * (v2 + 1)) / 2;
	    countmat[j][tmp1] += temp;
	    if(sexsp) { /* this cannot be right */
	      /* rf2[j] += nrec_bcsftb(v+1,v2+1, cur_rf[j], cross_scheme) * exp(gamma[v][v2] - s); */
	      countmat2[j][tmp1] += temp;
	    }
	  }
	}
     }

    } /* loop over individuals */

    /* rescale */
    for(j=0; j<*n_mar-1; j++) {
      /* rf[j] /= (double)*n_ind; */
      /*** this has to change ***/
      rf[j] = golden(countmat[j], n_gen, *maxit, *tol, cross_scheme);

      if(rf[j] < *tol/1000.0) rf[j] = *tol/1000.0;
      else if(rf[j] > 0.5-*tol/1000.0) rf[j] = 0.5-*tol/1000.0;
      
      if(sexsp) {
	rf2[j] /= (double)*n_ind;
	if(rf2[j] < *tol/1000.0) rf2[j] = *tol/1000.0;
	else if(rf2[j] > 0.5-*tol/1000.0) rf2[j] = 0.5-*tol/1000.0;
      }
      else rf2[j] = rf[j];
    }

    if(*verbose>1) {
      /* print estimates as we go along*/
      Rprintf("   %4d ", it+1);
      maxdif=0.0;
      for(j=0; j<*n_mar-1; j++) {
	temp = fabs(rf[j] - cur_rf[j])/(cur_rf[j]+*tol*100.0);
	if(maxdif < temp) maxdif = temp;
	if(sexsp) {
	  temp = fabs(rf2[j] - cur_rf2[j])/(cur_rf2[j]+*tol*100.0);
	  if(maxdif < temp) maxdif = temp;
	}
      }
      sprintf(text, "%s%s\n", "  max rel've change = ", pattern);
      Rprintf(text, maxdif);
    }

    /* check convergence */
    for(j=0, flag=0; j<*n_mar-1; j++) {
      if(fabs(rf[j] - cur_rf[j]) > *tol*(cur_rf[j]+*tol*100.0) || 
	 (sexsp && fabs(rf2[j] - cur_rf2[j]) > *tol*(cur_rf2[j]+*tol*100.0))) {
	flag = 1; 
	break;
      }
    }

    if(!flag) break;

  } /* end EM algorithm */
  
  if(flag) warning("Didn't converge!\n");

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<*n_ind; i++) { /* i = individual */
    /* initialize alpha */
    for(v=0; v<n_gen; v++) {
      alpha[v][0] = init_bcsftb(v+1, cross_scheme) + emit_bcsftb(Geno[0][i], v+1, *error_prob, cross_scheme);
    }
    /* forward equations */
    for(j=1; j<*n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + 
	  step_bcsftb(1, v+1, rf[j-1], rf2[j-1], cross_scheme);
	
	for(v2=1; v2<n_gen; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       step_bcsftb(v2+1,v+1,rf[j-1],rf2[j-1], cross_scheme));

	alpha[v][j] += emit_bcsftb(Geno[j][i],v+1,*error_prob, cross_scheme);
      }
    }

    curloglik = alpha[0][*n_mar-1];
    for(v=1; v<n_gen; v++)
      curloglik = addlog(curloglik, alpha[v][*n_mar-1]);

    *loglik += curloglik;
  }

  if(*verbose) {
    if(*verbose < 2) {
      /* print final estimates */
      Rprintf("  no. iterations = %d\n", it+1);
      maxdif=0.0;
      for(j=0; j<*n_mar-1; j++) {
	temp = fabs(rf[j] - cur_rf[j])/(cur_rf[j]+*tol*100.0);
	if(maxdif < temp) maxdif = temp;
	if(sexsp) {
	  temp = fabs(rf2[j] - cur_rf2[j])/(cur_rf2[j]+*tol*100.0);
	  if(maxdif < temp) maxdif = temp;
	}
      }
      sprintf(text, "%s%s\n", "  max rel've change at last step = ", pattern);
      Rprintf(text, maxdif);
    }
    
    Rprintf("  loglik: %10.4lf\n\n", *loglik);
  }

}

void argmax_geno_bcsft(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax)
{		    
  /* cross_scheme is hidden in argmax */
  int n_gen;
  n_gen = 2;
  if(argmax[1] > 0) n_gen = 3;

  argmax_geno(*n_ind, *n_pos, n_gen, geno, rf, rf, *error_prob,
	      argmax, init_bcsft, emit_bcsft, step_bcsft);
}

double myratio(double ct, double prob)
{
  if(prob > 0.0) prob = ct / prob;
  return(prob);
}

double nrec2_bcsft(int obs1, int obs2, double rf, int *cross_scheme)
{
  if((obs1 == 0) || (obs2 == 0)) return(0.0); /* shouldn't get here */

  static double oldrf = -1.0;
  static double transpr[10],transct[10];
  static int s = -1;
  static int t = -1;
  
  if(s != cross_scheme[0] || t != cross_scheme[1] || fabs(rf - oldrf) > TOL) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    oldrf = rf;
    if(rf < TOL) rf = TOL;

    prob_bcsft(rf, s, t, transpr);
    transpr[3] += transpr[4];
    count_bcsft(rf, s, t, transct);
    transct[3] += transct[4];
  }

  return(myratio(assign_bcsftc(obs1, obs2, transct), assign_bcsftc(obs1, obs2, transpr)));
}
double logprec_bcsft(int obs1, int obs2, double rf, int *cross_scheme)
{
  /* this routine is not correct yet */
  if((obs1 == 0) || (obs2 == 0)) return(log(-1.0)); /* shouldn't get here */

  static double transpr[10],margin[4];
  static double oldrf = -1.0;
  static int s = -1;
  static int t = -1;
  
  if(s != cross_scheme[0] || t != cross_scheme[1] || fabs(rf - oldrf) > TOL) {
    s = cross_scheme[0];
    t = cross_scheme[1];

    oldrf = rf;
    if(rf < TOL) rf = TOL;

    prob_bcsft(rf, s, t, transpr);
    transpr[3] += transpr[4];

    int k;
    for(k=1; k<4; k++)
      margin[k] = exp(transpr[6+k]);
  }

  double out1,out2,out3,out4;

  if(obs1 < 4) {
    if(obs2 < 4) /* known genotypes for both loci */
      return(log(assign_bcsft(obs1, obs2, transpr) / margin[obs1]));
    /* dominant genotype for locus 2 */
    out1 = assign_bcsft(obs1, obs2 - 3, transpr);
    out2 = assign_bcsft(obs1, obs2 - 2, transpr);
    return(log((out1 + out2) / margin[obs1]));
  }

  /* obs1 is 4 or 5 */
  double denom;
  denom =  margin[obs1 - 3] + margin[obs1 - 2];

  if(obs2 < 4) {
    /* dominant genotype for locus 1 */
    out1 = assign_bcsft(obs1 - 3, obs2, transpr);
    out2 = assign_bcsft(obs1 - 2, obs2, transpr);
    return(log((out1 + out2) / denom));
  }
  /* both loci have dominant genotypes */
  out1 = assign_bcsft(obs1 - 3, obs2 - 3, transpr);
  out2 = assign_bcsft(obs1 - 2, obs2 - 2, transpr);
  out3 = assign_bcsft(obs1 - 3, obs2 - 2, transpr);
  out4 = assign_bcsft(obs1 - 2, obs2 - 3, transpr);
  return(log((out1 + out2 + out3 + out4) / denom));
}

/* est_rfo_bcsft is ok for bc,f2, not for advanced crosses */
/* naive estimate based on number of recombinations in nrec2 does not generalize */
/* now use est_rf_bcsft with golden section search */

void est_rfo_bcsft(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol)
{
  int BC_gen,F_gen,meioses_per;
  BC_gen = rf[0];
  F_gen = rf[1];

  meioses_per = 2 * F_gen;
  if(BC_gen <= 0)
    meioses_per -= 2;
  if(BC_gen > 0) 
    meioses_per += BC_gen;

  est_rf(*n_ind, *n_mar, geno, rf, nrec2_bcsft, logprec_bcsft, 
	 *maxit, *tol, meioses_per);
}

void est_rf_bcsft(int *n_ind, int *n_mar, int *geno, double *rf, 
		  int *maxit, double *tol)
{
  int i, j1, j2, **Geno, n_mei=0, flag=0;
  double **Rf, next_rf=0.0;
  int cross_scheme[2];

  /* cross_scheme is hidden in rf */
  cross_scheme[0] = rf[0];
  cross_scheme[1] = rf[1];
  rf[0] = 0.0;
  rf[1] = 0.0;

  int meioses_per;
  meioses_per = 2 * cross_scheme[1];
  if(cross_scheme[0] <= 0)
    meioses_per -= 2;
  if(cross_scheme[0] > 0) 
    meioses_per += cross_scheme[0];

  /* reorganize geno and rf */
  reorg_geno(*n_ind, *n_mar, geno, &Geno);
  reorg_errlod(*n_mar, *n_mar, rf, &Rf);

  int obs1,obs2,n_gen,tmp1;
  double countmat[15];
  double temp,logprecval;

  n_gen = 2;
  if(cross_scheme[1] > 0)
    n_gen = 5;

  for(j1=0; j1<*n_mar; j1++) {

    /* count number of meioses */
    for(i=0, n_mei=0; i<*n_ind; i++) 
      if(Geno[j1][i] != 0) n_mei += meioses_per;
    Rf[j1][j1] = (double) n_mei;
    
    R_CheckUserInterrupt(); /* check for ^C */

    for(j2=j1+1; j2<*n_mar; j2++) {
      for(obs2=1; obs2<=n_gen; obs2++) {
	tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
	for(obs1=1; obs1<=obs2; obs1++) {
	  countmat[obs1 + tmp1] = 0.0;
	}
      }
      /* count meioses */
      n_mei = flag = 0;
      for(i=0; i<*n_ind; i++) {
	if(Geno[j1][i] != 0 && Geno[j2][i] != 0) {
	
	  /* make obs1 <= obs2 */
	  obs1 = Geno[j1][i];
	  obs2 = Geno[j2][i];
	  if(obs1 > obs2) {
	    tmp1 = obs2;
	    obs2 = obs1;
	    obs1 = tmp1;
	  }
	  tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
	  /* increment count by one */
	  countmat[obs1 + tmp1] += 1.0;
	}
      }
      /* check if marker pair is informative */
      for(obs2=1; obs2<=n_gen; obs2++) {
	tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
	for(obs1=1; obs1<=obs2; obs1++) {
	  temp = countmat[obs1 + tmp1];
	  if(temp > 0.0) {
	    logprecval = logprec_bcsft(obs1, obs2, 0.5, cross_scheme) -
	      logprec_bcsft(obs1, obs2, TOL, cross_scheme);
	    if(fabs(logprecval) > TOL) {
	      n_mei += (int) temp;
	      flag = 1;
	    }
	  }
	}
      }
      
      if(n_mei != 0 && flag == 1) {
	n_mei *= meioses_per;
	flag = 1;

	/* use golden section search of log likelihood instead of EM */
	next_rf = golden(countmat, n_gen, *maxit, *tol, cross_scheme);

	if(next_rf < 0.0) {
	  flag = 0;
	  next_rf = - next_rf;
	}
	if(!flag) warning("Markers (%d,%d) didn't converge\n", j1+1, j2+1);

	/* calculate LOD score */
	Rf[j1][j2] = next_rf;
	logprecval = 0.0;

	for(obs2=1; obs2<=n_gen; obs2++) {
	  tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
	  for(obs1=1; obs1<=obs2; obs1++) {
	    temp = countmat[obs1 + tmp1];
	    if(temp > 0.0)
	      logprecval += temp * (logprec_bcsft(obs1,obs2, next_rf, cross_scheme) -
				    logprec_bcsft(obs1,obs2, 0.5, cross_scheme));
	  }
	}
	Rf[j2][j1] = logprecval / log(10.0);
      }
      else { /* no informative meioses */
	Rf[j1][j2] = NA_REAL;
	Rf[j2][j1] = 0.0;
      }
	  
    } /* end loops over markers */
  }
}

void calc_pairprob_bcsft(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob) 
{
  /* cross_scheme is hidden in genoprob */
  int n_gen;
  n_gen = 2;
  if(genoprob[1] > 0) n_gen = 3;

  calc_pairprob(*n_ind, *n_mar, n_gen, geno, rf, rf, *error_prob, genoprob,
		pairprob, init_bcsft, emit_bcsft, step_bcsft);
}

void marker_loglik_bcsft(int *n_ind, int *geno,
		      double *error_prob, double *loglik)
{
  /* cross scheme hidden in loglik argument; used by hmm_bcsft */
  int n_gen,s,t;
  n_gen = 2;
  s = (int) ftrunc(*loglik / 1000.0);
  t = ((int) *loglik) - 1000 * s;
  if(t > 0) n_gen = 4;

  marker_loglik(*n_ind, n_gen, geno, *error_prob, 
		init_bcsftb, emit_bcsftb, loglik);
}

/*******************************************************************************************/
/* Reachable in BCs: */
/* D : AB.ab */
/* B1: Ab.ab, aB.ab */
/* A1: ab.ab */
/* Reachable only in Ft: */
/* A0: AB.AB */
/* B0: AB.Ab, AB.aB */
/* C : Ab.Ab, aB.aB */
/* E : Ab.aB */

void prob_bcs(double rf, int s, double *transpr)
{
  int k;
  double ws,s2,ss;
  /*#####BCs probabilities: */

  for(k=0; k<10; k++)
    transpr[k] = 0.0;
  transpr[3] = 1.0; /* PbD */

  if(s > 0) {
    ss = s;
    ws = R_pow(1.0 - rf, ss);
    s2 = R_pow(2.0, ss);

    /* PbA + 2*PbB + PbD = 1.0 regardless of s */
    transpr[0] = (s2 - 2.0 + ws)/ s2; /* PbA1 = PbDA1 + PbDBA1 */
    transpr[1] = (1.0 - ws)/ s2;      /* PbB1 = PbDB1 */
    transpr[3] = ws / s2;             /* PbD = PbDD */

    /* marginal probabilities for one marker */
    transpr[8] = -ss * M_LN2;             /* Aa */
    transpr[7] = log1p(-exp(transpr[8])); /* AA */
  }
}
void prob_ft(double rf, int t, double *transpr)
{
  int k;
  double t1,t2,t1m2,w,w2,r2,rw;
  double beta,gamma,beta1,sbeta1,sgamma1,SDt,SEt,sbetaBA,gamma1,beta2m1;
  /* compute transition probabilities to leave double het states */
  t1 = t - 1.0;
  t2 = R_pow(2.0, t); /* 2^t */
  t1m2 = 2.0 / t2;
  w = 1.0 - rf;
  w2 = w * w;
  r2 = rf * rf;
  rw = w * rf;

  for(k=0; k<10; k++)
    transpr[k] = 0.0;

  /* A = prob go from AB.ab to AB.AB at step t */
  /* D1 -> Dk or Ek -> Ak+1 -> At OR D1 -> Dk or Ek -> Bk+1 -> Ak+s -> At */
  beta = (w2 + r2) / 2.0; 
  gamma = (w2 - r2) / 2.0; 
  beta1 = R_pow(beta, t1);
  gamma1 = R_pow((w2 - r2) / 2.0, t1);
  /* beta1 = prob still in D or E at step t */
  /* sbeta1 = sum of beta1[k] from 1 to t-1 */
  /* Dt and Et depend on sbeta1 and sgamma1 */
  sbeta1 = (1.0 - beta1) / (1.0 - beta);
  sgamma1 = (1.0 - R_pow(gamma, t1)) / (1.0 - gamma);
  SDt = (sbeta1 + sgamma1) / 8.0;
  SEt = (sbeta1 - sgamma1) / 8.0;
  beta2m1 = 1.0 - 2.0 * beta;

  /* old code--make sure new code works first */
  /* w2pr2 = (w2 + r2); */
  /* sbetaBA = (rw / 2.0) * (sbeta1 - 2.0 * ((2.0 / t2) - beta1) / (1.0 - 2.0 * beta)); */
  /* transpr[1] = (2.0 * rw / t2) * ((1.0 - R_pow(w2pr2, t1)) / (1 - w2pr2)); */

  double s2beta1,sbeta2,s2beta2;
  s2beta1 = (t1m2 - beta1) / beta2m1;                   /* sum from 1 to t-1 of of (2*beta)^(k-1). */
  transpr[1] = rw * s2beta1;                            /* PfB1 = PfDB */
  transpr[6] = transpr[1];                              /* PfB0 = PfB1 */

  /* sbetaBA = sum beta1[k] * rw/2 * prob(B->A in remaining steps) */
  sbeta2 = 0.0;
  if(t > 2.0) sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); /* sum of beta^(k-1) from 1 to (t-1) */
  s2beta2 = (2.0 * t1m2 - (beta1 / beta)) / beta2m1;        /* sum of (2*beta)^(k-1) from 1 to t-2 */
  sbetaBA = 0.5 * rw * (sbeta2 - s2beta2);

  transpr[0] = SDt * w2 + SEt * r2 + sbetaBA;           /* PfA1 = PfDA + PfDEA + PfDBA */
  transpr[5] = transpr[0];                              /* PfA0 = PfA1 */

  transpr[2] = SDt * r2 + SEt * w2 + sbetaBA;           /* PfbC = PfDC + PfDEC + PfDBC */
  transpr[3] = (beta1 + gamma1) / 2.0;                  /* PfD = PfDD */
  transpr[4] = (beta1 - gamma1) / 2.0;                  /* PfE = PfDE */

  /* marginal probabilities for one marker */
  transpr[8] = -t1 * M_LN2;                             /* Aa */
  transpr[7] = log1p(-exp(transpr[8])) - M_LN2;         /* AA */
  transpr[9] = transpr[7];                              /* aa */

  return;
}
void prob_bcsft(double rf, int s, int t, double *transpr)
{ 
  /* BCsFt probabilities */
  /* D->x: PfDx using calculations from proc_ft, but need to multiple by PbD = BCs probability */
  /* PbDfx = PbD * PfDx (PfDx = probability of D->x in Ft) */

  if(s == 0) {
    prob_ft(rf, t, transpr);
    return;
  }

  if(t == 0) {
    prob_bcs(rf, s, transpr);
    return;
  }

  /* could silently pass transft via transpr with some care */
  double transbcs[10],transft[10];
  prob_bcs(rf, s, transbcs);
  prob_ft(rf, t+1, transft);

  double tm2,PbBfC;
  tm2 = R_pow(0.5, (double) t);
  PbBfC = transbcs[1] * 0.5 * (1.0 - tm2); /* PbBfC = PbB * sum from k=0 to t-1 (1/2)^k * (1/4) */

  transpr[5] = transbcs[3] * transft[0];               /* PbfA0 = PbD * PfDA1 */
  transpr[0] = transpr[5] + 2.0 * PbBfC + transbcs[0]; /* PbfA1 = PfA0 + 2 * PbBfC + PbA1 * 1 */

  transpr[6] = transbcs[3] * transft[1];               /* PbfB0 = PbD * PfDB1 */
  transpr[1] = transpr[6] + transbcs[1] * tm2;         /* PbfB1 = PfB0 + PbB * (1/2)^t */

  transpr[2] = transbcs[3] * transft[2] + PbBfC;       /* PbfC = PbD * PfDC + PbBfC */
  transpr[3] = transbcs[3] * transft[3];               /* PbfD = PbD * PfDD */
  transpr[4] = transbcs[3] * transft[4];               /* PbfE = PbD * PfDE */

  /* marginal probabilities for one marker */
  double s2,t2;
  s2 = -s * M_LN2;
  t2 = -t * M_LN2;
  transpr[8] = s2 + t2;                                /* Aa: log(2^-(s+t)) */
  transpr[9] = s2 + log1p(-exp(t2)) - M_LN2;           /* aa: log(2^-(s+1) * (1-2^-t)) */
  transpr[7] = addlog(log1p(-exp(s2)), transpr[9]);    /* AA: log(1-2^-s + 2^-(s+1) * (1-2^-t)) */

  return;
}
/*######################################################################################### */
double kptothek(double t, double p, double ptothet)
{
  /* sum of k * p^k from 0 to t-1 */
  double tmp;
  tmp = (1.0 - p);
  return((p - t * ptothet + (t - 1.0) * p * ptothet) / (tmp * tmp));
}
/*############################################################################################### */
void count_bcs(double rf, int s, double *transbcs, double *transct)
{
  int k;

  for(k=2; k<10; k++)
    transct[k] = 0.0; /* rest are 0 count for numerator */

  /* PbDA1 = (1 - rf) * (1 - PbD) / (1 + rf) */
  /* PbDBA1 = 1 - PbD + 2 * PbB1 + PbDA1 */
  /* NbA1 = 0 * PbDA1 + 1 * PbDBA1 */
  transct[0] = 1.0 - transbcs[3] - 2.0 * transbcs[1] - (1.0 - rf) * (1.0 - transbcs[3]) / (1.0 + rf);

  /* NbB1 = 1 * PbB1 */
  transct[1] = transbcs[1];
  return;
}

void count_ft(double rf, int t, double *transct)
{
  int k;

  if(t < 2) {
    for(k=0; k<10; k++)
      transct[k] = 0.0;
    return;
  }

  double t1,t1m2,r2,w2,rw,alpha;
  double beta,beta1,beta2,sbeta1,sbeta2,s2beta1,s2beta2;
  double gamma,gamma1,gamma2,sgamma1,sgamma2,s2gamma1,s2gamma2;
  double k1b,k2b,k1g,k2g,k1b2,k2b2,k1g2,k2g2;

  t1 = t - 1.0;
  t1m2 = R_pow(2.0, -t1);

  r2 = rf * rf;
  w2 = (1.0 - rf) * (1.0 - rf);
  rw = rf * (1.0 - rf);

  alpha = r2 / (w2 + r2);

  /* beta = probability of D or E at each step. */
  beta = (w2 + r2) / 2.0;
  beta1 = R_pow(beta, t1);
  beta2 = 1.0;
  if(t > 2.0) beta2 = beta1 / beta;
  /* Calculations for D->DE->B->AC */
  sbeta1 = (1.0 - beta1) / (1.0 - beta); /* SFt */
  sbeta2 = 0.0;
  if(t > 2.0) sbeta2 = (1.0 - beta1 / beta) / (1.0 - beta); /* SF(t-1) */
  s2beta1 = (t1m2 - beta1) / (1.0 - 2.0 * beta); /* sum from 1 to t-1 of of (2*beta)^(k-1). */
  s2beta2 = (2.0 * t1m2 - (beta1 / beta)) / (1.0 - 2.0 * beta); /* sum from 1 to t-2 of of (2*beta)^(k-1). */
  
  gamma = (w2 - r2) / 2.0; 
  gamma1 = 1.0;
  if(t > 1) gamma1 = R_pow(gamma, t1);
  gamma2 = 1.0;
  if(t > 2) gamma2 = R_pow(gamma, t1 - 1);
  sgamma1 = 1;
  sgamma2 = 0;
  s2gamma1 = t1m2;
  s2gamma2 = 0;
  if(t > 1) {
    sgamma2 = 1;
    s2gamma2 = t1m2 * 2.0;
  }
  if(gamma > 0) {
    sgamma1 = (1.0 - gamma1) / (1.0 - gamma); /* SGt */
    sgamma2 = (1.0 - gamma2) / (1.0 - gamma); /* SG(t-1) */
    s2gamma1 = (t1m2 - gamma1) / (1.0 - 2.0 * gamma); /* sum from 1 to t-1 of of (2*gamma)^(k-1). */
    s2gamma2 = (2.0 * t1m2 - (gamma1 / gamma)) / (1.0 - 2.0 * gamma); /* sum from 1 to t-2 of of (2*gamma)^(k-1). */
  }

  /* kptothek(t,p) = sum of k * p^k from 1 to t-1. */
  k1b = kptothek(t1, beta, beta1) / beta;
  k2b = t1m2 * kptothek(t1, 2.0 * beta, beta1 / t1m2) / (2 * beta);
  k1g = 0.0;
  k2g = 0.0;
  k1g2 = 0.0;
  k2g2 = 0.0;
  k1b2 = 0.0;
  k2b2 = 0.0;
  if(t > 2) {
    k1g = 1.0;
    k2g = t1m2;
    if(t > 3) {
      k1g2 = 1.0;
      k2g2 = 2.0 * t1m2;
    }
    k1b2 = kptothek(t1 - 1.0, beta, beta2) / beta;
    k2b2 = 2.0 * t1m2 * kptothek(t1 - 1.0, 2.0 * beta, beta2 / (2 * t1m2)) / (2 * beta);
  }
  if(gamma > 0) {
    /* Possible savings in doing the sum... */
    k1g = kptothek(t1, gamma, gamma1) / gamma;
    k2g = t1m2 * kptothek(t1, 2.0 * gamma, gamma1 / t1m2) / (2.0 * gamma);
    k1g2 = kptothek(t1 - 1.0, gamma, gamma2) / gamma;
    k2g2 = 2.0 * t1m2 * kptothek(t1 - 1.0, 2.0 * gamma, gamma2 / (2.0 * t1m2)) / (2.0 * gamma);
  }

  /* Below, x->y->z->w refers to generations 1, k, k+1, t
     In case of state at k+1 is B, there is a further transition to A or C at step k+s summed out */

  double ndda,ndea,nbab,nbag,NDDA,NDDC,NDEA,NDEC,NDBA,NEBA;
  ndda = (r2 / 2) * (k1b - k1g);
  NDDA = (w2 / 4) * ndda;
  NDDC = (r2 / 4) * (ndda + (sbeta1 + sgamma1));

  NDEA = 0.0;
  NDEC = 0.0;
  NDBA = 0.0;
  NEBA = 0.0;
  if(t > 2) {
    ndea = (r2 / 2.0) * (k1b + k1g);
    NDEA = (r2 / 4.0) * (ndea + (sbeta1 - sgamma1));
    NDEC = (w2 / 4.0) * ndea;
    
    nbab = rw * (0.25 * (sbeta2 - s2beta2) + (r2 / 2.0) * (0.5 * k1b2 - k2b2));
    nbag = rw * (0.25 * (sgamma2 - s2gamma2) - (r2 / 2.0) * (0.5 * k1g2 - k2g2));
    NDBA = nbab + nbag;
    if(t > 3) NEBA = nbab - nbag;
  }
  transct[0] = (NDDA + NDEA + NDBA + NEBA);      /* NfA1 = NfDDA + NfDEA + NfDBA + NfDEBA */
  transct[5] = transct[0];                    /* NfA0 = NfA1 */

  transct[1] = rw * (s2beta1 +  2.0 * r2 * k2b); /* NfB1 = NfDB + NfDEB */
  transct[6] = transct[1];                    /* NfB0 = NfB1 */

  transct[2] = (NDDC + NDEC + NDBA + NEBA);      /* NfC = NfDDC + NfDEC + NfDBA + NfDEBA */
  transct[3] = 0.5 * t1 * r2 * (beta2 - gamma2); /* nfD = NfDD */
  transct[4] = 0.5 * t1 * r2 * (beta2 + gamma2); /* NfE = NfDE */
  return;
}
void count_bcsft(double rf, int s, int t, double *transct)
{
  if(s == 0) {
    count_ft(rf, t, transct);
    return;
  }
  double transbcs[10];
  prob_bcs(rf, s, transbcs);
  if(t == 0) {
    count_bcs(rf, s, transbcs, transct);
    return;
  }

  double countbcs[10],countft[10];
  count_bcs(rf, s, transbcs, countbcs);
  count_ft(rf, t+1, countft);

  double tm2,PbBfC;
  tm2 = R_pow(0.5, (double) t);
  PbBfC = transbcs[1] * 0.5 * (1.0 - tm2); /* PbB1 * PfBC */

  /* numerator of counts */
  /* NbDfx = PbD * NfDx (NfDx = numerator of cound for D->x in Ft) */

  transct[5] = transbcs[3] * countft[0];                /* NbfA0 = PbD * NfA0 */
  transct[0] = transct[5] + countbcs[0] + 2.0 * PbBfC; /* NbfA1 = NfA0 + NbA1 * 1 + 2 * PbBfC */

  transct[6] = transbcs[3] * countft[1];                /* NbfB0 = PbD * NfDB */
  transct[1] = transct[6] + transbcs[1] * tm2;         /* NbfB1 = NfB0 + 1 * PbB1 * PfBB */

  transct[2] = transbcs[3] * countft[2] + PbBfC;        /* NbfC = PbD * NfDC + 1 * PbBfC */
  transct[3] = transbcs[3] * countft[3];                /* NbfD = PbD * NfDD */
  transct[4] = transbcs[3] * countft[4];                /* NbfE = PbD * NfDE */
  return;
}
void ratio_bcsft(double *transct, double *transexp)
{
  double tmp;
  int k;
  for(k=0; k<7; k++) {
    tmp = transexp[k];
    if(tmp > 0.0) tmp = transct[k] / tmp;
    transexp[k] = tmp;
  }
  return;
}
void expect_bcsft(double rf, int s, int t, double *transexp)
{
  double transct[10];
  prob_bcsft(rf, s, t, transexp);
  count_bcsft(rf, s, t, transct);

  ratio_bcsft(transct, transexp);
  return;
}
  
/* end of hmm_bcsft.c */
