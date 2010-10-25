/**********************************************************************
 * 
 * hmm_ft.c
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
 * Contains: init_ft, emit_ft, step_ft, init_ftb, emit_ftb, step_ftb,
 *           calc_genoprob_ft, calc_genoprob_special_ft, sim_geno_ft, est_map_ft, 
 *           argmax_geno_ft, errorlod_ft, calc_errorlod_ft, nrec2_ft,
 *           logprec_ft, est_rf_ft, calc_pairprob_ft, marker_loglik_ft
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
#include "hmm_ft.h"
#include "hmm_f2.h"

double init_ft(int true_gen, int *cross_scheme, double *cross_temps)
{
  if(true_gen==2) return(cross_temps[0]);
  else return(cross_temps[1]);
}

/* double emit_ft same as emit_f2? */
  
double step_ft(int gen1, int gen2, double rf, double junk, int *cross_scheme, double *cross_temps) 
{
  /* ref: Jiang and Zeng (1997 Genetics) */

  double t,t1, term1, term3, term4;

  t = cross_scheme[1];
  t1 = t - 1;
  term1 = exp(t1 * log(1.0 - 2.0 * rf * (1.0 - rf)));
  if(gen1 == 2) {
    if(gen2 == 2)
      return(log(term1));
    return(LN_05 + log(1.0 - term1));
  }
  if(gen2 == 2)
    return(cross_temps[3] + log(1.0 - term1));
  term3 = 1 + 2.0 * rf;
  term4 = exp(t * log(1 - 2 * rf)) / (2 * term3);
  if(gen1 == gen2)
    term3 = (cross_temps[2] / term3) - term4;
  else
    term3 = (2.0 * cross_temps[2] * rf / term3) + term4;
  return(cross_temps[3] + log(term3 - 1 + 0.5 * term1));
}

void calc_genoprob_ft(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob,
                      double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 3, geno, rf, rf, *error_prob,
                genoprob, init_ft, emit_f2, step_ft);
}
  
void sim_geno_ft(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws)
{
  sim_geno(*n_ind, *n_pos, 3, *n_draws, geno, rf, rf, *error_prob,
	   draws, init_ft, emit_f2, step_ft);
}

/* have to ask Karl about this. */

void est_map_ft(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose)
{
  est_map(*n_ind, *n_mar, 4, geno, rf, rf, *error_prob, 
	  init_ftb, emit_ftb, step_ftb, nrec_ftb, nrec_ftb,
	  loglik, *maxit, *tol, 0, *verbose);
}


void argmax_geno_ft(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax)
{		    
  argmax_geno(*n_ind, *n_pos, 3, geno, rf, rf, *error_prob,
	      argmax, init_ft, emit_ft, step_ft);
}


double errorlod_ft(int obs, double *prob, double error_prob)
{
  double p=0.0;

  switch(obs) {
  case 0: return(0.0);
  case 1: p=prob[0]; break;
  case 2: p=prob[1]; break;
  case 3: p=prob[2]; break;
  case 4: p=(prob[0]+prob[1]); break;
  case 5: p=(prob[1]+prob[2]); break;
  }
  
  p = (1.0-p)/p;

  if(obs==4 || obs==5) p *= (1.0-error_prob/2.0)/(error_prob/2.0);
  else p *= (1.0-error_prob)/error_prob;

  if(p < TOL) return(-12.0);
  else return(log10(p));
}


void calc_errorlod_ft(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 3, geno, *error_prob, genoprob,
		errlod, errorlod_ft);
}



double nrec2_ft(int obs1, int obs2, double rf)
{
  int temp;

  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    temp = obs2;
    obs2 = obs1;
    obs1 = temp;
  }

  switch(obs1) {
  case 1: 
    switch(obs2) {
    case 1: return(0.0);
    case 2: return(1.0);
    case 3: return(2.0);
    case 4: return(2.0*rf/(1.0+rf));
    case 5: return(2.0/(2.0-rf));
    }
  case 2:
    switch(obs2) {
    case 2: return(2*rf*rf/(rf*rf+(1.0-rf)*(1.0-rf)));
    case 3: return(1.0);
    case 4: case 5: return(rf*(1.0+rf)/(1.0-rf*(1.0-rf)));
    }
  case 3:
    switch(obs2) {
    case 3: return(0.0);
    case 4: return(2.0/(2.0-rf));
    case 5: return(2.0*rf/(1.0+rf));
    }
  case 4: case 5:
    if(obs1==obs2) return(4.0*rf/(3.0-2*rf+rf*rf));
    else return(2*rf*(2.0+rf)/(2+rf*rf));
  }
  return(log(-1.0)); /* shouldn't get here */
}


double logprec_ft(int obs1, int obs2, double rf)
{
  int temp;

  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    temp = obs2;
    obs2 = obs1;
    obs1 = temp;
  }

  switch(obs1) {
  case 1: 
    switch(obs2) {
    case 1: return(2.0*log(1.0-rf));
    case 2: return(LN_2 + log(rf) + log(1.0-rf));
    case 3: return(2.0*log(rf));
    case 4: return(log(1.0-rf*rf));
    case 5: return(log(1.0-(1.0-rf)*(1.0-rf)));
    }
  case 2:
    switch(obs2) {
    case 2: return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    case 3: return(log(rf*(1.0-rf)));
    case 4: case 5: return(log(1.0-rf*(1.0-rf)));
    }
  case 3:
    switch(obs2) {
    case 3: return(2.0*log(1.0-rf));
    case 4: return(log(1.0-(1.0-rf)*(1.0-rf)));
    case 5: return(log(1.0-rf*rf));
    }
  case 4: case 5:
    if(obs1==obs2) return(log(0.25*(1.0-rf)*(1.0-rf) + 0.5));
    else return(log(0.25*rf*rf+0.5));
  }
  return(log(-1.0)); /* shouldn't get here */
}


void est_rf_ft(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol)
{
  est_rf(*n_ind, *n_mar, geno, rf, nrec2_ft, logprec_ft, 
	 *maxit, *tol, 2);
}

void calc_pairprob_ft(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob) 
{
  calc_pairprob(*n_ind, *n_mar, 3, geno, rf, rf, *error_prob, genoprob,
		pairprob, init_ft, emit_ft, step_ft);
}

void marker_loglik_ft(int *n_ind, int *geno,
		      double *error_prob, double *loglik)
{
  marker_loglik(*n_ind, 4, geno, *error_prob, 
		init_ftb, emit_ftb, loglik);
}
  
/* end of hmm_ft.c */
