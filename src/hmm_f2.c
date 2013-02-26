/**********************************************************************
 * 
 * hmm_f2.c
 * 
 * copyright (c) 2001-2012, Karl W Broman
 *
 * last modified Apr, 2012
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
 * Contains: init_f2, emit_f2, step_f2, init_f2b, emit_f2b, step_f2b,
 *           calc_genoprob_f2, calc_genoprob_special_f2, sim_geno_f2, est_map_f2, 
 *           argmax_geno_f2, errorlod_f2, calc_errorlod_f2, nrec2_f2,
 *           logprec_f2, est_rf_f2, calc_pairprob_f2, marker_loglik_f2
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
#include "hmm_f2.h"

double init_f2(int true_gen, int *cross_scheme)
{
  if(true_gen==2) return(-M_LN2); /* ln(0.5) */
  else return(-2.0*M_LN2); /* ln(0.25) */
}

double emit_f2(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  switch(obs_gen) {
  case 0: return(0.0);
  case 1: case 2: case 3:
    if(obs_gen==true_gen) return(log(1.0-error_prob));
    else return(log(error_prob)-M_LN2);
  case 4: /* AA or AB (not BB) */
    if(true_gen != 3) return(log(1.0-error_prob/2.0));
    else return(log(error_prob));
  case 5: /* AB or BB (not AA) */
    if(true_gen != 1) return(log(1.0-error_prob/2.0));
    else return(log(error_prob));
  }
  return(0.0); /* shouldn't get here */
}
    
  
double step_f2(int gen1, int gen2, double rf, double junk, int *cross_scheme) 
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return(2.0*log(1.0-rf));
    case 2: return(M_LN2+log(1.0-rf)+log(rf));
    case 3: return(2.0*log(rf));
    }
  case 2:
    switch(gen2) {
    case 1: case 3: return(log(rf)+log(1.0-rf));
    case 2: return(log((1.0-rf)*(1.0-rf)+rf*rf));
    }
  case 3:
    switch(gen2) {
    case 1: return(2.0*log(rf));
    case 2: return(M_LN2+log(1.0-rf)+log(rf));
    case 3: return(2.0*log(1.0-rf));
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}






/* The following are the init, emit and step functions 
   when considering phase-known F2 genotypes 
   (i.e. the 4-state chain: AA, AB, BA, BB             */

double init_f2b(int true_gen, int *cross_scheme)
{
  return(-2.0*M_LN2);  /* ln(0.25) */
}

double emit_f2b(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  switch(obs_gen) {
  case 0: return(0.0);
  case 1: 
    switch(true_gen) {
    case 1: return(log(1.0-error_prob));
    case 2: case 3: case 4: return(log(error_prob)-M_LN2);
    }
  case 2: 
    switch(true_gen) {
    case 1: case 4: return(log(error_prob)-M_LN2);
    case 2: case 3: return(log(1.0-error_prob));
    }
  case 3:
    switch(true_gen) {
    case 4: return(log(1.0-error_prob));
    case 1: case 2: case 3: return(log(error_prob)-M_LN2);
    }
  case 4: /* AA or AB (not BB) */
    if(true_gen != 4) return(log(1.0-error_prob/2.0));
    else return(log(error_prob));
  case 5: /* AB or BB (not AA) */
    if(true_gen != 1) return(log(1.0-error_prob/2.0));
    else return(log(error_prob));
  }
  return(0.0); /* shouldn't get here */
}
    
  
double step_f2b(int gen1, int gen2, double rf, double junk, int *cross_scheme) 
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return(2.0*log(1.0-rf));
    case 2: case 3: return(log(1.0-rf)+log(rf));
    case 4: return(2.0*log(rf));
    }
  case 2: 
    switch(gen2) {
    case 1: case 4: return(log(rf)+log(1.0-rf));
    case 2: return(2.0*log(1.0-rf));
    case 3: return(2.0*log(rf));
    }
  case 3: 
    switch(gen2) {
    case 1: case 4: return(log(rf)+log(1.0-rf));
    case 3: return(2.0*log(1.0-rf));
    case 2: return(2.0*log(rf));
    }
  case 4:
    switch(gen2) {
    case 1: return(2.0*log(rf));
    case 2: case 3: return(log(1.0-rf)+log(rf));
    case 4: return(2.0*log(1.0-rf));
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double nrec_f2b(int gen1, int gen2, double rf, int *cross_scheme)
{
  switch(gen1) {
  case 1: 
    switch(gen2) {
    case 1: return(0.0);
    case 2: case 3: return(0.5);
    case 4: return(1.0);
    }
  case 2:
    switch(gen2) {
    case 1: case 4: return(0.5);
    case 2: return(0.0);
    case 3: return(1.0);
    }
  case 3:
    switch(gen2) {
    case 1: case 4: return(0.5);
    case 3: return(0.0);
    case 2: return(1.0);
    }
  case 4: 
    switch(gen2) {
    case 4: return(0.0);
    case 2: case 3: return(0.5);
    case 1: return(1.0);
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}




void calc_genoprob_f2(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 3, geno, rf, rf, *error_prob, genoprob,
		init_f2, emit_f2, step_f2);
}

  
void calc_genoprob_special_f2(int *n_ind, int *n_mar, int *geno, 
			      double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob_special(*n_ind, *n_mar, 3, geno, rf, rf, *error_prob, genoprob,
			init_f2, emit_f2, step_f2);
}

  
void sim_geno_f2(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws)
{
  sim_geno(*n_ind, *n_pos, 3, *n_draws, geno, rf, rf, *error_prob,
	   draws, init_f2, emit_f2, step_f2);
}

void est_map_f2(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose)
{
  est_map(*n_ind, *n_mar, 4, geno, rf, rf, *error_prob, 
	  init_f2b, emit_f2b, step_f2b, nrec_f2b, nrec_f2b,
	  loglik, *maxit, *tol, 0, *verbose);
}


void argmax_geno_f2(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax)
{		    
  argmax_geno(*n_ind, *n_pos, 3, geno, rf, rf, *error_prob,
	      argmax, init_f2, emit_f2, step_f2);
}


double errorlod_f2(int obs, double *prob, double error_prob)
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


void calc_errorlod_f2(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 3, geno, *error_prob, genoprob,
		errlod, errorlod_f2);
}



double nrec2_f2(int obs1, int obs2, double rf, int *cross_scheme)
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


double logprec_f2(int obs1, int obs2, double rf, int *cross_scheme)
{
  switch(obs1) {
  case 1: 
    switch(obs2) {
    case 1: return(2.0*log(1.0-rf));
    case 2: return(M_LN2 + log(rf) + log(1.0-rf));
    case 3: return(2.0*log(rf));
    case 4: return(log(1.0-rf*rf));
    case 5: return(log(1.0-(1.0-rf)*(1.0-rf)));
    }
  case 2:
    switch(obs2) {
    case 1: case 3: return(log(rf*(1.0-rf)));
    case 2: return(log(rf*rf+(1.0-rf)*(1.0-rf)));
    case 4: case 5: return(log(1.0-rf*(1.0-rf)));
    }
  case 3:
    switch(obs2) {
    case 1: return(2.0*log(rf));
    case 2: return(M_LN2 + log(rf) + log(1.0-rf));
    case 3: return(2.0*log(1.0-rf));
    case 4: return(log(1.0-(1.0-rf)*(1.0-rf)));
    case 5: return(log(1.0-rf*rf));
    }
  case 4:
    switch(obs2) {
    case 1: return(log((1.0-rf*rf) / 3.0));
    case 2: return(log((1.0-rf*(1.0-rf)) / 3.0) + M_LN2);
    case 3: return(log((1.0-(1.0-rf)*(1.0-rf)) / 3.0));
    case 4: return(log(((1.0-rf)*(1.0-rf) + 2.0) / 3.0));
    case 5: return(log((rf*rf+2.0) / 3.0));
    }
  case 5:
    switch(obs2) {
    case 1: return(log((1.0-(1.0-rf)*(1.0-rf)) / 3.0));
    case 2: return(log((1.0-rf*(1.0-rf)) / 3.0) + M_LN2);
    case 3: return(log((1.0-rf*rf) / 3.0));
    case 4: return(log((rf*rf+2.0) / 3.0));
    case 5: return(log(((1.0-rf)*(1.0-rf) + 2.0) / 3.0));
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}


void est_rf_f2(int *n_ind, int *n_mar, int *geno, double *rf, 
	       int *maxit, double *tol)
{
  est_rf(*n_ind, *n_mar, geno, rf, nrec2_f2, logprec_f2, 
	 *maxit, *tol, 2);
}

void calc_pairprob_f2(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob) 
{
  calc_pairprob(*n_ind, *n_mar, 3, geno, rf, rf, *error_prob, genoprob,
		pairprob, init_f2, emit_f2, step_f2);
}

void marker_loglik_f2(int *n_ind, int *geno,
		      double *error_prob, double *loglik)
{
  marker_loglik(*n_ind, 4, geno, *error_prob, 
		init_f2b, emit_f2b, loglik);
}
  
/* end of hmm_f2.c */
