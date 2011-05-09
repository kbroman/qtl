/**********************************************************************
 * 
 * hmm_bc.c
 * 
 * copyright (c) 2001-2010, Karl W Broman
 *
 * last modified Jul, 2010
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
 * Contains: init_bc, emit_bc, step_bc, nrec_bc, calc_genoprob_bc,
 *           calc_genoprob_special_bc, 
 *           sim_geno_bc, est_map_bc, argmax_geno_bc, errorlod_bc,
 *           calc_errorlod_bc, est_rf_bc, calc_pairprob_bc, 
 *           marker_loglik_bc
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the backcross.
 *
 * Genotype codes:  0=AA; 1=AB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "util.h"
#include "hmm_main.h"
#include "hmm_bc.h"

double init_bc(int true_gen, int *cross_scheme)
{
  return(-M_LN2); /* ln(0.5) */
}

double emit_bc(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  switch(obs_gen) {
  case 0: return(0.0);
  case 1: case 2:
    if(obs_gen==true_gen) return(log(1.0-error_prob));
    else return(log(error_prob));
  }
  return(0.0); /* shouldn't get here */
}

double step_bc(int gen1, int gen2, double rf, double junk, int *cross_scheme)
{
  if(gen1==gen2) return(log(1.0-rf));
  else return(log(rf));
  return(log(-1.0)); /* shouldn't get here */
}

double nrec_bc(int gen1, int gen2, double rf, int *cross_scheme) 
{
  if(gen1==gen2) return(0.0);
  else return(1.0);
}

void calc_genoprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 2, geno, rf, rf, *error_prob, genoprob,
		init_bc, emit_bc, step_bc);
}

void calc_genoprob_special_bc(int *n_ind, int *n_mar, int *geno, 
			      double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob_special(*n_ind, *n_mar, 2, geno, rf, rf, *error_prob, genoprob,
			init_bc, emit_bc, step_bc);
}

void sim_geno_bc(int *n_ind, int *n_pos, int *n_draws, int *geno,
		 double *rf, double *error_prob, int *draws)
{
  sim_geno(*n_ind, *n_pos, 2, *n_draws, geno, rf, rf, *error_prob,
	   draws, init_bc, emit_bc, step_bc);
}

void est_map_bc(int *n_ind, int *n_mar, int *geno, double *rf, 
		double *error_prob, double *loglik, int *maxit, 
		double *tol, int *verbose)
{
  est_map(*n_ind, *n_mar, 2, geno, rf, rf, *error_prob, 
	  init_bc, emit_bc, step_bc, nrec_bc, nrec_bc,
	  loglik, *maxit, *tol, 0, *verbose);
}

void argmax_geno_bc(int *n_ind, int *n_pos, int *geno, 
		   double *rf, double *error_prob, int *argmax)
{		    
  argmax_geno(*n_ind, *n_pos, 2, geno, rf, rf, *error_prob,
	      argmax, init_bc, emit_bc, step_bc);
}

double errorlod_bc(int obs, double *prob, double error_prob)
{
  double p=0.0;

  switch(obs) {
  case 0: return(0.0);
  case 1: p=prob[0]; break;
  case 2: p=prob[1]; break;
  }
  
  p = (1.0-p)/p*(1.0-error_prob)/error_prob;
  if(p < TOL) return(-12.0);
  else return(log10(p));
}

void calc_errorlod_bc(int *n_ind, int *n_mar, int *geno, 
		      double *error_prob, double *genoprob, 
		      double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 2, geno, *error_prob, genoprob,
		errlod, errorlod_bc);
}

void est_rf_bc(int *n_ind, int *n_mar, int *geno, double *rf) 
{
  int i, j1, j2, **Geno, a, b;
  double **Rf;

  /* reorganize geno and rf */
  reorg_geno(*n_ind, *n_mar, geno, &Geno);
  reorg_errlod(*n_mar, *n_mar, rf, &Rf);

  for(j1=0; j1< *n_mar; j1++) {

    /* count meioses */
    a = 0;
    for(i=0; i < *n_ind; i++) {
      if(Geno[j1][i] != 0) a++;
    }
    Rf[j1][j1] = (double) a;

    for(j2=j1+1; j2< *n_mar; j2++) {
      a=b=0;
      for(i=0; i< *n_ind; i++) {
	if(Geno[j1][i] != 0 && Geno[j2][i] != 0) {
	  a++;
	  if(Geno[j1][i] != Geno[j2][i]) b++;
	}
      }
      
      if(a != 0) { /* at least one informative meiosis */ 
	/*	if(b > a/2) b = a/2; */

	Rf[j1][j2] = (double)b/(double)a;
	if(b==0) /* no recombinations */
	  Rf[j2][j1] = (double)a*log10(1.0-Rf[j1][j2]);
	else 
	  Rf[j2][j1] = (double)b*log10(Rf[j1][j2]) +
	    (double)(a-b)*log10(1.0-Rf[j1][j2]);
	
	Rf[j2][j1] += (double)a*log10(2.0);

      }
      else {
	Rf[j1][j2] = NA_REAL;
	Rf[j2][j1] = 0.0;
      }

    } /* end loops over markers */
  }
}
 
void calc_pairprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob,
		      double *pairprob) 
{
  calc_pairprob(*n_ind, *n_mar, 2, geno, rf, rf, *error_prob, genoprob,
		pairprob, init_bc, emit_bc, step_bc);
}


void marker_loglik_bc(int *n_ind, int *geno,
		      double *error_prob, double *loglik)
{
  marker_loglik(*n_ind, 2, geno, *error_prob, init_bc, emit_bc,
		loglik);
}

/* end of hmm_bc.c */
