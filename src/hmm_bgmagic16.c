/**********************************************************************
 * 
 * hmm_bgmagic16.c
 * 
 * copyright (c) 2011-2012, Karl W Broman
 *
 * last modified Jul, 2012
 * first written Dec, 2011
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
 * Contains: init_bgmagic16, emit_bgmagic16, step_bgmagic16, 
 *           calc_genoprob_bgmagic16, calc_genoprob_special_bgmagic16,
 *           argmax_geno_bgmagic16, sim_geno_bgmagic16,
 *           est_map_bgmagic16, 
 *           marker_loglik_bgmagic16, calc_pairprob_bgmagic16, 
 *           errorlod_bgmagic16, calc_errorlod_bgmagic16
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for 8-way RIL by selfing.
 *
 * Genotype codes:    1-16
 * "Phenotype" codes: 0=missing; otherwise binary 1-65536, with bit i
 *                    indicating SNP compatible with parent i
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "hmm_bgmagic16.h"
#include "hmm_bc.h"
#include "util.h"

double init_bgmagic16(int true_gen, int *ignored)
{
  return(-4.0*M_LN2); /* log(1/16) */
}

double emit_bgmagic16(int obs_gen, int true_gen, double error_prob, int *ignored)
{
  if(obs_gen==0) return(0.0);
  if(obs_gen & (1 << (true_gen-1))) return(log(1.0-error_prob));
  else return(log(error_prob)); 
}
    
double step_bgmagic16(int gen1, int gen2, double rf, double junk, int *ignored) 
{
  int tempi;
  double p0, tempd;

  if(gen1 == gen2) {
    tempd = 1.0-rf;
    tempd = tempd*tempd;
    p0 = tempd*tempd;
  }
  else {
    if(gen1 > gen2) { /* order gen1 and gen2 */
      tempi = gen1; 
      gen1 = gen2; 
      gen2 = tempi;
    }
    if((gen1 == gen2 - 1) && (gen2 % 2 == 0)) { /* 1:2 case */
      p0 = rf*(1.0-rf)*(1.0-rf)*(1.0-rf);
    }
    else if((gen2 - gen1 <= 4) && ((gen2 % 4 == 3) || (gen2 % 4 == 0))) { /* 1:3 case */
      p0 = rf*(1.0-rf)*(1.0-rf)/2.0;
    }
    else if(gen2 <= 8 || gen1 > 8) { /* 1:5 case */
      p0 = rf*(1.0-rf)/4.0;
    }
    else { /* 1:9 case */
      p0 = rf/8.0;
    }
  
  }

  return( log( (1.0-rf)*(1.0-rf)*(1.0-rf) * (p0 - 1.0/16.0) + 1.0/16.0 ) );
}


void calc_genoprob_bgmagic16(int *n_ind, int *n_mar, int *geno, 
			   double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 16, geno, rf, rf, *error_prob, genoprob,
		init_bgmagic16, emit_bgmagic16, step_bgmagic16);
}

void calc_genoprob_special_bgmagic16(int *n_ind, int *n_mar, int *geno, 
				   double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob_special(*n_ind, *n_mar, 16, geno, rf, rf, *error_prob, genoprob,
			init_bgmagic16, emit_bgmagic16, step_bgmagic16);
}

void argmax_geno_bgmagic16(int *n_ind, int *n_pos, int *geno,
			 double *rf, double *error_prob, int *argmax)
{
  argmax_geno(*n_ind, *n_pos, 16, geno, rf, rf, *error_prob,
	      argmax, init_bgmagic16, emit_bgmagic16, step_bgmagic16);
}

void sim_geno_bgmagic16(int *n_ind, int *n_pos, int *n_draws, int *geno, 
		      double *rf, double *error_prob, int *draws) 
{
  sim_geno(*n_ind, *n_pos, 16, *n_draws, geno, rf, rf, *error_prob, 
	   draws, init_bgmagic16, emit_bgmagic16, step_bgmagic16);
}

void est_map_bgmagic16(int *n_ind, int *n_mar, int *geno, double *rf, 
		     double *error_prob, double *loglik, int *maxit, 
		     double *tol, int *verbose)
{
  est_map(*n_ind, *n_mar, 16, geno, rf, rf, *error_prob, 
	  init_bgmagic16, emit_bgmagic16, step_bgmagic16, nrec_bc, nrec_bc,
	  loglik, *maxit, *tol, 0, *verbose);
}



void marker_loglik_bgmagic16(int *n_ind, int *geno,
			   double *error_prob, double *loglik)
{
  marker_loglik(*n_ind, 16, geno, *error_prob, init_bgmagic16, emit_bgmagic16,
		loglik);
}

void calc_pairprob_bgmagic16(int *n_ind, int *n_mar, int *geno, 
			   double *rf, double *error_prob, 
			   double *genoprob, double *pairprob) 
{
  calc_pairprob(*n_ind, *n_mar, 16, geno, rf, rf, *error_prob, genoprob,
		pairprob, init_bgmagic16, emit_bgmagic16, step_bgmagic16);
}

double errorlod_bgmagic16(int obs, double *prob, double error_prob)
{
  double p=0.0, temp;
  int n=0, i;

  if(obs==0 || (obs == (1<<16) - 1)) return(0.0);
  for(i=0; i<16; i++) {
    if(obs & 1<<i) p += prob[i];
    else n++;
  }
  if(n==0 || n==16) return(0.0); /* shouldn't happen */
  
  p = (1.0-p)/p;
  temp = (double)n*error_prob/15.0;

  p *= (1.0 - temp)/temp;

  if(p < TOL) return(-12.0);
  else return(log10(p));
}

void calc_errorlod_bgmagic16(int *n_ind, int *n_mar, int *geno, 
			   double *error_prob, double *genoprob, 
			   double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 16, geno, *error_prob, genoprob,
		errlod, errorlod_bgmagic16);
}

/* end of hmm_bgmagic16.c */
