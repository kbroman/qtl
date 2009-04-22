/**********************************************************************
 * 
 * hmm_ri8self.c
 * 
 * copyright (c) 2009, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Apr, 2009
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
 * Contains: init_ri8self, emit_ri8self, step_ri8self, step_special_ri8self,
 *           calc_genoprob_ri8self, calc_genoprob_special_ri8self,
 *           argmax_geno_ri8self, sim_geno_ri8self,
 *           est_map_ri8self,
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the Collaborative Cross
 *
 * Genotype codes:    1-8
 * "Phenotype" codes: 0=missing; otherwise binary 1-255, with bit i
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
#include "hmm_ri8self.h"
#include "hmm_bc.h"

double init_ri8self(int true_gen)
{
  return(LN_0125);
}

double emit_ri8self(int obs_gen, int true_gen, double error_prob)
{
  if(obs_gen==0) return(0.0);
  if(obs_gen & (1 << (true_gen-1))) return(log(1.0-error_prob));
  else return(log(error_prob)); 
}
    
double step_ri8self(int gen1, int gen2, double rf, double junk) 
{
  if(gen1 == gen2) 
    return(2.0*log(1.0-rf)-log(1.0+2.0*rf));
  else if(gen1 == gen2 - 1 || gen2 == gen1 - 1) 
    return(log(rf)+log(1.0-rf)-log(1.0+2.0*rf));
  else
    return(log(rf) - LN_2 - log(1.0+2.0*rf));
}

/* this is need for est.map; estimated recombination fractions on the RIL scale */
double step_special_ri8self(int gen1, int gen2, double rf, double junk) 
{
  double RF;

  if(gen1 == gen2) 
    return(log(1.0-rf));
  else if(gen1 == gen2 - 1 || gen2 == gen1 - 1) {
    RF = 2.0-rf-sqrt(rf*rf-5.0*rf+4.0);
    return(log(RF)+log(1.0-RF)-log(1.0+2.0*RF));
  }
  else {
    RF = 2.0-rf-sqrt(rf*rf-5.0*rf+4.0);
    return(log(RF) - LN_2 - log(1.0+2.0*RF));
  }
}

void calc_genoprob_ri8self(int *n_ind, int *n_mar, int *geno, 
			   double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
		init_ri8self, emit_ri8self, step_ri8self);
}

void calc_genoprob_special_ri8self(int *n_ind, int *n_mar, int *geno, 
				   double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob_special(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
			init_ri8self, emit_ri8self, step_ri8self);
}

void argmax_geno_ri8self(int *n_ind, int *n_pos, int *geno,
			 double *rf, double *error_prob, int *argmax)
{
  argmax_geno(*n_ind, *n_pos, 8, geno, rf, rf, *error_prob,
	      argmax, init_ri8self, emit_ri8self, step_ri8self);
}

void sim_geno_ri8self(int *n_ind, int *n_pos, int *n_draws, int *geno, 
		      double *rf, double *error_prob, int *draws) 
{
  sim_geno(*n_ind, *n_pos, 8, *n_draws, geno, rf, rf, *error_prob, 
	   draws, init_ri8self, emit_ri8self, step_ri8self);
}

/* for estimating map, must do things with recombination fractions on the RIL scale */
void est_map_ri8self(int *n_ind, int *n_mar, int *geno, double *rf, 
		     double *error_prob, double *loglik, int *maxit, 
		     double *tol, int *verbose)
{
  int i;

  /* expand rf */
  for(i=0; i< *n_mar-1; i++) rf[i] = rf[i]*(1.0-rf[i])/(1.0+2.0*rf[i]);

  est_map(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, 
	  init_ri8self, emit_ri8self, step_special_ri8self, nrec_bc, nrec_bc,
	  loglik, *maxit, *tol, 0, *verbose);

  /* contract rf */
  for(i=0; i< *n_mar-1; i++) rf[i] = 2.0-rf[i]-sqrt(rf[i]*rf[i]-5.0*rf[i]+4.0);
}

/* end of hmm_ri8self.c */
