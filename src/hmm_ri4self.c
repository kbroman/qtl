/**********************************************************************
 * 
 * hmm_ri4self.c
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
 * Contains: init_ri4self, emit_ri4self, step_ri4self, step_special_ri4self,
 *           calc_genoprob_ri4self, calc_genoprob_special_ri4self, 
 *           argmax_geno_ri4self, sim_geno_ri4self,
 *           est_map_ri4self,
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the Collaborative Cross
 *
 * Genotype codes:    1-4
 * "Phenotype" codes: 0=missing; otherwise binary 1-15, with bit i
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
#include "hmm_ri4self.h"
#include "hmm_bc.h"

double init_ri4self(int true_gen)
{
  return(LN_025);
}

double emit_ri4self(int obs_gen, int true_gen, double error_prob)
{
  if(obs_gen==0) return(0.0);
  if(obs_gen & (1 << (true_gen-1))) return(log(1.0-error_prob));
  else return(log(error_prob)); 
}
    
  
double step_ri4self(int gen1, int gen2, double rf, double junk) 
{
  if(gen1 == gen2) 
    return(log(1.0-rf)-log(1.0+2.0*rf));
  else 
    return(log(rf)-log(1.0+2.0*rf));
}


/* this is needed for est.map; estimated recombination fractions on the RIL scale */
double step_special_ri4self(int gen1, int gen2, double rf, double junk) 
{
  if(gen1 == gen2) 
    return(log(1.0-rf));
  else 
    return(log(rf) - log(3.0));
}


void calc_genoprob_ri4self(int *n_ind, int *n_mar, int *geno, 
			  double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 4, geno, rf, rf, *error_prob, genoprob,
		init_ri4self, emit_ri4self, step_ri4self);
}

  
void calc_genoprob_special_ri4self(int *n_ind, int *n_mar, int *geno, 
				   double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob_special(*n_ind, *n_mar, 4, geno, rf, rf, *error_prob, genoprob,
			init_ri4self, emit_ri4self, step_ri4self);
}

  
void argmax_geno_ri4self(int *n_ind, int *n_pos, int *geno,
			double *rf, double *error_prob, int *argmax)
{
  argmax_geno(*n_ind, *n_pos, 4, geno, rf, rf, *error_prob,
	      argmax, init_ri4self, emit_ri4self, step_ri4self);
}

void sim_geno_ri4self(int *n_ind, int *n_pos, int *n_draws, int *geno, 
		      double *rf, double *error_prob, int *draws) 
{
  sim_geno(*n_ind, *n_pos, 4, *n_draws, geno, rf, rf, *error_prob, 
	   draws, init_ri4self, emit_ri4self, step_ri4self);
}


/* for estimating map, must do things with recombination fractions on the RIL scale */
void est_map_ri4self(int *n_ind, int *n_mar, int *geno, double *rf, 
		     double *error_prob, double *loglik, int *maxit, 
		     double *tol, int *verbose)
{
  int i;

  /* expand rf */
  for(i=0; i< *n_mar-1; i++) rf[i] = 3.0*rf[i]/(1.0+2.0*rf[i]);

  est_map(*n_ind, *n_mar, 4, geno, rf, rf, *error_prob, 
	  init_ri4self, emit_ri4self, step_special_ri4self, nrec_bc, nrec_bc,
	  loglik, *maxit, *tol, 0, *verbose);

  /* contract rf */
  for(i=0; i< *n_mar-1; i++) rf[i] = rf[i]/(3.0-2.0*rf[i]);
}

/* end of hmm_ri4self.c */
