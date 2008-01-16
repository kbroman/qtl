/**********************************************************************
 * 
 * hmm_cc.c
 * 
 * copyright (c) 2005, Karl W Broman
 *
 * last modified Mar, 2005
 * first written Mar, 2005
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 *
 * Contains: init_cc, emit_cc, step_cc, 
 *           calc_genoprob_cc, argmax_geno_cc
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the Collaborative Cross
 *
 * Genotype codes:  1-8
 * Phenotype codes: 0=missing; otherwise binary 1-256 which bit i
 *                  indicating SNP compatible with parent i
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "hmm_cc.h"

double init_cc(int true_gen)
{
  return(LN_0125);
}

double emit_cc(int obs_gen, int true_gen, double error_prob)
{
  if(obs_gen==0) return(0.0);
  if(obs_gen & (1 << (true_gen-1))) return(log(1.0-error_prob));
  else return(log(error_prob)); 
}
    
  
double step_cc(int gen1, int gen2, double rf, double junk) 
{
  if(gen1 == gen2) 
    return(log(1.0-rf)-log(1.0+6.0*rf));
  else 
    return(log(rf)-log(1.0+6.0*rf));
}


void calc_genoprob_cc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob) 
{
  calc_genoprob(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
		init_cc, emit_cc, step_cc);
}

  
void argmax_geno_cc(int *n_ind, int *n_pos, int *geno,
		    double *rf, double *error_prob, int *argmax)
{
  argmax_geno(*n_ind, *n_pos, 8, geno, rf, rf, *error_prob,
	      argmax, init_cc, emit_cc, step_cc);
}

/* end of hmm_cc.c */
