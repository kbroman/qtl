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

/* ref: Jiang and Zeng (1997 Genetics) */

void prob_bcsft(double rf, int s, int t, double *transpr);
void count_bcsft(double rf, int s, int t, double *transct);
void expect_bcsft(double rf, int s, int t, double *transexp);

void step_wrap(int *gen1, int *gen2, double *rf, int *cross_scheme, double *ret, double *transpr)
{
  prob_bcsft(*rf, cross_scheme[0], cross_scheme[1], transpr);

  ret[0] = step_bcsftb(*gen1, *gen2, *rf, *rf, cross_scheme);
  ret[1] = step_bcsft(*gen1, *gen2, *rf, *rf, cross_scheme);
  ret[2] = step_bc(*gen1, *gen2, *rf, *rf, cross_scheme);
  ret[3] = step_f2b(*gen1, *gen2, *rf, *rf, cross_scheme);
  ret[4] = step_f2(*gen1, *gen2, *rf, *rf, cross_scheme);
  return;
}
void init_wrap(int *gen, int *cross_scheme, double *ret)
{
  ret[0] = init_bcsftb(*gen, cross_scheme);
  ret[1] = init_f2b(*gen, cross_scheme);
  if(*gen < 4) {
    if(*gen < 4 || cross_scheme[1] > 0)
      ret[2] = init_bcsft(*gen,  cross_scheme);
    ret[3] = init_f2(*gen, cross_scheme);
    if(*gen < 3)
      ret[4] = init_bc(*gen,  cross_scheme);
  }
  return;
}
void nrec_wrap(int *gen1, int *gen2, double *rf, int *cross_scheme, double *ret)
{
  ret[0] = nrec_bcsftb(*gen1, *gen2, *rf, cross_scheme);
  ret[1] = nrec_f2b(*gen1, *gen2, *rf, cross_scheme);
  if(*gen1 < 3 && *gen2 < 3)
    ret[2] = nrec_bc(*gen1, *gen2, *rf, cross_scheme);
  return;
}
void rf_wrap(int *gen1, int *gen2, double *rf, int *cross_scheme, double *ret)
{
  ret[0] = nrec2_bcsft(*gen1, *gen2, *rf, cross_scheme);
  ret[1] = nrec2_f2(*gen1, *gen2, *rf, cross_scheme);
  ret[2] = logprec_bcsft(*gen1, *gen2, *rf, cross_scheme);
  ret[3] = logprec_f2(*gen1, *gen2, *rf, cross_scheme);
  return;
}
void bcsft_wrap(double *rf, int *cross_scheme, double *init, double *emit, double *step,
		double *stepb, double *nrec, double *transpr, double *transexp)
{
  int gen1,gen2;
  static double error_prob = 0.0001;

  prob_bcsft(*rf, cross_scheme[0], cross_scheme[1], transpr);
  expect_bcsft(*rf, cross_scheme[0], cross_scheme[1], transexp);

  for(gen1=0; gen1<4; gen1++) {
    if(gen1 < 3) {
      init[gen1] = init_bcsft(gen1+1, cross_scheme);
      init[gen1+3] = init_bc(gen1+1, cross_scheme);
    }
    for(gen2=0; gen2<3; gen2++) {
      if((gen1 < 3) && (gen2 < 3)) {
	emit[gen1 + 3 * gen2] = emit_bcsft(gen1+1, gen2+1, error_prob, cross_scheme);
	emit[gen1 + 3 * gen2 + 9] = emit_bc(gen1+1, gen2+1, error_prob, cross_scheme);
	step[gen1 + 3 * gen2] = step_bcsft(gen1+1, gen2+1, *rf, *rf, cross_scheme);
	step[gen1 + 3 * gen2 + 9] = step_bc(gen1+1, gen2+1, *rf, *rf, cross_scheme);
      }
      nrec[gen1 + 4 * gen2] = nrec_bcsftb(gen1+1, gen2+1, *rf, cross_scheme);
      nrec[gen1 + 4 * gen2 + 16] = nrec_bc(gen1+1, gen2+1, *rf, cross_scheme);
      stepb[gen1 + 4 * gen2] = step_bcsftb(gen1+1, gen2+1, *rf, *rf, cross_scheme);
      stepb[gen1 + 4 * gen2 + 16] = step_bc(gen1+1, gen2+1, *rf, *rf, cross_scheme);
    }
  }

  return;
}
  
/* end of test_bcsft.c */
