/**********************************************************************
 *
 * hmm_ri8selfIRIP.c
 *
 * Rohan Shah
 *
 * October, 2014
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
 * Contains: init_ri8IRIP1, emit_ri8IRIP1, step_ri8IRIP1, step_special_ri8IRIP1,
 *           calc_genoprob_ri8IRIP1, calc_genoprob_special_ri8IRIP1,
 *           argmax_geno_ri8IRIP1, sim_geno_ri8IRIP1,
 *           est_map_ri8IRIP1, nrec2_ri8IRIP1, logprec_ri8IRIP1, est_rf_ri8IRIP1,
 *           marker_loglik_ri8IRIP1, calc_pairprob_ri8IRIP1,
 *           errorlod_ri8IRIP1, calc_errorlod_ri8IRIP1
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for 8-way RIL by selfing, with 1 generation of .
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
#include "hmm_ri8selfIRIP1.h"
#include "hmm_bc.h"
#include "util.h"

#define M_LN56 4.0253516907351492333570491078177094338635851326626269
#define M_LN7 1.9459101490553133051053527434431797296370847295818611
double init_ri8selfIRIP1(int true_gen, int *cross_scheme)
{
  return(-3.0*M_LN2); /* log(1/8) */
}

double emit_ri8selfIRIP1(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  if(obs_gen==0) return(0.0);
  if(obs_gen & (1 << (true_gen-1))) return(log(1.0-error_prob));
  else return(log(error_prob));
}

double step_ri8selfIRIP1(int gen1, int gen2, double rf, double junk, int *cross_scheme)
{
  //See Teuscher and Broman 2007, "Haplotype Probabilities for Multiple-Strain Recombinant Inbred Lines", Equation (3), s=1
  double equalProb = (1-rf)*(1-rf)*(1-rf)+(2*rf)/8;

  if(gen1 == gen2)
    return(log(equalProb)-log(1+2*rf));
  else
    return(log(1-equalProb/(1+2*rf)) - M_LN7);
}

void calc_genoprob_ri8selfIRIP1(int *n_ind, int *n_mar, int *geno,
               double *rf, double *error_prob, double *genoprob)
{
  calc_genoprob(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
            init_ri8selfIRIP1, emit_ri8selfIRIP1, step_ri8selfIRIP1);
}

void calc_genoprob_special_ri8selfIRIP1(int *n_ind, int *n_mar, int *geno,
                   double *rf, double *error_prob, double *genoprob)
{
  calc_genoprob_special(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
            init_ri8selfIRIP1, emit_ri8selfIRIP1, step_ri8selfIRIP1);
}

void argmax_geno_ri8selfIRIP1(int *n_ind, int *n_pos, int *geno,
             double *rf, double *error_prob, int *argmax)
{
  argmax_geno(*n_ind, *n_pos, 8, geno, rf, rf, *error_prob,
          argmax, init_ri8selfIRIP1, emit_ri8selfIRIP1, step_ri8selfIRIP1);
}

void sim_geno_ri8IRIP1(int *n_ind, int *n_pos, int *n_draws, int *geno,
              double *rf, double *error_prob, int *draws)
{
  sim_geno(*n_ind, *n_pos, 8, *n_draws, geno, rf, rf, *error_prob,
       draws, init_ri8selfIRIP1, emit_ri8selfIRIP1, step_ri8selfIRIP1);
}

/* for estimating map, must do things with recombination fractions on the RIL scale */
void est_map_ri8IRIP1(int *n_ind, int *n_mar, int *geno, double *rf,
             double *error_prob, double *loglik, int *maxit,
             double *tol, int *verbose)
{
  est_map(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob,
      init_ri8selfIRIP1, emit_ri8selfIRIP1, step_ri8selfIRIP1, nrec_bc, nrec_bc,
      loglik, *maxit, *tol, 0, *verbose);
}

/* expected no. recombinants */
double nrec2_ri8selfIRIP1(int obs1, int obs2, double rf, int *cross_scheme)
{
  int n1, n2, n12, nr, and, i, nstr=8;
  double rf0, rf1, num;


  if(obs1==0 || obs2==0) return(-999.0); /* this shouldn't happen */

  n1=n2=n12=0;
  and = obs1 & obs2;

  /* count bits */
  for(i=0; i<nstr; i++) {
    n1 += ((obs1 & 1<<i)!=0);
    n2 += ((obs2 & 1<<i)!=0);
    n12 += ((and  & 1<<i)!=0);
  }

  nr = n1*n2-n12;

  /* joint prob for non-rec'ts */
  rf0 = ((1-rf)*(1-rf)*(1-rf)+(2*rf+1)/8)/ (8*(1+2 * rf));
  /* joint prob for recombinants */
  rf1 = (1 - 8 *rf0) / (8*7);

  num = (double)nr * rf1;

  return( num / (num + (double)n12*rf0));
}

/* log [joint probability * 8] */
double logprec_ri8selfIRIP1(int obs1, int obs2, double rf, int *cross_scheme)
{
  int n1, n2, n12, nr, and, i, nstr=8;
  double rf0, rf1;


  if(obs1==0 || obs2==0) return(-999.0); /* this shouldn't happen */

  n1=n2=n12=0;
  and = obs1 & obs2;

  /* count bits */
  for(i=0; i<nstr; i++) {
    n1 += ((obs1 & 1<<i)!=0);
    n2 += ((obs2 & 1<<i)!=0);
    n12 += ((and  & 1<<i)!=0);
  }

  nr = n1*n2-n12;

  /* joint prob for non-rec'ts */
  rf0 = ((1-rf)*(1-rf)*(1-rf)+(2*rf+1)/8)/ (8*(1+2 * rf));
  /* joint prob for recombinants */
  rf1 = (1 - 8 * rf0) / (8*7);

  return(log(8.0*((double)nr * rf1 + (double)n12 * rf0)));
}

void est_rf_ri8selfIRIP1(int *n_ind, int *n_mar, int *geno, double *rf,
           int *maxit, double *tol)
{
  est_rf(*n_ind, *n_mar, geno, rf, nrec2_ri8selfIRIP1, logprec_ri8selfIRIP1,
     *maxit, *tol, 1);
}

void marker_loglik_ri8selfIRIP1(int *n_ind, int *geno,
               double *error_prob, double *loglik)
{
  marker_loglik(*n_ind, 8, geno, *error_prob, init_ri8selfIRIP1, emit_ri8selfIRIP1,
        loglik);
}

void calc_pairprob_ri8selfIRIP1(int *n_ind, int *n_mar, int *geno,
               double *rf, double *error_prob,
               double *genoprob, double *pairprob)
{
  calc_pairprob(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
        pairprob, init_ri8selfIRIP1, emit_ri8selfIRIP1, step_ri8selfIRIP1);
}

double errorlod_ri8selfIRIP1(int obs, double *prob, double error_prob)
{
  double p=0.0, temp;
  int n=0, i;

  if(obs==0 || (obs == (1<<8) - 1)) return(0.0);
  for(i=0; i<8; i++) {
    if(obs & 1<<i) p += prob[i];
    else n++;
  }
  if(n==0 || n==8) return(0.0); /* shouldn't happen */

  p = (1.0-p)/p;
  temp = (double)n*error_prob/7.0;

  p *= (1.0 - temp)/temp;

  if(p < TOL) return(-12.0);
  else return(log10(p));
}

void calc_errorlod_ri8selfIRIP1(int *n_ind, int *n_mar, int *geno,
               double *error_prob, double *genoprob,
               double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 8, geno, *error_prob, genoprob,
        errlod, errorlod_ri8selfIRIP1);
}

/* end of hmm_ri8selfIRIP1.c */
