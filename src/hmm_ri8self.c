/**********************************************************************
 * 
 * hmm_ri8self.c
 * 
 * copyright (c) 2009-2010, Karl W Broman
 *
 * last modified Jul, 2010
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
 *           est_map_ri8self, nrec2_ri8self, logprec_ri8self, est_rf_ri8self,
 *           marker_loglik_ri8self, calc_pairprob_ri8self, 
 *           errorlod_ri8self, calc_errorlod_ri8self
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for 8-way RIL by selfing.
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

double init_ri8self(int true_gen, int *cross_scheme)
{
  return(-3.0*M_LN2); /* log(1/8) */
}

double emit_ri8self(int obs_gen, int true_gen, double error_prob, int *cross_scheme)
{
  if(obs_gen==0) return(0.0);
  if(obs_gen & (1 << (true_gen-1))) return(log(1.0-error_prob));
  else return(log(error_prob)); 
}
    
double step_ri8self(int gen1, int gen2, double rf, double junk, int *cross_scheme) 
{
  int temp;
  if(gen1 > gen2) { temp = gen1; gen1 = gen2; gen2=temp; }
  
  if(gen1 == gen2) 
    return(2.0*log(1.0-rf)-log(1.0+2.0*rf));
  else if((gen1==1 || gen1==3 || gen1==5 || gen1==7) && gen2==gen1+1)
    return(log(rf)+log(1.0-rf)-log(1.0+2.0*rf));
  else
    return(log(rf) - M_LN2 - log(1.0+2.0*rf));
}

/* this is needed for est.map; estimated recombination fractions on the RIL scale */
double step_special_ri8self(int gen1, int gen2, double rf, double junk, int *cross_scheme) 
{
  double RF;
  int temp;
  if(gen1 > gen2) { temp = gen1; gen1 = gen2; gen2=temp; }
  
  if(gen1 == gen2) 
    return(log(1.0-rf));
  else if((gen1==1 || gen1==3 || gen1==5 || gen1==7) && gen2==gen1+1) {
    RF = 2.0-rf-sqrt(rf*rf-5.0*rf+4.0);
    return(log(RF)+log(1.0-RF)-log(1.0+2.0*RF));
  }
  else {
    RF = 2.0-rf-sqrt(rf*rf-5.0*rf+4.0);
    return(log(RF) - M_LN2 - log(1.0+2.0*RF));
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

/* expected no. recombinants */
double nrec2_ri8self(int obs1, int obs2, double rf, int *cross_scheme)
{
  int n1, n2, n12, nr, and, i, nstr=8;
  int offby1;
  double rfcontract, rf0, rf1, rf2, num;
  

  if(obs1==0 || obs2==0) return(-999.0); /* this shouldn't happen */

  n1=n2=n12=0;
  and = obs1 & obs2;
  
  /* count bits */
  for(i=0; i<nstr; i++) {
    n1 += ((obs1 & 1<<i)!=0);
    n2 += ((obs2 & 1<<i)!=0);
    n12 += ((and  & 1<<i)!=0);
  }

  /* need to count the cases (1,2), (3,4), (5,6), (7,8) */
  offby1 = (((obs1 & 1) && (obs2 & 2))) +
           (((obs1 & 2) && (obs2 & 1))) +
           (((obs1 & 4) && (obs2 & 8))) +
           (((obs1 & 8) && (obs2 & 4))) +
           (((obs1 & 16) && (obs2 & 32))) +
           (((obs1 & 32) && (obs2 & 16))) +
           (((obs1 & 64) && (obs2 & 128))) +
           (((obs1 & 128) && (obs2 & 64)));

  nr = n1*n2-n12-offby1;

  /* rf on meiosis scale */
  rfcontract = 2.0 - rf - sqrt(rf*rf - 5.0*rf + 4.0);

  /* joint prob for non-rec'ts */
  rf0 = 1.0 - rf;
  /* joint prob for off-by-one rec'ts */
  rf1 = rfcontract * (1.0 - rfcontract) / (1.0 + 2.0 * rfcontract);
  /* joint prob for off-by-more-than-one rec'ts */
  rf2 = rfcontract / 2.0 / (1.0 + 2.0*rfcontract);

  num = (double)offby1 * rf1 + (double)nr * rf2;

  return( num / (num + (double)n12*rf0));
}

/* log [joint probability * 8] */
double logprec_ri8self(int obs1, int obs2, double rf, int *cross_scheme)
{
  int n1, n2, n12, nr, and, i, nstr=8;
  int offby1;
  double rfcontract, rf0, rf1, rf2;

  if(obs1==0 || obs2==0) return(-999.0); /* this shouldn't happen */

  n1=n2=n12=0;
  and = obs1 & obs2;
  
  /* count bits */
  for(i=0; i<nstr; i++) {
    n1 += ((obs1 & 1<<i)!=0);
    n2 += ((obs2 & 1<<i)!=0);
    n12 += ((and  & 1<<i)!=0);
  }

  /* need to count the cases (1,2), (3,4), (5,6), (7,8) */
  offby1 = (((obs1 & 1) && (obs2 & 2))) +
           (((obs1 & 2) && (obs2 & 1))) +
           (((obs1 & 4) && (obs2 & 8))) +
           (((obs1 & 8) && (obs2 & 4))) +
           (((obs1 & 16) && (obs2 & 32))) +
           (((obs1 & 32) && (obs2 & 16))) +
           (((obs1 & 64) && (obs2 & 128))) +
           (((obs1 & 128) && (obs2 & 64)));

  nr = n1*n2-n12-offby1;

  /* rf on meiosis scale */
  rfcontract = 2.0 - rf - sqrt(rf*rf - 5.0*rf + 4.0);

  /* joint prob for non-rec'ts */
  rf0 = 1.0 - rf;
  /* joint prob for off-by-one rec'ts */
  rf1 = rfcontract * (1.0 - rfcontract) / (1.0 + 2.0 * rfcontract);
  /* joint prob for off-by-more-than-one rec'ts */
  rf2 = rfcontract / 2.0 / (1.0 + 2.0*rfcontract);

  return( log( (double)offby1 * rf1 + (double)nr * rf2 + (double)n12*rf0 ) );
}

void est_rf_ri8self(int *n_ind, int *n_mar, int *geno, double *rf, 
		   int *maxit, double *tol)
{
  est_rf(*n_ind, *n_mar, geno, rf, nrec2_ri8self, logprec_ri8self, 
	 *maxit, *tol, 1);
}

void marker_loglik_ri8self(int *n_ind, int *geno,
			   double *error_prob, double *loglik)
{
  marker_loglik(*n_ind, 8, geno, *error_prob, init_ri8self, emit_ri8self,
		loglik);
}

void calc_pairprob_ri8self(int *n_ind, int *n_mar, int *geno, 
			   double *rf, double *error_prob, 
			   double *genoprob, double *pairprob) 
{
  calc_pairprob(*n_ind, *n_mar, 8, geno, rf, rf, *error_prob, genoprob,
		pairprob, init_ri8self, emit_ri8self, step_ri8self);
}

double errorlod_ri8self(int obs, double *prob, double error_prob)
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

void calc_errorlod_ri8self(int *n_ind, int *n_mar, int *geno, 
			   double *error_prob, double *genoprob, 
			   double *errlod)
{
  calc_errorlod(*n_ind, *n_mar, 8, geno, *error_prob, genoprob,
		errlod, errorlod_ri8self);
}

/* end of hmm_ri8self.c */
