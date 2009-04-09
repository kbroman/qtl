/**********************************************************************
 * 
 * vbscan.c
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written May, 2001
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
 * 
 * C function for performing QTL mapping with the model
 *     p_g = Pr(pheno unobserved | genotype g)
 *     y | pheno observed, g ~ N(mu_g, sigma) 
 *
 * Contains: vbscan, R_vbscan
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "R.h"
#include "Rmath.h"
#include <R_ext/Utils.h>
#include "vbscan.h"
#include "util.h"

void vbscan(int n_pos, int n_ind, int n_gen, double *genoprob, double *pheno,
	    int *survived, double *lod, int maxit, double tol)
{
  double ***Genoprob, **w;
  int i, j, k, p, flag, ncol;
  double *mu, sig=0.0, *pi, *mu_prev, sig_prev, *pi_prev;
  double mu0, sig0, pi0, q0, loglik0;
  double s1, s2, s3;

  ncol = 3;

  /* set up matrix for the genoprob */
  reorg_genoprob(n_ind, n_pos, n_gen, genoprob, &Genoprob);

  /* allocate space for things */
  allocate_dmatrix(n_gen, n_ind, &w);
  allocate_double(n_gen, &mu);
  allocate_double(n_gen, &pi);
  allocate_double(n_gen, &mu_prev);
  allocate_double(n_gen, &pi_prev);
  
  /* estimates under null model */
  pi0 = mu0 = sig0 = 0.0;
  for(i=0; i<n_ind; i++) {
    if(survived[i]) pi0 += 1.0;
    else { 
      mu0 += pheno[i];
      sig0 += (pheno[i]*pheno[i]);
    }
  }
  q0 = (double)n_ind - pi0;
  mu0 /= q0;
  sig0 = sqrt((sig0-mu0*mu0*q0)/(q0-1.0));
  pi0 /= (double)(n_ind);

  loglik0 = 0.0;
  for(i=0; i<n_ind; i++) 
    if(survived[i]) loglik0 += log(pi0);
    else loglik0 += (log(1.0-pi0)+dnorm(pheno[i],mu0,sig0,1));

  /* begin genome scan */
  for(p=0; p<n_pos; p++) { /* loop over positions */

    /* initial values for EM */
    sig_prev = 0.0;
    for(j=0; j<n_gen; j++) {
      mu_prev[j] = pi_prev[j] = s1=s2=s3=0.0;
      for(i=0; i<n_ind; i++) {
	if(survived[i]) pi_prev[j] += Genoprob[j][p][i];
	else {
	  mu_prev[j] += pheno[i]*Genoprob[j][p][i];
	  s2 += Genoprob[j][p][i];
	  s3 += pheno[i]*pheno[i]*Genoprob[j][p][i];
	}
	s1 += Genoprob[j][p][i];
      }
      pi_prev[j] /= s1;
      mu_prev[j] /= s2;
      sig_prev += (s3 - mu_prev[j]*mu_prev[j]*s2);
    }
    sig_prev = sqrt(sig_prev/q0);
      
    for(k=0; k<maxit; k++) { /* loop over iterations */

      R_CheckUserInterrupt(); /* check for ^C */

      /* E step */
      for(i=0; i<n_ind; i++) {
	s1=0.0;
	if(survived[i])
	  for(j=0; j<n_gen; j++) 
	    s1 += (w[j][i] = Genoprob[j][p][i]*pi_prev[j]);
	else 
	  for(j=0; j<n_gen; j++) 
	    s1 += (w[j][i] = Genoprob[j][p][i]*(1-pi_prev[j])*
		   dnorm(pheno[i],mu_prev[j],sig_prev,0));
	for(j=0; j<n_gen; j++) 
	  w[j][i] /= s1;
      }

      /* M step */
      sig = 0.0;
      for(j=0; j<n_gen; j++) {
	mu[j] = pi[j] = s1=s2=s3=0.0;
	for(i=0; i<n_ind; i++) {
	  if(survived[i]) pi[j] += w[j][i];
	  else {
	    mu[j] += pheno[i]*w[j][i];
	    s2 += w[j][i];
	    s3 += pheno[i]*pheno[i]*w[j][i];
	  }
	  s1 += w[j][i];
	}
	pi[j] /= s1;
	mu[j] /= s2;
	sig += (s3 - mu[j]*mu[j]*s2);
      }
      sig = sqrt(sig/q0);

      flag = 0;
      if(fabs(sig-sig_prev)>tol*(sig_prev+tol*100.0)) flag=1;
      for(j=0; j<n_gen; j++) {
	if(fabs(mu[j]-mu_prev[j])>tol*(fabs(mu_prev[j])+tol*100.0)) { flag=1; break; }
	if(fabs(pi[j]-pi_prev[j])>tol*(pi_prev[j]+tol*100.0)) { flag=1; break; }
      }
      if(flag==0) break;
      
      sig_prev = sig;
      for(j=0; j<n_gen; j++) {
	mu_prev[j] = mu[j];
	pi_prev[j] = pi[j];
      }
    } /* iterations */

    /* calculate likelihood */
    lod[p*ncol] = 0.0;
    for(i=0; i<n_ind; i++) {
      s1=0.0;
      if(survived[i]) 
	for(j=0; j<n_gen; j++) 
	  s1 += pi[j]*Genoprob[j][p][i];
      else 
	for(j=0; j<n_gen; j++) 
	  s1 += (1-pi[j])*Genoprob[j][p][i]*dnorm(pheno[i],mu[j],sig,0);
      lod[p*ncol] += log(s1);
    }

    /* repeat with pi's constant */ 
    for(j=0; j<n_gen; j++) mu_prev[j] = mu[j];
    sig_prev = sig;
    for(k=0; k<maxit; k++) { /* loop over iterations */
      R_CheckUserInterrupt(); /* check for ^C */

      /* E step */
      for(i=0; i<n_ind; i++) {
	s1=0.0;
	if(!survived[i]) {
	  for(j=0; j<n_gen; j++) 
	    s1 += (w[j][i] = Genoprob[j][p][i]*(1-pi0)*
		   dnorm(pheno[i],mu_prev[j],sig_prev,0));
	  for(j=0; j<n_gen; j++) 
	    w[j][i] /= s1;
	}
      }
      
      /* M step */
      sig = 0.0;
      for(j=0; j<n_gen; j++) {
	mu[j] = s2=s3=0.0;
	for(i=0; i<n_ind; i++) {
	  if(!survived[i]) {
	    mu[j] += pheno[i]*w[j][i];
	    s2 += w[j][i];
	    s3 += pheno[i]*pheno[i]*w[j][i];
	  }
	}
	mu[j] /= s2;
	sig += (s3 - mu[j]*mu[j]*s2);
      }
      sig = sqrt(sig/q0);

      flag = 0;
      if(fabs(sig-sig_prev)>tol*(sig_prev+tol*100.0)) flag=1;
      for(j=0; j<n_gen; j++) 
	if(fabs(mu[j]-mu_prev[j])>tol*(fabs(mu_prev[j])+tol*100.0)) { flag=1; break; }
      if(flag==0) break;
      
      sig_prev = sig;
      for(j=0; j<n_gen; j++)
	mu_prev[j] = mu[j];
    } /* iterations */

    /* calculate likelihood */
    lod[p*ncol+1] = 0.0;
    for(i=0; i<n_ind; i++) {
      s1=0.0;
      if(survived[i]) 
	for(j=0; j<n_gen; j++) 
	  s1 += pi0*Genoprob[j][p][i];
      else 
	for(j=0; j<n_gen; j++) 
	  s1 += (1-pi0)*Genoprob[j][p][i]*dnorm(pheno[i],mu[j],sig,0);
      lod[p*ncol+1] += log(s1);
    }

    /* repeat with mu's constant */
    for(k=0; k<maxit; k++) { /* loop over iterations */

      R_CheckUserInterrupt(); /* check for ^C */

      for(j=0; j<n_gen; j++) pi_prev[j] = pi[j];

      /* E step */
      for(i=0; i<n_ind; i++) {
	s1=0.0;
	if(survived[i]) 
	  for(j=0; j<n_gen; j++) 
	    s1 += (w[j][i] = Genoprob[j][p][i]*pi_prev[j]);
	else 
	  for(j=0; j<n_gen; j++) 
	    s1 += (w[j][i] = Genoprob[j][p][i]*(1-pi_prev[j]));

	for(j=0; j<n_gen; j++) 
	  w[j][i] /= s1;
      }
      
      /* M step */
      for(j=0; j<n_gen; j++) {
	pi[j] = 0.0;
	s1=0.0;
	for(i=0; i<n_ind; i++) {
	  if(survived[i]) pi[j] += w[j][i];
	  s1 += w[j][i];
	}
	pi[j] /= s1;
      }

      flag = 0;
      for(j=0; j<n_gen; j++) 
	if(fabs(pi[j]-pi_prev[j])>tol*(pi_prev[j]+tol*100.0)) { flag=1; break; }
      if(flag==0) break;
      
      for(j=0; j<n_gen; j++) 
	pi_prev[j] = pi[j];

    } /* iterations */

    /* calculate likelihood */
    lod[p*ncol+2] = 0.0;
    for(i=0; i<n_ind; i++) {
      s1=0.0;
      if(survived[i]) 
	for(j=0; j<n_gen; j++) 
	  s1 += pi[j]*Genoprob[j][p][i];
      else 
	for(j=0; j<n_gen; j++) 
	  s1 += (1-pi[j])*Genoprob[j][p][i]*dnorm(pheno[i],mu0,sig0,0);
      lod[p*ncol+2] += log(s1);
    }

    /* convert likelihoods to LOD scores */
    lod[p*ncol+2] = (lod[p*ncol]-lod[p*ncol+2])/log(10.0);
    lod[p*ncol+1] = (lod[p*ncol]-lod[p*ncol+1])/log(10.0);
    lod[p*ncol]   = (lod[p*ncol]-loglik0)/log(10.0);

  } /* end loop over positions */

} 
  

void R_vbscan(int *n_pos, int *n_ind, int *n_gen, double *genoprob, 
              double *pheno, 
	      int *survived, double *lod, int *maxit, double *tol)
{
  vbscan(*n_pos, *n_ind, *n_gen, genoprob, pheno, survived, 
         lod, *maxit, *tol);
}


/* end of vbscan.c */
