/**********************************************************************
 * 
 * effectscan.c
 *
 * copyright (c) 2007-8, Karl W Broman 
 *
 * last modified Jan, 2008
 * first written Sep, 2007
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
 * These functions are for calculating the estimated effects, by multiple
 * imputation, in a single-QTL scan along a chromosome.
 *
 * Contains: R_effectscan, effectscan
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include "effectscan.h"
#include "util.h"

/* R_effectscan: wrapper for effectscan */

void R_effectscan(int *nind, int *ngen, int *ndraws, int *npos,
		  int *draws, double *pheno, double *effectmapping,
		  double *beta, double *se, int *getse)
{
  int ***Draws;
  double **Beta, **SE;

  reorg_errlod(*ngen, *npos, beta, &Beta);
  reorg_errlod(*ngen, *npos, se, &SE);
  reorg_draws(*nind, *npos, *ndraws, draws, &Draws);

  effectscan(*nind, *ngen, *ndraws, *npos, Draws, pheno, effectmapping,
	     Beta, SE, *getse);
}
  

/**********************************************************************
 * effectscan
 *
 * nind   Number of individuals
 * ngen   Number of genotypes
 * ndraws Number of imputations
 * npos   Number of positions
 * Draws  The imputed genotypes (dim nind x npos x ndraws)
 * pheno  Phenotypes (length nind)
 * effectmapping  Matrix of size ngen x ngen, giving the design matrix 
 *                for each possible genotype
 * Beta   On exit, the estimated coefficients (dim npos x ngen)
 * SE     On exit, the estimated standard errors (dim npos x ngen)
 * getse  If 1, calculate SEs; if 0, don't
 *
 **********************************************************************/

void effectscan(int nind, int ngen, int ndraws, int npos,
		int ***Draws, double *pheno, double *mapping,
		double **Beta, double **SE, int getse)
{
  int nphe=1, lwork, info, i, j, s, k, ntrim, *index, *ng, *flag;
  double *resid, *dwork, *var, sigmasq=0.0, *x;
  double *wbeta, *wvar, *wrss, *wts, totwt=0.0;
  double lrss0, mpheno;
  
  lwork = 4*nind;
  ntrim = (int) floor( 0.5*log(ndraws)/log(2.0) );

  allocate_double(nind, &resid);
  allocate_double(lwork, &dwork);
  allocate_double(ngen*ndraws, &var);
  allocate_double(ngen*ndraws, &wbeta);
  allocate_double(ngen*ndraws, &wvar);
  allocate_double(ndraws, &wrss);
  allocate_double(ndraws, &wts);
  allocate_double(ngen*nind, &x);
  allocate_int(ndraws, &index);
  allocate_int(ngen, &ng);
  allocate_int(ndraws, &flag);
  
  /* log null RSS */
  lrss0 = mpheno = 0.0;
  for(i=0; i<nind; i++) mpheno += pheno[i];
  mpheno /= (double)nind;
  for(i=0; i<nind; i++) lrss0 += (pheno[i] - mpheno)*(pheno[i] - mpheno);
  lrss0 = log(lrss0);

  for(s=0; s<npos; s++) { /* loop over positions */

    for(i=0; i<ndraws; i++) { /* loop over imputations */

      /* form X matrix */
      for(k=0; k<ngen; k++) 
	for(j=0; j<nind; j++) 
	  x[j+k*nind] = mapping[Draws[i][s][j]+k*ngen];
	  
      memcpy(resid, pheno, nind*sizeof(double));


      /* look for 0's in the genotype matrix */
      flag[i] = 0;
      for(j=0; j<ngen; j++) ng[j] = 0;
      for(j=0; j<nind; j++) ng[Draws[i][s][j]]++;
      for(j=0; j<ngen; j++) if(ng[j]==0) flag[i] = 1;

      if(!flag[i]) {
	/* linear regression */
	F77_CALL(dgels)("N", &nind, &ngen, &nphe, x, &nind, resid, &nind,
			dwork, &lwork, &info);

        /* coefficient estimates */
        for(j=0; j<ngen; j++) wbeta[j+i*ngen] = resid[j];

        /* residual sum of squares and residual variance */
        wrss[i] = 0.0;
        for(j=ngen; j<nind; j++) wrss[i] += (resid[j]*resid[j]);
        if(getse) sigmasq = wrss[i] / (double)(nind-ngen);

        wts[i] = ((double)nind / 2.0) * (lrss0 - log(wrss[i]));
      }

      if(i==0) totwt = wts[i];
      else totwt = addlog(totwt, wts[i]);
    
      if(getse && !flag[i]) {
	for(j=0; j<ngen; j++)
	  memcpy(var+j*ngen, x+j*nind, ngen*sizeof(double));
  
	/* (X'X)^-1 */
	F77_CALL(dpotri)("U", &ngen, var, &ngen, &info);

	/* estimated variances of estimated coefficients */
	for(j=0; j<ngen; j++) wvar[j + i*ngen] = sigmasq * var[j+j*ngen];
      }

    } /* end loop over imputations */


    /* trim weights */
    for(i=0; i<ndraws; i++) index[i] = i;
    rsort_with_index(wrss, index, ndraws);
    for(i=0; i<ndraws; i++) wts[i] = exp(wts[i] - totwt);
    
    for(i=0; i<ntrim; i++) 
      wts[index[i]] = wts[index[ndraws-i-1]] = 0.0;

    /* now calculated overall estimates, 
       with var(beta) = E[ var(beta|imp) ] + var[ E(beta|imp) ] */
    for(j=0; j<ngen; j++) {
      Beta[s][j] = 0.0;
      if(getse) {
	SE[s][j] = 0.0;
	var[j] = 0.0;
      }
    }
    
    totwt = 0.0;
    for(i=0; i<ndraws; i++) {
      if(!flag[i]) {
	totwt += wts[i];

	for(j=0; j<ngen; j++) {
	  Beta[s][j] += wts[i]*wbeta[j+i*ngen];
      
	  if(getse) SE[s][j] += wts[i]*wvar[j+i*ngen];
	}
      }
    }

    for(j=0; j<ngen; j++) {
      Beta[s][j] /= totwt;
      if(getse) SE[s][j] /= totwt;
    }

    if(getse) {
      for(j=0; j<ngen; j++) {
	for(i=0; i<ndraws; i++) 
	  if(!flag[i]) 
	    var[j] += wts[i] * (wbeta[j+i*ngen] - Beta[s][j])*(wbeta[j+i*ngen] - Beta[s][j]);

	SE[s][j] += var[j] / ( (double)(ndraws-ntrim-1) * totwt / (double)(ndraws-ntrim) );

	SE[s][j] = sqrt(SE[s][j]);
      }
    }

  } /* end loop over positions */
}
 
/* end of effectscan.c */
