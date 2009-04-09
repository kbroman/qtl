/**********************************************************************
 * 
 * forwsel.c
 * 
 * copyright (c) 2007-8, Karl W Broman
 *
 * last modified Jan, 2008
 * first written Jan, 2007
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
 * This is a simple routine to do forward selection in regression
 * to a fixed number of covariates
 *
 * Contains: R_markerforwsel, markerforwsel, 
 *           R_markerforwself2, markerforwself2
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include "forwsel.h"


/* R wrappers */
void R_markerforwsel(int *n, int *m, double *x, double *y,
		     int *maxsize, int *chosen, double *rss)
{
  double **X;
  int i;

  /* reorganize x matrix */
  X = (double **)R_alloc(*m, sizeof(double *));
  X[0] = x;
  for(i=1; i< *m; i++) X[i] = X[i-1] + *n;

  markerforwsel(*n, *m, X, y, *maxsize, chosen, rss);
}

void R_markerforwself2(int *n, int *m, int *x, double *y,
		       int *maxsize, int *chosen, double *rss)
{
  double **X;
  int i, j;

  /* Create X matrix with 2 columns for each marker */
  X = (double **)R_alloc(*m * 2, sizeof(double *));
  X[0] = (double *)R_alloc(*m * 2 * *n, sizeof(double));
  for(i=1; i< *m * 2; i++) X[i] = X[i-1] + *n;

  /* fill up X matrix */
  for(i=0; i< *m; i++) {
    for(j=0; j< *n; j++) {
      if(x[i * *n + j] == 1) {
	X[i*2][j] = 1.0;
	X[i*2+1][j] = 0.0;
      }
      else if(x[i* *n + j] == 2) {
	X[i*2][j] = 0.0;
	X[i*2+1][j] = 1.0;
      }
      else 
	X[i*2][j] = X[i*2+1][j] = 0.0;
    }
  }

  markerforwself2(*n, *m, X, y, *maxsize, chosen, rss);
}

/**********************************************************************
 * markerforwsel 
 * 
 * n = number of individuals
 * m = number of covariates (not including intercept)
 *
 * X = covariate matrix, indexed as X[covariate][individual]
 * y = outcome
 *
 * maxsize = maximum number of covariates
 *
 * chosen = on output, index [0, 1, ..., (m-1)] of chosen covariates
 * rss = on output, rss for those models
 *
 **********************************************************************/
void markerforwsel(int n, int m, double **X, double *y,
		   int maxsize, int *chosen, double *rss)
{
  double minrss, *work, sxx, sxy, syy; 
  double minsxx=0.0, minsxy=0.0, currss, ym;
  int *ignore;
  int i, j, k;

  /* allocate space */
  work = (double *)R_alloc(m, sizeof(double));
  ignore = (int *)R_alloc(m, sizeof(int));
  
  for(j=0; j<m; j++) {
    ignore[j] = 0;
    work[j] = 0.0;
  }

  /* recenter everything */
  ym = 0.0;
  for(i=0; i<n; i++) {
    ym += y[i];
    for(j=0; j<m; j++) work[j] += X[j][i];
  }

  ym /= (double)n;
  for(j=0; j<m; j++) work[j] /= (double)n;

  minrss = 0.0;
  for(i=0; i<n; i++) {
    y[i] -= ym;
    minrss += (y[i]*y[i]);
    for(j=0; j<m; j++) X[j][i] -= work[j];
  }


  for(k=0; k<maxsize; k++) { 
    chosen[k] = -1;

    syy = minrss;

    /* loop over markers */
    for(j=0; j<m; j++) {
      if(!ignore[j]) {

	/* calculate RSS */
	sxx = sxy = 0.0;
	for(i=0; i<n; i++) {
	  sxx += (X[j][i]*X[j][i]);
	  sxy += (X[j][i]*y[i]);
	}
	
	currss = syy - sxy*sxy/sxx;

	if(currss < minrss) { /* best so far */
	  rss[k] = minrss = currss;
	  minsxx = sxx;
	  minsxy = sxy;
	  chosen[k] = j;
	}
      }

    }

    if(k==maxsize) break;

    ignore[chosen[k]] = 1;

    /* recenter y */
    for(i=0; i<n; i++) 
      y[i] -= (X[chosen[k]][i]*minsxy/minsxx);

    /* recenter other x's */
    for(j=0; j<m; j++) {
	
      if(!ignore[j]) {
      
	sxy = 0.0;
	for(i=0; i<n; i++) 
	  sxy += (X[j][i]*X[chosen[k]][i]);
	for(i=0; i<n; i++) 
	  X[j][i] -= (X[chosen[k]][i] * sxy/minsxx);
      }
    }

  }
}


/**********************************************************************
 * markerforwself2
 * 
 * the same as markerforwsel, but for an intercross, in which each
 * column must be expanded to two, and we must select on the pairs of
 * columns.
 *
 **********************************************************************/
void markerforwself2(int n, int m, double **X, double *y,
		     int maxsize, int *chosen, double *rss)
{
  double minrss, *work, sxx, sxy, syy, tsyy; 
  double currss, *tempy, ym;
  int *ignore;
  int i, j, k, s;

  /* allocate space */
  work = (double *)R_alloc(m*2, sizeof(double));
  tempy = (double *)R_alloc(n, sizeof(double));
  ignore = (int *)R_alloc(m, sizeof(int));
  
  for(j=0; j<m; j++) {
    ignore[j] = 0;
    work[j] = 0.0;
  }

  /* recenter everything */
  ym = 0.0;
  for(i=0; i<n; i++) {
    ym += y[i];
    for(j=0; j<m*2; j++) work[j] += X[j][i];
  }

  ym /= (double)n;
  for(j=0; j<m*2; j++) work[j] /= (double)n;

  minrss = 0.0;
  for(i=0; i<n; i++) {
    y[i] -= ym;
    minrss += (y[i]*y[i]);
    for(j=0; j<m*2; j++) X[j][i] -= work[j];
  }

  /* recenter second X relative to first X for each marker */
  for(j=0; j<m; j++) {
    sxx = sxy = 0.0;
    for(i=0; i<n; i++) {
      sxx += (X[j*2][i]*X[j*2][i]);
      sxy += (X[j*2+1][i]*X[j*2][i]);
    }
    for(i=0; i<n; i++) 
      X[j*2+1][i] -= (X[j*2][i] * sxy / sxx);
  }

  for(k=0; k<maxsize; k++) { 
    chosen[k] = -1;

    syy = minrss;

    /* loop over markers */
    for(j=0; j<m; j++) {
      if(!ignore[j]) {

	/* center y with first column */
	sxx = sxy = 0.0;
	for(i=0; i<n; i++) {
	  sxx += (X[j*2][i]*X[j*2][i]);
	  sxy += (X[j*2][i]*y[i]);
	}
	
	tsyy = 0.0;
	for(i=0; i<n; i++) {
	  tempy[i] = y[i] - X[j*2][i]*sxy/sxx;
	  tsyy += (tempy[i] * tempy[i]);
	}

	/* get rss from regr on second column */
	sxx = sxy = 0.0;
	for(i=0; i<n; i++) {
	  sxx += (X[j*2+1][i]*X[j*2+1][i]);
	  sxy += (X[j*2+1][i]*y[i]);
	}
	
	currss = tsyy - sxy*sxy/sxx;

	if(currss < minrss) { /* best so far */
	  rss[k] = minrss = currss;
	  chosen[k] = j;
	}
      }

    }

    if(k==maxsize) break;

    ignore[chosen[k]] = 1;

    /* recenter y and all of the x's with chosen columns */
    sxx = sxy = 0.0;
    for(i=0; i<n; i++) {
      sxx += (X[chosen[k]*2][i] * X[chosen[k]*2][i]);
      sxy += (X[chosen[k]*2][i] * y[i]);
    }
    for(i=0; i<n; i++)
      y[i] -= (X[chosen[k]*2][i] * sxy / sxx);

    sxx = sxy = 0.0;
    for(i=0; i<n; i++) {
      sxx += (X[chosen[k]*2+1][i] * X[chosen[k]*2+1][i]);
      sxy += (X[chosen[k]*2+1][i] * y[i]);
    }
    for(i=0; i<n; i++)
      y[i] -= (X[chosen[k]*2+1][i] * sxy / sxx);

    /* recenter other x's */
    for(j=0; j<m; j++) {
	
      if(!ignore[j]) {
      
	for(s=0; s<1; s++) {
	  sxx = sxy = 0.0;
	  for(i=0; i<n; i++) {
	    sxx += (X[chosen[k]*2][i] * X[chosen[k]*2][i]);
	    sxy += (X[chosen[k]*2][i] * X[j*2+s][i]);
	  }
	  for(i=0; i<n; i++)
	    X[j*2+s][i] -= (X[chosen[k]*2][i] * sxy / sxx);

	  sxx = sxy = 0.0;
	  for(i=0; i<n; i++) {
	    sxx += (X[chosen[k]*2+1][i] * X[chosen[k]*2+1][i]);
	    sxy += (X[chosen[k]*2+1][i] * X[j*2+s][i]);
	  }
	  for(i=0; i<n; i++)
	    X[j*2+s][i] -= (X[chosen[k]*2+1][i] * sxy / sxx);
	}
      }
    }

  }
}

/* end of forwsel.c */
