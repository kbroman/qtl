/**********************************************************************
 * 
 * hmm_util.c
 * 
 * copyright (c) 2001-9, Karl W Broman
 * modified from hmm_main.c by Brian S Yandell and Laura M Shannon (c) 2011
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
 * Contains: init_stepf, stepfc, forward, backward, golden
 *
 * These are used in hmm_bcsft to simplify calculations.
 * They could be used in hmm_main.c
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "util.h"
   
void init_stepf(double *rf, double *rf2, int n_gen, int n_mar, int *cross_scheme, 
		double stepf(int, int, double, double, int *),
		double **probmat)
{
  int j,obs1,obs2,tmp1;
  
  for(j=0; j<n_mar-1; j++) {
    for(obs2=1; obs2<=n_gen; obs2++) {
      tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
      for(obs1=1; obs1<=obs2; obs1++)
	probmat[j][obs1 + tmp1] = stepf(obs1, obs2, rf[j], rf2[j], cross_scheme);
    }
  }
}

double stepfc(int obs1, int obs2, int mar, double **probmat)
{
  int tmp1;
  
  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    tmp1 = obs2;
    obs2 = obs1;
    obs1 = tmp1;
  }
  tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
  return(probmat[mar][obs1 + tmp1]);
}

void forward_prob(int i, int n_mar, int n_gen, int curpos, int *cross_scheme, double error_prob,
	     int **Geno, double **probmat, double **alpha,
	     double initf(int, int *), 
	     double emitf(int, int, double, int *))
{
  /* forward equations */

  /* Note: true genotypes coded as 1, 2, ...
     but in the alpha's and beta's, we use 0, 1, ... */

  int j,v,v2;
  double errortol,salpha;

  /* initialize alpha */
  /* curpos = -1: use error_prob always */
  /* curpos >= 0: use TOL except when j == curpos, then use error_prob */

  errortol = error_prob;
  if(curpos > 0) errortol = TOL;
  for(v=0; v<n_gen; v++)
    alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, errortol, cross_scheme);
  if(curpos == 0) errortol = TOL;

  for(j=1; j<n_mar; j++) {
    if(curpos == j) errortol = error_prob;
    
    for(v=0; v<n_gen; v++) {
      salpha = alpha[0][j-1] + stepfc(1, v+1, j-1, probmat);
      
      for(v2=1; v2<n_gen; v2++)
	salpha = addlog(salpha, alpha[v2][j-1] + stepfc(v2+1, v+1, j-1, probmat));
      
      alpha[v][j] = salpha + emitf(Geno[j][i], v+1, errortol, cross_scheme);
    }
    if(curpos == j) errortol = TOL;
  }
}
void backward_prob(int i, int n_mar, int n_gen, int curpos, int *cross_scheme, double error_prob,
		   int **Geno, double **probmat, double **beta,
		   double initf(int, int *), 
		   double emitf(int, int, double, int *))
{
  /* backward equations */

  /* Note: true genotypes coded as 1, 2, ...
     but in the alpha's and beta's, we use 0, 1, ... */

  /* could divide this into forward and backward, then use forward second time
     in est_map */

  int j2,v,v2;
  double errortol,sbeta;

  /* initialize alpha and beta */
  for(v=0; v<n_gen; v++)
    beta[v][n_mar-1] = 0.0;

  /* curpos = -1: use error_prob always */
  /* curpos >= 0: use TOL except when j2+1 == curpos, then use error_prob */
  errortol = error_prob;
  if(curpos >= 0) errortol = TOL;

  for(j2=n_mar-2; j2>=0; j2--) {
    if(curpos == j2+1) errortol = error_prob;
    
    for(v=0; v<n_gen; v++) {
      sbeta = beta[0][j2+1] + stepfc(v+1, 1, j2, probmat) +
	emitf(Geno[j2+1][i], 1, errortol, cross_scheme);
      
      for(v2=1; v2<n_gen; v2++) {
	sbeta = addlog(sbeta, beta[v2][j2+1] + stepfc(v+1, v2+1, j2, probmat) +
		       emitf(Geno[j2+1][i], v2+1, errortol, cross_scheme));
      }
      beta[v][j2] = sbeta;
    }
    if(curpos == j2+1) errortol = TOL;
  }
}

void calc_probfb(int i, int n_mar, int n_gen, int curpos, double **alpha, double **beta,
		 double ***Genoprob)
{
  int j,v,j0,jmax;
  double s;

  j0 = 0;
  jmax = n_mar;
  if(curpos >= 0) {
    j0 = curpos;
    jmax = j0 + 1;
  }

  /* calculate genotype probabilities */
  for(j=j0; j<jmax; j++) {
    s = Genoprob[0][j][i] = alpha[0][j] + beta[0][j];
    for(v=1; v<n_gen; v++) {
      Genoprob[v][j][i] = alpha[v][j] + beta[v][j];
      s = addlog(s, Genoprob[v][j][i]);
    }
    for(v=0; v<n_gen; v++) 
      Genoprob[v][j][i] = exp(Genoprob[v][j][i] - s);
  }
}

double golden_search(double *countmat, int n_gen, int maxit, double tol, int *cross_scheme,
	      double comploglik(double, int, double *, int *))
{
  /* Golden section search. */
  /* en.wikipedia.org/wiki/Golden_section_search */

  static double resphi = 0.0;
  double x[4],y[4];
  int iter;

  if(resphi == 0.0)
    resphi = 1.5 - sqrt(5.0) / 2.0;

  x[0] = 0.0;
  x[2] = 1.0;
  y[0] = comploglik(0.0, n_gen, countmat, cross_scheme);
  y[2] = comploglik(0.5, n_gen, countmat, cross_scheme);

  if(y[2] < y[0]) {
    x[1] = x[0];
    x[0] = x[2];
    x[2] = x[1];
    y[1] = y[0];
    y[0] = y[2];
    y[2] = y[1];
  }
  
  x[1] = x[0] + resphi * (x[2] - x[0]);
  y[1] = comploglik(x[1], n_gen, countmat, cross_scheme);

  /* x0 and x2 are the current bounds; the minimum is between them.
   * x1 is the center point, which is closer to x0 than to x2. */

  for(iter=0; iter<maxit; iter++) {
    /* Create a new possible center in the area between x1 and x2, closer to x1. */
    x[3] = x[1] + resphi * (x[2] - x[1]);

    /* Evaluate termination criterion */
    if(fabs(x[2] - x[0]) < tol)
      break;
 
    y[3] = comploglik(x[3], n_gen, countmat, cross_scheme);

    if(y[3] >= y[1]) {
      x[0] = x[1];
      x[1] = x[3];
      y[0] = y[1];
      y[1] = y[3];
    }
    else {
      x[2] = x[0];
      x[0] = x[3];
      y[2] = y[0];
      y[0] = y[3];
    }
  }
  /* handle boundary situations cleanly */
  if((x[0] == 0.0 && y[0] >= y[1]) || (x[2] == 0.0 && y[2] >= y[1])) return(0.0);
  if((x[0] == 1.0 && y[0] >= y[1]) || (x[2] == 1.0 && y[2] >= y[1])) return(1.0);

  x[1] = (x[2] + x[0]) / 2.0;
  /* make negative if does not converge */
  if(iter >= maxit)
    x[1] = - x[1];
  return(x[1]);
}
  
/* end of hmm_util.c */
