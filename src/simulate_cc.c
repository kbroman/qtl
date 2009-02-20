/**********************************************************************
 * 
 * simulate_cc.c
 *
 * copyright (c) 2005-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Mar, 2005
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License, as
 *     published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version. 
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the
 *     GNU General Public License for more details.
 * 
 *     A copy of the GNU General Public License is available at
 *     http://www.r-project.org/Licenses/
 *
 * C functions for the R/qtl package
 *
 * These functions are for simulating experimental cross data;
 * I start with the simulation of RILs
 *
 * Contains: R_sim_ril, sim_ril, 
 *           allocate_individual, reallocate_individual, 
 *           copy_individual, docross, meiosis, sim_cc, R_sim_cc
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "simulate_cc.h"

/* wrapper for sim_ril, to be called from R */
void R_sim_ril(int *n_chr, int *n_mar, int *n_ril, double *map,
	       int *n_str, int *m, double *p, int *include_x, 
	       int *random_cross, int *selfing, int *cross, int *ril)
{
  GetRNGstate();

  sim_ril(*n_chr, n_mar, *n_ril, map, *n_str, *m, *p, *include_x, 
	  *random_cross, *selfing, cross, ril);

  PutRNGstate();
}
	  
/**********************************************************************
 * 
 * sim_ril
 * 
 * n_chr   Number of chromosomes
 * n_mar   Number of markers on each chromosome (vector of length n_chr)
 * n_ril   Number of RILs to simulate
 * 
 * map     Vector of marker locations, of length sum(n_mar)
 *         First marker on each chromosome should be at 0.
 *
 * n_str   Number of parental strains (either 2, 4, or 8)
 *
 * m       Interference parameter (0 is no interference)
 * p       For Stahl model, proportion of chiasmata from the NI model
 *
 * include_x   Whether the last chromosome is the X chromosome
 *
 * random_cross  Indicates whether the order of the strains in the cross
 *               should be randomized.
 *
 * selfing If 1, use selfing; if 0, use sib mating
 *
 * cross   On output, the cross used for each line 
 *         (vector of length n_ril x n_str)
 *
 * ril     On output, the simulated data 
 *         (vector of length sum(n_mar) x n_ril)
 *
 **********************************************************************/
void sim_ril(int n_chr, int *n_mar, int n_ril, double *map, 
	     int n_str, int m, double p, int include_x, 
	     int random_cross, int selfing, int *cross, int *ril) 
{
  int i, j, k, ngen, tot_mar, curseg;
  struct individual par1, par2, kid1, kid2;
  int **Ril, **Cross, maxwork, isX, flag, max_xo, *firstmarker;
  double **Map, maxlen, chrlen, *work;

 /* count total number of markers */
  for(i=0, tot_mar=0; i<n_chr; i++) 
    tot_mar += n_mar[i];

  reorg_geno(tot_mar, n_ril, ril, &Ril);
  reorg_geno(n_str, n_ril, cross, &Cross);

  /* allocate space */
  Map = (double **)R_alloc(n_chr, sizeof(double *));
  Map[0] = map;
  for(i=1; i<n_chr; i++)
    Map[i] = Map[i-1] + n_mar[i-1];

  /* location of first marker */
  firstmarker = (int *)R_alloc(n_chr, sizeof(int));
  firstmarker[0] = 0;
  for(i=1; i<n_chr; i++) 
    firstmarker[i] = firstmarker[i-1] + n_mar[i-1];

  /* maximum chromosome length (in cM) */
  maxlen = Map[0][n_mar[0]-1];
  for(i=1; i<n_chr; i++)
    if(maxlen < Map[i][n_mar[i]-1])
      maxlen =  Map[i][n_mar[i]-1];

  /* allocate space for individuals */
  max_xo = (int)qpois(1e-10, maxlen/100.0, 0, 0)*3;
  allocate_individual(&par1, max_xo);
  allocate_individual(&par2, max_xo);
  allocate_individual(&kid1, max_xo);
  allocate_individual(&kid2, max_xo);
  maxwork = (int)qpois(1e-10, (m+1)*maxlen/50.0, 0, 0)*3;
  work = (double *)R_alloc(maxwork, sizeof(double));

  for(i=0; i<n_ril; i++) {

    /* set up cross */
    for(j=0; j<n_str; j++) Cross[i][j] = j+1;
    if(random_cross) int_permute(Cross[i], n_str);

    for(j=0; j<n_chr; j++) {
      isX = include_x && j==n_chr-1;

      chrlen = Map[j][n_mar[j]-1];

      par1.n_xo[0] = par1.n_xo[1] = par2.n_xo[0] = par2.n_xo[1] = 0;

      /* initial generations */
      if(n_str==2) {
	par1.allele[0][0] = par2.allele[0][0] = 1;
	par1.allele[1][0] = par2.allele[1][0] = 2;
      }
      else if(n_str==4) {
	par1.allele[0][0] = 1;
	par1.allele[0][1] = 2;
	par2.allele[0][0] = 3;
	par2.allele[0][1] = 4;
      }
      else { /* 8 strain case */
	par1.allele[0][0] = 1;
	par1.allele[0][1] = 2;
	par2.allele[0][0] = 3;
	par2.allele[0][1] = 4;

	docross(par1, par2, &kid1, chrlen, m, p, 0, 
	      &maxwork, &work);

	par1.allele[0][0] = 5;
	par1.allele[0][1] = 6;
	par2.allele[0][0] = 7;
	par2.allele[0][1] = 8;

	docross(par1, par2, &kid2, chrlen, m, p, isX,
	      &maxwork, &work);

	copy_individual(&kid1, &par1);
	copy_individual(&kid2, &par2);
      }

      /* start inbreeding */
      ngen=1;
      while(1) {
	R_CheckUserInterrupt(); /* check for ^C */

	docross(par1, par2, &kid1, chrlen, m, p, 0,
		&maxwork, &work);
	if(!selfing) 
	  docross(par1, par2, &kid2, chrlen, m, p, isX,
		  &maxwork, &work);

	/* are we done? */
	flag = 0;
	if(selfing) {
	  if(kid1.n_xo[0] == kid1.n_xo[1]) {
	    for(k=0; k<kid1.n_xo[0]; k++) {
	      if(kid1.allele[0][k] != kid1.allele[1][k] ||
		 fabs(kid1.xoloc[0][k] - kid1.xoloc[1][k]) > 1e-6) {
		flag = 1;
		break;
	      }
	    }
	    if(kid1.allele[0][kid1.n_xo[0]] != kid1.allele[1][kid1.n_xo[0]])
	      flag = 1;
	  }
	  else flag = 1;
	}
	else {
	  if(kid1.n_xo[0] == kid1.n_xo[1] && 
	     kid1.n_xo[0] == kid2.n_xo[0] && 
	     kid1.n_xo[0] == kid2.n_xo[1]) {
	    for(k=0; k<kid1.n_xo[0]; k++) {
	      if(kid1.allele[0][k] != kid1.allele[1][k] ||
		 kid1.allele[0][k] != kid2.allele[0][k] ||
		 kid1.allele[0][k] != kid2.allele[1][k] ||
		 fabs(kid1.xoloc[0][k] - kid1.xoloc[1][k]) > 1e-6 ||
		 fabs(kid1.xoloc[0][k] - kid2.xoloc[0][k]) > 1e-6 ||
		 fabs(kid1.xoloc[0][k] - kid2.xoloc[1][k]) > 1e-6) {
		flag = 1;
		break;
	      }
	    }
	    if(kid1.allele[0][kid1.n_xo[0]] != kid1.allele[1][kid1.n_xo[0]] ||
	       kid1.allele[0][kid1.n_xo[0]] != kid2.allele[0][kid1.n_xo[0]] ||
	       kid1.allele[0][kid1.n_xo[0]] != kid2.allele[1][kid1.n_xo[0]]) 
	      flag = 1;
	  }
	  else flag = 1;
	}

	if(!flag) break; /* done inbreeding */

	/* go to next generation */
	copy_individual(&kid1, &par1);
	if(selfing) copy_individual(&kid1, &par2);
	else copy_individual(&kid2, &par2);

      } /* end with inbreeding of this chromosome */

      /* fill in alleles */
      curseg = 0;
      for(k=0; k<n_mar[j]; k++) { /* loop over markers */
	while(curseg < kid1.n_xo[0] && Map[j][k] > kid1.xoloc[0][curseg]) 
	  curseg++;
	  
	Ril[i][k+firstmarker[j]] = Cross[i][kid1.allele[0][curseg]-1];
      }

    } /* loop over chromosomes */

  } /* loop over lines */

}

/**********************************************************************
 * allocate_individual
 **********************************************************************/
void allocate_individual(struct individual *ind, int max_seg)
{
  (*ind).max_segments = max_seg;
  (*ind).n_xo[0] = (*ind).n_xo[1] = 0;
  (*ind).allele = (int **)R_alloc(2, sizeof(int *));
  (*ind).allele[0] = (int *)R_alloc(2*max_seg, sizeof(int));
  (*ind).allele[1] = (*ind).allele[0] + max_seg;
  (*ind).xoloc = (double **)R_alloc(2, sizeof(double *));
  (*ind).xoloc[0] = (double *)R_alloc(2*(max_seg-1), sizeof(double));
  (*ind).xoloc[1] = (*ind).xoloc[0] + (max_seg-1);
}


/**********************************************************************
 * reallocate_individual
 **********************************************************************/
void reallocate_individual(struct individual *ind, int old_max_seg, 
			   int new_max_seg)
{
  (*ind).max_segments = new_max_seg;
  (*ind).allele[0] = (int *)S_realloc((char *)(*ind).allele[0], 2*new_max_seg, 
				     2*old_max_seg, sizeof(int));
  (*ind).allele[1] = (*ind).allele[0] + new_max_seg;
  (*ind).xoloc[0] = (double *)S_realloc((char *)(*ind).xoloc[0], 2*(new_max_seg-1), 
				      2*(old_max_seg-1), sizeof(double));
  (*ind).xoloc[0] = (*ind).xoloc[1] + new_max_seg-1;
}

/**********************************************************************
 * copy_individual
 **********************************************************************/
void copy_individual(struct individual *from, struct individual *to)
{
  int i, j, n_xo;

  if((*from).max_segments > (*to).max_segments) 
    reallocate_individual(to, (*to).max_segments, (*from).max_segments);

  for(j=0; j<2; j++) {
    n_xo = (*to).n_xo[j] = (*from).n_xo[j];
    for(i=0; i<n_xo; i++) {
      (*to).allele[j][i] = (*from).allele[j][i];
      (*to).xoloc[j][i] = (*from).xoloc[j][i];
    }
    (*to).allele[j][n_xo] = (*from).allele[j][n_xo];
  }
}

/**********************************************************************
 * docross
 *
 * Cross two individuals and get a kid
 **********************************************************************/
void docross(struct individual par1, struct individual par2,
	     struct individual *kid, double L, int m,
	     double p, int isX, int *maxwork, double **work)
{
  int i, j, n_xo, curstrand, new_n_xo, curloc[2];

  meiosis(L, m, p, maxwork, work, &n_xo);

  /* need to reallocate? */
  i = par1.n_xo[0]+par1.n_xo[1]+n_xo;
  if((*kid).max_segments < i) 
    reallocate_individual(kid, (*kid).max_segments, i*2);

  if(unif_rand() < 0.5) curstrand = 0;
  else curstrand = 1;

  (*kid).allele[0][0] = par1.allele[curstrand][0];

  curloc[0] = curloc[1] = 0;
  new_n_xo = 0;
  for(i=0; i<n_xo; i++) {
    for(; curloc[curstrand] < par1.n_xo[curstrand] && 
	  par1.xoloc[curstrand][curloc[curstrand]] < (*work)[i]; curloc[curstrand]++) {
      (*kid).allele[0][new_n_xo+1] = par1.allele[curstrand][curloc[curstrand]+1];
      (*kid).xoloc[0][new_n_xo] = par1.xoloc[curstrand][curloc[curstrand]];
      new_n_xo++;
    }
    curstrand = 1-curstrand;

    /* skip ahead */
    for(; curloc[curstrand] < par1.n_xo[curstrand] &&
	  par1.xoloc[curstrand][curloc[curstrand]] < (*work)[i]; curloc[curstrand]++);

    (*kid).allele[0][new_n_xo+1] = par1.allele[curstrand][curloc[curstrand]];
    (*kid).xoloc[0][new_n_xo] = (*work)[i];

    new_n_xo++;
  } 

  for(; curloc[curstrand] < par1.n_xo[curstrand] && 
	par1.xoloc[curstrand][curloc[curstrand]] < L; curloc[curstrand]++) {
    (*kid).allele[0][new_n_xo+1] = par1.allele[curstrand][curloc[curstrand]+1];
    (*kid).xoloc[0][new_n_xo] = par1.xoloc[curstrand][curloc[curstrand]];
    new_n_xo++;
  }

  for(j=0, i=0; j<new_n_xo; j++) {
    if((*kid).allele[0][i] != (*kid).allele[0][j+1]) { 
      (*kid).allele[0][i+1] = (*kid).allele[0][j+1];
      (*kid).xoloc[0][i] = (*kid).xoloc[0][j];
      i++;
    }
  }
  (*kid).n_xo[0] = i;

  if(isX) { /* X chromosome */
    (*kid).n_xo[1] = par2.n_xo[0];
    for(j=0; j<par2.n_xo[0]; j++) {
      (*kid).allele[1][j] = par2.allele[0][j];
      (*kid).xoloc[1][j] = par2.xoloc[0][j];
    }
    (*kid).allele[1][par2.n_xo[0]] = 
      par2.allele[0][par2.n_xo[0]];
  }
  else {

    meiosis(L, m, p, maxwork, work, &n_xo);

    /* need to reallocate? */
    i = par2.n_xo[0]+par2.n_xo[1]+n_xo;
    if((*kid).max_segments < i) 
      reallocate_individual(kid, (*kid).max_segments, i*2);
    
    if(unif_rand() < 0.5) curstrand = 0;
    else curstrand = 1;

    (*kid).allele[1][0] = par2.allele[curstrand][0];

    curloc[0] = curloc[1] = 0;
    new_n_xo = 0;
    for(i=0; i<n_xo; i++) {
      for(; curloc[curstrand] < par2.n_xo[curstrand] && 
	    par2.xoloc[curstrand][curloc[curstrand]] < (*work)[i]; curloc[curstrand]++) {
	(*kid).allele[1][new_n_xo+1] = par2.allele[curstrand][curloc[curstrand]+1];
	(*kid).xoloc[1][new_n_xo] = par2.xoloc[curstrand][curloc[curstrand]];
	new_n_xo++;
      }
      curstrand = 1-curstrand;

      /* skip ahead */
      for(; curloc[curstrand] < par2.n_xo[curstrand] &&
	    par2.xoloc[curstrand][curloc[curstrand]] < (*work)[i]; curloc[curstrand]++);

      (*kid).allele[1][new_n_xo+1] = par2.allele[curstrand][curloc[curstrand]];
      (*kid).xoloc[1][new_n_xo] = (*work)[i];
      new_n_xo++;
    } 
    
    for(; curloc[curstrand] < par2.n_xo[curstrand] && 
	  par2.xoloc[curstrand][curloc[curstrand]] < L; curloc[curstrand]++) {
      (*kid).allele[1][new_n_xo+1] = par2.allele[curstrand][curloc[curstrand]+1];
      (*kid).xoloc[1][new_n_xo] = par2.xoloc[curstrand][curloc[curstrand]];
      new_n_xo++;
    }
    
    for(j=0, i=0; j<new_n_xo; j++) {
      if((*kid).allele[1][i] != (*kid).allele[1][j+1]) { 
	(*kid).allele[1][i+1] = (*kid).allele[1][j+1];
	(*kid).xoloc[1][i] = (*kid).xoloc[1][j];
	i++;
      }
    }
    (*kid).n_xo[1] = i;
  }
}


/**********************************************************************
 * 
 * meiosis
 *
 * chrlen Chromosome length (in cM) 
 *
 * m      interference parameter (0 corresponds to no interference)
 *
 * p      for stahl model, proportion of chiasmata from NI mechanism
 *
 * maxwork
 * work
 * 
 * n_xo
 *
 **********************************************************************/
void meiosis(double L, int m, double p, int *maxwork, double **work,
	     int *n_xo)
{
  int i, n, nn, j, first;

  if(m > 0 && p < 1.0) { /* crossover interference */

    /* simulate number of XOs and intermediates */
    n = (int)rpois(L*(double)(m+1)/50.0*(1.0-p));

    if(n > *maxwork) { /* need a bigger workspace */
      *work = (double *)S_realloc((char *)*work, n*2, *maxwork, sizeof(double));
      *maxwork = n*2;
    }

    for(i=0; i<n; i++) 
      (*work)[i] = L*unif_rand();
    /* sort them */
    R_rsort(*work, n);
    
    /* which is the first crossover? */
    first = random_int(0,m);

    for(i=first, j=0; i<n; i += (m+1), j++) 
      (*work)[j] = (*work)[i];
    n = j;
  
    /* thin with probability 1/2 */
    for(i=0, j=0; i<n; i++) {
      if(unif_rand() < 0.5) {
	(*work)[j] = (*work)[i]; 
	j++;
      }
    }
    n = j;

    nn = (int) rpois(L*p/100.0);
    if(n +nn > *maxwork) { /* need a bigger workspace */
      *work = (double *)S_realloc((char *)*work, (n+nn)*2, *maxwork, sizeof(double));
      *maxwork = (n+nn)*2;
    }
    
    for(i=0; i<nn; i++) 
      (*work)[i+n] = L*unif_rand();
    R_rsort(*work, n+nn);

    *n_xo = n+nn;
  }

  else { /* no crossover interference */
    n = (int) rpois(L/100.0);

    if(n > *maxwork) { /* need a bigger workspace */
      *work = (double *)S_realloc((char *)*work, n*2, *maxwork, sizeof(double));
      *maxwork = n*2;
    }

    for(i=0; i<n; i++) 
      (*work)[i] = L*unif_rand();
    /* sort them */
    R_rsort(*work, n);

    *n_xo = n;
  }
}


void R_meiosis(double *L, int *m, double *p, int *maxwork, double *work,
	       int *n_xo) 
{
  GetRNGstate();
  meiosis(*L, *m, *p, maxwork, &work, n_xo);
  PutRNGstate();
}


/**********************************************************************
 * 
 * sim_cc    Use the result of sim_all_ril with n_str=8 plus data on
 *           the SNP genotypes of the 8 parental strains to create 
 *           real SNP data for the Collaborative Cross
 *
 * n_ril     Number of RILs to simulate
 * tot_mar   Total number of markers
 *
 * Parents   SNP data for the 8 parental lines [dim tot_mar x 8]
 * 
 * Geno      On entry, the detailed genotype data; on exit, the 
 *           SNP data written bitwise.
 * 
 * error_prob  Probability of genotyping error
 * missing_prob  Probability a genotype will be missing
 *
 **********************************************************************/
void sim_cc(int n_ril, int tot_mar, int **Parents, int **Geno,
	    double error_prob, double missing_prob)
{
  int i, j, k, temp;

  for(i=0; i<n_ril; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    for(j=0; j<tot_mar; j++) {
      temp = Parents[Geno[j][i]-1][j];
      if(unif_rand() < error_prob)  /* switch the SNP genotype */
	temp = 1-temp;

      Geno[j][i] = 0;
      if(unif_rand() > missing_prob) {/* no error; convert to bit string */
	for(k=0; k<8; k++) 
	  if(temp == Parents[k][j]) Geno[j][i] += (1<<k);
      }
    }
  }
}

/* wrapper for calling sim_cc from R */
void R_sim_cc(int *n_ril, int *tot_mar, int *parents, int *geno,
	      double *error_prob, double *missing_prob)
{
  int **Parents, **Geno;

  reorg_geno(*tot_mar, 8, parents, &Parents);
  reorg_geno(*n_ril, *tot_mar, geno, &Geno);

  GetRNGstate();

  sim_cc(*n_ril, *tot_mar, Parents, Geno, *error_prob, *missing_prob);

  PutRNGstate();
}

/* end of simulate_cc.c */

