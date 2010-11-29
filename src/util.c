/**********************************************************************
 * 
 * util.c
 *
 * copyright (c) 2001-2010, Karl W Broman and Hao Wu
 *
 * This file written mostly by Karl Broman with some additions
 * from Hao Wu.
 *
 * last modified Nov, 2010
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
 * These are utility functions, mostly for the HMM engine.
 *
 * Other functions: addlog, subtrlog, reorg_geno, reorg_genoprob,
 *                  reorg_pairprob, allocate_int
 *                  allocate_alpha, reorg_draws, allocate_double,
 *                  sample_int, allocate_imatrix, allocate_dmatrix
 *                  reorg_errlod, double_permute, int_permute, 
 *                  random_int
 *                  wtaverage, comparegeno, R_comparegeno,
 *                  R_locate_xo, locate_xo, matmult, expand_col2drop
 *                  dropcol_xpx, dropcol_xpy, dropcol_x,
 *                  reviseMWril, R_reviseMWril, R_calcPermPval,
 *                  calcPermPval
 *
 **********************************************************************/

#include <R.h>
/* #include <R_ext/BLAS.h> */
#include "util.h"

#define THRESH 200.0

/**********************************************************************
 * 
 * addlog
 *
 * Calculate addlog(a,b) = log[exp(a) + exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.   
 *
 **********************************************************************/
double addlog(double a, double b)
{
  if(b > a + THRESH) return(b);
  else if(a > b + THRESH) return(a);
  else return(a + log1p(exp(b-a)));
}
		       
/**********************************************************************
 * 
 * subtrlog
 *
 * Calculate subtrlog(a,b) = log[exp(a) - exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.  
 *
 **********************************************************************/
double subtrlog(double a, double b)
{
  if(a > b + THRESH) return(a);
  else return(a + log1p(-exp(b-a)));
}

/**********************************************************************
 * 
 * reorg_geno
 *
 * Reorganize the genotype data so that it is a doubly indexed array
 * rather than a single long vector
 *
 * Afterwards, geno indexed like Geno[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
  int i;

  *Geno = (int **)R_alloc(n_pos, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_pos; i++) 
    (*Geno)[i] = (*Geno)[i-1] + n_ind;

}

/**********************************************************************
 * 
 * reorg_genoprob
 *
 * Reorganize the genotype probability data so that it is a triply 
 * indexed array rather than a single long vector
 *
 * Afterwards, genoprob indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_genoprob(int n_ind, int n_pos, int n_gen, 
		    double *genoprob, double ****Genoprob)
{
  int i, j;
  double **a;

  *Genoprob = (double ***)R_alloc(n_gen, sizeof(double **));

  a = (double **)R_alloc(n_pos*n_gen, sizeof(double *));

  (*Genoprob)[0] = a;
  for(i=1; i< n_gen; i++) 
    (*Genoprob)[i] = (*Genoprob)[i-1]+n_pos;
  
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_pos; j++) 
      (*Genoprob)[i][j] = genoprob + i*n_ind*n_pos + j*n_ind;
}

/**********************************************************************
 * 
 * reorg_pairprob
 *
 * Reorganize the joint genotype probabilities so that they form a 
 * quintuply indexed array rather than a single long vector
 *
 * Afterwards, pairprob indexed like 
 *    Pairprob[gen1][gen2][pos1][pos2][ind] with pos2 > pos1
 * 
 * You *must* refer to cases with pos2 > pos1, as cases with
 * pos2 <= pos1 point off into the ether.
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_pairprob(int n_ind, int n_pos, int n_gen, 
		    double *pairprob, double ******Pairprob)
{
  int i, j, k, s, n_pairs;
  double ****ptr1, ***ptr2, **ptr3;

  /* note: n_pos must be at least 2 */
  n_pairs = n_pos*(n_pos-1)/2;

  *Pairprob = (double *****)R_alloc(n_gen, sizeof(double ****));

  ptr1 = (double ****)R_alloc(n_gen*n_gen, sizeof(double ***));
  (*Pairprob)[0] = ptr1;
  for(i=1; i<n_gen; i++) 
    (*Pairprob)[i] = ptr1 + i*n_gen;

  ptr2 = (double ***)R_alloc(n_gen*n_gen*n_pos, sizeof(double **));
  for(i=0; i<n_gen; i++)
    for(j=0; j<n_gen; j++) 
      (*Pairprob)[i][j] = ptr2 + (i*n_gen+j)*n_pos;

  ptr3 = (double **)R_alloc(n_gen*n_gen*n_pos*n_pos, sizeof(double **));
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_gen; j++)
      for(k=0; k<n_pos; k++)
	(*Pairprob)[i][j][k] = ptr3 + ((i*n_gen+j)*n_pos + k)*n_pos;
  
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_gen; j++)
      for(k=0; k<n_pos-1; k++)
	for(s=(k+1); s<n_pos; s++) 
	  (*Pairprob)[i][j][k][s] = pairprob + (i*n_gen+j)*n_ind*n_pairs +
	    + n_ind*(2*n_pos-1-k)*k/2 + (s-k-1)*n_ind;
}

/**********************************************************************
 * 
 * allocate_alpha
 *
 * Allocate space for alpha and beta matrices
 *
 * Afterwards, indexed like alpha[gen][mar]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_alpha(int n_pos, int n_gen, double ***alpha)
{
  int i;

  *alpha = (double **)R_alloc(n_gen, sizeof(double *));

  (*alpha)[0] = (double *)R_alloc(n_gen*n_pos, sizeof(double));

  for(i=1; i< n_gen; i++) 
    (*alpha)[i] = (*alpha)[i-1] + n_pos;
}

/**********************************************************************
 * 
 * reorg_draws
 *
 * Reorganize the simulated genotypes so that it is a triply 
 * indexed array rather than a single long vector
 *
 * Afterwards, draws indexed like Draws[repl][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_draws(int n_ind, int n_pos, int n_draws, 
		 int *draws, int ****Draws)
{
  int i, j;
  int **a;

  *Draws = (int ***)R_alloc(n_draws, sizeof(int **));

  a = (int **)R_alloc(n_pos*n_draws, sizeof(int *));
  (*Draws)[0] = a;
  for(i=1; i<n_draws; i++) 
    (*Draws)[i] = (*Draws)[i-1]+n_pos;
  
  for(i=0; i<n_draws; i++) 
    for(j=0; j<n_pos; j++) 
      (*Draws)[i][j] = draws + (i*n_pos+j)*n_ind;
}

/**********************************************************************
 * 
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector)
{
  *vector = (double *)R_alloc(n, sizeof(double));
}

/**********************************************************************
 * 
 * allocate_int
 *
 * Allocate space for a vector of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_int(int n, int **vector)
{
  *vector = (int *)R_alloc(n, sizeof(int));
}

/**********************************************************************
 * 
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_dmatrix(int n_row, int n_col, double ***matrix)
{
  int i;

  *matrix = (double **)R_alloc(n_row, sizeof(double *));
  
  (*matrix)[0] = (double *)R_alloc(n_col*n_row, sizeof(double));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

/**********************************************************************
 * 
 * allocate_imatrix
 *
 * Allocate space for a matrix of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_imatrix(int n_row, int n_col, int ***matrix)
{
  int i;

  *matrix = (int **)R_alloc(n_row, sizeof(int *));
  
  (*matrix)[0] = (int *)R_alloc(n_col*n_row, sizeof(int));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

/**********************************************************************
 * 
 * sample_int
 *
 * Make a single draw from (1, ..., n) with probs (p_0, ..., p_(n-1))
 *
 **********************************************************************/
int sample_int(int n, double *p)
{
  int i;
  double r;

  /* R's random number generator */
  r = unif_rand();

  for(i=0; i<n; i++) {
    if(r < p[i]) return(i+1);
    else r -= p[i];
  }
  return(n); /* this shouldn't happen */
}

/**********************************************************************
 * 
 * reorg_errlod
 *
 * Just like reorg_geno(), only for a matrix of doubles.
 *
 * Afterwards, errlod indexed like Errlod[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod)
{
  int i;

  *Errlod = (double **)R_alloc(n_mar, sizeof(double *));

  (*Errlod)[0] = errlod;
  for(i=1; i< n_mar; i++) 
    (*Errlod)[i] = (*Errlod)[i-1] + n_ind;
}

/**********************************************************************
 * 
 * double_permute
 *
 *   This function randomly permutes a vector of doubles
 *   
 * Input:
 * 
 *   array = vector of doubles; on output, it contains a random 
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/
void double_permute(double *array, int len)
{
  int i, which;
  double tmp;
  
  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}

/**********************************************************************
 * 
 * int_permute
 *
 *   This function randomly permutes a vector of int
 *   
 * Input:
 * 
 *   array = vector of int; on output, it contains a random 
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/
void int_permute(int *array, int len)
{
  int i, which;
  int tmp;
  
  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}

/**********************************************************************
 * 
 * random_int
 *   
 * Generates a random int integer between "low" and "high", inclusive.
 *
 *  Input:
 * 
 *    low
 *
 *    high
 *
 **********************************************************************/
int random_int(int low, int high)
{
  return((int)(unif_rand()*(double)(high - low + 1)) + low);
}

/**********************************************************************
 * wtaverage
 * calculate the weight average of the LOD scores
 *********************************************************************/
double wtaverage(double *LOD, int n_draws)
{
  int k, idx, nnewLOD;
  double sum, sums, meanLOD, varLOD, *newLOD;

  /* calculate the number of LOD needs to be thrown */
  idx = (int) floor( 0.5*log(n_draws)/log(2) );
  nnewLOD = n_draws-2*idx; /* number of items in newLOD vector */
  /* allocate memory for newLOD */  
  newLOD = (double *)R_alloc( nnewLOD, sizeof(double) );

  /* sort the LOD scores in ascending order */
  R_rsort(LOD, n_draws);

  /* get a new list of LOD scores, throwing the biggest 
     and smallest idx LOD scores. */
  for(k=idx, sum=0.0; k<n_draws-idx; k++) {
    newLOD[k-idx] = LOD[k];
    sum += LOD[k]; /* calculate the sum of newLOD in the same loop */
  }

  /* calculate the mean of newLOD */
  meanLOD = sum / nnewLOD; 
  /* calculate the variance of newLOD */
  if(nnewLOD > 1) {
    for(k=0,sums=0.0; k<nnewLOD; k++) 
      sums += (newLOD[k]-meanLOD) * (newLOD[k]-meanLOD);
    varLOD = sums/(nnewLOD-1);
  }
  else varLOD = 0.0;

  /* return the weight average */
  return( meanLOD+0.5*log(10.0)*varLOD );

}

/**********************************************************************
 * comparegeno
 * 
 * Count number of matches in the genotypes for all pairs of
 * individuals.
 *
 * Input:
 *   
 **********************************************************************/
void comparegeno(int **Geno, int n_ind, int n_mar, 
		 int **N_Match, int **N_Missing)
{
  int i, j, k;

  for(i=0; i<n_ind;i++) {
    for(k=0; k<n_mar; k++) {
      if(Geno[k][i]==0) 
	(N_Missing[i][i])++;
      else
	(N_Match[i][i])++;
    }

    for(j=(i+1); j<n_ind; j++) {
      R_CheckUserInterrupt(); /* check for ^C */

      for(k=0; k<n_mar; k++) {
	if(Geno[k][i]==0 || Geno[k][j]==0) (N_Missing[i][j])++;
	else if(Geno[k][i] == Geno[k][j]) (N_Match[i][j])++;
      }
      N_Missing[j][i] = N_Missing[i][j];
      N_Match[j][i] = N_Match[i][j];
    }
  }

}

/**********************************************************************
 * R_comparegeno: wrapper for R
 **********************************************************************/
void R_comparegeno(int *geno, int *n_ind, int *n_mar, 
		   int *n_match, int *n_missing)
{
  int **Geno, **N_Match, **N_Missing;
  int i;

  /* allocate space */
  Geno = (int **)R_alloc(*n_mar, sizeof(int *));
  N_Match = (int **)R_alloc(*n_ind, sizeof(int *));
  N_Missing = (int **)R_alloc(*n_ind, sizeof(int *));

  Geno[0] = geno;
  N_Match[0] = n_match;
  N_Missing[0] = n_missing;
  for(i=1; i< *n_mar; i++) 
    Geno[i] = Geno[i-1] + *n_ind;
  for(i=1; i< *n_ind; i++) {
    N_Match[i] = N_Match[i-1] + *n_ind;
    N_Missing[i] = N_Missing[i-1] + *n_ind;
  }

  comparegeno(Geno, *n_ind, *n_mar, N_Match, N_Missing);
}
  
void R_locate_xo(int *n_ind, int *n_mar, int *type,
		 int *geno, double *map, 
		 double *location, int *nseen,
		 int *ileft, int *iright, double *left, double *right,
		 int *ntyped, int *full_info)
{
  int **Geno, **iLeft, **iRight, **nTyped;
  double **Location, **Left, **Right;

  reorg_geno(*n_ind, *n_mar, geno, &Geno);
  reorg_errlod(*n_ind, (*type+1)*(*n_mar-1), location, &Location);
  if(*full_info) {
    reorg_errlod(*n_ind, (*type+1)*(*n_mar-1), left, &Left);
    reorg_errlod(*n_ind, (*type+1)*(*n_mar-1), right, &Right);
    reorg_geno(*n_ind, (*type+1)*(*n_mar-1), ileft, &iLeft);
    reorg_geno(*n_ind, (*type+1)*(*n_mar-1), iright, &iRight);
    reorg_geno(*n_ind, (*type+1)*(*n_mar-1), ntyped, &nTyped);
  }

  locate_xo(*n_ind, *n_mar, *type, Geno, map, Location,
	    nseen, iLeft, iRight, Left, Right, nTyped, *full_info);
}

/* Note: type ==0 for backcross and ==1 for intercross */
void locate_xo(int n_ind, int n_mar, int type, int **Geno,
	       double *map, double **Location, int *nseen,
	       int **iLeft, int **iRight, double **Left, double **Right,
	       int **nTyped, int full_info)
{
  int i, j, k, curgen, number, icurpos;
  double curpos;

  for(i=0; i<n_ind; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    curgen = Geno[0][i];
    curpos = map[0];
    icurpos = 0;
    nseen[i]=0;
    for(j=1; j<n_mar; j++) {
      if(curgen==0) { /* haven't yet seen a genotype */
	curgen = Geno[j][i];
	curpos = map[j];
	icurpos = j;
      }
      else {
	if(Geno[j][i] == 0) { /* not typed */
	}
	else {
	  if(Geno[j][i] == curgen) {
	    curpos = map[j];
	    icurpos = j;
	  }
	  else {
	    if(type==0) {
	      Location[nseen[i]][i] = (map[j]+curpos)/2.0;

	      if(full_info) {
		Left[nseen[i]][i] = curpos;
		Right[nseen[i]][i] = map[j];
		iLeft[nseen[i]][i] = icurpos+1;
		iRight[nseen[i]][i] = j+1;
	      }

	      curgen = Geno[j][i];
	      curpos = map[j];
	      icurpos=j;
	      nseen[i]++;
	    }
	    else {
	      number = 0; /* number of XOs; indicates to set Location[] */
	      switch(Geno[j][i]) {
	      case 1:
		switch(curgen) {
		case 2: curgen=1; number=1; break;
		case 3: curgen=1; number=2; break;
		case 4: curgen=1; break;
		case 5: curgen=1; number=1; break;
		} break;
	      case 2:
		switch(curgen) {
		case 1: curgen=2; number=1; break;
		case 3: curgen=2; number=1; break;
		case 4: curgen=2; break;
		case 5: curgen=2; break;
		} break;
	      case 3:
		switch(curgen) {
		case 1: curgen=3; number=2; break;
		case 2: curgen=3; number=1; break;
		case 4: curgen=3; number=1; break;
		case 5: curgen=3; break;
		} break;
	      case 4:
		switch(curgen) {
		case 1: break;
		case 2: break;
		case 3: curgen=2; number=1; break;
		case 5: curgen=2; break;
		} break;
	      case 5:
		switch(curgen) {
		case 1: curgen=2; number=1; break;
		case 2: break;
		case 3: break;
		case 4: curgen=2; break;
		} break;
	      }
	      
	      if(number==1) {
		Location[nseen[i]][i] = (curpos+map[j])/2.0;
		if(full_info) {
		  Left[nseen[i]][i] = curpos;
		  Right[nseen[i]][i] = map[j];
		  iLeft[nseen[i]][i] = icurpos+1;
		  iRight[nseen[i]][i] = j+1;
		}
		nseen[i]++;
	      }
	      else if(number==2) { /* two crossovers in interval: place 1/3 and 2/3 along */
		Location[nseen[i]][i] = (curpos+2.0*map[j])/3.0;
		Location[nseen[i]+1][i] = (2.0*curpos+map[j])/3.0;
		if(full_info) {
		  Left[nseen[i]][i] = Left[nseen[i]+1][i] = curpos; 
		  Right[nseen[i]][i] = Right[nseen[i]+1][i] = map[j];
		  iLeft[nseen[i]][i] = iLeft[nseen[i]+1][i] = icurpos+1; 
		  iRight[nseen[i]][i] = iRight[nseen[i]+1][i] = j+1;
		}
		nseen[i] += 2;
	      }
	      curpos = map[j];
	      icurpos = j;
	    }
	  }
	}
      }
    } /* end loop over markers */

    /* count number of typed markers between adjacent crossovers */
    if(full_info) {
      for(j=0; j<nseen[i]-1; j++) {
	nTyped[j][i] = 0;
	for(k=iRight[j][i]-1; k<=iLeft[j+1][i]-1; k++) 
	  if(Geno[k][i] != 0) nTyped[j][i]++;
      }
    }

  } /* end loop over individuals */
}
	  
/* multiply two matrices - I'm using dgemm from lapack here */
/*
void matmult2(double *result, double *a, int nrowa,
             int ncola, double *b, int ncolb)

{ 
  double alpha=1.0, beta=1.0;
  F77_CALL(dgemm)("N", "N", &nrowa, &ncolb, &ncola, &alpha, a, &nrowa,
             b, &ncola, &beta, result, &nrowa);
}
*/

void matmult(double *result, double *a, int nrowa,
             int ncola, double *b, int ncolb)

{
  int i, j, k;

  for(i=0; i<nrowa; i++)
    for(j=0; j<ncolb; j++) {
      /* clear the content of result */
      result[j*nrowa+i] = 0.0;
      /*result[i*ncolb+j] = 0.0;*/
      for(k=0; k<ncola; k++)
        result[j*nrowa+i] += a[k*nrowa+i]*b[j*ncola+k];
    }

}


/**********************************************************************
 * 
 * expandcol2drop
 *
 * Used in scantwo_1chr_em for the X chromosome, to figure out 
 * what columns to drop in the presence of covariates when certain
 * genotype columns must be dropped
 *
 **********************************************************************/

void expand_col2drop(int n_gen, int n_addcov, int n_intcov, 
		     int *col2drop, int *allcol2drop)
{
  int k1, k2, s, ss, j;

  for(k1=0, s=0, ss=0; k1<n_gen; k1++, ss++, s++) 
    allcol2drop[s] = col2drop[ss];


  for(k2=0; k2<n_gen-1; k2++, s++, ss++) 
    allcol2drop[s] = col2drop[ss];

  for(j=0; j<n_addcov; j++, s++) 
    allcol2drop[s]=0;

  for(j=0; j<n_intcov; j++) {
    for(k1=0, ss=0; k1<n_gen-1; k1++, s++, ss++)
      allcol2drop[s] = col2drop[ss];

    ss++;

    for(k2=0; k2<n_gen-1; k2++, s++, ss++)
      allcol2drop[s] = col2drop[ss];
  }

  for(k1=0, ss=2*n_gen-1; k1<n_gen-1; k1++) 
    for(k2=0; k2<n_gen-1; k2++, s++, ss++) 
      allcol2drop[s] = col2drop[ss];

  for(j=0; j<n_intcov; j++)
    for(k1=0, ss=2*n_gen-1; k1<n_gen-1; k1++)
      for(k2=0; k2<n_gen-1; k2++, s++, ss++)
	allcol2drop[s] = col2drop[ss];
}


void dropcol_xpx(int *n_col, int *col2drop, double *xpx)
{
  int i, j, s, n;

  n=0;
  for(i=0, s=0; i< *n_col; i++) {
    if(!col2drop[i]) n++;
    for(j=0; j < *n_col; j++) {
      if(!col2drop[i] && !col2drop[j]) {
	xpx[s] = xpx[j+i*(*n_col)];
	s++;
      }
    }
  }
	
  *n_col = n;
}


void dropcol_xpy(int n_col, int *col2drop, double *xpy)
{
  int i, s;

  for(i=0, s=0; i<n_col; i++) {
    if(!col2drop[i]) {
      xpy[s] = xpy[i];
      s++;
    }
  }
}

void dropcol_x(int *n_col, int n_row, int *col2drop, double *x)
{
  int i, j, n, s;
  
  n=0;
  for(i=0, s=0; i<*n_col; i++) {
    if(!col2drop[i]) {
      n++;
      for(j=0; j<n_row; j++)
	x[j+s*n_row] = x[j+i*n_row];
      s++;
    }
  }
  *n_col = n;
}


/**********************************************************************
 * 
 * reviseMWril    Revise genotypes for 4- or 8-way RIL 
 *                to form encoding the founders' genotypes
 *
 * n_ril     Number of RILs to simulate
 * n_mar     Number of markers
 * n_str     Number of founder strains
 *
 * Parents   SNP data for the founder strains [dim n_str x n_mar]
 * 
 * Geno      On entry, the detailed genotype data; on exit, the 
 *           SNP data written bitwise. [dim n_ril x n_mar]
 * 
 * Crosses   The crosses [n_ril x n_str]
 *
 * missingval  Integer indicating missing value
 *
 **********************************************************************/
void reviseMWril(int n_ril, int n_mar, int n_str, 
		 int **Parents, int **Geno, int **Crosses,
		 int missingval)
{
  int i, j, k, temp;

  for(i=0; i<n_ril; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    for(j=0; j<n_mar; j++) {
      temp = 0;
      for(k=0; k<n_str; k++) {
	if(Geno[j][i] == missingval) Geno[j][i] = 0;
	else if(Parents[j][Crosses[k][i]-1]==missingval ||
		Geno[j][i] == Parents[j][Crosses[k][i]-1])
	  temp += (1 << k);
      }
      Geno[j][i] = temp;
    }
  }
}

/* wrapper for calling reviseMWril from R */
void R_reviseMWril(int *n_ril, int *n_mar, int *n_str, 
		   int *parents, int *geno, int *crosses,
		   int *missingval)
{
  int **Parents, **Geno, **Crosses;

  reorg_geno(*n_str, *n_mar, parents, &Parents);
  reorg_geno(*n_ril, *n_mar, geno, &Geno);
  reorg_geno(*n_ril, *n_str, crosses, &Crosses);

  reviseMWril(*n_ril, *n_mar, *n_str, Parents, Geno, Crosses,
	      *missingval);
}


/* wrapper for calcPermPval */
void R_calcPermPval(double *peaks, int *nc_peaks, int *nr_peaks,
		    double *perms, int *n_perms, double *pval)
{
  double **Peaks, **Perms, **Pval;
  
  reorg_errlod(*nr_peaks, *nc_peaks, peaks, &Peaks);
  reorg_errlod(*n_perms, *nc_peaks, perms, &Perms);
  reorg_errlod(*nr_peaks, *nc_peaks, pval, &Pval);

  calcPermPval(Peaks, *nc_peaks, *nr_peaks, Perms, *n_perms, Pval);
}

/* calculate permutation p-values for summary.scanone() */
void calcPermPval(double **Peaks, int nc_peaks, int nr_peaks,
		  double **Perms, int n_perms, double **Pval)
{
  int i, j, k, count;

  for(i=0; i<nc_peaks; i++) {
    for(j=0; j<nr_peaks; j++) {
      count = 0;
      for(k=0; k<n_perms; k++) 
	if(Perms[i][k] >= Peaks[i][j]) count++;
      Pval[i][j] = (double)count/(double)n_perms;
    }
  }
}


/* end of util.c */
