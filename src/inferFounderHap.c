/**********************************************************************
 * 
 * inferFounderHap.c
 * 
 * copyright (c) 2011, Karl W Broman
 *
 * last modified Dec, 2011
 * first written Dec, 2011
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
 * Contains: constructFounderHap, whichUnique
 *
 * These are for reconstructing the founder haplotypes in inbred lines
 * by a crude method using groups of adjacent SNPs
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "inferFounderHap.h"
#include "util.h"

void R_inferFounderHap(int *n_snp, int *n_founders, int *n_ind,
		       int *foundergen, int *indgen, int *max_snp,
		       int *hap, int *verbose)
{
  int **founderGen, **indGen, **Hap;
  reorg_geno(*n_founders, *n_snp, foundergen, &founderGen);
  reorg_geno(*n_ind, *n_snp, indgen, &indGen);
  reorg_geno(*n_ind, *n_snp, hap, &Hap);

  inferFounderHap(*n_snp, *n_founders, *n_ind, founderGen, indGen, *max_snp,
		  Hap, *verbose);
}
  
  

void inferFounderHap(int n_snp, int n_founders, int n_ind, int **founderGen,
		     int **indGen, int max_snp, int **Hap, int verbose)
{
  int i, j, k, flag, left, offset, n_unique;
  int *fhap, *fhapunique, *indhap;

  allocate_int(n_founders, &fhap);
  allocate_int(n_founders, &fhapunique);
  allocate_int(n_ind, &indhap);
  
  for(left=0; left<n_snp; left++) {
    if(verbose) Rprintf("marker %d of %d\n", left+1, n_snp);
    for(i=0; i<n_founders; i++) fhap[i] = 0;
    for(i=0; i<n_ind; i++) indhap[i] = 0;

    for(offset=0; offset<max_snp && left+offset<n_snp; offset++) {
      R_CheckUserInterrupt(); /* check for ^C */

      /* founder haplotypes as integers */
      for(i=0; i<n_founders; i++) {
	if(founderGen[left+offset][i])
	  fhap[i] += (1 << offset);
      }

      /* individual haplotypes as integers */
      for(i=0; i<n_ind; i++) {
	if(Hap[left][i] == 0) { /* haven't figured this one out yet */
	  if(indGen[left+offset][i] < 0)  /* missing genotype */
	    Hap[left][i] = -1;
	  else if(indGen[left+offset][i])
	    indhap[i] += (1 << offset);
	}
#ifdef UNDEFINED	
	if(verbose && offset > 8 && left < 20 && i < 5) {
	  Rprintf("    ind=%d  geno=", i+1);
	  for(k=offset; k>=0; k--) {
	    if(indGen[left+k][i] < 0)
	      Rprintf("9");
	    else
	      Rprintf("%d", indGen[left+k][i]);
	  }
	  Rprintf("  hap=%d  Hap=%d\n", indhap[i], Hap[left][i]);
	}
#endif

      }

      /* which founder haplotypes are unique */
      whichUnique(fhap, n_founders, fhapunique, &n_unique);
      if(verbose) Rprintf("    offset=%d   n_unique=%d\n", offset, n_unique);

      if(n_unique>0) { /* at least one unique founder haplotype */
	flag = 0;
	for(i=0; i<n_ind; i++) {
	  if(Hap[left][i] == 0) { /* haven't figured this one out yet */
	    for(j=0; j<n_founders; j++) {
	      if(fhapunique[j] && fhap[j]==indhap[i]) {
		Hap[left][i] = j+1;

		if(verbose) {
		  if(!flag) {
		    Rprintf("    ");
		    for(k=0; k<n_founders; k++) {
		      Rprintf(" %d=%d", k+1, fhap[k]);
		    }
		    Rprintf("\n");
		    flag = 1;
		  }
		  Rprintf("    ind=%d  gen=", i+1);
		  for(k=offset; k>=0; k--) {
		    if(indGen[left+k][i] < 0)
		      Rprintf("9");
		    else
		      Rprintf("%d", indGen[left+k][i]);
		  }
		  Rprintf("  hap=%d  Hap=%d\n", indhap[i], Hap[left][i]);
		}

	      }
	    }
	  }
	}
      }

      if(n_unique == n_founders) { 
	if(verbose) Rprintf("    break\n");
	break; /* stop extending; go to next SNP */
      }
    }
  }
}


/* for vector x, identify which elements are unique (is_unique -> 1) 
   and then count the unique ones */
void whichUnique(int *x, int n_x, int *is_unique, int *n_unique)
{
  int i, j;

  for(i=0; i<n_x; i++)
    is_unique[i] = 1;

  for(i=0; i<n_x-1; i++) {
    if(is_unique[i]) {
      for(j=i+1; j<n_x; j++) {
	if(is_unique[j]) {
	  if(x[i] == x[j]) {
	    is_unique[i] = is_unique[j] = 0;
	  }
	}
      }
    }
  }
  
  *n_unique = 0;
  for(i=0; i<n_x; i++) 
    *n_unique += is_unique[i];
}

/* end of inferFounderHap.c */
