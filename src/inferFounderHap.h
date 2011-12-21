/**********************************************************************
 * 
 * inferFounderHap.h
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
 * Contains: inferFounderHap, whichUnique
 *
 * These are for reconstructing the founder haplotypes in inbred lines
 * by a crude method using groups of adjacent SNPs
 *
 **********************************************************************/

void R_inferFounderHap(int *n_snp, int *n_founders, int *n_ind,
		       int *foundergen, int *indgen, int *max_snp,
		       int *hap, int *verbose);

void inferFounderHap(int n_snp, int n_founders, int n_ind, int **founderGen,
		     int **indGen, int max_snp, int **Hap, int verbose);

void whichUnique(int *x, int n_x, int *is_unique, int *n_unique);

/* end of inferFounderHap.h */
