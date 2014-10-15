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
 * Contains: inferFounderHap, R_inferFounderHap, whichUnique,
 *           restoreMWrilGeno, R_restoreMWrilGeno
 *
 * These are for reconstructing the founder haplotypes in inbred lines
 * by a crude method using groups of adjacent SNPs
 *
 **********************************************************************/

void R_inferFounderHap(int *n_snp, int *n_founders, int *n_ind,
                       int *foundergen, int *indgen, int *max_offset,
                       int *hap);

void inferFounderHap(int n_snp, int n_founders, int n_ind, int **founderGen,
                     int **indGen, int max_offset, int **Hap);

void whichUnique(unsigned int *x, int n_x, int *is_unique, int *n_unique);

/**********************************************************************
 *
 * restoreMWrilGeno    Do the reverse of reviseMWril, to get
 *                     genotypes back
 *
 * n_ril     Number of RILs to simulate
 * n_mar     Number of markers
 * n_str     Number of founder strains
 *
 * Parents   SNP data for the founder strains [dim n_str x n_mar]
 *
 * Geno      On exit, the detailed genotype data; on entry, the
 *           SNP data written bitwise. [dim n_ril x n_mar]
 *
 * Crosses   The crosses [n_ril x n_str]
 *
 * missingval  Integer indicating missing value
 *
 **********************************************************************/
void restoreMWrilGeno(int n_ril, int n_mar, int n_str,
                      int **Parents, int **Geno, int **Crosses,
                      int missingval);

/* wrapper for calling restoreMWrilGeno from R */
void R_restoreMWrilGeno(int *n_ril, int *n_mar, int *n_str,
                        int *parents, int *geno, int *crosses,
                        int *missingval);

/* end of inferFounderHap.h */
