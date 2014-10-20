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
                       int *foundergen, int *indgen, int *max_offset,
                       int *hap)
{
    int **founderGen, **indGen, **Hap;
    reorg_geno(*n_founders, *n_snp, foundergen, &founderGen);
    reorg_geno(*n_ind, *n_snp, indgen, &indGen);
    reorg_geno(*n_ind, *n_snp, hap, &Hap);

    inferFounderHap(*n_snp, *n_founders, *n_ind, founderGen, indGen, *max_offset,
                    Hap);
}



void inferFounderHap(int n_snp, int n_founders, int n_ind, int **founderGen,
                     int **indGen, int max_offset, int **Hap)
{
    int i, j, left, offset, n_unique;
    unsigned int *fhap, *indhap;
    int *fhapunique;

    allocate_uint(n_founders, &fhap);
    allocate_int(n_founders, &fhapunique);
    allocate_uint(n_ind, &indhap);

    for(left=0; left<n_snp; left++) {
        for(i=0; i<n_founders; i++) fhap[i] = 0;
        for(i=0; i<n_ind; i++) indhap[i] = 0;

        for(offset=0; offset<max_offset && left+offset<n_snp && left-offset>=0; offset++) {
            R_CheckUserInterrupt(); /* check for ^C */

            /* founder haplotypes as integers */
            for(i=0; i<n_founders; i++) {
                if(founderGen[left+offset][i])
                    fhap[i] += (1 << (offset*2));
                if(offset > 0 && founderGen[left-offset][i])
                    fhap[i] += (1 << (offset*2+1));
            }

            /* individual haplotypes as integers */
            for(i=0; i<n_ind; i++) {
                if(Hap[left][i] == 0) { /* haven't figured this one out yet */
                    if(indGen[left+offset][i] < 0 ||
                       (offset > 0 && indGen[left-offset][i] < 0))  /* missing genotype */
                        Hap[left][i] = -1;
                    else {
                        if(indGen[left+offset][i])
                            indhap[i] += (1 << (offset*2));
                        if(offset > 0 && indGen[left-offset][i])
                            indhap[i] += (1 << (offset*2+1));
                    }
                }
            }

            /* which founder haplotypes are unique */
            whichUnique(fhap, n_founders, fhapunique, &n_unique);

            if(n_unique>0) { /* at least one unique founder haplotype */
                for(i=0; i<n_ind; i++) {
                    if(Hap[left][i] == 0) { /* haven't figured this one out yet */
                        for(j=0; j<n_founders; j++) {
                            if(fhapunique[j] && fhap[j]==indhap[i]) {
                                Hap[left][i] = j+1;
                            }
                        }
                    }
                }
            }

            if(n_unique == n_founders)
                break; /* stop extending; go to next SNP */

        }
    }
}


/* for vector x, identify which elements are unique (is_unique -> 1)
   and then count the unique ones */
void whichUnique(unsigned int *x, int n_x, int *is_unique, int *n_unique)
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
                      int missingval)
{
    int i, j, k;

    for(i=0; i<n_ril; i++) {
        R_CheckUserInterrupt(); /* check for ^C */

        for(j=0; j<n_mar; j++) {
            if(Geno[j][i] == 0) Geno[j][i] = missingval;
            else {
                for(k=0; k<n_str; k++) {
                    if(Parents[j][Crosses[k][i]-1]!=missingval) {
                        if(Geno[j][i] & (1 << k))
                            Geno[j][i] = Parents[j][Crosses[k][i]-1];
                        else
                            Geno[j][i] = 1-Parents[j][Crosses[k][i]-1];
                        break;
                    }
                }
            }
        }
    }
}

/* wrapper for calling restoreMWrilGeno from R */
void R_restoreMWrilGeno(int *n_ril, int *n_mar, int *n_str,
                        int *parents, int *geno, int *crosses,
                        int *missingval)
{
    int **Parents, **Geno, **Crosses;

    reorg_geno(*n_str, *n_mar, parents, &Parents);
    reorg_geno(*n_ril, *n_mar, geno, &Geno);
    reorg_geno(*n_ril, *n_str, crosses, &Crosses);

    restoreMWrilGeno(*n_ril, *n_mar, *n_str, Parents, Geno, Crosses,
                     *missingval);
}

/* end of inferFounderHap.c */
