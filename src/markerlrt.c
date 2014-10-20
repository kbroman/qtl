/**********************************************************************
 *
 * markerlrt.c
 *
 * copyright (c) 2010, Karl W Broman
 *
 * last modified Jul, 2010
 * first written Jul, 2010
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
 * These functions are for performing a general likelihood ratio test for
 * each pair of markers, to assess their association.
 *
 * Contains: R_markerlrt, markerlrt
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "markerlrt.h"
#include "util.h"

void R_markerlrt(int *n_ind, int *n_mar, int *geno, int *maxg,
                 double *lod)
{
    int **Geno;
    double **Lod;

    reorg_geno(*n_ind, *n_mar, geno, &Geno);
    reorg_errlod(*n_mar, *n_mar, lod, &Lod);

    markerlrt(*n_ind, *n_mar, Geno, *maxg, Lod);
}


void markerlrt(int n_ind, int n_mar, int **Geno, int maxg, double **Lod)
{
    int i, j, k, k2, *n1, *n2, **n12, n;

    allocate_imatrix(maxg, maxg, &n12);
    allocate_int(maxg, &n1);
    allocate_int(maxg, &n2);

    for(i=0; i<n_mar; i++) {

        /* count no. typed individuals */
        for(k=0, n=0; k<n_ind; k++)
            if(Geno[i][k]) n++;
        Lod[i][i] = n;

        for(j=(i+1); j<n_mar; j++) {

            /* zero out */
            for(k=0; k<maxg; k++) {
                n1[k] = n2[k] = 0;
                for(k2=0; k2<maxg; k2++) n12[k][k2] = 0;
            }
            n=0;

            /* get counts */
            for(k=0; k<n_ind; k++) {
                if(Geno[i][k] && Geno[j][k]) {
                    n1[Geno[i][k]-1]++;
                    n2[Geno[j][k]-1]++;
                    n12[Geno[i][k]-1][Geno[j][k]-1]++;
                    n++;
                }
            }

            /* calculate LOD score */
            Lod[i][j] = 0;
            for(k=0; k<maxg; k++) {
                for(k2=0; k2<maxg; k2++) {
                    if(n12[k][k2]) {
                        Lod[i][j] += n12[k][k2]*(log10((double)n12[k][k2]) + log10((double)n) -
                                                 log10((double)n1[k]) - log10((double)n2[k2]));
                    }
                }
            }
            Lod[j][i] = Lod[i][j];

        } /* loop over marker pairs */
    }
}


/* end of markerlrt.c */
