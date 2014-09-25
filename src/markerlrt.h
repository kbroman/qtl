/**********************************************************************
 *
 * markerlrt.h
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

void R_markerlrt(int *n_ind, int *n_mar, int *geno, int *maxg,
                 double *lod);


void markerlrt(int n_ind, int n_mar, int **Geno, int maxg, double **Lod);

/* end of markerlrt.h */
