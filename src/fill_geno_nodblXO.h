/**********************************************************************
 * 
 * fill_geno_nodblXO.h
 *
 * copyright (c) 2010, Karl W Broman
 *
 * last modified May, 2010
 * first written May, 2010
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
 * This function fills in missing genotype data only between markers with 
 * exactly the same genotype.
 *
 **********************************************************************/

void R_fill_geno_nodblXO(int *n_ind, int *n_mar, int *geno);

void fill_geno_nodblXO(int n_ind, int n_mar, int **Geno);

/* end of fill_geno_nodblXO.h */
