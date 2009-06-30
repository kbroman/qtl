/**********************************************************************
 * 
 * findDupMarkers_notexact.h
 *
 * copyright (c) 2009, Karl W Broman
 *
 * last modified Jun, 2009
 * first written Jun, 2009
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
 * This function is for identifying duplicate markers 
 * where the observed genotypes for one marker match those of another marker, 
 * with no observed genotypes for which the other is missing
 *
 * Contains: R_findDupMarkers_notexact, findDupMarkers_notexact
 *  
 **********************************************************************/

void R_findDupMarkers_notexact(int *nind, int *nmar, int *geno,
			       int *order, int *markerloc, 
			       int *adjacent_only, int *result);

void findDupMarkers_notexact(int nind, int nmar, int **Geno,
			     int *order, int *markerloc, 
			     int adjacent_only, int *result);

/* end of findDupMarkers_notexact.h */

