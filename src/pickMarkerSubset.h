/**********************************************************************
 *
 * pickMarkerSubset.h
 *
 * copyright (c) 2011, Karl W Broman
 *
 * last modified Nov, 2011
 * first written Nov, 2011
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
 * These functions is for selecting the largest set of markers for
 * which adjacent markers are min_distance apart.
 *
 * Contains: R_pickMarkerSubset, pickMarkerSubset
 *
 **********************************************************************/

void R_pickMarkerSubset(double *locations, int *n_locations, double *weights,
                        double *min_distance, int *path, int *n_path);

void pickMarkerSubset(double *locations, int n_locations, double *weights,
                      double min_distance, int *path, int *n_path);

/* end of pickMarkerSubset.h */
