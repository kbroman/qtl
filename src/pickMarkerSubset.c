/**********************************************************************
 * 
 * pickMarkerSubset.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "pickMarkerSubset.h"

void R_pickMarkerSubset(double *locations, int *n_locations, double *weights,
			double *min_distance, int *path, int *n_path)
{
  GetRNGstate();
  pickMarkerSubset(locations, *n_locations, weights, *min_distance, path, n_path);
  PutRNGstate();
}


void pickMarkerSubset(double *locations, int n_locations, double *weights, 
		      double min_distance, int *path, int *n_path)
{
  int i, j;
  double *total_weights, themax;
  int *prev_marker, *max_to_choose, n_max_to_choose;

  total_weights = (double *)R_alloc(n_locations, sizeof(double));
  prev_marker = (int *)R_alloc(n_locations, sizeof(int));
  max_to_choose = (int *)R_alloc(n_locations, sizeof(int));


  /* first location */
  prev_marker[0] = -1;
  total_weights[0] = weights[0];

  for(i=1; i<n_locations; i++) {
    if(locations[i] < locations[0] + min_distance) { 
      /* no markers to left of i that are > min_distance away */
      total_weights[i] = weights[i];
      prev_marker[i] = -1;
    }
    else {

      /* look for maxima */
      n_max_to_choose = 1;
      max_to_choose[0] = 0;
      themax = total_weights[0];
      for(j=1; j<i; j++) {

        R_CheckUserInterrupt(); /* check for ^C */

        if(locations[i] < locations[j] + min_distance) break;
       
        if(total_weights[j] > themax) {
	  n_max_to_choose = 1;
	  max_to_choose[0] = j;
	  themax = total_weights[j];
	}
	else if(total_weights[j] == themax) {
	  max_to_choose[n_max_to_choose] = j;
	  n_max_to_choose++;
	}
      }
     
      /* now choose among the maxima at random */
      total_weights[i] = themax + weights[i];
      if(n_max_to_choose == 1) prev_marker[i] = max_to_choose[0];
      else /* pick random */
	prev_marker[i] = max_to_choose[(int)(unif_rand()*(double)n_max_to_choose)]; 
    }
  }

  /* now find global max */
  themax = total_weights[0];
  n_max_to_choose = 1;
  max_to_choose[0] = 0;
  
  for(i=1; i<n_locations; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    if(total_weights[j] > themax) {
      themax = total_weights[i];
      n_max_to_choose = 1;
      max_to_choose[0] = i;
    }
    else if(total_weights[i] == themax) {
      max_to_choose[n_max_to_choose] = i;
      n_max_to_choose++;
    }
  }

  /* right-most marker at global maximum */
  if(n_max_to_choose == 1) path[0] = max_to_choose[0];
  else /* pick random */
    path[0] = max_to_choose[(int)(unif_rand()*(double)n_max_to_choose)]; 

  *n_path=1;

  /* trace back */
  while(prev_marker[path[*n_path-1]] > -1) {
    R_CheckUserInterrupt(); /* check for ^C */

    path[*n_path] = prev_marker[path[*n_path-1]];
    (*n_path)++;
  }

  /* note: result is backwards and has indices starting at 0 */
}

/* end of pickMarkerSubset.c */
