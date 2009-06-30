/**********************************************************************
 * 
 * findDupMarkers_notexact.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "findDupMarkers_notexact.h"

void R_findDupMarkers_notexact(int *nind, int *nmar, int *geno,
			       int *order, int *markerloc, 
			       int *adjacent_only, int *result)
{
  int **Geno;

  reorg_geno(*nind, *nmar, geno, &Geno);

  findDupMarkers_notexact(*nind, *nmar, Geno, order, markerloc,
			  *adjacent_only, result);

}

void findDupMarkers_notexact(int nind, int nmar, int **Geno,
			     int *order, int *markerloc, 
			     int adjacent_only, int *result)
{
  int i, j, oi, oj, k, nna, flag;

  for(i=0; i<nmar-1; i++) {
    oi = order[i]-1;
    for(j=(i+1); j<nmar; j++) {
      oj = order[j]-1;

      if(result[oj] != 0 || 
	 (adjacent_only && abs(markerloc[oi] - markerloc[oj]) > 1)) {
	/* skip */
      }
      else {
	flag = 0;
	for(k=0; k<nind; k++) {
	  if((Geno[oi][k]==0 && Geno[oj][k]!=0) || 
	     (Geno[oi][k]!=0 && Geno[oj][k]!=0 && Geno[oi][k] != Geno[oj][k])) {
	    flag = 1;
	    break;
	  }
	}
	if(!flag) { /* it worked */
	  if(result[oi] != 0) result[oj] = result[oi];
	  else result[oj] = oi+1;
	}
      }
    }
  }
}


/* end of findDupMarkers_notexact.c */

