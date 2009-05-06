/**********************************************************************
 * 
 * ripple.c
 *
 * copyright (c) 2002-9, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Mar, 2002
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
 * These functions are for comparing marker orders by counts of 
 * obligate crossovers
 *
 * Contains: R_ripple_bc, R_ripple_f2, R_ripple_4way, ripple, 
 *           countxo_bc, countxo_f2, countxo_4way
 *           R_ripple_ril48, countxo_ril48
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
#include "ripple.h"

/**********************************************************************
 * 
 * ripple
 *
 * This function inspects each of a set of marker orders and counts 
 * the number of obligate crossovers for each order.
 *
 * Input:
 *
 * n_ind    = no. individuals
 *
 * n_mar    = no. markers
 *
 * geno     = genotype data [n_ind x n_mar]
 *
 * n_orders = no. different orders
 *
 * orders   = matrix of marker orders [n_orders x n_mar]
 *            (Note: each row contains {0, 1, ..., n_mar-1}
 *
 * nxo      = the output; the number of obligate crossovers for each
 *            order (a vector of length n_orders)
 *
 * print_by = How often to print out information?
 *
 * countxo  = function to count the number of obligate crossovers in
 *            an interval and to update the current inferred genotype
 *            (specific for backcross, intercross, and four-way cross)
 *
 **********************************************************************/

void ripple(int n_ind, int n_mar, int *geno,
	    int n_orders, int *orders, int *nxo, 
	    int print_by, int countxo(int *curgen, int nextgen))
{
  int **Geno, **Orders;
  int i, j, k, curgen;

  /* reorganize genotype data and marker order matrix */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  reorg_geno(n_orders, n_mar, orders, &Orders);
  
  for(i=0; i<n_orders; i++) { /* loop over possible orders */

    R_CheckUserInterrupt(); /* check for ^C */

    /*    if(((int)((i+1)/print_by))*print_by == i+1) 
	  printf("    --Order %d\n", i+1); */

    nxo[i] = 0;
    for(j=0; j<n_ind; j++) { /* loop over individuals */
      /* genotype at first marker */
      curgen = Geno[Orders[0][i]][j];
      for(k=1; k<n_mar; k++) /* loop over markers */
	/* count no obligate crossovers and update current genotype */
	nxo[i] += countxo(&curgen, Geno[Orders[k][i]][j]);
    }
  }

}


/**********************************************************************
 * 
 * R_ripple_bc
 *
 * Wrapper for call from R for a backcross
 * 
 **********************************************************************/

void R_ripple_bc(int *n_ind, int *n_mar, int *geno, 
		 int *n_orders, int *orders,
		 int *nxo, int *print_by)
{
  ripple(*n_ind, *n_mar, geno, *n_orders, orders, nxo,
	 *print_by, countxo_bc);
}


/**********************************************************************
 * 
 * countxo_bc
 * 
 * count no. obligate crossovers in a backcross
 *
 **********************************************************************/

int countxo_bc(int *curgen, int nextgen) 
{
  if(*curgen == 0) { /* missing data */
    *curgen=nextgen; 
    return(0);
  }
  else { 
    if(nextgen == 0) return(0); /* missing data */
    else {
      if(nextgen != *curgen) {
	*curgen = nextgen;
	return(1); /* crossover */
      }
      else return(0); /* no crossover */
    }
  }
}


/**********************************************************************
 * 
 * R_ripple_f2
 *
 * Wrapper for call from R for an intercross
 * 
 **********************************************************************/

void R_ripple_f2(int *n_ind, int *n_mar, int *geno, 
		 int *n_orders, int *orders,
		 int *nxo, int *print_by)
{
  ripple(*n_ind, *n_mar, geno, *n_orders, orders, nxo,
	 *print_by, countxo_f2);
}


/**********************************************************************
 * 
 * countxo_f2
 * 
 * count no. obligate crossovers in an intecross
 *
 **********************************************************************/

int countxo_f2(int *curgen, int nextgen) 
{
  if(nextgen == 0) return(0);

  switch(*curgen) {
  case 0: /* missing data */
    *curgen=nextgen; 
    return(0);
  case 1: 
    switch(nextgen) {
    case 1: return(0);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(2);
    case 4: return(0);
    case 5: *curgen=2; return(1);
    }
  case 2:
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: return(0);
    case 3: *curgen=3; return(1);
    case 4: return(0);
    case 5: return(0);
    }
  case 3:
    switch(nextgen) {
    case 1: *curgen=1; return(2);
    case 2: *curgen=2; return(1);
    case 3: return(0);
    case 4: *curgen=2; return(1);
    case 5: return(0);
    }
  case 4:
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(1);
    case 4: return(0);
    case 5: *curgen=2; return(0);
    }
  case 5:
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(0);
    case 4: *curgen=2; return(0);
    case 5: return(0);
    }
  default: return(0); /* shouldn't get here! */
  }
  
}


/**********************************************************************
 * 
 * R_ripple_4way
 *
 * Wrapper for call from R for a four-way cross
 * 
 **********************************************************************/

void R_ripple_4way(int *n_ind, int *n_mar, int *geno, 
		   int *n_orders, int *orders,
		   int *nxo, int *print_by)
{
  ripple(*n_ind, *n_mar, geno, *n_orders, orders, nxo,
	 *print_by, countxo_4way);
}


/**********************************************************************
 * 
 * countxo_4way
 * 
 * count no. obligate crossovers in a four-way cross
 *
 **********************************************************************/

int countxo_4way(int *curgen, int nextgen) 
{
  if(nextgen == 0) return(0);

  switch(*curgen) {
  case 0: /* missing data */
    *curgen=nextgen; 
    return(0);
  case 1: 
    switch(nextgen) {
    case 1: return(0);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(1);
    case 4: *curgen=4; return(2);

    case 5: return(0);
    case 6: *curgen=2; return(1);
    case 7: return(0);
    case 8: *curgen=3; return(1);

    case 9: return(0);
    case 10: *curgen=10; return(1);

    case 11: *curgen=10; return(1);
    case 12: return(0);
    case 13: return(0);
    case 14: return(0);
    }
  case 2: 
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: return(0);
    case 3: *curgen=3; return(2);
    case 4: *curgen=4; return(1);

    case 5: *curgen=1; return(1);
    case 6: return(0);
    case 7: return(0);
    case 8: *curgen=4; return(1);

    case 9: *curgen=9; return(1);
    case 10: return(0);

    case 11: return(0);
    case 12: *curgen=9; return(1);
    case 13: return(0);
    case 14: return(0);
    }
  case 3: 
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: *curgen=2; return(2);
    case 3: return(0);
    case 4: *curgen=4; return(1);

    case 5: return(0);
    case 6: *curgen=4; return(1);
    case 7: *curgen=1; return(1);
    case 8: return(0);

    case 9: *curgen=9; return(1);
    case 10: return(0);

    case 11: return(0);
    case 12: return(0);
    case 13: *curgen=9; return(1);
    case 14: return(0);
    }
  case 4: 
    switch(nextgen) {
    case 1: *curgen=1; return(2);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(1);
    case 4: return(0);

    case 5: *curgen=3; return(1);
    case 6: return(0);
    case 7: *curgen=2; return(1);
    case 8: return(0);

    case 9: return(0);
    case 10: *curgen=10; return(1);

    case 11: return(0);
    case 12: return(0);
    case 13: return(0);
    case 14: *curgen=10; return(1);
    }
  case 5: 
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(0);
    case 4: *curgen=4; return(1);

    case 5: return(0);
    case 6: *curgen=6; return(1);
    case 7: *curgen=1; return(0);
    case 8: *curgen=3; return(0);

    case 9: *curgen=1; return(0);
    case 10: *curgen=3; return(0);

    case 11: *curgen=3; return(0);
    case 12: return(0);
    case 13: *curgen=1; return(0);
    case 14: return(0);
    }
  case 6: 
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(1);
    case 4: *curgen=4; return(0);

    case 5: *curgen=5; return(1);
    case 6: return(0);
    case 7: *curgen=2; return(0);
    case 8: *curgen=4; return(0);

    case 9: *curgen=4; return(0);
    case 10: *curgen=2; return(0);

    case 11: return(0);
    case 12: *curgen=4; return(0);
    case 13: return(0);
    case 14: *curgen=2; return(0);
    }
  case 7: 
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(1);
    case 4: *curgen=4; return(1);

    case 5: *curgen=1; return(0);
    case 6: *curgen=2; return(0);
    case 7: return(0);
    case 8: *curgen=8; return(1);

    case 9: *curgen=1; return(0);
    case 10: *curgen=2; return(0);

    case 11: *curgen=2; return(0);
    case 12: *curgen=1; return(0);
    case 13: return(0);
    case 14: return(0);
    }
  case 8: 
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(0);
    case 4: *curgen=4; return(0);

    case 5: *curgen=3; return(0);
    case 6: *curgen=4; return(0);
    case 7: *curgen=7; return(1);
    case 8: return(0);

    case 9: *curgen=4; return(0);
    case 10: *curgen=3; return(0);

    case 11: return(0);
    case 12: return(0);
    case 13: *curgen=4; return(0);
    case 14: *curgen=3; return(0);
    }
  case 9: 
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(1);
    case 4: *curgen=4; return(0);

    case 5: *curgen=1; return(0);
    case 6: *curgen=4; return(0);
    case 7: *curgen=1; return(0);
    case 8: *curgen=4; return(0);

    case 9: *curgen=9; return(0);
    case 10: *curgen=10; return(1);

    case 11: *curgen=4; return(0);
    case 12: return(0);
    case 13: return(0);
    case 14: *curgen=1; return(0);
    }
  case 10: 
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(0);
    case 4: *curgen=4; return(1);

    case 5: *curgen=3; return(0);
    case 6: *curgen=2; return(0);
    case 7: *curgen=2; return(0);
    case 8: *curgen=3; return(0);

    case 9: *curgen=9; return(1);
    case 10: return(0);

    case 11: return(0);
    case 12: *curgen=3; return(0);
    case 13: *curgen=2; return(0);
    case 14: return(0);
    }
  case 11: 
    switch(nextgen) {
    case 1: *curgen=1; return(1);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(0);
    case 4: *curgen=4; return(0);

    case 5: *curgen=3; return(0);
    case 6: *curgen=6; return(0);
    case 7: *curgen=2; return(0);
    case 8: *curgen=8; return(0);

    case 9: *curgen=4; return(0);
    case 10: *curgen=10; return(0);

    case 11: return(0);
    case 12: *curgen=8; return(0);
    case 13: *curgen=6; return(0);
    case 14: *curgen=10; return(0);
    }
  case 12: 
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(1);
    case 3: *curgen=3; return(0);
    case 4: *curgen=4; return(0);

    case 5: *curgen=5; return(0);
    case 6: *curgen=4; return(0);
    case 7: *curgen=1; return(0);
    case 8: *curgen=8; return(0);

    case 9: *curgen=9; return(0);
    case 10: *curgen=3; return(0);

    case 11: *curgen=8; return(0);
    case 12: return(0);
    case 13: *curgen=9; return(0);
    case 14: *curgen=5; return(0);
    }
  case 13: 
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(1);
    case 4: *curgen=4; return(0);

    case 5: *curgen=1; return(0);
    case 6: *curgen=6; return(0);
    case 7: *curgen=7; return(0);
    case 8: *curgen=4; return(0);

    case 9: *curgen=9; return(0);
    case 10: *curgen=2; return(0);

    case 11: *curgen=6; return(0);
    case 12: *curgen=9; return(0);
    case 13: return(0);
    case 14: *curgen=7; return(0);
    }
  case 14: 
    switch(nextgen) {
    case 1: *curgen=1; return(0);
    case 2: *curgen=2; return(0);
    case 3: *curgen=3; return(0);
    case 4: *curgen=4; return(1);

    case 5: *curgen=5; return(0);
    case 6: *curgen=2; return(0);
    case 7: *curgen=7; return(0);
    case 8: *curgen=3; return(0);

    case 9: *curgen=1; return(0);
    case 10: *curgen=10; return(0);

    case 11: *curgen=10; return(0);
    case 12: *curgen=5; return(0);
    case 13: *curgen=7; return(0);
    case 14: return(0);
    }
  default: return(0); /* shouldn't get here */
  }
}


/**********************************************************************
 * 
 * R_ripple_ril48
 *
 * Wrapper for call from R for 4- or 8-way RIL
 * 
 **********************************************************************/

void R_ripple_ril48(int *n_ind, int *n_mar, int *geno, 
		    int *n_orders, int *orders,
		    int *nxo, int *print_by)
{
  ripple(*n_ind, *n_mar, geno, *n_orders, orders, nxo,
	 *print_by, countxo_ril48);
}


/**********************************************************************
 * 
 * countxo_ril48
 * 
 * count no. obligate crossovers in 4- or 8-way RIL
 *
 **********************************************************************/

int countxo_ril48(int *curgen, int nextgen) 
{
  int and;

  if(nextgen == 0) return(0);

  and = *curgen & nextgen;
  if(and == 0) {
    *curgen = nextgen; 
    return(1);
  }
  else {
    *curgen = and;
    return(0);
  }
}

/* end of ripple.c */

