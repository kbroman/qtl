/**********************************************************************
 *
 * mqmprob.cpp
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Apr, 2009
 * first written Feb, 2009
 *
 * Original version R.C Jansen
 * first written before 2000
 *
 * Much modified by Danny Arends and Pjotr Prins (c) 2009
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
 * Probability functions used by the MQM algorithm
 *
 **********************************************************************/

#include "mqm.h"

/* Chooses the starting probability (when a marker is the first, or unlinked)
 * based on the experimental cross type; used by the augmentation and mixture
 * methods
 *
 * The start probabilities are:
 *
 *          MAA    MBB   MAB
 *   F2     1/4    1/4   1/2
 *   RIL    1/2    1/2     x
 *   BC     1/2      x   1/2
 */

double start_prob(const MQMCrossType crosstype, const char markertype) {
  switch (crosstype) {
    case CF2:
      return (markertype==MH ? 0.5 : 0.25);
      break;
    case CRIL:
      if (markertype==MH) warning("Strange: start_prob function trying to find H in RIL");
      return 0.5;
      break;
    case CBC:
      if (markertype==MBB) warning("Strange: start_prob function trying to find BB in BC");
      return (markertype==MH ? 0.5 : 0.5);
      break;
    default:
      warning("Strange: unknown crosstype in start_prob");
  }
  error("Should get here");
  return NAN;
}

/*
 * Return probability comparing loci[j][i] versus loci[j+1][i],
 * OR if markertype is set, loci[j][i] versus comparemarkertype. Probability
 * is simply based on the markertypes and known recombination rate.
 *
 * Specify an ADJ to adjust loci[j][i] to a specific location in the r[j+ADJ]
 */

double prob(const cmatrix loci, const vector rs, const int i, const int j, const
char checkmarker, const MQMCrossType crosstype, const int ADJ) {
  char compareto;

  const double r = rs[j+ADJ];
  const double r2 = r*r;
  const double rr = 1.0-r; // right side recombination frequency
  const double rr2 = rr*rr;
  const char markertype = loci[j][i];
  // From this point on we don't use loci, r, i, j, ADJ(!)

  if (checkmarker != MUNUSED) {
    compareto = checkmarker;
  } else {
    error("We never get here, all calls pass in the markertype");
    compareto = loci[j+1][i];
  }
  // number of recombinations recombinations
  const double recombinations = fabs((double)markertype-(double)compareto);
  double probj = rr;  // default to no recombinations (1-r)
  switch (crosstype) {
    case CF2:
      if ((markertype==MH)&&(compareto==MH)) {
        probj = r2 + rr2;
      } else if (recombinations==0) {
        probj = rr2;
      } else if (recombinations==1) {
        if (ADJ!=0) {  // FIXME: now this is not clear to me
          probj = ((markertype==MH) ? 2.0*r*rr : r*rr);
        } else {
          probj = ((compareto==MH) ? 2.0*r*rr : r*rr);
        }
      } else {
        probj = r2;  // two recombinations
      }
      break;
    case CRIL:
      if (compareto==MH) {
        error("Strange: prob function trying to find H in RIL");
        return 0.0; // No chance finding a 1 or H in an RIL
      }
      if (recombinations) probj = r;
      break;
    case CBC:
      if (compareto==MBB) {
        error("Strange: prob function trying to find BB in BC");
        return 0.0; // No chance finding a 2/BB in a BC
      }
      if (recombinations) probj = r;
      break;
    default:
      warning("Strange: unknown crosstype in start_prob");
  }
  return probj;
}

/*
 * Return the probability of a marker being of markertype (at j), using the
 * information from the right neighbouring marker and the known recombination
 * frequencies. This function is used by augmentation.
 */

double probright(const char markertype, const int j, const cvector imarker, const vector rs, const cvector position, const MQMCrossType crosstype) {
  double prob0, prob1, prob2;
  if ((position[j]==MRIGHT)||(position[j]==MUNLINKED)) {
    //We're at the end of a chromosome or an unlinked marker
    return 1.0;
  }
  const double r = rs[j];
  const double r2 = r*r;
  const double rr = 1.0-r; // right side recombination frequency
  const double rr2 = rr*rr;
  const char rightmarker = imarker[j+1];

  // markertype markerr diff recombinations
  //   AA        AA      0     0      1-r
  //   AA        BB     -2     2       r
  //   BB        BB      0     0      1-r
  const double recombinations = fabs(markertype-rightmarker);
  switch (crosstype) {
    case CF2:
      if ((rightmarker==MAA)||(rightmarker==MH)||(rightmarker==MBB)) {
        //NEXT marker is known
        if ((markertype==MH)&&(rightmarker==MH)) {
          //special case in which we observe a H after an H then we can't know if we recombinated or not
          return r2+rr2;
        } else {
          //The number of recombinations between observed marker and the next marker
          if (recombinations==0) {
            //No recombination
            return rr2;
          } else if (recombinations==1) {
            if (rightmarker==MH) {
              //the chances of having a H after 1 recombination are 2 times the chance of being either A or B
              return 2.0*r*rr;
            } else {
              //Chance of 1 recombination
              return r*rr;
            }
          } else {
            //Both markers could have recombinated which has a very low chance
            return r2;
          }
        }
      } else if (rightmarker==MNOTAA) {
        //SEMI unknown next marker known is it is not an A
        if (markertype==MAA) {
          //Observed marker is an A
          prob1= 2.0*r*rr;
          prob2= r2;
        } else if (markertype==MH) {
          //Observed marker is an H
          prob1= r2+rr2;
          prob2= r*rr;
        } else {
          //Observed marker is an B
          prob1= 2.0*r*rr;
          prob2= rr2;
        }
        return prob1*probright(MH, j+1, imarker, rs, position, crosstype) + prob2*probright(MBB, j+1, imarker, rs, position, crosstype);
      } else if (rightmarker==MNOTBB) {
        //SEMI unknown next marker known is it is not a B
        if (markertype==MAA) {
          //Observed marker is an A
          prob0= rr2;
          prob1= 2.0*r*rr;
        } else if (markertype==MH) {
          //Observed marker is an H
          prob0= r*rr;
          prob1= r2+rr2;
        } else {
          //Observed marker is an B
          prob0= r2;
          prob1= 2.0*r*rr;
        }
        return prob0*probright(MAA, j+1, imarker, rs, position, crosstype) + prob1*probright(MH, j+1, imarker, rs, position, crosstype);
      } else {
        // Unknown next marker so estimate all posibilities (imarker[j+1]==MMISSING)
        if (markertype==MAA) {
          //Observed marker is an A
          prob0= rr2;
          prob1= 2.0*r*rr;
          prob2= r2;
        } else if (markertype==MH) {
          //Observed marker is an H
          prob0= r*rr;
          prob1= r2+rr2;
          prob2= r*rr;
        } else {
          //Observed marker is an B
          prob0= r2;
          prob1= 2.0*r*rr;
          prob2= rr2;
        }
        return prob0*probright(MAA, j+1, imarker, rs, position, crosstype) + prob1*probright(MH, j+1, imarker, rs, position, crosstype) + prob2*probright(MBB, j+1, imarker, rs, position, crosstype);
      }
      break;
    case CRIL:
      if (markertype==MH) {
        error("Strange: encountered heterozygous genotype in RIL");
        return 0.0;
      }
      if ((rightmarker==MAA)||(rightmarker==MBB)) {
        if (recombinations==0) {
          return rr;  // no recombination (1-r) 
        } else {
          return r;   // recombination (r)
        }
      } else {
        // [pjotr:] I think this code is never reached (FIXME)
        // Both markers could have recombinated which has a very low chance
        warning("Unreachable code");
        if (markertype==MAA) {
          prob0= rr;
          prob2= r;
        } else { // MBB
          prob0= r;
          prob2= rr;
        }
        return prob0*probright(MAA, j+1, imarker, rs, position, crosstype) + prob2*probright(MBB, j+1, imarker, rs, position, crosstype);
      }
      break;
    case CBC:
      if (markertype==MBB) {
        error("Strange: encountered BB genotype in BC");
        return 0.0;
      }
      if ((rightmarker==MAA)||(rightmarker==MH)) {
        if (recombinations==0) {
          return rr;
        } else {
          return r;
        }
      } else {
        // [pjotr:] I think this code is never reached (FIXME)
        // Both markers could have recombinated which has a very low chance
        warning("Unreachable code");
        if (markertype==MAA) {
          prob0= rr;
          prob2= r;
        } else {
          prob0= r;
          prob2= rr;
        }
        return prob0*probright(MAA, j+1, imarker, rs, position, crosstype) + prob2*probright(MH, j+1, imarker, rs, position, crosstype);
      }
      break;
    default:
      warning("Strange: unknown crosstype in start_prob");
  }
  return 1.0;
}

