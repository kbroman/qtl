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
 * first written <2000 (unknown)
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

/* Chooses the starting probability based on the experimental cross type;
 * used by the augmentation and mixture methods
 */

double start_prob(const char crosstype, const char markertype) {
  switch (crosstype) {
    case CF2:
      return (markertype==MH ? 0.5 : 0.25);
      break;
    case CRIL:
      return (markertype==MH ? 0.0 : 0.5);
      break;
    case CBC:
      return (markertype==MBB ? 0.0 : 0.5);
      break;
    default:
      warning("Strange: unknown crosstype in start_prob");
  }
  return 0.0;
}

/*
 * Return probability comparing loci[j][i] versus loci[j+1][i],
 * OR if markertype is set, loci[j][i] versus comparemarkertype. Probability
 * is simply based on the markertypes and known recombination rate.
 *
 * Specify an ADJ to adjust loci[j][i] to a specific location in the r[j+ADJ]
 */

double prob(const cmatrix loci, const vector r, const int i, const int j, const
char markertype, const char crosstype, const int ADJ, const int start) {
  double probj =0.0;
  double Nrecom;
  char compareto;
  // Rprintf("Prob called: values:\n(i, j, ADJ)=(%d, %d, %d)\nR[j+ADJ] value: %f Loci[j][i]=%c\n", i, j, ADJ, r[j+ADJ], loci[j][i]);

  if (markertype != MUNUSED) {
    //Rprintf("C %d %d\n", i, j);
    compareto = markertype;
  } else {
    error("We never get here, all calls pass in the markertype");
    //Rprintf("loci[j+1][i] %d\n", j);
    compareto = loci[j+1][i];
  }
  switch (crosstype) {
    case CF2:
      if (start) {
        return (loci[j][i]==MH ? 0.5 : 0.25);
      }
      Nrecom= fabs((double)loci[j][i]-(double)compareto);
      if ((loci[j][i]==MH)&&(compareto==MH)) {
        //Rprintf("SCase %c <-> %c:\n", compareto, loci[j][i]);
        probj = (r[j+ADJ]*r[j+ADJ]+(1.0-r[j+ADJ])*(1.0-r[j+ADJ]));
      } else if (Nrecom==0) {
        //Rprintf("Nrecom=0 %c <-> %c:\n", compareto, loci[j][i]);
        probj = (1.0-r[j+ADJ])*(1.0-r[j+ADJ]);
      } else if (Nrecom==1) {
        //Rprintf("Nrecom=1 %c <-> %c:\n", compareto, loci[j][i]);
        if (ADJ!=0) {
          probj = ((loci[j][i]==MH) ? 2.0*r[j+ADJ]*(1.0-r[j+ADJ]) : r[j+ADJ]*(1.0-r[j+ADJ]));
        } else {
          probj = ((compareto==MH) ? 2.0*r[j+ADJ]*(1.0-r[j+ADJ]) : r[j+ADJ]*(1.0-r[j+ADJ]));
        }
      } else {
        //Rprintf("Nrecom=2 %c <-> %c:\n", compareto, loci[j][i]);
        probj = r[j+ADJ]*r[j+ADJ];
      }
      //Rprintf("after IF\n", j);
      break;
    case CRIL:
      if (start) {
        return 0.5;
      }
      if (compareto==MH) {
        warning("Strange: prob function trying to find H in RIL");
        return 0.0; // No chance finding a 1 or H in an RIL
      }
      Nrecom = fabs((double)loci[j][i]-(double)compareto);
      if (Nrecom==0) {
        //No recombination has a chance of r[j]
        probj = 1.0-r[j+ADJ];
      } else {
        // Recombination between markers has a chance of r[j-1]
        probj = r[j+ADJ];
      }
      break;
    case CBC:
      if (start) {
        return 0.5;
      }
      if (compareto==MBB) {
        warning("Strange: prob function trying to find BB in BC");
        return 0.0; // No chance finding a 2/BB in a BC
      }
      Nrecom= fabs((double)loci[j][i]-(double)compareto);
      if (Nrecom==0) {
        //No recombination has a chance of r[j]
        probj =  (1.0-r[j+ADJ]);
      } else {
        // Recombination between markers has a chance of r[j-1]
        probj = r[j+ADJ];
      }
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

double probright(const char markertype, const int j, const cvector imarker, const vector r, const cvector position, const char crosstype) {
  double nrecom, prob0, prob1, prob2;
  if ((position[j]==MRIGHT)||(position[j]==MUNLINKED)) {
    //We're at the end of a chromosome or an unlinked marker
    return 1.0;
  }
  const double rj = r[j];
  const char rightmarker = imarker[j+1];

  switch (crosstype) {
    case CF2:
      if ((rightmarker==MAA)||(rightmarker==MH)||(rightmarker==MBB)) {
        //NEXT marker is known
        if ((markertype==MH)&&(rightmarker==MH)) {
          //special case in which we observe a H after an H then we can't know if we recombinated or not
          return rj*rj+(1.0-rj)*(1.0-rj);
        } else {
          //The number of recombinations between observed marker and the next marker
          nrecom = fabs(markertype-rightmarker);
          if (nrecom==0) {
            //No recombination
            return (1.0-rj)*(1.0-rj);
          } else if (nrecom==1) {
            if (rightmarker==MH) {
              //the chances of having a H after 1 recombination are 2 times the chance of being either A or B
              return 2.0*rj*(1.0-rj);
            } else {
              //Chance of 1 recombination
              return rj*(1.0-rj);
            }
          } else {
            //Both markers could have recombinated which has a very low chance
            return rj*rj;
          }
        }
      } else if (rightmarker==MNOTAA) {
        //SEMI unknown next marker known is it is not an A
        if (markertype==MAA) {
          //Observed marker is an A
          prob1= 2.0*rj*(1.0-rj);
          prob2= rj*rj;
        } else if (markertype==MH) {
          //Observed marker is an H
          prob1= rj*rj+(1.0-rj)*(1.0-rj);
          prob2= rj*(1.0-rj);
        } else {
          //Observed marker is an B
          prob1= 2.0*rj*(1.0-rj);
          prob2= (1.0-rj)*(1-rj);
        }
        return prob1*probright(MH, j+1, imarker, r, position, crosstype) + prob2*probright(MBB, j+1, imarker, r, position, crosstype);
      } else if (rightmarker==MNOTBB) {
        //SEMI unknown next marker known is it is not a B
        if (markertype==MAA) {
          //Observed marker is an A
          prob0= (1.0-rj)*(1.0-rj);
          prob1= 2.0*rj*(1.0-rj);
        } else if (markertype==MH) {
          //Observed marker is an H
          prob0= rj*(1.0-rj);
          prob1= rj*rj+(1.0-rj)*(1.0-rj);
        } else {
          //Observed marker is an B
          prob0= rj*rj;
          prob1= 2.0*rj*(1.0-rj);
        }
        return prob0*probright(MAA, j+1, imarker, r, position, crosstype) + prob1*probright(MH, j+1, imarker, r, position, crosstype);
      } else {
        // Unknown next marker so estimate all posibilities (imarker[j+1]==MMISSING)
        if (markertype==MAA) {
          //Observed marker is an A
          prob0= (1.0-rj)*(1.0-rj);
          prob1= 2.0*rj*(1.0-rj);
          prob2= rj*rj;
        } else if (markertype==MH) {
          //Observed marker is an H
          prob0= rj*(1.0-rj);
          prob1= rj*rj+(1.0-rj)*(1.0-rj);
          prob2= rj*(1.0-rj);
        } else {
          //Observed marker is an B
          prob0= rj*rj;
          prob1= 2.0*rj*(1.0-rj);
          prob2= (1.0-rj)*(1.0-rj);
        }
        return prob0*probright(MAA, j+1, imarker, r, position, crosstype) + prob1*probright(MH, j+1, imarker, r, position, crosstype) + prob2*probright(MBB, j+1, imarker, r, position, crosstype);
      }
      break;
    case CRIL:
      if (markertype==MH) {
        error("Strange: encountered heterozygous genotype in RIL");
        return 0.0;
      }
      if ((rightmarker==MAA)||(rightmarker==MBB)) {
        // markertype markerr diff nrecom
        //   AA        AA      0     0      1-r
        //   AA        BB     -2     2       r
        //   BB        BB      0     0      1-r
        nrecom = fabs(markertype-rightmarker);
        if (nrecom==0) {
          return (1.0-rj);
        } else {
          return rj;
        }
      } else {
        // [pjotr:] I think this code is never reached (FIXME)
        // Both markers could have recombinated which has a very low chance
        warning("Unreachable code");
        if (markertype==MAA) {
          prob0= (1.0-rj);
          prob2= rj;
        } else { // MBB
          prob0= rj;
          prob2= (1.0-rj);
        }
        return prob0*probright(MAA, j+1, imarker, r, position, crosstype) + prob2*probright(MBB, j+1, imarker, r, position, crosstype);
      }
      break;
    case CBC:
      if (markertype==MBB) {
        error("Strange: encountered BB genotype in BC");
        return 0.0;
      }
      if ((rightmarker==MAA)||(rightmarker==MH)) {
        nrecom = fabs(markertype-rightmarker);
        if (nrecom==0) {
          return 1.0-rj;
        } else {
          return rj;
        }
      } else {
        // [pjotr:] I think this code is never reached (FIXME)
        // Both markers could have recombinated which has a very low chance
        warning("Unreachable code");
        if (markertype==MAA) {
          prob0= 1.0-rj;
          prob2= rj;
        } else {
          prob0= rj;
          prob2= 1.0-rj;
        }
        return prob0*probright(MAA, j+1, imarker, r, position, crosstype) + prob2*probright(MH, j+1, imarker, r, position, crosstype);
      }
      break;
    default:
      warning("Strange: unknown crosstype in start_prob");
  }
  return 1.0;
}

