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
  }
  return 0.0;
}

/*
 * Return probability comparing loci[j][i] versus loci[j+1][i],
 * OR if JorC is set to 1 loci[j][i] versus comparemarkertype.
 *
 * Specify an ADJ to adjust loci[j][i] to a specific location in the r[j+ADJ]
 */

double prob(const cmatrix loci, const vector r, const int i, const int j, const char comparemarkertype, const char crosstype, const int JorC, const int ADJ, const int start) {
  double calc_i=0.0;
  double Nrecom;
  char compareto;
  // Rprintf("Prob called: values:\n(i, j, ADJ)=(%d, %d, %d)\nR[j+ADJ] value: %f Loci[j][i]=%c\n", i, j, ADJ, r[j+ADJ], loci[j][i]);

  if (JorC==1) {
    //Rprintf("C %d %d\n", i, j);
    compareto = comparemarkertype;
  } else {
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
      calc_i= (r[j+ADJ]*r[j+ADJ]+(1.0-r[j+ADJ])*(1.0-r[j+ADJ]));
    } else if (Nrecom==0) {
      //Rprintf("Nrecom=0 %c <-> %c:\n", compareto, loci[j][i]);
      calc_i= (1.0-r[j+ADJ])*(1.0-r[j+ADJ]);
    } else if (Nrecom==1) {
      //Rprintf("Nrecom=1 %c <-> %c:\n", compareto, loci[j][i]);
      if (ADJ!=0) {
        calc_i= ((loci[j][i]==MH) ? 2.0*r[j+ADJ]*(1.0-r[j+ADJ]) : r[j+ADJ]*(1.0-r[j+ADJ]));
      } else {
        calc_i= ((compareto==MH) ? 2.0*r[j+ADJ]*(1.0-r[j+ADJ]) : r[j+ADJ]*(1.0-r[j+ADJ]));
      }
    } else {
      //Rprintf("Nrecom=2 %c <-> %c:\n", compareto, loci[j][i]);
      calc_i= r[j+ADJ]*r[j+ADJ];
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
      calc_i = 1.0-r[j+ADJ];
    } else {
      // Recombination between markers has a chance of r[j-1]
      calc_i = r[j+ADJ];
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
      calc_i =  (1.0-r[j+ADJ]);
    } else {
      // Recombination between markers has a chance of r[j-1]
      calc_i = r[j+ADJ];
    }
    break;
  }
  return calc_i;
}

/*
 * Return the probability of a marker being of markertype (at jloc), using the
 * information from the right neighbouring marker and the known recombination
 * frequencies. This function is used by augmentation.
 */

double probright(const char markertype, const int jloc, const cvector imarker, const vector r, const cvector position, const char crosstype) {
  double nrecom, prob0, prob1, prob2;
  if ((position[jloc]==MRIGHT)||(position[jloc]==MUNLINKED)) {
    //We're at the end of a chromosome or an unlinked marker
    return 1.0;
  }
  switch (crosstype) {
    case CF2:
      if ((imarker[jloc+1]==MAA)||(imarker[jloc+1]==MH)||(imarker[jloc+1]==MBB)) {
        //NEXT marker is known
        if ((markertype==MH)&&(imarker[jloc+1]==MH)) {
          //special case in which we observe a H after an H then we can't know if we recombinated or not
          return r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
        } else {
          //The number of recombinations between observed marker and the next marker
          nrecom = fabs(markertype-imarker[jloc+1]);
          if (nrecom==0) {
            //No recombination
            return (1.0-r[jloc])*(1.0-r[jloc]);
          } else if (nrecom==1) {
            if (imarker[jloc+1]==MH) {
              //the chances of having a H after 1 recombination are 2 times the chance of being either A or B
              return 2.0*r[jloc]*(1.0-r[jloc]);
            } else {
              //Chance of 1 recombination
              return r[jloc]*(1.0-r[jloc]);
            }
          } else {
            //Both markers could have recombinated which has a very low chance
            return r[jloc]*r[jloc];
          }
        }
      } else if (imarker[jloc+1]==MNOTAA) {
        //SEMI unknown next marker known is it is not an A
        if (markertype==MAA) {
          //Observed marker is an A
          prob1= 2.0*r[jloc]*(1.0-r[jloc]);
          prob2= r[jloc]*r[jloc];
        } else if (markertype==MH) {
          //Observed marker is an H
          prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
          prob2= r[jloc]*(1.0-r[jloc]);
        } else {
          //Observed marker is an B
          prob1= 2.0*r[jloc]*(1.0-r[jloc]);
          prob2= (1.0-r[jloc])*(1-r[jloc]);
        }
        return prob1*probright(MH, jloc+1, imarker, r, position, crosstype) + prob2*probright(MBB, jloc+1, imarker, r, position, crosstype);
      } else if (imarker[jloc+1]==MNOTBB) {
        //SEMI unknown next marker known is it is not a B
        if (markertype==MAA) {
          //Observed marker is an A
          prob0= (1.0-r[jloc])*(1.0-r[jloc]);
          prob1= 2.0*r[jloc]*(1.0-r[jloc]);
        } else if (markertype==MH) {
          //Observed marker is an H
          prob0= r[jloc]*(1.0-r[jloc]);
          prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
        } else {
          //Observed marker is an B
          prob0= r[jloc]*r[jloc];
          prob1= 2.0*r[jloc]*(1.0-r[jloc]);
        }
        return prob0*probright(MAA, jloc+1, imarker, r, position, crosstype) + prob1*probright(MH, jloc+1, imarker, r, position, crosstype);
      } else {
        // Unknown next marker so estimate all posibilities (imarker[j+1]==MMISSING)
        if (markertype==MAA) {
          //Observed marker is an A
          prob0= (1.0-r[jloc])*(1.0-r[jloc]);
          prob1= 2.0*r[jloc]*(1.0-r[jloc]);
          prob2= r[jloc]*r[jloc];
        } else if (markertype==MH) {
          //Observed marker is an H
          prob0= r[jloc]*(1.0-r[jloc]);
          prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
          prob2= r[jloc]*(1.0-r[jloc]);
        } else {
          //Observed marker is an B
          prob0= r[jloc]*r[jloc];
          prob1= 2.0*r[jloc]*(1.0-r[jloc]);
          prob2= (1.0-r[jloc])*(1.0-r[jloc]);
        }
        return prob0*probright(MAA, jloc+1, imarker, r, position, crosstype) + prob1*probright(MH, jloc+1, imarker, r, position, crosstype) + prob2*probright(MBB, jloc+1, imarker, r, position, crosstype);
      }
      break;
    case CRIL:
      if (markertype==MH) {
        error("Strange: encountered heterozygous genotype in RIL");
        return 0.0;
      }
      if ((imarker[jloc+1]==MAA)||(imarker[jloc+1]==MBB)) {
        // markertype markerr diff nrecom
        //   AA        AA      0     0      1-r
        //   AA        BB     -2     2       r
        //   BB        BB      0     0      1-r
        nrecom = fabs(markertype-imarker[jloc+1]);
        if (nrecom==0) {
          return (1.0-r[jloc]);
        } else {
          return r[jloc];
        }
      } else {
        // [pjotr:] I think this code is never reached (FIXME)
        // Both markers could have recombinated which has a very low chance
        warning("Unreachable code");
        if (markertype==MAA) {
          prob0= (1.0-r[jloc]);
          prob2= r[jloc];
        } else { // MBB
          prob0= r[jloc];
          prob2= (1.0-r[jloc]);
        }
        return prob0*probright(MAA, jloc+1, imarker, r, position, crosstype) + prob2*probright(MBB, jloc+1, imarker, r, position, crosstype);
      }
      break;
    case CBC:
      if (markertype==MBB) {
        error("Strange: encountered BB genotype in BC");
        return 0.0;
      }
      if ((imarker[jloc+1]==MAA)||(imarker[jloc+1]==MH)) {
        nrecom = fabs(markertype-imarker[jloc+1]);
        if (nrecom==0) {
          return 1.0-r[jloc];
        } else {
          return r[jloc];
        }
      } else {
        // [pjotr:] I think this code is never reached (FIXME)
        // Both markers could have recombinated which has a very low chance
        warning("Unreachable code");
        if (markertype==MAA) {
          prob0= 1.0-r[jloc];
          prob2= r[jloc];
        } else {
          prob0= r[jloc];
          prob2= 1.0-r[jloc];
        }
        return prob0*probright(MAA, jloc+1, imarker, r, position, crosstype) + prob2*probright(MH, jloc+1, imarker, r, position, crosstype);
      }
      break;
  }
  return 1.0;
}

