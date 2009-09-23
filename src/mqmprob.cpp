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

/* This function walks the marker list and determins for every position whether the marker is in the Middle, Left,Right of the chromosome
 When there is only 1 marker on a chromosome it is defined Unlinked*/
 
// FIXME Naming of function (set_relative_marker_position / findLeftRightMarkerpos)
cvector locate_markers(const int nmark,const ivector chr) 
{
  cvector position = newcvector(nmark);
  info("Calculating relative genomepositions of the markers");
  for (int j=0; j<nmark; j++) {
    if (j==0) {
      if (chr[j]==chr[j+1]) 
        position[j]=MLEFT;
      else 
        position[j]=MUNLINKED;
    } else if (j==nmark-1) {
      if (chr[j]==chr[j-1]) 
        position[j]=MRIGHT;
      else 
        position[j]=MUNLINKED;
    } else if (chr[j]==chr[j-1]) {
      if (chr[j]==chr[j+1]) 
        position[j]=MMIDDLE;
      else 
        position[j]=MRIGHT;
    } else {
      if (chr[j]==chr[j+1]) 
        position[j]=MLEFT;
      else 
        position[j]=MUNLINKED;
    }
  }
  return position;
}

/*Using haldane we calculate Recombination frequencies. Using relative marker positions and distances Return array of rec. frequencies (unknown = 999.0)*/
//FIXME 999.0 to RFUNKNOWN
//FIXME unsigned int als array index (uint)
//NOTE checking for r[j] <0 (marker ordering) can ahppen at contract
vector recombination_frequencies(const int nmark, const cvector position, const vector mapdistance) 
{
  // contract: if DEBUG is_valid(positionarray)
  
  info("Estimating recombinant frequencies");
  vector r = newvector(nmark);
  for (int j=0; j<nmark; j++) {
    r[j]= 999.0;
    if ((position[j]==MLEFT)||(position[j]==MMIDDLE)) {
      r[j]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-mapdistance[j])));
      if (r[j]<0) {
        Rprintf("ERROR: Position=%d r[j]=%f\n", position[j], r[j]);
        fatal("Recombination frequency is negative, (Marker ordering problem ?)");
        return NULL;
      }
    }
    //RRprintf("recomfreq:%d, %f\n", j, r[j]);
  }
  return r;
}



/* 
 * Ascertain a marker is valid (AA, BB or AB) for the cross type
 */

void validate_markertype(const MQMCrossType crosstype, const char markertype)
{
  if (markertype==MNOTAA || markertype==MNOTBB || markertype==MUNKNOWN)
    fatal("validate_markertype: Undecided markertype");
  if (crosstype==CRIL && markertype==MH) 
    fatal("validate_markertype: Found markertype H (AB) in RIL");
  if (crosstype==CBC && markertype==MBB) 
    fatal("validate_markertype: Found markertype BB in back cross (BC)");
}

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
 
 in the case of marketypes that are invalid for a cross (here: x) a probability of 0.0 is returned
 This is because to make the code in augmentation & Regression more generic.
 */

 //FIXME: Special function for special case (start_prob for regression)
 //FIXME: MQMmarker ipv cont char en dan als index (dus in datatypes -> '0'-> 0)
double start_prob(const MQMCrossType crosstype, const char markertype) {
  //validate_markertype(crosstype,markertype);
  switch (crosstype) {
    case CF2:
      switch (markertype) {
      case MAA:
        return 0.25;
      case MH:
        return 0.5;
      case MBB:
        return 0.25;
      default:
        Rprintf("Strange: Probability requested for invalid markertype: %c",markertype);
        return 0.0;
      }
    case CRIL:
      switch (markertype) {
      case MAA:
        return 0.5;
      case MH:
        return 0.0;
      case MBB:
        return 0.5;
      default:
        Rprintf("Strange: Probability requested for invalid markertype: %c",markertype);
        return 0.0;
      }
    case CBC:
      switch (markertype) {
      case MAA:
        return 0.5;
      case MH:
        return 0.5;
      case MBB:
        return 0.0;
      default:
        Rprintf("Strange: Probability requested for invalid markertype: %c",markertype);
        return 0.0;
      }
      //return (markertype==MH ? 0.5 : 0.5);
    default:
      fatal("Strange: unknown crosstype in start_prob");
  }
  fatal("Should not get here");
  return NAN;
}


/*
 * Return probability comparing loci[j][i] versus loci[j+1][i],
 * OR if markertype is set, loci[j][i] versus comparemarkertype. Probability
 * is simply based on the markertypes and known recombination rate.
 *
 * Specify an ADJ to adjust loci[j][i] to a specific location in the r[j+ADJ]
 
 FIXME zie start_prob for the same issue with markertypes not bellong to the current cross
 */
// IPV rs geven we de 'single Rec Freq' dan vervalt ADJ en kan de berekening korter en mooier
//refactor name left_prob


double prob_new(const double r, const char markerL,const char markerR,const MQMCrossType crosstype){
  const double r2 = r*r;  //Double recombination
  const double rr = 1.0-r; // No recombination
  const double rr2 = rr*rr; //Double Norecombination
  //Calculate the number of recombinations
  const int recombinations = fabs((double)markerL-(double)markerR);

  switch (crosstype) {
    case CF2:
      if ((markerL==MH)&&(markerR==MH)) {
        return r2 + rr2; //Special case H after H  So either double recombination or double Norecombination
      }else if (recombinations==0) {
        return rr2;
      }else if (recombinations==1) {
        return ((markerR==MH) ? 2.0*r*rr : r*rr); //If the marker was a H then we have 2 * recombination*No recombination chance otherwise just single chance
      }else{
        return r2; // two recombinations
      }
      return rr; // Not Recombinated
      break;
    case CRIL:
      if (markerR==MH) {
        return 0.0; // No chance finding a 1 or H in an RIL
      }
      if (recombinations){
        return r; //Recombinated
      }
      return rr; // Not Recombinated
      break;
    case CBC:
      if (markerR==MBB) {
        return 0.0; // No chance finding a 2/BB in a BC
      }
      if (recombinations){
        return r; //Recombinated
      }
      return rr; // Not Recombinated
      break;
    default:
      fatal("Strange: unknown crosstype in prob");
  }
  fatal("Should not get here");
  return NAN;
}


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
    fatal("We never get here, all calls happen to pass in the markertype");
    compareto = loci[j+1][i]; // FIXME
  }

  //validate_markertype(crosstype,compareto);
  //validate_markertype(crosstype,markertype);

  // number of recombinations recombinations
  const int recombinations = fabs((double)markertype-(double)compareto);
  double prob = rr;  // default to no recombinations (1-r)
  switch (crosstype) {
    case CF2:
      if ((markertype==MH)&&(compareto==MH)) {
        prob = r2 + rr2;
      } else if (recombinations==0) {
        prob = rr2;
      } else if (recombinations==1) {
        if (ADJ!=0) {  // FIXME: now this is not clear to me
          prob = ((markertype==MH) ? 2.0*r*rr : r*rr);
        } else {
          prob = ((compareto==MH) ? 2.0*r*rr : r*rr);
        }
      } else {
        prob = r2;  // two recombinations
      }
      break;
    case CRIL:
      if (compareto==MH) {
        return 0.0; // No chance finding a 1 or H in an RIL
      }
      if (recombinations) prob = r;
      break;
    case CBC:
      if (compareto==MBB) {
        return 0.0; // No chance finding a 2/BB in a BC
      }
      if (recombinations) prob = r;
      break;
    default:
      fatal("Strange: unknown crosstype in prob");
  }
  return prob;
}

/*
 * Return the probability of a marker being of markertype (at j), using the
 * information from the right neighbouring marker and the known recombination
 * frequencies. This function is used by augmentation.
 */
//fixme IF to switch case
//Split into three function !
//markkertype -> Marker because it is a marker not a Type of marker
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
  //validate_markertype(crosstype,markertype);
  //validate_markertype(crosstype,rightmarker);

  // markertype markerr diff recombinations
  //   AA        AA      0     0      1-r
  //   AA        BB     -2     2       r
  //   BB        BB      0     0      1-r
  const int recombinations = fabs(markertype-rightmarker);
  switch (crosstype) {
    case CF2:
      //fixme smallfunction: bool is_knownMarker() ?
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
        //FIXME If (semi) unknown marker we need to step untill we reach a known one
        //TESTME: What happens if a chromosome ends with unknown markers AAAAAHHHHHUUUUUUUU
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
    //FIXME: RIL case is same as BC case if we use a function is_known_marker_for_cross()
      if (markertype==MH) {
        //info("Strange: encountered heterozygous genotype in RIL");
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
        //info("Unreachable code");
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
        //info("Strange: encountered BB genotype in BC");
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
        //[Danny] Code is reached in MQMaugment, We could have an unknown or semi known next marker in the cross e.g.  A A A H U U A A
        // Both markers could have recombinated which has a very low chance
        //info("Unreachable code");
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
      info("Strange: unknown crosstype in start_prob");
  }
  return 1.0;
}

