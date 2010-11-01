/**********************************************************************
 *
 * mqmprob.h
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Danny Arends and Pjotr Prins
 * last modified September 2009
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
 **********************************************************************/


#ifdef __cplusplus
  extern "C" {
#endif



cvector relative_marker_position(const unsigned int nmark,const ivector chr);
vector recombination_frequencies(const unsigned int nmark, const cvector position, const vector mapdistance);
double recombination_frequentie(const double cmdistance);
void validate_markertype(const MQMCrossType crosstype, const char markertype);
//double probright(const char c, const int j, const cvector imarker, const vector rs, const cvector position,const MQMCrossType crosstype);

double left_prob(const double r, const MQMMarker markerL,const MQMMarker markerR,const MQMCrossType crosstype);
double right_prob_F2(const char markerL, const int j, const MQMMarkerVector imarker, const vector rs, const cvector position);
double right_prob_BC(const char markerL, const int j, const MQMMarkerVector imarker, const vector rs, const cvector position);
double right_prob_RIL(const char markerL, const int j, const MQMMarkerVector imarker, const vector rs, const cvector position);
//double prob(const cmatrix loci, const vector rs, const int i, const int j, const char markertype, const MQMCrossType crosstype, const int ADJ);

double start_prob(const MQMCrossType crosstype,MQMMarker markertype);
bool is_knownMarker(const char marker,const MQMCrossType crosstype);

#ifdef __cplusplus
  }
#endif
/* end of mqmprob.h */
