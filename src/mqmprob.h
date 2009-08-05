/**********************************************************************
 *
 * mqmprob.h
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Feb, 2009
 * first written Feb, 2009
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
 * C external functions used by the MQM algorithm
 * Contains (stabile): prob, start_prob, probright
 *
 **********************************************************************/


#ifdef __cplusplus
  extern "C" {
#endif


cvector locate_markers(const int nmark,const ivector chr);
vector recombination_frequencies(const int nmark, const cvector position, const vector mapdistance);

double probright(const char c, const int j, const cvector imarker, const vector rs, const cvector position,const MQMCrossType crosstype);

double prob(const cmatrix loci, const vector rs, const int i, const int j, const char markertype, const MQMCrossType crosstype, const int ADJ);

double start_prob(const MQMCrossType crosstype,const char markertype);

#ifdef __cplusplus
  }
#endif
/* end of mqmprob.h */
