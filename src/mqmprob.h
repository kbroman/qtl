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

double probright(const char c, const int jloc, const cvector imarker, const vector r, const cvector position,const char crosstype);

double prob(const cmatrix loci, const vector r, const int i, const int j, const char markertype, const char crosstype, const int ADJ);

double start_prob(const char crosstype,const char markertype);

#ifdef __cplusplus
  }
#endif
/* end of mqmprob.h */
