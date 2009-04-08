/**********************************************************************
 * 
 * MQMprob.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License, as
 *     published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version. 
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the
 *     GNU General Public License for more details.
 * 
 *     A copy of the GNU General Public License is available at
 *     http://www.r-project.org/Licenses/
 *
 * C external functions used by the MQM algorithm
 * Contains (stabile): prob, start_prob, probright
 * Contains (unstabile): create_lookup_table, probnew
 *
 **********************************************************************/

void create_lookup_table(double ***MendelM,int Nmark,vector r,char crosstype);

double probnew(double ***MendelM,cmatrix loci, vector r, int i, int j,char c,char crosstype,int JorC,int ADJ,int start);

double probright(char c, int jloc, cvector imarker, vector r, cvector position,char crosstype);

double prob(cmatrix loci, vector r,int i,int j,char c,char crosstype,int JorC,int ADJ,int start);

double start_prob(char crosstype,char c);

/* end of MQMprob.h */
