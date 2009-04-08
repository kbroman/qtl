/**********************************************************************
 * 
 * MQMprob.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/

void create_lookup_table(double ***MendelM,int Nmark,vector r,char crosstype);

double probnew(double ***MendelM,cmatrix loci, vector r, int i, int j,char c,char crosstype,int JorC,int ADJ,int start);

double probright(char c, int jloc, cvector imarker, vector r, cvector position,char crosstype);

double prob(cmatrix loci, vector r,int i,int j,char c,char crosstype,int JorC,int ADJ,int start);

double start_prob(char crosstype,char c);

/* end of MQMprob.h */
