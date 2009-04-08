/**********************************************************************
 * 
 * MQMmixture.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/
/* ML estimation of recombination frequencies via EM;
    calculation of multilocus genotype probabilities;
    ignorance of unlikely genotypes*/
double rmixture(cmatrix marker, vector weight, vector r,
              cvector position, ivector ind,
              int Nind, int Naug, int Nmark,vector *mapdistance, char reestimate,char crosstype,Mmatrix MendelM,int verbose);


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
              vector y, ivector ind, int Nind, int Naug,
              int Nloci,
              double *variance, int em, vector *weight,char REMLorML,char fitQTL,char dominance,char crosstype,Mmatrix MendelM,int verbose);

/* end of MQMmixture.c */
