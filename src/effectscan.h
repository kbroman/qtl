/**********************************************************************
 * 
 * effectscan.h
 *
 * copyright (c) 2007, Karl W Broman 
 *
 * last modified Sep, 2007
 * first written Sep, 2007
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/qtl package
 *
 * These functions are for calculating the estimated effects, by multiple
 * imputation, in a single-QTL scan along a chromosome.
 *
 * Contains: R_effectscan, effectscan
 *
 **********************************************************************/

/* R_effectscan: wrapper for effectscan */
void R_effectscan(int *nind, int *ngen, int *ndraws, int *npos,
		  int *draws, double *pheno, double *effectmapping,
		  double *beta, double *se, int *getse);
  

/**********************************************************************
 * effectscan
 *
 * nind   Number of individuals
 * ngen   Number of genotypes
 * ndraws Number of imputations
 * npos   Number of positions
 * Draws  The imputed genotypes (dim nind x npos x ndraws)
 * pheno  Phenotypes (length nind)
 * effectmapping  Matrix of size ngen x ngen, giving the design matrix 
 *                for each possible genotype
 * Beta   On exit, the estimated coefficients (dim npos x ngen)
 * SE     On exit, the estimated standard errors (dim npos x ngen)
 * getse  If 1, calculate SEs; if 0, don't
 *
 **********************************************************************/
void effectscan(int nind, int ngen, int ndraws, int npos,
		int ***Draws, double *pheno, double *mapping,
		double **Beta, double **SE, int getse);
 
/* end of effectscan.h */
