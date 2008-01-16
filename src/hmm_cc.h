/**********************************************************************
 * 
 * hmm_cc.h
 * 
 * copyright (c) 2005, Karl W Broman
 *
 * last modified Mar, 2005
 * first written Mar, 2005
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * 
 * C functions for the R/qtl package
 *
 * Contains: init_cc, emit_cc, step_cc, 
 *           calc_genoprob_cc, argmax_geno_cc
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for the Collaborative Cross
 *
 * Genotype codes:  1-8
 * Phenotype codes: 0=missing; otherwise binary 1-256 which bit i
 *                  indicating SNP compatible with parent i
 *
 **********************************************************************/

double init_cc(int true_gen);

double emit_cc(int obs_gen, int true_gen, double error_prob);
  
double step_cc(int gen1, int gen2, double rf, double junk);

void calc_genoprob_cc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *genoprob);
  
void argmax_geno_cc(int *n_ind, int *n_pos, int *geno,
		    double *rf, double *error_prob, int *argmax);

/* end of hmm_cc.h */
