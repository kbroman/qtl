/**********************************************************************
 * 
 * simulate_ril.h
 *
 * copyright (c) 2005-9, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Mar, 2005
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
 * These functions are for simulating experimental cross data;
 * I start with the simulation of RILs
 *
 * Contains: R_sim_ril, sim_ril, 
 *           allocate_individual, reallocate_individual, 
 *           copy_individual, cross, meiosis, 
 *           convertMWril, R_convertMWril
 *  
 **********************************************************************/

/* struct for individual */
struct individual {
  int max_segments;
  int n_xo[2];
  int **allele;
  double **xoloc;
};

/* wrapper for sim_ril, to be called from R */
void R_sim_ril(int *n_chr, int *n_mar, int *n_ril, double *map,
	       int *n_str, int *m, double *p, int *include_x, 
	       int *random_cross, int *selfing, int *cross, int *ril,
	       int *origgeno, double *error_prob, double *missing_prob, 
	       int *errors);

/**********************************************************************
 * 
 * sim_ril
 * 
 * n_chr   Number of chromosomes
 * n_mar   Number of markers on each chromosome (vector of length n_chr)
 * n_ril   Number of RILs to simulate
 * 
 * map     Vector of marker locations, of length sum(n_mar)
 *         First marker on each chromosome should be at 0.
 *
 * n_str   Number of parental strains (either 2, 4, or 8)
 *
 * m       Interference parameter (0 is no interference)
 * p       For Stahl model, proportion of chiasmata from the NI model
 *
 * include_x   Whether the last chromosome is the X chromosome
 *
 * random_cross  Indicates whether the order of the strains in the cross
 *               should be randomized.
 *
 * selfing If 1, use selfing; if 0, use sib mating
 *
 * cross   On output, the cross used for each line 
 *         (vector of length n_ril x n_str)
 *
 * ril     On output, the simulated data 
 *         (vector of length sum(n_mar) x n_ril)
 *
 * origgeno       Like ril, but with no missing data
 *
 * error_prob     Genotyping error probability (used only with n_str==2)
 *
 * missing_prob   Rate of missing genotypes
 *
 * errors         Error indicators (n_mar x n_ril)
 *
 **********************************************************************/
void sim_ril(int n_chr, int *n_mar, int n_ril, double *map, 
	     int n_str, int m, double p, int include_x, 
	     int random_cross, int selfing, int *cross, int *ril,
	     int *origgeno, double error_prob, double missing_prob, 
	     int *errors);

/**********************************************************************
 * allocate_individual
 **********************************************************************/
void allocate_individual(struct individual *ind, int max_seg);

/**********************************************************************
 * reallocate_individual
 **********************************************************************/
void reallocate_individual(struct individual *ind, int old_max_seg, 
			   int new_max_seg);

/**********************************************************************
 * copy_individual
 **********************************************************************/
void copy_individual(struct individual *from, struct individual *to);

/* void print_ind(struct individual ind); */

void docross(struct individual par1, struct individual par2,
	     struct individual *kid, double L, int m,
	     double p, int isX, int *maxwork, double **work);

/**********************************************************************
 * 
 * meiosis
 *
 * chrlen Chromosome length (in cM) 
 *
 * m      interference parameter (0 corresponds to no interference)
 *
 * p      for stahl model, proportion of chiasmata from NI mechanism
 *
 * maxwork
 * work
 * 
 * n_xo
 *
 **********************************************************************/
void meiosis(double L, int m, double p, int *maxwork, double **work,
	     int *n_xo);

/**********************************************************************
 * 
 * convertMWril    Convert RIL genotypes using genotypes in founders
 *                 (and the cross types).  [for a single chr]
 *
 * n_ril     Number of RILs to simulate
 * n_mar     Number of markers
 * n_str     Number of founder strains
 *
 * Parents   SNP data for the founder strains [dim n_mar x n_str]
 * 
 * Geno      On entry, the detailed genotype data; on exit, the 
 *           SNP data written bitwise. [dim n_ril x n_mar]
 * 
 * Cross     The crosses [n_ril x n_str]
 *
 * all_snps  0/1 indicator of whether all parent genotypes are snps
 *
 * error_prob  Genotyping error probability (used only if all_snps==1)
 *
 * Errors      Error indicators
 *
 **********************************************************************/
void convertMWril(int n_ril, int n_mar, int n_str, 
		  int **Parents, int **Geno, int **Crosses, 
		  int all_snps, double error_prob, 
		  int **Errors);

/* wrapper for calling convertMWril from R */
void R_convertMWril(int *n_ril, int *n_mar, int *n_str, 
		    int *parents, int *geno, int *crosses,
		    int *all_snps, double *error_prob,
		    int *errors);

/* end of simulate_ril.h */

