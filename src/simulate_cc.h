/**********************************************************************
 * 
 * simulate_cc.h
 *
 * copyright (c) 2005-6, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Mar, 2005
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
 * C functions for the R/qtl package
 *
 * These functions are for simulating experimental cross data;
 * I start with the simulation of RILs
 *
 * Contains: R_sim_ril, sim_ril, 
 *           allocate_individual, reallocate_individual, 
 *           copy_individual, cross, meiosis, sim_cc, R_sim_cc
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
	       int *random_cross, int *selfing, int *cross, int *ril);

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
 **********************************************************************/
void sim_ril(int n_chr, int *n_mar, int n_ril, double *map, 
	     int n_str, int m, double p, int include_x, 
	     int random_cross, int selfing, int *cross, int *ril);

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
 * sim_cc    Use the result of sim_all_ril with n_str=8 plus data on
 *           the SNP genotypes of the 8 parental strains to create 
 *           real SNP data for the Collaborative Cross
 *
 * n_ril     Number of RILs to simulate
 * tot_mar   Total number of markers
 *
 * Parents   SNP data for the 8 parental lines [dim tot_mar x 8]
 * 
 * Geno      On entry, the detailed genotype data; on exit, the 
 *           SNP data written bitwise.
 * 
 * error_prob  Probability of genotyping error
 * missing_prob  Probability a genotype will be missing
 *
 **********************************************************************/
void sim_cc(int n_ril, int tot_mar, int **Parents, int **Geno,
	    double error_prob, double missing_prob);

/* wrapper for calling sim_cc from R */
void R_sim_cc(int *n_ril, int *tot_mar, int *parents, int *geno,
	      double *error_prob, double *missing_prob);

/* end of simulate_cc.h */

