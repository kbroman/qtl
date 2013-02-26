/**********************************************************************
 * 
 * hmm_main.h
 *
 * copyright (c) 2001-2010, Karl W Broman
 *
 * last modified Jul, 2010
 * first written Feb, 2001
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
 * These functions are for the main HMM engine
 *
 * Contains: calc_genoprob, calc_genoprob_special, sim_geno, est_map, argmax_geno,
 *           calc_errorlod, est_rf, calc_pairprob, calc_pairprob_condindep,
 *           R_calc_pairprob_condindep, marker_loglik
 *  
 **********************************************************************/

#define TOL       1.0e-12

/**********************************************************************
 * 
 * calc_genoprob
 *
 * This function uses the Lander-Green algorithm to calculate the 
 * genotype probabilities at each of marker and (optionally) at points
 * in-between markers, conditional on all marker data for a chromosome.
 * This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * genoprob     Genotype probabilities (the output); a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              genotype
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_genoprob(int n_ind, int n_pos, int n_gen, int *geno, 
		   double *rf, double *rf2, 
		   double error_prob, double *genoprob, 
		   double initf(int, int *), 
		   double emitf(int, int, double, int *),
		   double stepf(int, int, double, double, int *));


/**********************************************************************
 * 
 * calc_genoprob_special
 *
 * This is a special version of calc_genoprob, rerun specially for
 * each individual at each marker, assuming that that genotype may
 * be in error but others are without error
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_genoprob_special(int n_ind, int n_pos, int n_gen, int *geno, 
			   double *rf, double *rf2, 
			   double error_prob, double *genoprob, 
			   double initf(int, int *), 
			   double emitf(int, int, double, int *),
			   double stepf(int, int, double, double, int *));

/**********************************************************************
 * 
 * sim_geno
 *
 * This function simulates from the joint distribution Pr(g | O)
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              simulate genotypes)
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of simulation replicates
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps
 *
 * error_prob   Genotyping error probability
 *
 * draws        Simulated genotypes (the output), a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              draws
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void sim_geno(int n_ind, int n_pos, int n_gen, int n_draws,
	      int *geno, double *rf, double *rf2, 
	      double error_prob, int *draws,
	      double initf(int, int *), 
	      double emitf(int, int, double, int *),
	      double stepf(int, int, double, double, int *));


/**********************************************************************
 * 
 * est_map
 *
 * This function re-estimates the genetic map for a chromosome
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * n_gen        Number of different genotypes
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          Second set of recombination fractions (may not be needed)
 *
 * error_prob   Genotyping error probability
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 * nrecf1       Function returning number of recombinations associated
 *              with (g_1, g_2)
 *
 * nrecf2       Another such function, used only in the case of a sex-
 *              specific map
 *
 * loglik       Value of loglik at final estimates of rec fracs.
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 * sexsp        Indicates whether sex-specific maps should be estimated
 *
 * verbose         Indicates whether to print initial and final rec fracs
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void est_map(int n_ind, int n_mar, int n_gen, int *geno, double *rf, 
	     double *rf2, double error_prob, double initf(int, int *), 
	     double emitf(int, int, double, int *),
	     double stepf(int, int, double, double, int *), 
	     double nrecf1(int, int, double, int*), double nrecf2(int, int, double, int*), 
	     double *loglik, int maxit, double tol, int sexsp, 
	     int verbose);


/**********************************************************************
 * 
 * argmax_geno
 *
 * This function uses the Viterbi algorithm to calculate the most 
 * likely sequence of underlying genotypes, given the observed marker
 * data for a chromosome.
 * This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              find most likely genotypes
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * argmax       Matrix of most likely genotypes (the output); a single 
 *              vector stored by columns (ind moves fastest, then pos)
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void argmax_geno(int n_ind, int n_pos, int n_gen, int *geno, 
		 double *rf, double *rf2, 
		 double error_prob, int *argmax, 
		 double initf(int, int *), 
		 double emitf(int, int, double, int*),
		 double stepf(int, int, double, double, int *));

/**********************************************************************
 * 
 * calc_errorlod
 *
 * Uses the results of calc_genoprob to calculate a LOD score for 
 * each genotype, indicating whether it is likely to be in error.
 *
 * n_ind, n_mar, n_gen, geno        These are all as in the above funcs
 * error_prob, genoprob 
 *
 * errlod          The output, as a single vector stored by columns, 
 *                 of size n_ind x n_mar
 * 
 * errorlod        Function taking observed genotype, genotype probs,
 *                 and error probability, and returning the error LOD
 *
 **********************************************************************/

void calc_errorlod(int n_ind, int n_mar, int n_gen, int *geno, 
		   double error_prob, double *genoprob, double *errlod, 
		   double errorlod(int, double *, double));


/**********************************************************************
 * 
 * est_rf
 *
 * Estimate sex-averaged recombination fractions for all pairs of loci
 *
 * This is for f2 and 4way crosses; backcrosses don't need the EM 
 * algorithm, since there is no partially missing data.
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers
 *
 * geno         Matrix of genotype data (n_ind x n_mar), stored as a 
 *              single vector (by columns)
 * 
 * rf           The output: matrix of doubles (n_mar x n_mar), stored
 *              as a single vector (by columns).  The diagonal will 
 *              contain the number of meioses, the upper triangle will
 *              contain the est'd rec fracs, and the lower triangle
 *              will contain the LOD scores (testing rf=0.5)
 *
 * erec         Function returning the expected number of recombination
 *              events given observed marker genotypes
 * 
 * logprec      Function returning the log probability of a pair of 
 *              observed genotypes, given the recombination fraction
 *              (for calculating the LOD score)
 * 
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence of the EM 
 *
 * meioses_per  No. meioses per individual
 *
 **********************************************************************/

void est_rf(int n_ind, int n_mar, int *geno, double *rf, 
	    double erec(int, int, double, int *), 
	    double logprec(int, int, double, int *), 
	    int maxit, double tol, int meioses_per);

/**********************************************************************
 * 
 * calc_pairprob
 *
 * This function uses the hidden Markov model technology to calculate 
 * the joint genotype probabilities for all pairs of putative QTLs.
 * This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * genoprob     Genotype probabilities (the output); a single vector, 
 *              of length n_ind x n_pos x n_gen, stored by columns 
 *              (ind moves fastest, then mar, then genotype
 *
 * pairprob     Joint genotype probabilities for pairs of positions.
 *              A single vector of length n_ind x n_pos x (n_pos-1)/2 x
 *              n_gen^2.  We only calculate probabilities for 
 *              pairs (i,j) with i < j.
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_pairprob(int n_ind, int n_pos, int n_gen, int *geno, 
		   double *rf, double *rf2, 
		   double error_prob, double *genoprob, 
		   double *pairprob, 
		   double initf(int, int *), 
		   double emitf(int, int, double, int *),
		   double stepf(int, int, double, double, int *));

/**********************************************************************
 * 
 * calc_pairprob_condindep
 *
 * This function calculates the joint genotype probabilities assuming
 * conditional independence of QTL genotypes given the marker data
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * genoprob     QTL genotype probabilities given the marker data
 *
 * pairprob     Joint genotype probabilities for pairs of positions.
 *              A single vector of length n_ind x n_pos x (n_pos-1)/2 x
 *              n_gen^2.  We only calculate probabilities for 
 *              pairs (i,j) with i < j.
 *
 **********************************************************************/
void calc_pairprob_condindep(int n_ind, int n_pos, int n_gen, 
			     double ***Genoprob, double *****Pairprob);


/* wrapper for calc_pairprob_condindep */
void R_calc_pairprob_condindep(int *n_ind, int *n_pos, int *n_gen, 
			       double *genoprob, double *pairprob);



/**********************************************************************
 * 
 * marker_loglik
 *
 * This function calculates the log likelihood for a fixed marker
 *
 * n_ind        Number of individuals
 *
 * n_gen        Number of different genotypes
 *
 * geno         Genotype data, as a single vector
 *
 * error_prob   Genotyping error probability
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * loglik       Loglik at return
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2 */

void marker_loglik(int n_ind, int n_gen, int *geno, 
		   double error_prob, double initf(int, int *), 
		   double emitf(int, int, double, int *),
		   double *loglik);

/* end of hmm_main.h */
