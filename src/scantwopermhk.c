/**
 * scantwo_hk.c
 **/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "lapackutil.h"
#include "scantwo_hk.h"
#include "scantwopermhk.h"
#define TOL 1e-12

/**
 * R_scantwopermhk_1chr
 **/

void R_scantwopermhk_1chr(int *n_ind, int *n_pos, int *n_gen,
                          double *genoprob, double *pairprob, 
                          double *addcov, int *n_addcov, 
                          double *pheno, int* n_perm, double *weights, 
                          double *result, int *n_col2drop, int *col2drop)
{
  double ***Genoprob, **Result, **Addcov=0, *****Pairprob;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
  reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
  reorg_errlod(*n_perm, 5, result, &Result);

  /* reorganize addcov (if not empty) */
  if(*n_addcov > 0) {
    reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);

    scantwopermhk_1chr(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob, 
                       Addcov, *n_addcov,
                       pheno, *n_perm, weights, Result, *n_col2drop,
                       col2drop);
  }
  else {
    scantwopermhk_1chr_nocovar(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob, 
                               pheno, *n_perm, weights, Result, *n_col2drop,
                               col2drop);
  }
}

/**
 * scantwopermhk_1chr_nocovar: with no covariates, can do calculations in batch
 */

void scantwopermhk_1chr_nocovar(int n_ind, int n_pos, int n_gen,
                                double ***Genoprob, double *****Pairprob,
                                double *pheno, int n_perm, double *weights, 
                                double **Result, int n_col2drop, int *col2drop)
{

}



/**
 * scantwopermhk_1chr
 */

void scantwopermhk_1chr(int n_ind, int n_pos, int n_gen,
                        double ***Genoprob, double *****Pairprob,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, double *weights, double **Result,
                        int n_col2drop, int *col2drop)
{

}




/**
 * R_scantwopermhk_2chr
 **/

void R_scantwopermhk_2chr(int *n_ind, int *n_pos1, int *n_pos2, 
                          int *n_gen1, int *n_gen2,
                          double *genoprob1, double *genoprob2,
                          double *addcov, int *n_addcov, 
                          double *pheno, int *n_perm, double *weights,
                          double *result)
{
  double ***Genoprob1, ***Genoprob2, **Result, **Addcov=0;

  reorg_genoprob(*n_ind, *n_pos1, *n_gen1, genoprob1, &Genoprob1);
  reorg_genoprob(*n_ind, *n_pos2, *n_gen2, genoprob2, &Genoprob2);
  reorg_errlod(*n_perm, 5, result, &Result);

  /* reorganize addcov (if not empty) */
  if(*n_addcov > 0) {
    reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
    
    scantwopermhk_2chr(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
                       Genoprob1, Genoprob2, Addcov, *n_addcov,
                       pheno, *n_perm, weights, Result);
  }
  else {
    scantwopermhk_2chr_nocovar(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2, 
                               Genoprob1, Genoprob2,
                               pheno, *n_perm, weights, Result);
  }
}


/**
 * scantwo permhk_2chr_nocovar: with no covariates, can do calculations in batch
 **/

void scantwopermhk_2chr_nocovar(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                                int n_gen2, double ***Genoprob1, double ***Genoprob2,
                                double *pheno, int n_perm, double *weights, 
                                double **Result)
{

}


/**
 * scantwopermhk_2chr
 */

void scantwopermhk_2chr(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                        int n_gen2, double ***Genoprob1, double ***Genoprob2,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, double *weights, double **Result)
{

}
