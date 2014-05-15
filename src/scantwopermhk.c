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
#include "scanone_hk.h"
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

  GetRNGstate();
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

  PutRNGstate();
}

/**
 * scantwopermhk_1chr_nocovar: with no covariates, can do calculations in batch
 */

void scantwopermhk_1chr_nocovar(int n_ind, int n_pos, int n_gen,
                                double ***Genoprob, double *****Pairprob,
                                double *pheno, int n_perm, double *weights,
                                double **Result, int n_col2drop, int *col2drop)
{
  double *phematrix, *scanone_result, **scanone_Result;
  double *scantwo_result, ***scantwo_Result;
  int *ind_noqtl;

  /* setup */
  create_shuffled_phematrix(n_ind, n_perm, pheno, &phematrix);
  create_zero_vector(&ind_noqtl, n_ind);
  allocate_double(n_perm*n_pos, &scanone_result);
  reorg_errlod(n_pos, n_perm, scanone_result, &scanone_Result);
  allocate_double(n_perm*n_pos*n_pos, &scantwo_result);
  reorg_genoprob(n_pos, n_pos, n_perm, scantwo_result, &scantwo_Result);

  /* scanone */
  scanone_hk(n_ind, n_pos, n_gen, Genoprob,
             0, 0, 0, 0, /* null covariates */
             phematrix, n_perm, weights, scanone_Result,
             ind_noqtl);

  /* scantwo */
  scantwo_1chr_hk(n_ind, n_pos, n_gen, Genoprob, Pairprob,
                  0, 0, 0, 0, /* null covariates */
                  phematrix, n_perm, weights,
                  scantwo_Result, n_col2drop, col2drop);
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
  int *ind_index, i;
  double *dwork;
  double *scanone_result, **scanone_Result;
  double *scantwo_result, ***scantwo_Result;
  int *ind_noqtl;

  /* setup */
  allocate_int(n_ind, &ind_index);
  for(i=0; i<n_ind; i++) ind_index[i] = i;
  create_zero_vector(&ind_noqtl, n_ind);
  allocate_double(n_pos, &scanone_result);
  reorg_errlod(n_pos, 1, scanone_result, &scanone_Result);
  allocate_double(n_pos*n_pos, &scantwo_result);
  reorg_genoprob(n_pos, n_pos, 1, scantwo_result, &scantwo_Result);


  /* space for shuffles */
  allocate_double(n_ind, &dwork);

  for(i=0; i<n_perm; i++) {
    shuffle_covar_and_phe(n_ind, ind_index, pheno, Addcov, n_addcov, dwork);

    /* scanone */
    scanone_hk(n_ind, n_pos, n_gen, Genoprob,
               Addcov, n_addcov, 0, 0, /* null inter've covariates */
               pheno, 1, weights, scanone_Result,
               ind_noqtl);

    /* scantwo */
    scantwo_1chr_hk(n_ind, n_pos, n_gen, Genoprob, Pairprob,
                    Addcov, n_addcov, 0, 0, /* null inter've covariates */
                    pheno, 1, weights,
                    scantwo_Result, n_col2drop, col2drop);

  }
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

  GetRNGstate();

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
  PutRNGstate();
}


/**
 * scantwo permhk_2chr_nocovar: with no covariates, can do calculations in batch
 **/

void scantwopermhk_2chr_nocovar(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                                int n_gen2, double ***Genoprob1, double ***Genoprob2,
                                double *pheno, int n_perm, double *weights,
                                double **Result)
{
  double *phematrix;
  double *scanone_result1, **scanone_Result1;
  double *scanone_result2, **scanone_Result2;
  double *scantwo_result_full, ***scantwo_Result_Full;
  double *scantwo_result_add, ***scantwo_Result_Add;
  int *ind_noqtl;

  /* setup */
  create_shuffled_phematrix(n_ind, n_perm, pheno, &phematrix);
  create_zero_vector(&ind_noqtl, n_ind);
  allocate_double(n_perm*n_pos1, &scanone_result1);
  reorg_errlod(n_pos1, n_perm, scanone_result1, &scanone_Result1);
  allocate_double(n_perm*n_pos2, &scanone_result2);
  reorg_errlod(n_pos2, n_perm, scanone_result2, &scanone_Result2);
  allocate_double(n_perm*n_pos1*n_pos2, &scantwo_result_full);
  reorg_genoprob(n_pos2, n_pos1, n_perm, scantwo_result_full, &scantwo_Result_Full);
  allocate_double(n_perm*n_pos1*n_pos2, &scantwo_result_add);
  reorg_genoprob(n_pos1, n_pos2, n_perm, scantwo_result_add, &scantwo_Result_Add);

  /* scanone */
  scanone_hk(n_ind, n_pos1, n_gen1, Genoprob1,
             0, 0, 0, 0, /* null covariates */
             phematrix, n_perm, weights, scanone_Result1,
             ind_noqtl);

  scanone_hk(n_ind, n_pos2, n_gen2, Genoprob2,
             0, 0, 0, 0, /* null covariates */
             phematrix, n_perm, weights, scanone_Result2,
             ind_noqtl);

  /* scantwo */
  scantwo_2chr_hk(n_ind, n_pos1, n_pos2, n_gen1, n_gen2,
                  Genoprob1, Genoprob2,
                  0, 0, 0, 0, /* null covariates */
                  phematrix, n_perm, weights,
                  scantwo_Result_Full, scantwo_Result_Add);

}


/**
 * scantwopermhk_2chr
 */

void scantwopermhk_2chr(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                        int n_gen2, double ***Genoprob1, double ***Genoprob2,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, double *weights, double **Result)
{
  int *ind_index, i;
  double *dwork;
  double *scanone_result1, **scanone_Result1;
  double *scanone_result2, **scanone_Result2;
  double *scantwo_result_full, ***scantwo_Result_Full;
  double *scantwo_result_add, ***scantwo_Result_Add;
  int *ind_noqtl;

  /* setup */
  allocate_int(n_ind, &ind_index);
  for(i=0; i<n_ind; i++) ind_index[i] = i;
  create_zero_vector(&ind_noqtl, n_ind);
  allocate_double(n_pos1, &scanone_result1);
  reorg_errlod(n_pos1, 1, scanone_result1, &scanone_Result1);
  allocate_double(n_pos2, &scanone_result2);
  reorg_errlod(n_pos2, 1, scanone_result2, &scanone_Result2);
  allocate_double(n_pos1*n_pos2, &scantwo_result_full);
  reorg_genoprob(n_pos2, n_pos1, 1, scantwo_result_full, &scantwo_Result_Full);
  allocate_double(n_pos1*n_pos2, &scantwo_result_add);
  reorg_genoprob(n_pos1, n_pos2, 1, scantwo_result_add, &scantwo_Result_Add);

  /* space for shuffles */
  allocate_double(n_ind, &dwork);

  for(i=0; i<n_perm; i++) {
    shuffle_covar_and_phe(n_ind, ind_index, pheno, Addcov, n_addcov, dwork);

    /* scanone */
    scanone_hk(n_ind, n_pos1, n_gen1, Genoprob1,
               Addcov, n_addcov, 0, 0, /* null inter've covariates */
               pheno, 1, weights, scanone_Result1,
               ind_noqtl);

    scanone_hk(n_ind, n_pos2, n_gen2, Genoprob2,
               Addcov, n_addcov, 0, 0, /* null inter've covariates */
               pheno, 1, weights, scanone_Result2,
               ind_noqtl);

    /* scantwo */
    scantwo_2chr_hk(n_ind, n_pos1, n_pos2, n_gen1, n_gen2,
                    Genoprob1, Genoprob2,
                    Addcov, n_addcov, 0, 0, /* null inter've covariates */
                    pheno, 1, weights,
                    scantwo_Result_Full, scantwo_Result_Add);

  }

}


/* create a matrix of shuffled phenotypes */
void create_shuffled_phematrix(int n_ind, int n_perm, double *pheno, double **phematrix)
{
  int i, j;

  allocate_double(n_ind*n_perm, phematrix);

  for(i=0; i<n_perm; i++) {
    for(j=0; j<n_ind; j++)
      (*phematrix)[i*n_ind + j] = pheno[j];
    double_permute((*phematrix)+i*n_ind, n_ind);
  }
}

/* shuffle the rows in the phenotype data and covariate matrix */
void shuffle_covar_and_phe(int n_ind, int *ind_index, double *pheno,
                           double **Addcov, int n_addcov,
                           double *dwork)
{
  int i, j;

  /* shuffle index */
  int_permute(ind_index, n_ind);

  /* reorder phenotypes */
  memcpy(dwork, pheno, n_ind*sizeof(double));
  for(i=0; i<n_ind; i++)
    pheno[i] = dwork[ind_index[i]];

  /* reorder rows in Addcov */
  for(j=0; j<n_addcov; j++) {
    memcpy(dwork, Addcov[j], n_ind*sizeof(double));
    for(i=0; i<n_ind; i++)
      Addcov[j][i] = dwork[ind_index[i]];
  }
}

/* create a vector of 0's */
void create_zero_vector(int **vector, int n)
{
  int i;

  allocate_int(n, vector);

  for(i=0; i<n; i++) (*vector)[i] = 0;
}
