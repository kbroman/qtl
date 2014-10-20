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

/* R wrapper for function to perform scantwo permutations by
   Haley-Knott regression within a chromosome */
void R_scantwopermhk_1chr(int *n_ind, int *n_pos, int *n_gen,
                          double *genoprob, double *pairprob,
                          double *addcov, int *n_addcov,
                          double *pheno, int* n_perm, int *permindex,
                          double *weights, double *result,
                          int *n_col2drop, int *col2drop)
{
    double ***Genoprob, **Result, **Addcov=0, *****Pairprob;
    int **Permindex;

    reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
    reorg_pairprob(*n_ind, *n_pos, *n_gen, pairprob, &Pairprob);
    reorg_errlod(*n_perm, 6, result, &Result);
    reorg_geno(*n_ind, *n_perm, permindex, &Permindex);

    /* reorganize addcov (if not empty) */
    if(*n_addcov > 0) {
        reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);

        scantwopermhk_1chr(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob,
                           Addcov, *n_addcov,
                           pheno, *n_perm, Permindex, weights,
                           Result, *n_col2drop, col2drop);
    }
    else {
        scantwopermhk_1chr_nocovar(*n_ind, *n_pos, *n_gen, Genoprob, Pairprob,
                                   pheno, *n_perm, Permindex, weights,
                                   Result, *n_col2drop, col2drop);
    }
}


/* scantwo perms within chr; with no covariates, can do calculations in batch */
void scantwopermhk_1chr_nocovar(int n_ind, int n_pos, int n_gen,
                                double ***Genoprob, double *****Pairprob,
                                double *pheno, int n_perm, int **Permindex,
                                double *weights, double **Result,
                                int n_col2drop, int *col2drop)
{
    int i;
    double *phematrix, **Phematrix, *scanone_result, **scanone_Result;
    double *scantwo_result, ***scantwo_Result;
    int *ind_noqtl;

    /* setup */
    allocate_double(n_perm*n_ind, &phematrix);
    reorg_errlod(n_ind, n_perm, phematrix, &Phematrix);
    create_zero_vector(&ind_noqtl, n_ind);
    allocate_double(n_perm*n_pos, &scanone_result);
    reorg_errlod(n_pos, n_perm, scanone_result, &scanone_Result);
    allocate_double(n_perm*n_pos*n_pos, &scantwo_result);
    reorg_genoprob(n_pos, n_pos, n_perm, scantwo_result, &scantwo_Result);

    /* order phenotypes as shuffled*/
    fill_phematrix(n_ind, n_perm, pheno, Permindex, Phematrix);

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

    min3d_uppertri(n_pos, n_perm, scantwo_Result, Result[0]); /* full */
    min3d_lowertri(n_pos, n_perm, scantwo_Result, Result[3]); /* add */
    min2d(n_pos, n_perm, scanone_Result, Result[5]); /* scanone */

    for(i=0; i<n_perm; i++) {
        Result[1][i] = Result[0][i] - Result[5][i]; /* fv1 */
        Result[2][i] = Result[0][i] - Result[3][i]; /* int */
        Result[4][i] = Result[3][i] - Result[5][i]; /* av1 */
    }
}



/* scantwo perms within chr; with covariates, must do permutations singly */
void scantwopermhk_1chr(int n_ind, int n_pos, int n_gen,
                        double ***Genoprob, double *****Pairprob,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, int **Permindex,
                        double *weights, double **Result,
                        int n_col2drop, int *col2drop)
{
    int i;
    double *pheno_shuffled, *addcov_shuffled, **Addcov_shuffled;
    double *scanone_result, **scanone_Result;
    double *scantwo_result, ***scantwo_Result;
    int *ind_noqtl;

    /* setup */
    create_zero_vector(&ind_noqtl, n_ind);
    allocate_double(n_pos, &scanone_result);
    reorg_errlod(n_pos, 1, scanone_result, &scanone_Result);
    allocate_double(n_pos*n_pos, &scantwo_result);
    reorg_genoprob(n_pos, n_pos, 1, scantwo_result, &scantwo_Result);
    allocate_double(n_ind, &pheno_shuffled);
    allocate_double(n_ind*n_addcov, &addcov_shuffled);
    reorg_errlod(n_ind, n_addcov, addcov_shuffled, &Addcov_shuffled);

    for(i=0; i<n_perm; i++) {
        fill_covar_and_phe(n_ind, Permindex[i], pheno, Addcov, n_addcov, pheno_shuffled, Addcov_shuffled);

        /* scanone */
        scanone_hk(n_ind, n_pos, n_gen, Genoprob,
                   Addcov_shuffled, n_addcov, 0, 0, /* null inter've covariates */
                   pheno_shuffled, 1, weights, scanone_Result,
                   ind_noqtl);

        /* scantwo */
        scantwo_1chr_hk(n_ind, n_pos, n_gen, Genoprob, Pairprob,
                        Addcov_shuffled, n_addcov, 0, 0, /* null inter've covariates */
                        pheno_shuffled, 1, weights,
                        scantwo_Result, n_col2drop, col2drop);

        min3d_uppertri(n_pos, 1, scantwo_Result, Result[0]+i); /* full */
        min3d_lowertri(n_pos, 1, scantwo_Result, Result[3]+i); /* add */
        min2d(n_pos, 1, scanone_Result, Result[5]+i); /* scanone */
        Result[1][i] = Result[0][i] - Result[5][i]; /* fv1 */
        Result[2][i] = Result[0][i] - Result[3][i]; /* int */
        Result[4][i] = Result[3][i] - Result[5][i]; /* av1 */

    }
}



/* R wrapper for function to perform scantwo permutations by
   Haley-Knott regression for a pair of chromosomes */
void R_scantwopermhk_2chr(int *n_ind, int *n_pos1, int *n_pos2,
                          int *n_gen1, int *n_gen2,
                          double *genoprob1, double *genoprob2,
                          double *addcov, int *n_addcov,
                          double *pheno, int *n_perm, int *permindex,
                          double *weights, double *result)
{
    double ***Genoprob1, ***Genoprob2, **Result, **Addcov=0;
    int **Permindex;

    reorg_genoprob(*n_ind, *n_pos1, *n_gen1, genoprob1, &Genoprob1);
    reorg_genoprob(*n_ind, *n_pos2, *n_gen2, genoprob2, &Genoprob2);
    reorg_errlod(*n_perm, 6, result, &Result);
    reorg_geno(*n_ind, *n_perm, permindex, &Permindex);

    /* reorganize addcov (if not empty) */
    if(*n_addcov > 0) {
        reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);

        scantwopermhk_2chr(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2,
                           Genoprob1, Genoprob2, Addcov, *n_addcov,
                           pheno, *n_perm, Permindex, weights, Result);
    }
    else {
        scantwopermhk_2chr_nocovar(*n_ind, *n_pos1, *n_pos2, *n_gen1, *n_gen2,
                                   Genoprob1, Genoprob2,
                                   pheno, *n_perm, Permindex, weights, Result);
    }
}


/* scantwo perms for chr pair; with no covariates, can do calculations in batch */
void scantwopermhk_2chr_nocovar(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                                int n_gen2, double ***Genoprob1, double ***Genoprob2,
                                double *pheno, int n_perm, int **Permindex, double *weights,
                                double **Result)
{
    int i;
    double *phematrix, **Phematrix;
    double *scanone_result1, **scanone_Result1;
    double *scanone_result2, **scanone_Result2;
    double *scantwo_result_full, ***scantwo_Result_Full;
    double *scantwo_result_add, ***scantwo_Result_Add;
    int *ind_noqtl;

    /* setup */
    allocate_double(n_perm*n_ind, &phematrix);
    reorg_errlod(n_ind, n_perm, phematrix, &Phematrix);
    create_zero_vector(&ind_noqtl, n_ind);
    allocate_double(n_perm*n_pos1, &scanone_result1);
    reorg_errlod(n_pos1, n_perm, scanone_result1, &scanone_Result1);
    allocate_double(n_perm*n_pos2, &scanone_result2);
    reorg_errlod(n_pos2, n_perm, scanone_result2, &scanone_Result2);
    allocate_double(n_perm*n_pos1*n_pos2, &scantwo_result_full);
    reorg_genoprob(n_pos2, n_pos1, n_perm, scantwo_result_full, &scantwo_Result_Full);
    allocate_double(n_perm*n_pos1*n_pos2, &scantwo_result_add);
    reorg_genoprob(n_pos1, n_pos2, n_perm, scantwo_result_add, &scantwo_Result_Add);

    /* order phenotypes as shuffled*/
    fill_phematrix(n_ind, n_perm, pheno, Permindex, Phematrix);

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

    min2d(n_pos1, n_perm, scanone_Result1, Result[0]);
    min2d(n_pos2, n_perm, scanone_Result2, Result[5]);
    for(i=0; i<n_perm; i++)
        if(Result[0][i] < Result[5][i])
            Result[5][i] = Result[0][i];

    min3d(n_pos2, n_pos1, n_perm, scantwo_Result_Full, Result[0]);
    min3d(n_pos1, n_pos2, n_perm, scantwo_Result_Add, Result[3]);
    for(i=0; i<n_perm; i++) {
        Result[1][i] = Result[0][i] - Result[5][i]; /* fv1 */
        Result[2][i] = Result[0][i] - Result[3][i]; /* int */
        Result[4][i] = Result[3][i] - Result[5][i]; /* av1 */
    }
}


/* scantwo perms for chr pair; with covariates, must do perms singly */
void scantwopermhk_2chr(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                        int n_gen2, double ***Genoprob1, double ***Genoprob2,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, int **Permindex, double *weights, double **Result)
{
    int i;
    double *pheno_shuffled, *addcov_shuffled, **Addcov_shuffled;
    double *scanone_result1, **scanone_Result1;
    double *scanone_result2, **scanone_Result2;
    double *scantwo_result_full, ***scantwo_Result_Full;
    double *scantwo_result_add, ***scantwo_Result_Add;
    int *ind_noqtl;

    /* setup */
    create_zero_vector(&ind_noqtl, n_ind);
    allocate_double(n_pos1, &scanone_result1);
    reorg_errlod(n_pos1, 1, scanone_result1, &scanone_Result1);
    allocate_double(n_pos2, &scanone_result2);
    reorg_errlod(n_pos2, 1, scanone_result2, &scanone_Result2);
    allocate_double(n_pos1*n_pos2, &scantwo_result_full);
    reorg_genoprob(n_pos2, n_pos1, 1, scantwo_result_full, &scantwo_Result_Full);
    allocate_double(n_pos1*n_pos2, &scantwo_result_add);
    reorg_genoprob(n_pos1, n_pos2, 1, scantwo_result_add, &scantwo_Result_Add);
    allocate_double(n_ind, &pheno_shuffled);
    allocate_double(n_ind*n_addcov, &addcov_shuffled);
    reorg_errlod(n_ind, n_addcov, addcov_shuffled, &Addcov_shuffled);

    for(i=0; i<n_perm; i++) {
        fill_covar_and_phe(n_ind, Permindex[i], pheno, Addcov, n_addcov, pheno_shuffled, Addcov_shuffled);

        /* scanone */
        scanone_hk(n_ind, n_pos1, n_gen1, Genoprob1,
                   Addcov_shuffled, n_addcov, 0, 0, /* null inter've covariates */
                   pheno_shuffled, 1, weights, scanone_Result1,
                   ind_noqtl);

        scanone_hk(n_ind, n_pos2, n_gen2, Genoprob2,
                   Addcov_shuffled, n_addcov, 0, 0, /* null inter've covariates */
                   pheno_shuffled, 1, weights, scanone_Result2,
                   ind_noqtl);

        /* scantwo */
        scantwo_2chr_hk(n_ind, n_pos1, n_pos2, n_gen1, n_gen2,
                        Genoprob1, Genoprob2,
                        Addcov_shuffled, n_addcov, 0, 0, /* null inter've covariates */
                        pheno_shuffled, 1, weights,
                        scantwo_Result_Full, scantwo_Result_Add);


        min2d(n_pos1, 1, scanone_Result1, Result[0]+i);
        min2d(n_pos2, 1, scanone_Result2, Result[5]+i);
        if(Result[0][i] < Result[5][i])
            Result[5][i] = Result[0][i];

        min3d(n_pos2, n_pos1, 1, scantwo_Result_Full, Result[0]+i);
        min3d(n_pos1, n_pos2, 1, scantwo_Result_Add, Result[3]+i);
        Result[1][i] = Result[0][i] - Result[5][i]; /* fv1 */
        Result[2][i] = Result[0][i] - Result[3][i]; /* int */
        Result[4][i] = Result[3][i] - Result[5][i]; /* av1 */

    }

}



/* fill matrix with shuffled phenotypes */
void fill_phematrix(int n_ind, int n_perm, double *pheno, int **Permindex, double **Phematrix)
{
    int i,j;

    for(i=0; i<n_ind; i++)
        for(j=0; j<n_perm; j++)
            Phematrix[j][i] = pheno[Permindex[j][i]];
}

/* fill covar and phe with shuffled versions */
void fill_covar_and_phe(int n_ind, int *Permindex, double *pheno, double **Addcov,
                        int n_addcov, double *pheno_shuffled, double **Addcov_shuffled)
{
    int i,j;
    for(i=0; i<n_ind; i++) {
        pheno_shuffled[i] = pheno[Permindex[i]];
        for(j=0; j<n_addcov; j++)
            Addcov_shuffled[j][i] = Addcov[j][Permindex[i]];
    }
}

/* create a vector of 0's */
void create_zero_vector(int **vector, int n)
{
    int i;

    allocate_int(n, vector);

    for(i=0; i<n; i++) (*vector)[i] = 0;
}

/* minimize over first two dimensions for each value of 3rd and place in results */
void min3d(int d1, int d2, int d3, double ***Values, double *results)
{
    int i, j, k;

    for(k=0; k<d3; k++) {
        results[k] = Values[k][0][0];
        for(i=0; i<d1; i++) {
            for(j=0; j<d2; j++) {
                if(Values[k][j][i] < results[k])
                    results[k] = Values[k][j][i];
            }
        }
    }
}


/* minimize over first dimension for each value of 2nd and place in results */
void min2d(int d1, int d2, double **Values, double *results)
{
    int i, k;

    for(k=0; k<d2; k++) {
        results[k] = Values[k][0];
        for(i=0; i<d1; i++) {
            if(Values[k][i] < results[k])
                results[k] = Values[k][i];
        }
    }
}


/* minimize over upper triangle of square matrix d1xd1 for each value of d3 */
void min3d_uppertri(int d1, int d3, double ***Values, double *results)
{
    int i, j, k;

    for(k=0; k<d3; k++) {
        results[k] = R_PosInf;
        for(i=0; i<d1; i++) {
            for(j=i+1; j<d1; j++) {
                if(Values[k][i][j] < results[k])
                    results[k] = Values[k][i][j];
            }
        }
    }
}


/* minimize over lower triangle of square matrix d1xd1 for each value of d3 */
void min3d_lowertri(int d1, int d3, double ***Values, double *results)
{
    int i, j, k;

    for(k=0; k<d3; k++) {
        results[k] = R_PosInf;
        for(i=0; i<d1; i++) {
            for(j=i+1; j<d1; j++) {
                if(Values[k][j][i] < results[k])
                    results[k] = Values[k][j][i];
            }
        }
    }
}
