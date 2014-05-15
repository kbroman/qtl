/**
 * scantwo_hk.h
 **/

void R_scantwopermhk_1chr(int *n_ind, int *n_pos, int *n_gen,
                          double *genoprob, double *pairprob,
                          double *addcov, int *n_addcov,
                          double *pheno, int* n_perm, double *weights,
                          double *result, int *n_col2drop, int *col2drop);

void scantwopermhk_1chr_nocovar(int n_ind, int n_pos, int n_gen,
                                double ***Genoprob, double *****Pairprob,
                                double *pheno, int n_perm, double *weights,
                                double **Result, int n_col2drop, int *col2drop);

void scantwopermhk_1chr(int n_ind, int n_pos, int n_gen,
                        double ***Genoprob, double *****Pairprob,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, double *weights, double **Result,
                        int n_col2drop, int *col2drop);

void R_scantwopermhk_2chr(int *n_ind, int *n_pos1, int *n_pos2,
                          int *n_gen1, int *n_gen2,
                          double *genoprob1, double *genoprob2,
                          double *addcov, int *n_addcov,
                          double *pheno, int *n_perm, double *weights,
                          double *result);

void scantwopermhk_2chr_nocovar(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                                int n_gen2, double ***Genoprob1, double ***Genoprob2,
                                double *pheno, int n_perm, double *weights,
                                double **Result);

void scantwopermhk_2chr(int n_ind, int n_pos1, int n_pos2, int n_gen1,
                        int n_gen2, double ***Genoprob1, double ***Genoprob2,
                        double **Addcov, int n_addcov, double *pheno,
                        int n_perm, double *weights, double **Result);

void create_shuffled_phematrix(int n_ind, int n_perm, double *pheno, double **phematrix);

void shuffle_covar_and_phe(int n_ind, int *ind_index, double *pheno,
                           double **Addcov, int n_addcov,
                           double *dwork);

void create_zero_vector(int **vector, int n);
