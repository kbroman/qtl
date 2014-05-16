/**
 * perm_strat.h: code to do a stratified permutation
 */

void validate_strata(int n, int n_strata, int *strata);

void permute_by_strata_int(int n, int *index, int n_strata, int *strata);

void permute_by_strata_double(int n, double *index, int n_strata, int *strata);

void R_permute_by_strata_int(int *n, int *index, int *n_strata, int *strata);

void R_permute_by_strata_double(int *n, double *index, int *n_strata, int *strata);
