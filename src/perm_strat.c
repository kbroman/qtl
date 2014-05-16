/**
 * perm_strat.c: code to do a stratified permutation
 */

#include <R.h>
#include "util.h"
#include "perm_strat.h"

/* check that the strata vector is okay */
void validate_strata(int n, int n_strata, int *strata)
{
  int i;

  for(i=0; i<n; i++) {
    if(strata[i] < 0 || strata[i] > n_strata-1)
      error("strata values must be in 0, ..., n_strata-1");
  }
}


/**
 * shuffle indices within strata 
 *
 * n = length of index and strata
 * index = vector to be shuffled
 * n_strata = number of strata
 * strata = vector of stratum assignments; values in {0, ..., n_strata-1}
 **/

void permute_by_strata_int(int n, int *index, int n_strata, int *strata)
{
  int i, *n_per_strata, *curindex, **split_index, this_stratum;

  /* allocate space */
  allocate_int(n_strata, &n_per_strata);
  allocate_int(n_strata, &curindex);
  split_index = (int **)R_alloc(n_strata, sizeof(int *));
  split_index[0] = (int *)R_alloc(n, sizeof(int));
  
  /* count strata */
  for(i=0; i<n_strata; i++) 
    n_per_strata[i] = 0;
  for(i=0; i<n; i++)
    n_per_strata[strata[i]]++;

  /* split index pointers */
  for(i=1; i<n_strata; i++)
    split_index[i] = split_index[i-1] + n_per_strata[i-1];

  /* fill up indices */
  for(i=0; i<n_strata; i++)
    curindex[i] = 0;
  for(i=0; i<n; i++) {
    this_stratum = strata[i];
    split_index[this_stratum][curindex[this_stratum]] = index[i];
    curindex[this_stratum]++;
  }

  /* permute within strata */
  for(i=0; i<n_strata; i++)
    int_permute(split_index[i], n_per_strata[i]);

  /* unsplit */
  for(i=0; i<n_strata; i++)
    curindex[i] = 0;
  for(i=0; i<n; i++) {
    this_stratum = strata[i];
    index[i] = split_index[this_stratum][curindex[this_stratum]];
    curindex[this_stratum]++;
  }

}


void permute_by_strata_double(int n, double *index, int n_strata, int *strata)
{
  int i, *n_per_strata, *curindex, this_stratum;
  double **split_index;

  /* allocate space */
  allocate_int(n_strata, &n_per_strata);
  allocate_int(n_strata, &curindex);
  split_index = (double **)R_alloc(n_strata, sizeof(double *));
  split_index[0] = (double *)R_alloc(n, sizeof(double));
  
  /* count strata */
  for(i=0; i<n_strata; i++) 
    n_per_strata[i] = 0;
  for(i=0; i<n; i++)
    n_per_strata[strata[i]]++;

  /* split index pointers */
  for(i=1; i<n_strata; i++)
    split_index[i] = split_index[i-1] + n_per_strata[i-1];

  /* fill up indices */
  for(i=0; i<n_strata; i++)
    curindex[i] = 0;
  for(i=0; i<n; i++) {
    this_stratum = strata[i];
    split_index[this_stratum][curindex[this_stratum]] = index[i];
    curindex[this_stratum]++;
  }

  /* permute within strata */
  for(i=0; i<n_strata; i++)
    double_permute(split_index[i], n_per_strata[i]);

  /* unsplit */
  for(i=0; i<n_strata; i++)
    curindex[i] = 0;
  for(i=0; i<n; i++) {
    this_stratum = strata[i];
    index[i] = split_index[this_stratum][curindex[this_stratum]];
    curindex[this_stratum]++;
  }

}


/* for testing from R */
void R_permute_by_strata_int(int *n, int *index, int *n_strata, int *strata)
{
  GetRNGstate();
  
  /* check that the strata vector is ok */
  validate_strata(*n, *n_strata, strata);
  
  /* stratified permutation */
  permute_by_strata_int(*n, index, *n_strata, strata);

  PutRNGstate();
}

void R_permute_by_strata_double(int *n, double *index, int *n_strata, int *strata)
{
  GetRNGstate();
  
  /* check that the strata vector is ok */
  validate_strata(*n, *n_strata, strata);
  
  /* stratified permutation */
  permute_by_strata_double(*n, index, *n_strata, strata);

  PutRNGstate();
}
