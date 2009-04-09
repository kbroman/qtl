/**********************************************************************
 * 
 * summary_scantwo.c
 *
 * copyright (c) 2006, Karl W Broman
 *
 * last modified Dec, 2006
 * first written Oct, 2006
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
 * These constitute a subroutine for getting the maximum LOD scores
 * for each pair of chromosomes in scantwo output.
 *
 * Contains: R_summary_scantwo, summary_scantwo
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "summary_scantwo.h"

/**********************************************************************
 * 
 * R_summary_scantwo
 *
 * Wrapper for call from R
 * 
 **********************************************************************/

void R_summary_scantwo(int *n_pos, int *n_phe, double *lod, int *n_chr, 
		       int *chr, double *pos, int *xchr, double *scanoneX, 
		       int *n_chrpair, int *chr1, int *chr2, int *chrpair, 
		       double *pos1_jnt, double *pos2_jnt, 
		       double *pos1_add, double *pos2_add,
		       double *pos1_int, double *pos2_int, 
		       double *jnt_lod_full, double *jnt_lod_add, 
		       double *add_lod_full, double *add_lod_add, 
		       double *int_lod_full, double *int_lod_add, 
		       double *lod_1qtl)
{
  double ***Lod, **ScanoneX;
  double **Pos1_jnt, **Pos2_jnt;
  double **Pos1_add, **Pos2_add;
  double **Pos1_int, **Pos2_int;
  double **JNT_lod_full, **JNT_lod_add;
  double **ADD_lod_full, **ADD_lod_add;
  double **INT_lod_full, **INT_lod_add;
  double **LOD_1qtl;
  int i, j, k, **Chrpair;

  *n_chrpair = *n_chr*(*n_chr+1)/2;
  /* re-organize matrices */
  reorg_genoprob(*n_pos, *n_pos, *n_phe, lod, &Lod);
  reorg_errlod(*n_chrpair, *n_phe, pos1_jnt, &Pos1_jnt);
  reorg_errlod(*n_chrpair, *n_phe, pos2_jnt, &Pos2_jnt);
  reorg_errlod(*n_chrpair, *n_phe, pos1_add, &Pos1_add);
  reorg_errlod(*n_chrpair, *n_phe, pos2_add, &Pos2_add);
  reorg_errlod(*n_chrpair, *n_phe, pos1_int, &Pos1_int);
  reorg_errlod(*n_chrpair, *n_phe, pos2_int, &Pos2_int);
  reorg_errlod(*n_chrpair, *n_phe, jnt_lod_full, &JNT_lod_full);
  reorg_errlod(*n_chrpair, *n_phe, jnt_lod_add, &JNT_lod_add);
  reorg_errlod(*n_chrpair, *n_phe, add_lod_full, &ADD_lod_full);
  reorg_errlod(*n_chrpair, *n_phe, add_lod_add, &ADD_lod_add);
  reorg_errlod(*n_chrpair, *n_phe, int_lod_full, &INT_lod_full);
  reorg_errlod(*n_chrpair, *n_phe, int_lod_add, &INT_lod_add);
  reorg_errlod(*n_chrpair, *n_phe, lod_1qtl, &LOD_1qtl);
  reorg_errlod(*n_pos, *n_phe, scanoneX, &ScanoneX);
  reorg_geno(*n_chr, *n_chr, chrpair, &Chrpair);

  for(i=0, k=0; i<*n_chr; i++) {
    for(j=i; j<*n_chr; j++, k++) {
      chr1[k] = i;
      chr2[k] = j;
      Chrpair[j][i] = Chrpair[i][j] = k;
    }
  }

  summary_scantwo(*n_pos, *n_phe, Lod, *n_chr, chr, pos, xchr, 
		  ScanoneX, *n_chrpair, Chrpair, 
		  Pos1_jnt, Pos2_jnt,
		  Pos1_add, Pos2_add, 
		  Pos1_int, Pos2_int, 
		  JNT_lod_full, JNT_lod_add,
		  ADD_lod_full, ADD_lod_add,
		  INT_lod_full, INT_lod_add,
		  LOD_1qtl);
}


/**********************************************************************
 * 
 * summary_scantwo
 *
 * Function to pull out the major bits from scantwo output
 *
 * n_pos: Total number of positions
 * n_phe: Number of phenotype columns
 * Lod:   Array of LOD scores indexed as [phe][pos2][pos1]
 *        diagonal = scanone results; upper.tri = add've LOD; lower = full
 * n_chr  Number of distinct chromosomes
 * chr    Index of chromosomes; length n_pos, taking values in 1..n_chr
 * pos    cM positions; length n_pos
 * xchr   Index of xchr; length n_chr; 0=autosome, 1=X chromosome
 * ScanoneX   special X scanone; matrix indexed as [phe][pos]
 *
 * n_chrpair             Number of pairs of chromosomes
 * Chrpair               Matrix giving chrpair index for a pair of chromosomes
 * Pos1_jnt, Pos2_jnt    On output, positions of maximum joint LOD
 *                       Matrices indexed as [phe][chrpair]
 * Pos1_add, Pos2_add    On output, positions of maximum add've LOD
 *                       Matrices indexed as [phe][chrpair]
 * Pos1_int, Pos2_int    On output, positions of maximum int've LOD
 *                       Matrices indexed as [phe][chrpair]
 * JNT_lod_*             On output, joint and add've LOD at pos'ns with 
 *                       maximum joint LOD; matrices indexed as [phe][chrpair]
 * ADD_lod_*             On output, joint and add've LOD at pos'ns with 
 *                       maximum add've LOD; matrices indexed as [phe][chrpair]
 * INT_lod_*             On output, joint and add've LOD at pos'ns with 
 *                       maximum int've LOD; matrices indexed as [phe][chrpair]
 * LOD_1qtl              On output, maximum 1-QTL LOD for each pair of chr
 *                       (selected from either scanone or scanoneX
 * 
 **********************************************************************/
void summary_scantwo(int n_pos, int n_phe, double ***Lod, int n_chr, 
		     int *chr, double *pos, int *xchr, double **ScanoneX, 
		     int n_chrpair, int **Chrpair, 
		     double **Pos1_jnt, double **Pos2_jnt, 
		     double **Pos1_add, double **Pos2_add, 
		     double **Pos1_int, double **Pos2_int, 
		     double **JNT_lod_full, double **JNT_lod_add,
		     double **ADD_lod_full, double **ADD_lod_add,
		     double **INT_lod_full, double **INT_lod_add,
		     double **LOD_1qtl)
{
  int i, j, k, c1, c2, thepair;
  double *maxone, *maxoneX;

  allocate_double(n_chr, &maxone);
  allocate_double(n_chr, &maxoneX);

  for(i=0; i<n_phe; i++) { /* loop over phenotype columns */

    /* get maximum single-QTL results */
    for(j=0; j<n_chr; j++) 
      maxone[j] = maxoneX[j] = 0.0;
  
    for(j=0; j<n_pos; j++) {
      if(Lod[i][j][j] > maxone[chr[j]-1])
	maxone[chr[j]-1] = Lod[i][j][j];

      if(ScanoneX[i][j] > maxoneX[chr[j]-1])
	maxoneX[chr[j]-1] = ScanoneX[i][j];
    }

    /* zero out the matrices for the maximum LOD scores */
    for(j=0; j<n_chrpair; j++) {
      JNT_lod_full[i][j] = JNT_lod_add[i][j] = 
	ADD_lod_full[i][j] = ADD_lod_add[i][j] = 
	INT_lod_full[i][j] = INT_lod_add[i][j] = 
	Pos1_add[i][j] = Pos2_add[i][j] =
	Pos1_int[i][j] = Pos2_int[i][j] =
	Pos1_jnt[i][j] = Pos2_jnt[i][j] = 0.0; 
    }

    /* maximum joint, add've, and int've LOD scores */
    for(j=0; j<n_pos; j++) {
      for(k=j; k<n_pos; k++) {
	R_CheckUserInterrupt(); /* check for ^C */
	c1 = chr[j]-1;
	c2 = chr[k]-1;
	thepair = Chrpair[c1][c2];

	if(Lod[i][j][k] > JNT_lod_full[i][thepair]) {
	  JNT_lod_full[i][thepair] = Lod[i][j][k];
	  JNT_lod_add[i][thepair] = Lod[i][k][j];
	  Pos1_jnt[i][thepair] = pos[j];
	  Pos2_jnt[i][thepair] = pos[k];
	}

	if(Lod[i][k][j] > ADD_lod_add[i][thepair]) {
	  ADD_lod_add[i][thepair] = Lod[i][k][j];
	  ADD_lod_full[i][thepair] = Lod[i][j][k];
	  Pos1_add[i][thepair] = pos[j];
	  Pos2_add[i][thepair] = pos[k];
	}

	if(Lod[i][j][k] - Lod[i][k][j] > 
	   INT_lod_full[i][thepair] - INT_lod_add[i][thepair]) {
	  INT_lod_full[i][thepair] = Lod[i][j][k];
	  INT_lod_add[i][thepair] = Lod[i][k][j];
	  Pos1_int[i][thepair] = pos[j];
	  Pos2_int[i][thepair] = pos[k];
	}
      }
    }

    /* pull out biggest single-QTL LOD scores */
    for(j=0; j<n_chr; j++) {
      for(k=j; k<n_chr; k++) {
	R_CheckUserInterrupt(); /* check for ^C */

	thepair = Chrpair[j][k];

	/* calculate the two 2 v 1 conditional LOD scores */
	if(xchr[j] || xchr[k]) {
	  if(maxoneX[j] > maxoneX[k]) LOD_1qtl[i][thepair] = maxoneX[j];
	  else LOD_1qtl[i][thepair] = maxoneX[k];
	}
	else {
	  if(maxone[j] > maxone[k]) LOD_1qtl[i][thepair] = maxone[j];
	  else LOD_1qtl[i][thepair] = maxone[k];
	}
      }
    }

  } /* end loop over phenotype columns */
}





/* end of summary_scantwo.c */

