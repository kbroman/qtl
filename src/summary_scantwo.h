/**********************************************************************
 * 
 * summary_scantwo.h
 *
 * copyright (c) 2006, Karl W Broman
 *
 * last modified Oct, 2006
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
		       double *lod_1qtl);

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
 * Pos1_add, Pos2_add    On output, positions of maximum add've LOD
 *                       Matrices indexed as [phe][chrpair]
 * Pos1_full, Pos2_full  On output, positions of maximum joint LOD
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
		     double **LOD_1qtl);

/* end of summary_scantwo.h */

