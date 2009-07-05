/**********************************************************************
 *
 * MQMdata.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Feb, 2009
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
 * Several datastructures needed by the MQM algorithm are defined here
 * Contains:
 *
 **********************************************************************/


#ifdef __cplusplus
  extern "C" {
#endif

/*------------------------------------------------------------------------
Datastructures for matrix and vector calculus
------------------------------------------------------------------------ */
typedef double*** Mmatrix;
typedef double** matrix;
typedef double*  vector;
typedef char**   cmatrix;
typedef char*    cvector;
typedef int*  ivector;

char determin_cross(int *Nmark,int *Nind,int **Geno,int *crosstype);

void change_coding(int *Nmark,int *Nind,int **Geno,cmatrix markers, int crosstype);


/*------------------------------------------------------------------------
Basic routines for matrix and vector calculus
------------------------------------------------------------------------ */
vector newvector(int dim);
ivector newivector(int dim);
cvector newcvector(int dim);
matrix newmatrix(int rows, int cols);
Mmatrix newMmatrix(int rows, int cols,int depth);
void   printmatrix(matrix m, int rows, int cols);
void   printcmatrix(cmatrix m, int rows, int cols);
cmatrix newcmatrix(int rows, int cols);
void delmatrix(matrix m);
void delMmatrix(Mmatrix m);
void delcmatrix(cmatrix m);
void copyvector(vector vsource, vector vdestination, int dim);

#ifdef __cplusplus
  }
#endif

/* end of MQMdata.h */
