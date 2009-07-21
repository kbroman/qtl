/**********************************************************************
 *
 * mqmdatatypes.cpp
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Apr, 2009
 * first written Feb, 2009
 *
 * Original version R.C Jansen
 * first written <2000 (unknown)
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
 * Basic datatypes used by R/qtl-MQM
 *
 **********************************************************************/

#include "mqm.h"

/* 
 * Determine the experimental cross type from the R/qtl dataset. Returns the
 * type.
 */

char determine_cross(int *Nmark,int *Nind,int **Geno,int *crosstype) {
  for (int i=0; i< *Nmark; i++) {
    for (int j=0; j< *Nind; j++) {
      //Some checks to see if the cross really is the cross we got (So BC can't contain 3's (BB) and RILS can't contain 2's (AB)
      if (Geno[i][j] != 9 && Geno[i][j] > 3 && (*crosstype) != 1) {
        Rprintf("INFO: Unexpected genotype pattern, switching to F2\n");
        Rprintf("ind = %d marker = %d Geno = %d\n", i+1, j+1, Geno[i][j]);
        (*crosstype) = 1;
        break;
      }
      if (Geno[i][j] == 3 && (*crosstype) == 2) {
        Rprintf("INFO: Unexpected genotype pattern, switching from BC to F2\n");
        (*crosstype) = 1;
        break;
      }
      //IF we have a RIL and find AB then the set is messed up; we have a BC genotype
      if (Geno[i][j] == 2 && (*crosstype) == 3) {
        Rprintf("INFO: Unexpected genotype pattern, switching from RISELF to BC\n");
        (*crosstype) = 2;
        break;
      }

    }
    //Rprintf("\n");
  }

  unsigned char cross = 0;
  switch(*crosstype) {
    case 1: cross=CF2;
            break;
    case 2: cross=CBC;
            break;
    case 3: cross=CRIL;
            break;
    default:
            error("Unknown cross type %d",*crosstype);
  }

  return cross;
}

/*
 * Change all the genotypes from default R/qtl format to MQM internal.  R/qtl
 * uses internally { 'AA' => 1, 'H' => 2, 'BB' => 3, 'NOTBB' => 4, 'NOTAA' => 5,
 * '-' => 'NA' }
 *
 */

void change_coding(int *Nmark,int *Nind,int **Geno,cmatrix markers, int crosstype) {
  for (int i=0; i< *Nmark; i++) {
    for (int j=0; j< *Nind; j++) {
      switch (Geno[i][j]) {
        case 1: markers[i][j] = MAA;
                break;
        case 2: markers[i][j] = MH;
          // [karl:] I think this needs to be changed, but my fix doesn't work.
          //			  if(crosstype!=3) markers[i][j] = MH; // non-RIL
          //			  else markers[i][j] = MBB;  // RIL
                break;
        case 3: markers[i][j] = MBB;
                break;
        case 4: markers[i][j] = MNOTBB;
                break;
        case 5: markers[i][j] = MNOTAA;
                break;
        case 9: markers[i][j] = MMISSING;
                break;
        default:
                error("Can not convert R/qtl genotype with value %d",Geno[i][j]);
      }
    }
  }
}

/* 
 * Allocate a memory block using the 'safe' R method calloc_init, but with
 * guarantee all data has been zeroed
 */

void *calloc_init(size_t num, size_t size) {
  void *buf;
  buf = R_chk_calloc(num,size);
  if (buf) memset(buf,0,num*size);
  return buf;
}

vector newvector(int dim) {
  vector v;
  v = (double *)calloc_init(dim, sizeof(double));
  if (v==NULL) {
    warning("Not enough memory for new vector of dimension %d",(dim+1));
  }
  return v;
}

ivector newivector(int dim) {
  ivector v;
  v = (int *)calloc_init(dim, sizeof(int));
  if (v==NULL) {
    warning("Not enough memory for new vector of dimension %d",(dim+1));
  }
  return v;
}

cvector newcvector(int dim) {
  cvector v;
  v = (char *)calloc_init(dim, sizeof(char));
  if (v==NULL) {
    warning("Not enough memory for new vector of dimension %d",(dim+1));
  }
  return v;
}

matrix newmatrix(int rows, int cols) {
  matrix m;
  m = (double **)calloc_init(rows, sizeof(double*));
  if (m==NULL) {
    warning("Not enough memory for new double matrix");
  }
  for (int i=0; i<rows; i++) {
    m[i]= newvector(cols);
  }
  return m;
}

Mmatrix newMmatrix(int rows, int cols,int depth) {
  Mmatrix m;
  m = (double ***)calloc_init(rows, sizeof(double**));
  if (m==NULL) {
    warning("Not enough memory for new double matrix");
  }
  for (int i=0; i<rows; i++) {
    m[i]= newmatrix(cols,depth);
  }
  return m;
}

void printmatrix(matrix m, int rows, int cols) {

  for (int r=0; r<rows; r++) {
    for (int c=0; c<cols; c++) {
      Rprintf("%f\t",m[r][c]);
    }
    Rprintf("\n");
  }
}

void printcmatrix(cmatrix m, int rows, int cols) {

  for (int r=0; r<rows; r++) {
    for (int c=0; c<cols; c++) {
      Rprintf("%c\t",m[r][c]);
    }
    Rprintf("\n");
  }
}

cmatrix newcmatrix(int rows, int cols) {
  cmatrix m;
  m = (char **)calloc_init(rows, sizeof(char*));
  if (m==NULL) {
    warning("Not enough memory for new char matrix");
  }
  for (int i=0; i<rows; i++) {
    m[i]= newcvector(cols);
  }
  return m;
}

void freevector(void *v) {
  Free(v);
}

void freematrix(void **m, size_t rows) {
  for (int i=0; i<rows; i++) {
    Free(m[i]);
  }
  Free(m);
}

void delmatrix(matrix m, size_t rows) {
  freematrix((void**)m,rows);
}

void delMmatrix(Mmatrix m, size_t rows) {
  freematrix((void**)m,rows);
}

void delcmatrix(cmatrix m, size_t rows) {
  freematrix((void **)m,rows);
}

void copyvector(vector vsource, vector vdestination, int dim) {

  for (int i=0; i<dim; i++) {
    vdestination[i]= vsource[i];
  }
}
