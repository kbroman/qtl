/**********************************************************************
 *
 * MQMdatatypes.cpp
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

#include "MQM.h"

extern "C" {

  char determin_cross(int *Nmark,int *Nind,int **Geno,int *crosstype) {
    for (int i=0; i< *Nmark; i++) {
      for (int j=0; j< *Nind; j++) {
        //Some checks to see if the cross really is the cross we got (So BC can't contain 3's (BB) and RILS can't contain 2's (AB)
        if (Geno[i][j] != 9 && Geno[i][j] > 3 && (*crosstype) != 1) {
          Rprintf("INFO: Strange genotype pattern, switching to F2\n");
          Rprintf("ind = %d marker = %d Geno = %d\n", i+1, j+1, Geno[i][j]);
          (*crosstype) = 1;
          break;
        }
        if (Geno[i][j] == 3 && (*crosstype) == 2) {
          Rprintf("INFO: Strange genotype pattern, switching from BC to F2\n");
          (*crosstype) = 1;
          break;
        }
        //IF we have a RIL and find AB then the set is messed up; we have a BC genotype
        if (Geno[i][j] == 2 && (*crosstype) == 3) {
          Rprintf("INFO: Strange genotype pattern, switching from RISELF to BC\n");
          (*crosstype) = 2;
          break;
        }

      }
      //Rprintf("\n");
    }

    char cross = 'F';
    if ((*crosstype) == 1) {
      cross = 'F';
    }
    if ((*crosstype) == 2) {
      cross = 'B';
    }
    if ((*crosstype) == 3) {
      cross = 'R';
    }
    return cross;
  }


  void change_coding(int *Nmark,int *Nind,int **Geno,cmatrix markers, int crosstype) {
    // Change all the genotypes from default R/qtl format to MQM internal
    for (int i=0; i< *Nmark; i++) {
      for (int j=0; j< *Nind; j++) {
        markers[i][j] = '9';
        if (Geno[i][j] == 1) {				//AA
          markers[i][j] = '0';
        }
        if (Geno[i][j] == 2) {				//AB
          // [karl:] I think this needs to be changed, but my fix doesn't work.
          //			  if(crosstype!=3) markers[i][j] = '1'; // non-RIL
          //			  else markers[i][j] = '2';  // RIL
          markers[i][j] = '1';
        }
        if (Geno[i][j] == 3) {				//BB
          markers[i][j] = '2';
        }
        if (Geno[i][j] == 4) {				//AA of AB
          markers[i][j] = '4';
        }
        if (Geno[i][j] == 5) {				//BB of AB
          markers[i][j] = '3';
        }
      }
    }
  }

  vector newvector(int dim) {
    vector v;
    v = (double *)Calloc(dim, double);
    if (v==NULL) {
      warning("Not enough memory for new vector of dimension %d",(dim+1));
    }
    return v;
  }

  ivector newivector(int dim) {
    ivector v;
    v = (int *)Calloc(dim, int);
    if (v==NULL) {
      warning("Not enough memory for new vector of dimension %d",(dim+1));
    }
    return v;
  }

  cvector newcvector(int dim) {
    cvector v;
    v = (char *)Calloc(dim, char);
    if (v==NULL) {
      warning("Not enough memory for new vector of dimension %d",(dim+1));
    }
    return v;
  }

  matrix newmatrix(int rows, int cols) {
    matrix m;
    m = (double **)Calloc(rows, double*);
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
    m = (double ***)Calloc(rows, double**);
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
    m = (char **)Calloc(rows, char*);
    if (m==NULL) {
      warning("Not enough memory for new char matrix");
    }
    for (int i=0; i<rows; i++) {
      m[i]= newcvector(cols);
    }
    return m;
  }

  void delmatrix(matrix m) {

    Free(m);
  }

  void delMmatrix(Mmatrix m) {

    Free(m);
  }

  void delcmatrix(cmatrix m) {

    Free(m);
  }

  void copyvector(vector vsource, vector vdestination, int dim) {

    for (int i=0; i<dim; i++) {
      vdestination[i]= vsource[i];
    }
  }

}

/* end of MQMdata.c */
