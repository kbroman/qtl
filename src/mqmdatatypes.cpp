/**********************************************************************
 *
 * mqmdatatypes.cpp
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Danny Arends and Pjotr Prins
 * last modified November 2009
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
 **********************************************************************/

#include "mqm.h"

/*
 * Determine the experimental cross type from the R/qtl dataset. Returns the
 * type.
 */

MQMCrossType determine_MQMCross(const int Nmark, const int Nind, const int **Geno, const RqtlCrossType rqtlcrosstype) {
  MQMCrossType crosstype = CUNKNOWN;
  switch(rqtlcrosstype) {
    case RC_F2: crosstype = CF2;
                break;
    case RC_RIL: crosstype = CRIL;
                break;
    case RC_BC: crosstype = CBC;
                break;
  }
  for (int j=0; j < Nmark; j++) {
    for (int i=0; i < Nind; i++) {
      //Some checks to see if the cross really is the cross we got (So BC can't contain 3's (BB) and RILS can't contain 2's (AB)
      if (Geno[j][i] != 9 && Geno[j][i] > 3 && (rqtlcrosstype) != 1) {
        Rprintf("INFO: ind = %d marker = %d Geno = %d\n", i+1, j+1, Geno[j][i]);
        Rprintf("INFO: Unexpected genotype pattern, switching to F2\n");
        crosstype = CF2;
        break;
      }
      if (Geno[j][i] == 3 && (rqtlcrosstype) == 2) {
        Rprintf("INFO: Unexpected genotype pattern, switching from BC to F2\n");
        crosstype = CF2;
        break;
      }
      //IF we have a RIL and find AB then the set is messed up; we have a BC genotype
      if (Geno[j][i] == 2 && (rqtlcrosstype) == 3) {
        Rprintf("INFO: Unexpected genotype pattern, switching from RIL to BC\n");
        crosstype = CBC;
        break;
      }
    }
  }
  return crosstype;
}

/*
 * Change all the genotypes from default R/qtl format to MQM internal.  R/qtl
 * uses internally { 'AA' => 1, 'H' => 2, 'BB' => 3, 'NOTBB' => 4, 'NOTAA' => 5,
 * '-' => 'NA' }
 *
 */

void change_coding(int *Nmark, int *Nind, int **Geno, MQMMarkerMatrix markers, const MQMCrossType crosstype) {
  //info("Convert codes R/qtl -> MQM");
  for (int j=0; j < *Nmark; j++) {
    for (int i=0; i < *Nind; i++) {
      switch (Geno[j][i]) {
        case 1: markers[j][i] = MAA;
                break;
        case 2: markers[j][i] = MH;
                if (crosstype == CRIL) markers[j][i] = MBB; // FIXME test
                break;
        case 3: markers[j][i] = MBB;
                break;
        case 4: markers[j][i] = MNOTBB;
                break;
        case 5: markers[j][i] = MNOTAA;
                break;
        case 9: markers[j][i] = MMISSING;
                break;
        default:
                error("Can not convert R/qtl genotype with value %d",Geno[j][i]);
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
  buf = S_alloc(num,(int)size);
  return buf;
}

vector newvector(int dim) {
  vector v = (double *)calloc_init(dim, sizeof(double));
  if(!v){ warning("Not enough memory for new vector of dimension %d", (dim+1)); }
  return v;
}

ivector newivector(int dim) {
  ivector v = (int *)calloc_init(dim, sizeof(int));
  if(!v){ warning("Not enough memory for new vector of dimension %d", (dim+1)); }
  return v;
}

cvector newcvector(int dim) {
  cvector v = (char *)calloc_init(dim, sizeof(char));
  if(!v){ warning("Not enough memory for new vector of dimension %d", (dim+1)); }
  return v;
}

MQMMarkerVector newMQMMarkerVector(int dim) {
  MQMMarkerVector v = (MQMMarker *)calloc_init(dim, sizeof(MQMMarker));
  if(!v){ warning("Not enough memory for new mqm marker vector of dimension %d", (dim+1)); }
  return v;
}

relmarkerarray newRelMarkerPos(int dim){
  relmarkerarray v = (MQMRelMarkerPos *)calloc_init(dim, sizeof(char));
  if(!v){ warning("Not enough memory for the relative marker position vector with dimension %d",(dim+1)); }
  return v;
}

matrix newmatrix(int rows, int cols) {
  matrix m = (double **)calloc_init(rows, sizeof(double*));
  if(!m){ warning("Not enough memory for new double matrix"); }
  for (int i=0; i<rows; i++) {
    m[i]= newvector(cols);
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
  cmatrix m = (char **)calloc_init(rows, sizeof(char*));
  if(!m){ warning("Not enough memory for new char matrix"); }
  for (int i=0; i<rows; i++) {
    m[i]= newcvector(cols);
  }
  return m;
}


MQMMarkerMatrix newMQMMarkerMatrix(int rows, int cols) {
  MQMMarkerMatrix m = (MQMMarkerMatrix)calloc_init(rows, sizeof(MQMMarkerVector));
  if(!m){ warning("Not enough memory for new markermatrix"); }
  for (int i=0; i<rows; i++) {
    m[i]= newMQMMarkerVector(cols);
  }
  return m;
}

void copyvector(vector vsource, vector vdestination, int dim) {
  for(int i=0; i<dim; i++) {
    vdestination[i]= vsource[i];
  }
}

#ifndef ENABLE_C99_MACROS
  void fatal(const char* s, ...){
      Rprintf("%s", s);
      Rf_error("Fatal error\n");
  }
  void debug_trace(const char*, ...){}
  void info(const char* s, ...){
      Rprintf("INFO: %s\n", s);
  }
#endif
