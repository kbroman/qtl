/**********************************************************************
 *
 * mqmdatatypes.h
 *
 * Copyright (c) 1996-2009 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Pjotr Prins and Danny Arends
 * last modified December 2009
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

#ifdef __cplusplus
  extern "C" {
#endif

// Cross types
typedef unsigned int RqtlCrossType;
enum MQMCrossType { CUNKNOWN = 'U', CF2 = 'F', CBC = 'B', CRIL = 'R' };
enum MQMCofactorType { MNOCOF = '0', MCOF ='1', MSEX = '2', MQTL = '3' };
enum MQMRelMarkerPos { MLEFT = 'L', MRIGHT = 'R', MMIDDLE = 'M', MUNKNOWN = 'U', MUNLINKED = '-' };
enum MQMMarker { MAA = '0', MH = '1', MBB = '2', MNOTAA    = '3', MNOTBB = '4',  MMISSING = '9', MUNUSED = '-'};

const RqtlCrossType RC_F2  = 1;
const RqtlCrossType RC_BC  = 2;
const RqtlCrossType RC_RIL = 3;
const double RFUNKNOWN = 999.0;
const double TRAITUNKNOWN = 999.0;
const double POSITIONUNKNOWN = 999.0;

// Marker locations/relations

//FIXME : Enum of relmarkerposition type
//FIXME : Typedef relmarkerposition* relmarkerarray to replace cvector
//FIXME : duplicate allocation of new type
//const unsigned char MLEFT     = 'L';
//const unsigned char MRIGHT    = 'R';
//const unsigned char MMIDDLE   = 'M';
//const unsigned char MUNLINKED = 'U';
//const unsigned char MUNKNOWN  = 0;

// Marker genotypes (scored at marker)
/*
const unsigned char MAA       = '0';  // Homozygous AA
const unsigned char MH        = '1';  // Heterozygous AB 
const unsigned char MBB       = '2';  // Homozygous BB
const unsigned char MNOTAA    = '3';  // Not AA
const unsigned char MNOTBB    = '4';  // Not BB 
const unsigned char MMISSING  = '9';  // Unknown (marker genotype missing)
const unsigned char MUNUSED   = '-';  // Unused parameter
*/


/*------------------------------------------------------------------------
Datastructures for matrix and vector calculus
FIXME : CamelCase for TYPEDEFS classes / structs! So ivector should be IVector
------------------------------------------------------------------------ */
typedef double** matrix;
typedef double*  vector;
typedef char**   cmatrix;
typedef MQMMarker **MQMMarkerMatrix;
typedef MQMMarker *MQMMarkerVector;
typedef char*    cvector;
typedef MQMRelMarkerPos* relmarkerarray;
typedef int*  ivector;


MQMCrossType determine_MQMCross(const int Nmark, const int Nind, const int **Geno, const RqtlCrossType rqtlcrosstype);

void change_coding(int *Nmark,int *Nind,int **Geno,MQMMarkerMatrix markers, const MQMCrossType crosstype);


/*------------------------------------------------------------------------
Basic routines for matrix and vector calculus
------------------------------------------------------------------------ */
vector newvector(int dim);
ivector newivector(int dim);
cvector newcvector(int dim);
MQMMarkerVector newMQMMarkerVector(int dim);
relmarkerarray newRelMarkerPos(int dim);
matrix newmatrix(int rows, int cols);
void printmatrix(matrix m, int rows, int cols);
void printcmatrix(cmatrix m, int rows, int cols);
cmatrix newcmatrix(int rows, int cols);
MQMMarkerMatrix newMQMMarkerMatrix(int rows, int cols);
void freematrix(void **m, size_t rows);
void freevector(void *v);
void delmatrix(matrix m, size_t rows);
void delcmatrix(cmatrix m, size_t rows);
void delMQMMarkerMatrix(MQMMarkerMatrix m,size_t rows);
void copyvector(vector vsource, vector vdestination, int dim);
void *calloc_init(size_t num, size_t size);

#ifndef ENABLE_C99_MACROS
  void info(const char *,...);

  void verbose(const char*,...);

  void debug_trace(const char *, ...);
#endif

#ifdef __cplusplus
  }
#endif

/* end of mqmdata.h */
