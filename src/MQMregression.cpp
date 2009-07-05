/**********************************************************************
 *
 * MQMRegression.cpp
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
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
 * C external functions used by the MQM algorithm
 * Contains: regression
 *
 **********************************************************************/
#include "MQM.h"
#include <Rmath.h>


/* regression of trait on multiple cofactors  y=xb+e with weight w
*							(xtwx)b=(xtw)y
*							b=inv(xtwx)(xtw)y
*/

double regression(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y,
                  vector *weight, ivector ind, int Naug, double *variance,
                  vector Fy, char biasadj, char fitQTL, char dominance) {
  // Rprintf("regression IN\n");
  /*
  cofactor[j] at locus j:
  '0': no cofactor at locus j
  '1': cofactor at locus j
  '2': QTL at locus j, but QTL effect is not included in the model
  '3': QTL at locu j and QTL effect is included in the model
  */
//	for (int j=0; j<Naug; j++){
//	   Rprintf("J:%d, COF:%d, VAR:%f, WEIGHT:%f, Trait:%f, IND[j]:%d\n", j, cofactor[j], *variance, (*weight)[j], y[j], ind[j]);
  //  }

  matrix XtWX;
  cmatrix Xt;
  vector XtWY;
  int dimx=1, j, jj;
  for (int j=0; j<Nmark; j++)
    if (cofactor[j]=='1') dimx+= (dominance=='n' ? 1 : 2);  // per QTL only additivity !!
    else if (cofactor[j]=='2') {
      dimx+=1;
    }

  XtWX= newmatrix(dimx+2, dimx+2);
  Xt= newcmatrix(dimx+2, Naug);
  XtWY= newvector(dimx+2);
  dimx=1;
  for (j=0; j<Nmark; j++)
    if ((cofactor[j]=='1')||(cofactor[j]=='3')) dimx+= (dominance=='y' ? 2 : 1);
  cvector xtQTL; // '0'=mu; '1'=cofactor; '2'=QTL (additive); '3'= QTL (dominance);
  xtQTL= newcvector(dimx);
  int jx=0;
  for (int i=0; i<Naug; i++) Xt[jx][i]= '1';
  xtQTL[jx]= '0';

  for (j=0; j<Nmark; j++)
    if (cofactor[j]=='1') { // cofactor (not a QTL moving along the chromosome)
      jx++;
      xtQTL[jx]= '1';
      if (dominance=='y') {
        for (int i=0; i<Naug; i++)
          if      (marker[j][i]=='1') {
            Xt[jx][i]=48;  //ASCII code 47, 48 en 49 voor -1, 0, 1;
            Xt[jx+1][i]=49;
          } else if (marker[j][i]=='0') {
            Xt[jx][i]=47;  // '/' stands for -1
            Xt[jx+1][i]=48;
          } else                        {
            Xt[jx][i]=49;
            Xt[jx+1][i]=48;
          }
        jx++;
        xtQTL[jx]= '1';
      } else {
        for (int i=0; i<Naug; i++) {
          if      (marker[j][i]=='1') {
            Xt[jx][i]=48;  //ASCII code 47, 48 en 49 voor -1, 0, 1;
          } else if (marker[j][i]=='0') {
            Xt[jx][i]=47;  // '/' stands for -1
          } else                        {
            Xt[jx][i]=49;
          }
        }
      }
    } else if (cofactor[j]=='3') { // QTL
      jx++;
      xtQTL[jx]= '2';
      if (dominance=='y') {
        jx++;
        xtQTL[jx]= '3';
      }
    }

  //Rprintf("calculate xtwx and xtwy\n");
  /* calculate xtwx and xtwy */
  double xtwj, yi, wi, calc_i;
  for (j=0; j<dimx; j++) {
    XtWY[j]= 0.0;
    for (jj=0; jj<dimx; jj++) XtWX[j][jj]= 0.0;
  }


  if (fitQTL=='n')
    for (int i=0; i<Naug; i++) {
      yi= y[i];
      wi= (*weight)[i];
      //in the original version when we enable Dominance , we crash around here
      for (j=0; j<dimx; j++) {
        xtwj= ((double)Xt[j][i]-48.0)*wi;
        XtWY[j]+= xtwj*yi;
        for (jj=0; jj<=j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
      }
    }
  else // QTL is moving along the chromosomes
    for (int i=0; i<Naug; i++) {
      wi= (*weight)[i]+ (*weight)[i+Naug]+ (*weight)[i+2*Naug];
      yi= y[i];
      //Changed <= to < to prevent chrashes, this could make calculations a tad different then before
      for (j=0; j<dimx; j++)
        if (xtQTL[j]<='1') {
          xtwj= ((double)Xt[j][i]-48.0)*wi;
          XtWY[j]+= xtwj*yi;
          for (jj=0; jj<=j; jj++)
            if (xtQTL[jj]<='1') XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            else if (xtQTL[jj]=='2') // QTL: additive effect if QTL='0' or '2'
            {  // QTL=='0'
              XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*(*weight)[i]*(47.0-48.0);
              // QTL=='2'
              XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*(*weight)[i+2*Naug]*(49.0-48.0);
            } else // (xtQTL[jj]=='3')  QTL: dominance effect only if QTL='1'
            {  // QTL=='1'
              XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*(*weight)[i+Naug]*(49.0-48.0);
            }
        } else if (xtQTL[j]=='2') { // QTL: additive effect if QTL='0' or '2'
          xtwj= -1.0*(*weight)[i]; // QTL=='0'
          XtWY[j]+= xtwj*yi;
          for (jj=0; jj<j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
          XtWX[j][j]+= xtwj*-1.0;
          xtwj= 1.0*(*weight)[i+2*Naug]; // QTL=='2'
          XtWY[j]+= xtwj*yi;
          for (jj=0; jj<j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
          XtWX[j][j]+= xtwj*1.0;
        } else { // (xtQTL[j]=='3') QTL: dominance effect only if QTL='1'
          xtwj= 1.0*(*weight)[i+Naug]; // QTL=='1'
          XtWY[j]+= xtwj*yi;
          // j-1 is for additive effect, which is orthogonal to dominance effect
          for (jj=0; jj<j-1; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
          XtWX[j][j]+= xtwj*1.0;
        }
    }
  for (j=0; j<dimx; j++)
    for (jj=j+1; jj<dimx; jj++) XtWX[j][jj]= XtWX[jj][j];


  int d;
  ivector indx;
  indx= newivector(dimx);
  /* solve equations */
  // Rprintf("LUcmp equations\nPrintinf matrix XiWX\n");
  // printmatrix(XtWX, dimx, dimx);
//	 for (jj=0; jj<dimx; jj++){
//		Rprintf("%f", indx);
//	 }
//	 Rprintf("\n");
  ludcmp(XtWX, dimx, indx, &d);

//	 Rprintf("LUsolve equations\nPrintinf indX\n");
  //for (jj=0; jj<dimx; jj++){
  //	Rprintf("%f", indx);
  //}
  //Rprintf("\n");

  lusolve(XtWX, dimx, indx, XtWY);
  // luinvert(xtwx, inv, dimx, indx);
  //  Rprintf("Parameter Estimates\n");
  //for (jj=0; jj<dimx; jj++){
  //	Rprintf("%d %f\n", jj, XtWY[jj]);
  // }

  long double *indL;
  // long double indL[Nind];
  indL = (long double *)Calloc(Nind, long double);
  int newNaug;
  newNaug= (fitQTL=='n' ? Naug : 3*Naug);
  vector fit, resi;
  fit= newvector(newNaug);
  resi= newvector(newNaug);
  // cout << "Calculate residuals" << endl;
  if (*variance<0) {
    *variance= 0.0;
    if (fitQTL=='n')
      for (int i=0; i<Naug; i++) {
        fit[i]= 0.0;
        for (j=0; j<dimx; j++)
          fit[i]+=((double)Xt[j][i]-48.0)*XtWY[j];
        resi[i]= y[i]-fit[i];
        *variance += (*weight)[i]*pow(resi[i], 2.0);
      }
    else
      for (int i=0; i<Naug; i++) {
        fit[i]= 0.0;
        fit[i+Naug]= 0.0;
        fit[i+2*Naug]= 0.0;
        for (j=0; j<dimx; j++)
          if (xtQTL[j]<='1') {
            calc_i =((double)Xt[j][i]-48.0)*XtWY[j];
            fit[i]+= calc_i;
            fit[i+Naug]+= calc_i;
            fit[i+2*Naug]+= calc_i;
          } else if (xtQTL[j]=='2') {
            fit[i]+=-1.0*XtWY[j];
            fit[i+2*Naug]+=1.0*XtWY[j];
          } else
            fit[i+Naug]+=1.0*XtWY[j];
        resi[i]= y[i]-fit[i];
        resi[i+Naug]= y[i]-fit[i+Naug];
        resi[i+2*Naug]= y[i]-fit[i+2*Naug];
        *variance +=(*weight)[i]*pow(resi[i], 2.0);
        *variance +=(*weight)[i+Naug]*pow(resi[i+Naug], 2.0);
        *variance +=(*weight)[i+2*Naug]*pow(resi[i+2*Naug], 2.0);
      }
    *variance/= (biasadj=='n' ? Nind : Nind-dimx); // to compare results with Johan; variance/=Nind;
    if (fitQTL=='n')
      for (int i=0; i<Naug; i++) Fy[i]= Lnormal(resi[i], *variance);
    else
      for (int i=0; i<Naug; i++) {
        Fy[i]       = Lnormal(resi[i], *variance);
        Fy[i+Naug]  = Lnormal(resi[i+Naug], *variance);
        Fy[i+2*Naug]= Lnormal(resi[i+2*Naug], *variance);
      }
  } else {
    if (fitQTL=='n')
      for (int i=0; i<Naug; i++) {
        fit[i]= 0.0;
        for (j=0; j<dimx; j++)
          fit[i]+=((double)Xt[j][i]-48.0)*XtWY[j];
        resi[i]= y[i]-fit[i];
        Fy[i]  = Lnormal(resi[i], *variance); // ????
      }
    else
      for (int i=0; i<Naug; i++) {
        fit[i]= 0.0;
        fit[i+Naug]= 0.0;
        fit[i+2*Naug]= 0.0;
        for (j=0; j<dimx; j++)
          if (xtQTL[j]<='1') {
            calc_i =((double)Xt[j][i]-48.0)*XtWY[j];
            fit[i]+= calc_i;
            fit[i+Naug]+= calc_i;
            fit[i+2*Naug]+= calc_i;
          } else if (xtQTL[j]=='2') {
            fit[i]+=-1.0*XtWY[j];
            fit[i+2*Naug]+=1.0*XtWY[j];
          } else
            fit[i+Naug]+=1.0*XtWY[j];
        resi[i]= y[i]-fit[i];
        resi[i+Naug]= y[i]-fit[i+Naug];
        resi[i+2*Naug]= y[i]-fit[i+2*Naug];
        Fy[i]       = Lnormal(resi[i], *variance);
        Fy[i+Naug]  = Lnormal(resi[i+Naug], *variance);
        Fy[i+2*Naug]= Lnormal(resi[i+2*Naug], *variance);
      }
  }
  /* calculation of logL */
  // cout << "calculate logL" << endl;
  long double logL=0.0;
  for (int i=0; i<Nind; i++) {
    indL[i]= 0.0;
  }
  if (fitQTL=='n') {
    for (int i=0; i<Naug; i++) indL[ind[i]]+=(*weight)[i]*Fy[i];
  } else {
    for (int i=0; i<Naug; i++) {
      indL[ind[i]]+=(*weight)[i]*       Fy[i];
      indL[ind[i]]+=(*weight)[i+Naug]*  Fy[i+Naug];
      indL[ind[i]]+=(*weight)[i+2*Naug]*Fy[i+2*Naug];
    }

  }

  for (int i=0; i<Nind; i++) {
    //Sum up log likelyhoods for each individual

    logL+= log(indL[i]);
  }
  //RRprintf("LLhood: %f\n", logL);

  Free(indL);
  Free(indx);
  Free(xtQTL);
  delmatrix(XtWX);
  delcmatrix(Xt);
  Free(fit);
  Free(resi);
  Free(XtWY);
  return (double)logL;
}

/* LU decomposition (from Johan via Numerical Recipes in C) Given an n x n matrix a[1..n][1..n], this routine replaces it by the LU
  decomposition of a rowwise permutation of itself.   A and n are input.  a is output. indx[1..n] is an output vector which records the row
  permutation effected by the partial pivoting; d is output as +-1 depending on whether the number of row interchanges was even or odd,
  respectively. This routine is used in combination with lusolve to solve
  linear equations or to invert a matrix.
*/
void ludcmp(matrix m, int dim, ivector ndx, int *d) {
  int r, c, rowmax, i;
  double max, temp, sum;
  vector scale, swap;
  scale= newvector(dim);
  *d=1;
  for (r=0; r<dim; r++) {
    for (max=0.0, c=0; c<dim; c++) if ((temp=fabs(m[r][c])) > max) max=temp;
    if (max==0.0) {
      Rprintf("ERROR: Singular matrix\n");
    }
    scale[r]=1.0/max;
  }
  for (c=0; c<dim; c++) {
    for (r=0; r<c; r++) {
      for (sum=m[r][c], i=0; i<r; i++) sum-= m[r][i]*m[i][c];
      m[r][c]=sum;
    }
    for (max=0.0, rowmax=c, r=c; r<dim; r++) {
      for (sum=m[r][c], i=0; i<c; i++) sum-= m[r][i]*m[i][c];
      m[r][c]=sum;
      if ((temp=scale[r]*fabs(sum)) > max) {
        max=temp;
        rowmax=r;
      }
    }
    if (max==0.0) {
      Rprintf("ERROR: Singular matrix\n");
    }
    if (rowmax!=c) {
      swap=m[rowmax];
      m[rowmax]=m[c];
      m[c]=swap;
      scale[rowmax]=scale[c];
      (*d)= -(*d);
    }
    ndx[c]=rowmax;
    temp=1.0/m[c][c];
    for (r=c+1; r<dim; r++) m[r][c]*=temp;
  }
  Free(scale);
  //Free(swap);
}

/* Solve the set of n linear equations AX=B.
Here a[1..n][1..n] is input as the LU decomposition of A.
b[1..n] is input as the right hand side vector B, and returns
with the solution vector X.
a, n and indx are not modified by this routine and can be left
for successive calls with different right-hand sides b.
*/
void lusolve(matrix lu, int dim, ivector ndx, vector b) {
  int r, c;
  double sum;
  for (r=0; r<dim; r++) {
    sum=b[ndx[r]];
    b[ndx[r]]=b[r];
    for (c=0; c<r; c++) sum-= lu[r][c]*b[c];
    b[r]=sum;
  }
  for (r=dim-1; r>-1; r--) {
    sum=b[r];
    for (c=r+1; c<dim; c++) sum-= lu[r][c]*b[c];
    b[r]=sum/lu[r][r];
  }
}



double inverseF(int df1, int df2, double alfa, int verbose) {
  double prob=0.0, prob2=0.0, minF=0.0, maxF=100.0, halfway=50.0, absdiff=1.0;
  int count=0;
  while ((absdiff>0.001)&&(count<100)) {
    //Rprintf("INFO df1:%d df2:%d alpha:%f\n", df1, df2, alfa);
    count++;
    halfway= (maxF+minF)/2.0;
    //prob= betai(df2/2.0, df1/2.0, df2/(df2+df1*halfway));
    //USE R FUNCTIONALITY TO REPLACE OLD C ROUTINES
    prob = pbeta(df2/(df2+df1*halfway), df2/2.0, df1/2.0, 1, 0);
    //Rprintf("->(%f, %f, %f) %f %f\n", df2/(df2+df1*halfway), df2/2.0, df1/2.0, prob, prob2);
    if (prob<alfa) maxF= halfway;
    else minF= halfway;
    absdiff= fabs(prob-alfa);
  }
  if (verbose==1) {
    Rprintf("INFO: Prob=%f Alfa=%f\n", prob, alfa);
  }
  return halfway;
}

