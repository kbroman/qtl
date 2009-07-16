/**********************************************************************
 *
 * MQMextra.cpp
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Apr, 2009
 * first written Apr, 2009
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
 * Contains: Extra functions for accessing internals
 *
 **********************************************************************/

#include "mqm.h"

void R_Lnorm(double *a,double *b) {
  Rprintf("Lnormal with parameters: (%f,%f)\n",*a,*b);
  double *ans;
  *ans = Lnormal(*a,*b);
  Rprintf("Result: %f\n",*ans);
  b = ans;
}

