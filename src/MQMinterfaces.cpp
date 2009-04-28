/**********************************************************************
 * 
 * MQMinterfaces.cpp
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Apr, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License, as
 *     published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version. 
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the
 *     GNU General Public License for more details.
 * 
 *     A copy of the GNU General Public License is available at
 *     http://www.r-project.org/Licenses/
 *
 * C external functions used by the MQM algorithm
 * Contains: Several functions to expose underlying elements from the C-code to R
 *
 **********************************************************************/
 
extern "C"
{
#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include <math.h>
#include "MQMData.h"
#include "MQMscan.h"


void R_Lnorm(double *a,double *b){
	Rprintf("Lnormal with parameters: (%f,%f)\n",*a,*b);
	double *ans;
	*ans = Lnormal(*a,*b);
	Rprintf("Result: %f\n",*ans);
	b = ans;
}

}
