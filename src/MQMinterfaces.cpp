/**********************************************************************
 * 
 * MQMinterfaces.cpp
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Apr, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
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
