
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "teststatic.h"

/* R wrapper */
void R_teststatic(int *start, int *num)
{
  teststatic(*start, *num);
}

/* outer function (along the lines of an HMM main function */
void teststatic(int start, int num)
{
  int i;

  for(i=0; i<num; i++) 
    teststatic_sub(start, i);
}

/* function using static variables (like out step or init function */
void teststatic_sub(int start, int val)
{
  static int didcalc=0, staticvec[2]; /* storage of these should be maintained until R halts */
  int val2, val3;

  if(!didcalc || start != staticvec[0]) { /* check to see whether you need to re-initialize */
    didcalc = 1;
    Rprintf(" --Doing function initialization\n");
    staticvec[0] = start;
    staticvec[1] = start*2;
  }

  val2 = staticvec[0]+val;
  val3 = staticvec[1]+val;
  Rprintf(" %d %d %d\n", val, val2, val3);
}
