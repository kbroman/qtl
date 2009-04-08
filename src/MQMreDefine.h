/**********************************************************************
 * 
 * reDefine.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Mrt, 2009
 * first written Mrt, 2009
 *
 * Defines if we build a Rpackage or a standalone app, and any other functions etc we want to have changed AFTER reading all the libraries & dependancies
 * Contains: 
 *
 **********************************************************************/
 
//#define ALONE

#ifdef ALONE
       #undef Rprintf
       #define Rprintf(args...) printf(args)
#endif
