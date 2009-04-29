/**********************************************************************
 * 
 * reDefine.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Mrt, 2009
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
 * Defines if we build a Rpackage or a standalone app, and any other functions etc we want to have changed AFTER reading all the libraries & dependancies
 * Contains: 
 *
 **********************************************************************/
 
//#define ALONE

#ifdef ALONE
       #undef Rprintf
       #define Rprintf(args...) printf(args)
#endif
