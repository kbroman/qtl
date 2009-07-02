######################################################################
#
# test_qtl.R
#
# copyright (c) 2009, Karl W Broman, Pjotr Prins
# first written July 2009
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Some basic regression/integration testing for some of the QTL mapping routines
#
######################################################################

library(qtl)

data(listeria)
if (nind(listeria)!=121) stop("Number of individuals incorrect")

mr = scanone(listeria,method='mr')
if (mr[15,]$lod != 0.9663805) stop("scanone_mr gives an incorrect result")

