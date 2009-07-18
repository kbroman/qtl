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
# You can run it with:
#
#   R --no-save --no-restore --no-readline --slave < ./tests/test_qtl.R 

######################################################################

library(qtl)

data(listeria)
if (nind(listeria)!=120) stop("Number of individuals incorrect")

mr = scanone(listeria,method='mr')
test = round(mr[15,]$lod*1000)
cat(mr[15,]$lod,test)
if (test != 966) stop("scanone_mr gives an incorrect result")

augmentedcross <- MQMaugment(listeria,neglect=1)
if (nind(augmentedcross)!=116) stop("Number of individuals incorrect")

result <- scanMQM(augmentedcross)
if (round(result[5,5]*1000) != 153) stop("MQM gives an unexpected result (1)")
if (round(max(result[,5])*1000) != 5467) stop("MQM gives an unexpected result (2)")

cat("test_qtl.R tests succesfully run!")
