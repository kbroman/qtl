######################################################################
#
# Regression test
#
# copyright (c) 2009 Pjotr Prins
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

script='mqm_listeria1'

library(qtl)

data(listeria)

augmentedcross <- mqmaugment(listeria, minprob=1.0)
result <- mqmscan(augmentedcross)
sink(paste('regression/',script,'.rnew',sep=''))
result
sink()

cat(script,'successful')

