#####################################################################
#
# summary.scanonevar.R
#
# copyright (c) 2001-2014, Karl W Broman
# modified by Robert Corty in March 2015
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
# Part of the R/qtl package
# Contains: summary.scanonevar,
#
######################################################################

######################################################################
#
# summary.scanonevar:
# Summarize the key features of a scanonevar
#
######################################################################

summary.scanonevar <- function(s) {

	print('Effect estimates when no marker present')
	print(paste('mean effects:', round(s$null.effects[[1]], 3)))
	print(paste('variance effects:', round(s$null.effects[[2]], 3)))
}
