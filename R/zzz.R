######################################################################
#
# zzz.R
#
# copyright (c) 2001, Karl W Broman
# written Feb, 2001
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
#
# .First.lib is run when the package is loaded with library(qtl)
#
######################################################################

.First.lib <- function(lib, pkg) library.dynam("qtl", pkg, lib)

# end of zzz.R
