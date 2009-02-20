######################################################################
#
# zzz.R
#
# copyright (c) 2001, Karl W Broman
# written Feb, 2001
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License, as
#     published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version. 
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the
#     GNU General Public License for more details.
# 
#     A copy of the GNU General Public License is available at
#     http://www.r-project.org/Licenses/
#
# Part of the R/qtl package
#
# .First.lib is run when the package is loaded with library(qtl)
#
######################################################################

.First.lib <- function(lib, pkg) library.dynam("qtl", pkg, lib)

# end of zzz.R
