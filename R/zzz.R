######################################################################
#
# zzz.R
#
# copyright (c) 2001, Karl W Broman
# written Feb, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/qtl package
#
# .First.lib is run when the package is loaded with library(qtl)
#
######################################################################

.First.lib <- function(lib, pkg) library.dynam("qtl", pkg, lib)

# end of zzz.R
