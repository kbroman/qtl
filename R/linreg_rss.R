######################################################################
# linreg_rss.R
#
# copyright (c) 2012, Karl W Broman
# last modified Apr, 2012
# first written Apr, 2012
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
# calculate residual sum of squares from linear regression
# [really just to test my linreg_rss code]
######################################################################

linreg.rss <-
function(x, y, tol=1e-8)
{
  if(!is.matrix(y)) y <- as.matrix(y)

  if(nrow(x) != nrow(y))
    stop("nrow(x) != nrow(y): ", nrow(x), " != ", nrow(y))

  .C("R_linreg_rss",
     as.integer(nrow(x)),
     as.integer(ncol(x)),
     as.double(x),
     as.integer(ncol(y)),
     as.double(y),
     rss=as.double(rep(0,ncol(y))),
     as.double(tol),
     PACKAGE="qtl")$rss
}

# end of linreg_rss.R
