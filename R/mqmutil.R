#####################################################################
#
# MQMutil.R
#
# copyright (c) 2009, Danny Arends
# last modified Apr, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: ourstop, ourline
#
######################################################################

ourstop <- function(...){
	ourline()
	stop(...)
#	ourline() [commented out, as this would never be called anyway]
}

ourline <- function(){
	cat("------------------------------------------------------------------\n")		
}

# end of MQMutil.R

