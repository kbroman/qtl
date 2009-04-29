#####################################################################
#
# MQMutil.R
#
# copyright (c) 2009, Danny Arends
# last modified Apr, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: ourcat, ourstop, ourline
#
######################################################################

ourcat <- function(...,a=TRUE){
	if(a){
		cat(...)
	}else{
		return
	}
}

ourstop <- function(...){
	ourline()
	stop(...)
	ourline()
}

ourline <- function(){
	cat("------------------------------------------------------------------\n")		
}

# end of MQMutil.R

