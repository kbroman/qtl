#####################################################################
#
# prepareMQM.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
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
# Contains: MQMLoadanalysis, MQMfind.marker
#
######################################################################

MQMLoadanalysis <- function(file="MQM_output.txt"){
	data <- read.table(file)
	class(data) <- c("scanone",class(data))
	plot(data)
	data
}


MQMfind.marker <- function(cross,mqmscan=NULL,perm=NULL,alpha=0.05,verbose=FALSE){
	
	chr <- summary(mqmscan,alpha=alpha,perms=perm,pvalues=FALSE)$'chr'
	pos <- summary(mqmscan,alpha=alpha,perms=perm,pvalues=FALSE)$'pos (Cm)'
	if(verbose) cat("INFO: Found",length(chr),"markers with alpha <",alpha,".\n")
	ret <- NULL
	for(i in 1:length(chr)){
		ret <- rbind(ret,cbind(find.marker(cross,chr=chr[i],pos=pos[i]),as.integer(chr[i]),as.double(pos[i])))
	}
	colnames(ret) <- c("marker","chr","pos (Cm)")
	ret
}



# end of prepareMQM.R
