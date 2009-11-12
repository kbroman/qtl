#####################################################################
#
# mqmutil.R
#
# copyright (c) 2009, Danny Arends
# last modified Nov, 2009
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

testnormal <- function(cross,pheno.col=1){
	if(any(rownames(installed.packages())=="nortest")){
		library(nortest)
		if(pearson.test(cross$pheno[[pheno.col]])$p.value < 0.05){
			cat("Trait distribution not normal\n")
		}else{
			cat("Trait distribution normal\n")
		}
	}else{
		cat("Please install package: nortest\n")
	}
}

# end of mqmutil.R

