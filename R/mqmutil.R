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
	stop(...)
}

ourline <- function(){
	cat("------------------------------------------------------------------\n")		
}

mqmtestnormal <- function(cross, pheno.col=1){
	if(pheno.col <0 || pheno.col > nphe(cross)){
		stop("No such phenotype (pheno.col = ",pheno.col,")")
	}
	if(!is.numeric(cross$pheno[[pheno.col]])){
		stop("Please supply a numeric trait (pheno.col = ",pheno.col," is not numeric)")
	}
	if(any(rownames(installed.packages())=="nortest")){
		library(nortest)
		if(pearson.test(cross$pheno[[pheno.col]])$p.value < 0.05){
			cat("Trait distribution not normal\n")
		}else{
			cat("Trait distribution normal\n")
		}
	}else{
		cat("Please install package: nortest to enable testing of normality\n")
	}
}

# end of mqmutil.R

