#'  @title Get Genotype Probabilities From a Cross Object By Marker (or Psuedomarker) Name
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{get.genoprobs.by.marker.name} is an accessor function that returns the 
#'    genotype probabilities of each individual in the cross at the given marker or pseudomarker.
#'  
#'  @param cross The cross from which the genetic information will be extracted.
#'  @param marker.name The name of the marker where we want to know each individuals most likely genotype.
#'    
#'  @return Probability of each genotype at the given locus for each individual.
#'
#'  
get.genoprobs.by.marker.name <- function(cross, marker.name) {

	idx.of.chr.with.marker <- which(sapply(X = cross$geno,
																				 FUN = function(x, mn) { mn %in% colnames(x$data)},
																				 marker.name))

	name.of.chr.with.marker <- names(cross$geno)[idx.of.chr.with.marker]

	if (length(name.of.chr.with.marker) == 0) {	stop('Marker not found.')	}
	if (length(name.of.chr.with.marker) >= 2) {	stop('Marker found 2+ times.')	}

	genoprobs <- cross$geno[[name.of.chr.with.marker]]$prob[,marker.name,]
	colnames(genoprobs) <- paste0(marker.name, '_', colnames(genoprobs))

	return(genoprobs)
}
