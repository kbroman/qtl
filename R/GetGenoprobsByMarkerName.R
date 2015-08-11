GetGenoprobsByMarkerName <- function(cross, marker.name) {

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
