convert.varscan.to.empirical.ps <- function(scan, null.scan.maxes) {

	lod.columns <- grep(pattern = 'lod', names(scan), value = TRUE)
	chr.types <- unique(scan$chrtype)
	scan.as.emp.ps <- scan %>%
		select(-matches('lod'))
# 		mutate(lod.full = NA, lod.mean = NA, lod.var = NA)

	for (chr.type in chr.types) {

		null.lods <- dplyr::filter(null.scan.maxes, chrtype == chr.type)
		obs.lods <- dplyr::filter(scan, chrtype == chr.type)

		for (lod.column in lod.columns) {

			evd <- fgev(null.lods[[paste0('max.', lod.column)]])

			emp.ps <- pgev(q = obs.lods[[lod.column]],
										 loc = fitted(evd)[1],
										 scale = fitted(evd)[2],
										 shape = fitted(evd)[3],
										 lower.tail = FALSE)

			scan.as.emp.ps[scan$chrtype == chr.type, lod.column] <- emp.ps
		}
	}

	class(scan.as.emp.ps) <- class(scan)
	mostattributes(scan.as.emp.ps) <- attributes(scan)
	attr(scan.as.emp.ps, 'names') <- c(names(scan)[-(5:7)],
																		 'emp.p.full',
																		 'emp.p.mean',
																		 'emp.p.var')
	attr(scan.as.emp.ps, 'units') <- 'emp.ps'

	return(scan.as.emp.ps)
}
