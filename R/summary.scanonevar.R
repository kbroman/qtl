#####################################################################
#
# summary.scanonevar.R
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
# Contains: summary.scanonevar,
#
######################################################################

######################################################################
#
# summary.scanonevar:
# Summarize the key features of a scanonevar
#
######################################################################

summary.scanonevar <- function(vs, digits = 3, thresh) {

	if (!any(class(vs) == "scanonevar")) {
		stop("Input should have class \"scanonevar\".")
	}

  if (missing(thresh)) {
    if (units(vs) == 'lods') { thresh <- 3 }
    if (units(vs) == 'emp.ps') { thresh <- 0.05 }    
  }
  
	# Print Effects Fitted by Null Model
	null.effects <- attr(vs, 'null.effects')
	message('Effect estimates when no marker present')
	message('mean covariate effects:')
	print(round(null.effects[[1]], digits))
	message('variance covariate effects:')
	print(round(null.effects[[2]], digits))

	peaks <- GetPeaksFromVarscan(vs, thresh)
	
	
	message('Full Model Peaks:')
	print(peaks %>% 
	        dplyr::filter(full.peak == TRUE, emp.p.full < thresh) %>%
	        select(-matches('chrtype|effect|baseline|peak')))
	message('Mean Model Peaks:')
	print(peaks %>% 
	        dplyr::filter(mean.peak == TRUE, emp.p.mean < thresh) %>%
	        select(-matches('chrtype|effect|baseline|peak')))
	message('Var Model Peaks:')
	print(peaks %>% 
	        dplyr::filter(var.peak == TRUE, emp.p.var < thresh) %>%
	        select(-matches('chrtype|effect|baseline|peak')))
}
