#'  @title Summary of Peaks in Scanonevar
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{summary.scanonevar} prints out the loci in a scanonevar object
#'    that exceed \code{thresh}.  It is an S3 generic for summary().  It handles scanonevar
#'    objects in both LOD units and empirical p value units.
#'    
#'  @param scanonevar the scanonevar object to be summarized
#'  @param thresh the threshold over which (for LODs) or under which (for emprirical p values)
#'    a locus will be printed.
#'  
#'  @return None.  Only prints results to screen.
#'  
#'  
summary.scanonevar <- function(scanonevar, digits = 3, thresh) {

	if (!any(class(scanonevar) == "scanonevar")) {
		stop("Input should have class \"scanonevar\".")
	}

  if (missing(thresh)) {
    if (units(scanonevar) == 'lods') { thresh <- 3 }
    if (units(scanonevar) == 'emp.ps') { thresh <- 0.05 }    
  }


	peaks <- get.peaks.from.scanonevar(scanonevar, thresh)
	
	if (units(scanonevar) == 'lods') {
	  message('Full Model Peaks:')
	  print(peaks %>% 
	          dplyr::filter(full.peak == TRUE, full.lod > thresh) %>%
	          select(-matches('chrtype|effect|baseline|peak')))
	  message('Mean Model Peaks:')
	  print(peaks %>% 
	          dplyr::filter(mean.peak == TRUE, mean.lod > thresh) %>%
	          select(-matches('chrtype|effect|baseline|peak')))
	  message('Var Model Peaks:')
	  print(peaks %>% 
	          dplyr::filter(var.peak == TRUE, var.lod > thresh) %>%
	          select(-matches('chrtype|effect|baseline|peak')))
	}

	if (units(scanonevar) == 'emp.ps') {
	  message('Full Model Peaks:')
	  print(peaks %>% 
	          dplyr::filter(full.peak == TRUE, emp.p.full.lod < thresh) %>%
	          select(-matches('chrtype|effect|baseline|peak')))
	  message('Mean Model Peaks:')
	  print(peaks %>% 
	          dplyr::filter(mean.peak == TRUE, emp.p.mean.lod < thresh) %>%
	          select(-matches('chrtype|effect|baseline|peak')))
	  message('Var Model Peaks:')
	  print(peaks %>% 
	          dplyr::filter(var.peak == TRUE, emp.p.var.lod < thresh) %>%
	          select(-matches('chrtype|effect|baseline|peak')))
	}	
}
