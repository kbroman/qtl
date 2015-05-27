units <- function(scan) {

	if (!class(scan)[1] == 'scanonevar') {
		stop("This 'units' function only applies to scanonevar class.")
	}

	return(attr(scan, 'units'))
}
