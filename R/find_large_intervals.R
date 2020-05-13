# find large intervals in a map
#
#
# input: map = list of chromosomes that are vectors of marker positions
#              [can also be a cross object, in which case pull.map() is used]
#        min_length = minimum distance between markers to be flagged
#
# example use:
#   library(qtl)
#   data(hyper)
#   find_large_intervals(hyper, 20)

find_large_intervals <-
    function(map, min_length=35)
{
    output <- NULL

    # if input is a cross, pull out the map
    if("cross" %in% class(map)) {
        map <- pull.map(map)
    }

    # make sure there are names
    if(is.null(names(map))) {
        names(map) <- seq_along(map)
    }

    for(i in names(map)) {

        d <- diff(map[[i]])
        big <- which(d > min_length)
        n_big <- length(big)
        mar <- names(map[[i]])

        if(n_big > 0) {
            this <- data.frame(chr=rep(i, n_big),
                               left=mar[big],
                               right=mar[big+1],
                               interval_size=d[big],
                               stringsAsFactors=FALSE)

            if(is.null(output)) output <- this
            else output <- rbind(output, this)
        }

    }

    if(!is.null(output)) rownames(output) <- 1:nrow(output)
    output
}
