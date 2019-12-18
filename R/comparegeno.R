######################################################################
# comparegeno
######################################################################
comparegeno <-
    function(cross, what=c("proportion","number", "both"))
{
    if(!inherits(cross, "cross"))
        stop("Input should have class \"cross\".")

    what <- match.arg(what)
    g <- pull.geno(cross)
    g[is.na(g)] <- 0
    n.ind <- nrow(g)
    n.mar <- ncol(g)
    z <- .C("R_comparegeno",
            as.integer(g),
            as.integer(n.ind),
            as.integer(n.mar),
            n.match=as.integer(rep(0,n.ind^2)),
            n.missing=as.integer(rep(0,n.ind^2)),
            PACKAGE="qtl")

    if(what=="number") {
        z <- matrix(z$n.match,n.ind,n.ind)
    }
    else {
        if(what=="proportion") {
            z <- matrix(z$n.match/(n.mar-z$n.missing),n.ind,n.ind)
            diag(z) <- NA
        }
        else {
            prop <- matrix(z$n.match/(n.mar-z$n.missing),n.ind,n.ind)
            z <- matrix(z$n.match,n.ind,n.ind)
            z[lower.tri(z)] <- prop[lower.tri(z)]
        }
    }

    id <- getid(cross)
    if(is.null(id)) id <- as.character(seq_len(nind(cross)))

    dimnames(z) <- list(id, id)
    class(z) <- c("comparegeno", class(z))
    attr(z, "what") <- what

    z
}

# report which pairs of individuals have nearly-matching genotypes
summary.comparegeno <-
    function(object, thresh=0.9, ...)
{
    what <- attr(object, "what")
    if(is.null(what)) what <- "proportion"

    if(what=="number" && thresh < 1) thresh <- thresh*max(diag(object), na.rm=TRUE)

    wh <- which(!is.na(object) & object >= thresh & row(object) > col(object), arr=TRUE)
    if(length(wh)==0) {
        result <- data.frame(ind1=character(0),
                             ind2=character(0),
                             prop_match=numeric(0),
                             stringsAsFactors=FALSE)
        class(result) <- c("summary.comparegeno", "data.frame")
        return(result)
    }

    id <- rownames(object)
    if(is.null(id)) id <- as.character(seq_len(nrow(object)))


    # create results object
    result <- data.frame(ind1=id[wh[,2]],
                         ind2=id[wh[,1]],
                         prop_match=rep(0, nrow(wh)),
                         stringsAsFactors=FALSE)

    # fill in values with lower triangle
    for(i in seq_len(nrow(wh))) {
        result[i,3] <- object[wh[i,1], wh[i,2]]
    }


    # if both count and proportion is provided, determine number of markers
    if(what=="both") {
        result <- cbind(result,
                        n_markers=rep(0, nrow(wh)))

        for(i in seq_len(nrow(wh))) {
            result[i,4] <- round(object[wh[i,2], wh[i,1]]/object[wh[i,1], wh[i,2]])
        }
    }

    if(what=="number") { # change column name
        colnames(result)[3] <- "number_match"
    }

    # sort by decreasing percent matching
    result <- result[order(result[,3], decreasing=TRUE),,drop=FALSE]
    rownames(result) <- 1:nrow(result)

    class(result) <- c("summary.comparegeno", "data.frame")

    result

}


print.summary.comparegeno <-
    function(x, ...)
{
    if(nrow(x)==0) {
        cat("No pairs above threshold.\n")
    } else {
        print.data.frame(x, digits=3)
    }
}


plot.comparegeno <-
    function(x, breaks=NULL, main="", xlab="Proportion matching genotypes", ...)
{
    vals <- x[lower.tri(x)]

    if(is.null(breaks)) breaks <- 2*sqrt(length(vals))

    if(attr(x, "what")=="number" && missing(xlab)) {
        xlab <- "Number matching genotypes"
    }

    hist(vals, breaks=breaks, main=main, xlab=xlab, ...)
    rug(vals)
    invisible()
}
