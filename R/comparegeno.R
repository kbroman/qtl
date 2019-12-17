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

    class(z) <- c("comparegeno", class(z))

    z
}

# report which pairs of lines have nearly-matching genotypes
summary.comparegeno <-
    function(object, line_name="line", thresh=0.8)
{
    wh <- which(!is.na(cg) & cg > thresh & row(cg) < col(cg), arr=TRUE)
    if(length(wh)==0) return(NULL)

    if(is.null(line_name))
        id <- 1:nrow(cg)
    else
        id <- cross$pheno[,line_name]

    g <- pull.geno(cross)

    match <- not_mis <- rep(0, nrow(wh))
    for(i in 1:nrow(wh)) {
        match[i] <- sum(g[wh[i,1],] == g[wh[i,2],], na.rm=TRUE)
        not_mis[i] <- sum(!is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]))
    }

    result <- data.frame(line1=id[wh[,1]],
                         line2=id[wh[,2]],
                         percent_match=match/not_mis*100,
                         match=match,
                         total=not_mis)

    # sort by decreasing percent matching
    result <- result[order(result$percent_match, decreasing=TRUE),,drop=FALSE]
    rownames(result) <- 1:nrow(result)

    result
}
