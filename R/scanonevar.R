# scanonevar
# single-QTL genome scan for QTL affecting variance
# with code from Lars Ronnegard

scanonevar <-
function(cross, pheno.col=1, mean_covar = NULL, var_covar = NULL,
         additive_alleles=TRUE, family = gaussian(), maxit = 25 ,
         fixed_gamma_disp = FALSE, quiet=TRUE)
{
    # check input
    crosstype <- class(cross)[1]
    if(!(crosstype %in% c("bc", "dh", "f2", "haploid", "risib", "riself")))
       stop('scanonevar not implemented for cross type "', crosstype, '"')
    chrtype <- sapply(cross$geno, class)
    if(any(chrtype=="X")) {
        warning("Analysis of X chromosome not implemented for scanonevar; omitted.")
        cross <- subset(cross, chr=(chrtype != "X"))
    }

    # grab phenotype
    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }
    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(any(is.na(num))) {
            if(sum(is.na(num)) > 1)
                stop("Couldn't identify phenotypes ", paste(paste("\"", pheno.col[is.na(num)], "\"", sep=""),
                                                            collapse=" "))
            else
                stop("Couldn't identify phenotype \"", pheno.col[is.na(num)], "\"")
        }
        pheno.col <- num
    }
    if(any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")
    pheno <- cross$pheno[,pheno.col]
    if(is.matrix(pheno) && ncol(pheno) > 1) {
        pheno <- pheno[,1]
        warning('scanonevar requires a single phenotype; all but "', phenames(cross)[pheno.col[1]], '" omitted.')
    }

    N <- length(pheno) # No. individuals
    n.chr <- nchr(cross) #No. chromosomes
    chr.names <- chrnames(cross)

    # need to run calc.genoprob?
    if(!("prob" %in% names(cross$geno[[1]]))) {
        warning("First running calc.genoprob")
        cross <- calc.genoprob(cross)
    }

    scan.logPm <- scan.logPd <- chr.names.out <- NULL

    # set up data and formulas
    X <- cbind(pheno=pheno, add=rep(0, length(pheno)))
    mean_formula <- var_formula <- "pheno ~ add"
    if(crosstype=="f2" && !additive_alleles) {
        X <- cbind(X, dom=rep(0, nrow(X)))
        mean_formula <- var_formula <- "pheno ~ add + dom"
    }
    if(!is.null(mean_covar)) {
        ncolX <- ncol(X)
        X <- cbind(X, mean_covar)
        meancovarnames <- paste0("meancov", 1:(ncol(X)-ncolX))
        colnames(X)[-(1:ncolX)] <- meancovarnames
        mean_formula <- paste(mean_formula, "+",
                              paste(meancovarnames, collapse="+"))
    }
    if(!is.null(var_covar)) {
        ncolX <- ncol(X)
        X <- cbind(X, var_covar)
        varcovarnames <- paste0("varcov", 1:(ncol(X)-ncolX))
        colnames(X)[-(1:ncolX)] <- varcovarnames
        var_formula <- paste(var_formula, "+",
                             paste(varcovarnames, collapse="+"))
    }
    X <- as.data.frame(X)
    mean_formula <- as.formula(mean_formula)
    var_formula <- as.formula(var_formula)

    result <- NULL
    for(j in seq(along=cross$geno)) { # loop over chromosomes
        if(!quiet) message("Chr ", chr.names[j])

        if (crosstype=="f2") {
            g11 <- cross$geno[[j]]$prob[,,1]
            g12 <- cross$geno[[j]]$prob[,,2]
            g13 <- cross$geno[[j]]$prob[,,3]
            a1  <- g11 + g12/2
            d1 <-  g12 - (g11+g13)/2
        }
        else {
            a1 <- cross$geno[[j]]$prob[,,1]
        }

        n.loci <- dim(a1)[2]

        logP.m <- logP.d <- numeric(n.loci)

        for(i in 1:n.loci) { # loop over positions within chromosome
            # fill in genotype probs for this locus
            X[,2] <- a1[,i]
            if(crosstype=="f2" && !additive_alleles)
                X[,3] <- d1[,i]

            d.fit<-try(dglm::dglm(formula=mean_formula, dformula=var_formula,
                                  data=X, family = family, control=dglm::dglm.control(maxit=maxit)), silent=TRUE)

            if ((class(d.fit)!="dglm")[1]) {
                warning("dglm did not converge on chr ", chr.names[j], " position ", i)
            }
            else {
                p.mean <- summary(d.fit)$coef[2,4]
                if (!fixed_gamma_disp)
                    p.disp<- summary(d.fit$dispersion)$coef[2,4]
                if (fixed_gamma_disp)
                    p.disp<-(summary(d.fit)$dispersion.summary$coeff[2,4])
                logP.m[i]<- -log10(p.mean)
                logP.d[i]<- -log10(p.disp)
                if (d.fit$iter==maxit) {
                    logP.d[i]=0
                    warning("dglm did not converge on chr", chr.names[j], " position ", i)
                }
            }
        }

        # set up the output
        map <- attr(cross$geno[[j]]$prob,"map")
        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",chr.names[j],".",w[o],sep="")
        thischr <- data.frame(chr=rep(chr.names[j], length(w)),
                              pos=map,
                              neglogP_mean=logP.m,
                              neglogP_disp=logP.d, stringsAsFactors=FALSE)
        rownames(thischr) <- w

        if(is.null(result)) result <- thischr
        else result <- rbind(result, thischr)
    }

    class(result) <- c("scanone", "data.frame")
    attr(result, "method") <- "scanonevar"

    result
}