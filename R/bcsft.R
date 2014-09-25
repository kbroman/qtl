convert2bcsft <- function(cross, BC.gen = 0, F.gen = 0, estimate.map = TRUE,
                          error.prob=0.0001, map.function=c("haldane","kosambi","c-f","morgan"),
                          verbose=FALSE)
{
    cross.class <- class(cross)[1]

    if((cross.class %in% c("bc","f2"))) {
        class(cross)[1] <- "bcsft"
        ## If BC.gen = 0 and F.gen = 0, then set to BC1F0 (bc) or BC0F2 (f2).
        if(cross.class == "bc" & F.gen > 0) {
            stop("input cross has only 2 genotypes--cannot have F.gen > 0")
            if(BC.gen == 0)
                BC.gen <- 1
        }
        if(cross.class == "f2") {
            if(F.gen == 0) {
                if(BC.gen == 0)
                    F.gen <- 2
                else
                    stop("input cross has 3 genotypes--cannot have F.gen = 0")
            }
        }
        if(BC.gen < 0 | F.gen < 0)
            stop("BC.gen and F.gen cannot be negative")

        attr(cross, "scheme") <- c(BC.gen, F.gen)
        cross
    }
    else stop("cross object has to be of class bc or f2 to be converted to bcsft")

    # re-estimate map?
    if(estimate.map) {
        cat(" --Estimating genetic map\n")
        newmap <- est.map(cross, error.prob=error.prob, map.function=map.function, verbose=verbose)
        cross <- replace.map(cross, newmap)
    }

    cross
}

read.cross.bcsft <- function(..., BC.gen = 0, F.gen = 0, cross = NULL, force.bcsft = FALSE,
                             estimate.map=TRUE)
{
    ## Must specify s = BC.gen and t = F.gen.
    ## Later: Could import in clever way from qtlcart? See qtlcart_io.R and their software.

    ## Make sure we only estimate map once!
    if(is.null(cross)) # Estimate map at end of this routine (called read.cross.bcsft directly).
        cross <- read.cross(..., estimate.map = FALSE)
    else # Estimate map in parent read.cross() call (read.cross.bcsft is pass-through from read.cross).
        estimate.map <- FALSE

    force.bcsft <- force.bcsft | (BC.gen > 0 | F.gen > 0)
    if((class(cross)[1] %in% c("bc","f2")) & force.bcsft)
        cross <- convert2bcsft(cross, BC.gen, F.gen, estimate.map = estimate.map, ...)

    cross
}

######################################################################

sim.cross.bcsft <- function(map,model,n.ind,error.prob,missing.prob,
                            partial.missing.prob,keep.errorind,
                            m,p,map.function,
                            cross.scheme)
{
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    if(any(sapply(map,is.matrix)))
        stop("Map must not be sex-specific.")

    ## cross.scheme = c(s,t) for bcsft.
    if(missing(cross.scheme))
        stop("must specify cross.scheme for bcsft")
    if(length(cross.scheme) != 2)
        stop("cross.scheme for bcsft must have 2 values")
    cross.scheme <- round(cross.scheme)
    if(min(cross.scheme) < 0)
        stop("cross.scheme for bcsft must have 2 non-negative integers")
    n.eff <- 3 + (cross.scheme[2] > 0)

    ## chromosome types
    chr.type <- sapply(map,function(a)
                       if(is.null(class(a))) return("A")
                       else return(class(a)))

    n.chr <- length(map)
    if(is.null(model)) n.qtl <- 0
    else {
        if(!((!is.matrix(model) && length(model) == n.eff) ||
             (is.matrix(model) && ncol(model) == n.eff))) {
            stop(paste("Model must be a matrix with ", n.eff, " columns (chr, pos and effect",
                       ifelse(n.eff == 4, "s", ""), ").", sep = ""))
        }
        if(!is.matrix(model)) model <- rbind(model)
        n.qtl <- nrow(model)
        if(any(model[,1] < 0 | model[,1] > n.chr))
            stop("Chromosome indicators in model matrix out of range.")
        model[,2] <- model[,2]+1e-14 # so QTL not on top of marker
    }

    # if any QTLs, place qtls on map
    if(n.qtl > 0) {
        for(i in 1:n.qtl) {
            temp <- map[[model[i,1]]]
            if(model[i,2] < min(temp)) {
                temp <- c(model[i,2],temp)
                names(temp)[1] <- paste("QTL",i,sep="")
            }
            else if(model[i,2] > max(temp)) {
                temp <- c(temp,model[i,2])
                names(temp)[length(temp)] <- paste("QTL",i,sep="")
            }
            else {
                j <- max((seq(along=temp))[temp < model[i,2]])
                temp <- c(temp[1:j],model[i,2],temp[(j+1):length(temp)])
                names(temp)[j+1] <- paste("QTL",i,sep="")
            }
            map[[model[i,1]]] <- temp
        }
    }

    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    n.mar <- sapply(map,length)
    mar.names <- lapply(map,names)
    BC.gen <- cross.scheme[1]
    F.gen <- cross.scheme[2] - (BC.gen == 0)

    for(i in 1:n.chr) {

        # simulate genotype data
        bcallele1 <- sim.bcg(n.ind, map[[i]], m, p, map.function) - 1
        ## BCs: multiply independent instances of meiosis together.
        if(BC.gen > 0) {
            if(BC.gen > 1) for(j in seq(2, BC.gen))
                bcallele1 <- bcallele1 * (sim.bcg(n.ind, map[[i]], m, p, map.function) - 1)
        }

        if(F.gen == 0) ## BCs only.
            thedata <- bcallele1 + 1
        else {
            if(chr.type[i] == "X") {
                if(F.gen > 1) for(j in seq(F.gen))
                    bcallele1 <- bcallele1 * (sim.bcg(n.ind, map[[i]], m, p, map.function) - 1)

                thedata <- bcallele1 + 1
            }
            else { ## chr.type[i] != "X"
                if(BC.gen > 0) { ## Two unique alleles from BC(s).
                    bcallele2 <- bcallele1 * (sim.bcg(n.ind, map[[i]], m, p, map.function) - 1)
                    bcallele1 <- bcallele1 * (sim.bcg(n.ind, map[[i]], m, p, map.function) - 1)
                }
                else ## Starting from F(1) with two unique alleles.
                    bcallele2 <- sim.bcg(n.ind, map[[i]], m, p, map.function) - 1

                if(F.gen > 1) for(j in seq(2, F.gen)) {
                    ## need two instances.
                    allelemask1 <- sim.bcg(n.ind, map[[i]], m, p, map.function) - 1
                    allelemask2 <- sim.bcg(n.ind, map[[i]], m, p, map.function) - 1
                    bcallele1 <- bcallele1 * allelemask1 + bcallele2 * (1 - allelemask1)
                    bcallele2 <- bcallele2 * allelemask2 + bcallele1 * (1 - allelemask2)
                }
                thedata <- bcallele1 + bcallele2 + 1
            }
        }

        dimnames(thedata) <- list(NULL,mar.names[[i]])

        geno[[i]] <- list(data = thedata, map = map[[i]])
        class(geno[[i]]) <- chr.type[i]
        class(geno[[i]]$map) <- NULL

    } # end loop over chromosomes

    # simulate phenotypes
    pheno <- rnorm(n.ind,0,1)

    if(n.qtl > 0) {
        # find QTL positions in genotype data
        QTL.chr <- QTL.loc <- NULL
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0) {
                QTL.chr <- c(QTL.chr,rep(i,length(o)))
                QTL.loc <- c(QTL.loc,o)
            }
        }

        # incorporate QTL effects
        for(i in 1:n.qtl) {
            QTL.geno <- geno[[QTL.chr[i]]]$data[,QTL.loc[i]]
            pheno[QTL.geno==2] <- pheno[QTL.geno==2] + model[i,n.eff]
            if(n.eff == 4) {
                pheno[QTL.geno==1] <- pheno[QTL.geno==1] - model[i,3]
                pheno[QTL.geno==3] <- pheno[QTL.geno==3] + model[i,3]
            }
        }

    } # end simulate phenotype

    n.mar <- sapply(geno, function(a) length(a$map))

    # add errors
    if(error.prob > 0) {
        for(i in 1:n.chr) {
            if(chr.type[i]=="X") {
                a <- sample(0:1,n.mar[i]*n.ind,replace=TRUE,
                            prob=c(1-error.prob,error.prob))
                geno[[i]]$data[a == 1] <- 3 - geno[[i]]$data[a == 1]
            }
            else {
                a <- sample(0:2,n.mar[i]*n.ind,replace=TRUE,
                            prob=c(1-error.prob,error.prob/2,error.prob/2))
                if(any(a>0 & geno[[i]]$data==1))
                    geno[[i]]$data[a>0 & geno[[i]]$data==1] <-
                        (geno[[i]]$data+a)[a>0 & geno[[i]]$data==1]
                if(any(a>0 & geno[[i]]$data==2)) {
                    geno[[i]]$data[a>0 & geno[[i]]$data==2] <-
                        (geno[[i]]$data+a)[a>0 & geno[[i]]$data==2]
                    geno[[i]]$data[geno[[i]]$data>3] <- 1
                }
                if(any(a>0 & geno[[i]]$data==3))
                    geno[[i]]$data[a>0 & geno[[i]]$data==3] <-
                        (geno[[i]]$data-a)[a>0 & geno[[i]]$data==3]
            }

            if(keep.errorind) {
                errors <- matrix(0,n.ind,n.mar[i])
                errors[a>0] <- 1
                colnames(errors) <- colnames(geno[[i]]$data)
                geno[[i]]$errors <- errors
            }

        } # end loop over chromosomes
    } # end simulate genotyping errors

    ## add partial missing
    if(partial.missing.prob > 0) {
        for(i in 1:n.chr) {
            if(chr.type[i] != "X") {
                o <- sample(c(TRUE,FALSE),n.mar[i],replace=TRUE,
                            prob=c(partial.missing.prob,1-partial.missing.prob))
                if(any(o)) {
                    o2 <- grep("^QTL[0-9]+",mar.names[[i]])
                    if(length(o2)>0)
                        x <- geno[[i]]$data[,o2]
                    m <- (1:n.mar[i])[o]
                    for(j in m) {
                        if(runif(1) < 0.5)
                            geno[[i]]$data[geno[[i]]$data[,j]==1 | geno[[i]]$data[,j]==2,j] <- 4
                        else
                            geno[[i]]$data[geno[[i]]$data[,j]==3 | geno[[i]]$data[,j]==2,j] <- 5
                    }
                    if(length(o2)>0)
                        geno[[i]]$data[,o2] <- x
                }
            }

        } # end loop over chromosomes
    } # end simulate partially missing data

    # add missing
    if(missing.prob > 0) {
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0)
                x <- geno[[i]]$data[,o]
            geno[[i]]$data[sample(c(TRUE,FALSE),n.mar[i]*n.ind,replace=TRUE,
                                  prob=c(missing.prob,1-missing.prob))] <- NA
            if(length(o)>0)
                geno[[i]]$data[,o] <- x
        }
    }

    pheno <- data.frame(phenotype=pheno)

    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c("bcsft","cross")
    attr(cross, "scheme") <- cross.scheme

    cross
}

## End bcsft.R
