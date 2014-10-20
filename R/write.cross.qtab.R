######################################################################
# write.cross.qtab.R
#
# copyright (c) 2012, Karl W Broman and Danny Arends
# last modified Jul, 2012
# first written Jul, 2012
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
# Contains: write.cross.qtab
#       rqtl.to.qtab.* where * = symbols, location, genotypes, phenotypes, founder
#       get.qtlHD.ID, get.indID.for.qtab, getgenonames.for.qtab
#       get.qtab.geno.symbols, get.phenotype.type
#
######################################################################


# write cross in a set of qtab-format files
write.cross.qtab <-
    function(cross, filestem="data", descr, verbose=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")

    # for now, omit X chromosome
    chrtype <- sapply(cross$geno, class)
    if(any(chrtype == "X")) {
        cross <- subset(cross, chr=names(chrtype)[chrtype != "X"])
        warning("Omitting X chromosome.")
    }

    if(verbose) cat("Writing symbols\n")
    rqtl.to.qtab.symbols(cross, paste(filestem,"_symbols.qtab",sep=""), descr=descr)

    if(verbose) cat("Writing founder file\n")
    rqtl.to.qtab.founder(cross, paste(filestem,"_founder.qtab",sep=""), descr=descr)

    if(verbose) cat("Writing genetic map\n")
    rqtl.to.qtab.location(cross, paste(filestem,"_location.qtab",sep=""), descr=descr)

    if(verbose) cat("Writing genotypes\n")
    rqtl.to.qtab.genotypes(cross, paste(filestem,"_genotypes.qtab",sep=""), descr=descr)

    if(verbose) cat("Writing phenotypes\n")
    rqtl.to.qtab.phenotypes(cross, paste(filestem,"_phenotypes.qtab",sep=""), descr=descr)
}


# version number for qtlHD
get.qtlHD.ID <-
    function(){
        VER <- "0.1"
        ID <- paste("qtlHD-in-", VER, sep="")
        ID
    }

# individual IDs for qtab files
get.indID.for.qtab <-
    function(cross)
{
    id <- getid(cross)
    if(is.null(id)) id <- 1:nind(cross)
    paste("ID_", id, sep="")
}

# genotype codes for qtab files
getgenonames.for.qtab <-
    function(cross)
{
    gnames <- getgenonames(class(cross)[1], "A", "full", getsex(cross), attributes(cross))
    if(class(cross)[1] == "f2") {
        gnames <- c(gnames, paste(gnames[1], "or", gnames[2], sep=""), paste(gnames[2], "or", gnames[3], sep=""))
    }

    c("-", gnames)
}

# qtab genotypes symbols
get.qtab.geno.symbols <-
    function(cross)
{
    crtype <- class(cross)[1]
    if(crtype == "bc") {
        return(c("None", "0,0", "0,1"))
    }
    if(crtype == "riself" || crtype == "risib") {
        return(c("None", "0,0", "1,1"))
    }
    if(crtype == "f2") {
        return(c("None", "0,0", "0,1 1,0", "1,1", "0,0 0,1 1,0", "0,1 1,0 1,1"))
    }
    stop("cross type \"", crtype, "\" not yet supported for qtab.")
}

# write qtab symbols
rqtl.to.qtab.symbols <-
    function(cross, filename="symbols.qtab",descr)
{
    if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")

    cat(file=filename, "# --- ",get.qtlHD.ID()," Symbol ",descr, "\n",sep="");
    cat(file=filename, "# --- Genotype Symbol begin\n", append=TRUE)

    gnames <- getgenonames.for.qtab(cross)
    symbols <- get.qtab.geno.symbols(cross)
    for(i in seq(along=gnames)) {
        cat(file=filename, gnames[i], " as ", symbols[i], "\n", sep="", append=TRUE)
    }

    cat(file=filename, "# --- Genotype Symbol end\n", append=TRUE)
}

# write qtab marker map
rqtl.to.qtab.location <-
    function(cross, filename="locations.qtab", descr)
{
    if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")

    cat(file=filename, "# --- ",get.qtlHD.ID()," Location ", descr, "\n", sep="")
    cat(file=filename, "# --- Data Location begin\n", append=TRUE)
    cat(file=filename, "#\tChr\tPos\n", append=TRUE)
    map <- pull.map(cross, as.table=TRUE)
    map <- cbind(rownames(map), map)
    write.table(map, file=filename, append=TRUE, quote=FALSE, sep="\t", na="-", row.names=FALSE, col.names=FALSE)
    cat(file=filename, "# --- Data Location end\n", append=TRUE)
}

rqtl.to.qtab.genotypes <-
    function(cross, filename="genotypes.qtab", descr)
{
    if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")

    cat(file=filename, "# --- ",get.qtlHD.ID()," Genotype ", descr, "\n", sep="")
    cat(file=filename, "# --- Data Genotype begin\n", append=TRUE)

    # pull out genotypes data; convert to strings
    genotypes <- pull.geno(cross)
    genotypes[is.na(genotypes)] <- 0
    gnames <- getgenonames.for.qtab(cross)
    gstr <- matrix(rep("", prod(dim(genotypes))), ncol=ncol(genotypes))
    for(i in seq(along=gnames))
        gstr[genotypes == (i-1)] <- gnames[i]

    # add column with individual IDs
    id <- get.indID.for.qtab(cross)
    gstr <- cbind(id, gstr)

    # marker names
    cat(file=filename, "#", paste(colnames(genotypes), collapse="\t"), "\n", sep="", append=TRUE)

    # genotypes
    write.table(gstr, file=filename, append=TRUE, quote=FALSE, sep="\t", na="-", row.names=FALSE, col.names=FALSE)

    cat(file=filename, "# --- Data Genotype end\n", append=TRUE)
}

rqtl.to.qtab.phenotypes <-
    function(cross, filename="phenotypes.qtab", descr)
{
    if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")

    cat(file=filename, "# --- ",get.qtlHD.ID()," Phenotype ", descr, "\n", sep="")
    cat(file=filename, "# --- Type Phenotype begin\n", append=TRUE)

    for(phename in colnames(cross$pheno)) {
        cat(file=filename, phename, "\t", get.phenotype.type(cross,phename), "\n", sep="", append=TRUE)
    }

    cat(file=filename, "# --- Type Phenotype end", "\n", sep="", append=TRUE)

    cat(file=filename, "# --- Data Phenotype begin", "\n", sep="", append=TRUE)

    cat(file=filename, "#", paste(colnames(cross$pheno), collapse="\t"), "\n", sep="", append=TRUE)

    # add column with individual IDs
    id <- get.indID.for.qtab(cross)
    phe <- cbind(id, cross$pheno)

    write.table(phe, file=filename, append=TRUE, quote=FALSE, sep="\t", na="-", row.names=FALSE, col.names=FALSE)

    cat(file=filename, "# --- Data Phenotype end", "\n", sep="", append=TRUE)
}

get.phenotype.type <-
    function(cross, phenotype)
{
    if(class(cross$pheno[,phenotype])=="numeric") return("Float")
    if(class(cross$pheno[,phenotype])=="character") return("Char")
    if(class(cross$pheno[,phenotype])=="factor") return("Char")
    return("Float")
}

rqtl.to.qtab.founder <-
    function(cross, filename="founder.qtab", descr)
{
    if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")

    cat(file=filename, "# --- ",get.qtlHD.ID()," Founder ", descr, "\n", sep="")
    cat(file=filename, "# --- Set Founder begin\n", append=TRUE)
    cat(file=filename, "Cross\t", toupper(class(cross)[1]), "\n", sep="", append=TRUE)
    cat(file=filename, "# --- Set Founder end\n", append=TRUE)
}

# end of write.cross.qtab.R
