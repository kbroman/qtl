######################################################################
# write.cross.mq.R
#
# copyright (c) 2014, 2015, INRA (author: Timothee Flutre)
# last modified March, 2015
# first written May, 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License,
# version 3, as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  See the GNU
# General Public License, version 3, for more details.
#
# A copy of the GNU General Public License, version 3, is available
# at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/qtl package
# Contains: write.cross.mq, write.cross.mq.loc, write.cross.mq.map,
# and write.cross.mq.qua
#
######################################################################

######################################################################
#
# write.cross.mq: write data from an experimental cross in MapQTL (and
# JoinMap) format.
#
# Three files are written: a "loc" file containing the genotype data,
# a "map" file containing the linkage group assignments and map
# positions, and a "qua" file containing the phenotypes.
#
# File formats are described in the MapQTL manual available online at
# http://www.kyazma.nl/docs/MQ6Manual.pdf
# Only 4-way crosses are supported ("CP" type in MapQTL/JoinMap).
#
######################################################################

write.cross.mq <-
    function(cross, filestem, digits)
{
    if(class(cross)[1] != "4way"){
        msg <- paste("population type", class(cross)[1],
                     "is not supported for writing in MapQTL format (yet)")
        stop(msg, call.=FALSE)
    }
    write.cross.mq.loc(cross, paste0(filestem, ".loc"))
    write.cross.mq.map(cross, paste0(filestem, ".map"), digits)
    write.cross.mq.qua(cross, paste0(filestem, ".qua"), digits)
}

write.cross.mq.loc <-
    function(cross, locfile)
{
    sink(locfile)

    cat(paste0("name = cross", "\n"))
    if(class(cross)[1] == "4way"){
        cat(paste0("popt = ", "CP", "\n"))
    } else{
        sink()
        msg <- paste("population type", class(cross)[1],
                     "is not supported (yet)")
        stop(msg, call.=FALSE)
    }
    cat(paste0("nloc = ", totmar(cross), "\n"))
    cat(paste0("nind = ", nind(cross), "\n"))
    cat("\n")

    for(chr in names(cross$geno)){
        for(marker in colnames(cross$geno[[chr]]$data)){
            cat(marker)
            tmp <- sort(names(table(cross$geno[[chr]]$data[,marker])))
            if(length(tmp) == 4 && all(tmp == c("1","2","3","4"))){
                cat("\t<abxcd>") # arbitrary choice
                cat("\t{00}") # arbitrary choice, could also be {01}, {10}, {11}
                for(i in 1:nind(cross)){
                    if(is.na(cross$geno[[chr]]$data[i, marker])){
                        cat("\t--")
                    } else if(cross$geno[[chr]]$data[i, marker] == 1){
                        cat("\tac")
                    } else if(cross$geno[[chr]]$data[i, marker] == 3){
                        cat("\tad")
                    } else if(cross$geno[[chr]]$data[i, marker] == 2){
                        cat("\tbc")
                    } else if(cross$geno[[chr]]$data[i, marker] == 4){
                        cat("\tbd")
                    }
                }
                cat("\n")
            } else if((length(tmp) == 3 && all(tmp == c("1","10","4"))) ||
                      (length(tmp) == 5 && all(tmp == c("1","10","11","14","4")))){
                cat("\t<hkxhk>")
                cat("\t{00}") # arbitrary choice, could also be {11}
                for(i in 1:nind(cross)){
                    if(is.na(cross$geno[[chr]]$data[i, marker])){
                        cat("\t--")
                    } else if(cross$geno[[chr]]$data[i, marker] == 1){
                        cat("\thh")
                    } else if(cross$geno[[chr]]$data[i, marker] == 10){
                        cat("\thk")
                    } else if(cross$geno[[chr]]$data[i, marker] == 4){
                        cat("\tkk")
                    } else if(cross$geno[[chr]]$data[i, marker] == 14){
                        cat("\th-")
                    } else if(cross$geno[[chr]]$data[i, marker] == 11){
                        cat("\tk-")
                    }
                }
                cat("\n")
            } else if((length(tmp) == 3 && all(tmp == c("2","3","9"))) ||
                      (length(tmp) == 5 && all(tmp == c("12","13","2","3","9")))){
                cat("\t<hkxhk>")
                cat("\t{01}") # arbitrary choice, could also be {10}
                for(i in 1:nind(cross)){
                    if(is.na(cross$geno[[chr]]$data[i, marker])){
                        cat("\t--")
                    } else if(cross$geno[[chr]]$data[i, marker] == 9){
                        cat("\thk")
                    } else if(cross$geno[[chr]]$data[i, marker] == 3){
                        cat("\thh")
                    } else if(cross$geno[[chr]]$data[i, marker] == 2){
                        cat("\tkk")
                    } else if(cross$geno[[chr]]$data[i, marker] == 12){
                        cat("\th-")
                    } else if(cross$geno[[chr]]$data[i, marker] == 13){
                        cat("\tk-")
                    }
                }
                cat("\n")
            } else if(length(tmp) == 2 && all(tmp == c("5","6"))){
                cat("\t<lmxll>")
                cat("\t{0-}") # arbitrary choice, could also be {1-}
                for(i in 1:nind(cross)){
                    if(is.na(cross$geno[[chr]]$data[i, marker])){
                        cat("\t--")
                    } else if(cross$geno[[chr]]$data[i, marker] == 5){
                        cat("\tll")
                    } else if(cross$geno[[chr]]$data[i, marker] == 6){
                        cat("\tlm")
                    }
                }
                cat("\n")
            } else if(length(tmp) == 2 && all(tmp == c("7","8"))){
                cat("\t<nnxnp>")
                cat("\t{-0}") # arbitrary choice, could also be {-1}
                for(i in 1:nind(cross)){
                    if(is.na(cross$geno[[chr]]$data[i, marker])){
                        cat("\t--")
                    } else if(cross$geno[[chr]]$data[i, marker] == 7){
                        cat("\tnn")
                    } else if(cross$geno[[chr]]$data[i, marker] == 8){
                        cat("\tnp")
                    }
                }
                cat("\n")
            } else{
                sink()
                msg <- paste("unrecognized segregation type at marker",
                             marker, "on chromosome", chr)
                stop(msg, call.=FALSE)
            }
        }
    }

    sink()
}

write.cross.mq.map <-
    function(cross, mapfile, digits=NULL)
{
    if(is.matrix(cross$geno[[1]]$map)){
        mapfile.female <- sub(pattern="\\.map",
                              replacement="_female.map",
                              x=mapfile)
        sink(mapfile.female)
        for(chr in names(cross$geno)){
            cat(paste0("group ", chr, "\n"))
            for(m in 1:ncol(cross$geno[[chr]]$map)){
                cat(paste0(colnames(cross$geno[[chr]]$map)[m],
                           "\t",
                           ifelse(is.null(digits),
                                  cross$geno[[chr]]$map[1,m],
                                  round(cross$geno[[chr]]$map[1,m], digits)),
                           "\n"))
            }
            cat("\n")
        }
        sink()
        mapfile.male <- sub(pattern="\\.map",
                            replacement="_male.map",
                            x=mapfile)
        sink(mapfile.male)
        for(chr in names(cross$geno)){
            cat(paste0("group ", chr, "\n"))
            for(m in 1:ncol(cross$geno[[chr]]$map)){
                cat(paste0(colnames(cross$geno[[chr]]$map)[m],
                           "\t",
                           ifelse(is.null(digits),
                                  cross$geno[[chr]]$map[2,m],
                                  round(cross$geno[[chr]]$map[2,m], digits)),
                           "\n"))
            }
            cat("\n")
        }
        sink()
    } else{
        sink(mapfile)
        for(chr in names(cross$geno)){
            cat(paste0("group ", chr, "\n"))
            for(m in 1:length(cross$geno[[chr]]$map)){
                cat(paste0(names(cross$geno[[chr]]$map)[m],
                           "\t",
                           ifelse(is.null(digits),
                                  cross$geno[[chr]]$map[m],
                                  round(cross$geno[[chr]]$map[m], digits)),
                           "\n"))
            }
            cat("\n")
        }
        sink()
    }
}

write.cross.mq.qua <-
    function(cross, quafile, digits=NULL)
{
    sink(quafile)

    cat(paste0("ntrt = ", nphe(cross), "\n"))
    cat(paste0("nind = ", nind(cross), "\n"))
    cat(paste0("miss = ", NA, "\n"))
    cat("\n")

    colnames(cross$pheno) <- gsub(pattern=" ", replacement="_",
                                  x=colnames(cross$pheno))
    if(any(nchar(colnames(cross$pheno)) > 20)){
        sink()
        msg <- paste("phenotype name",
                     colnames(cross$pheno)[which(nchar(colnames(cross$pheno)) > 20)[1]],
                     "is longer than 20 characters")
        stop(msg, call.=FALSE)
    }
    cat(colnames(cross$pheno)[1])
    if(ncol(cross$pheno) > 1)
        for(j in 2:ncol(cross$pheno))
            cat(paste0("\t", colnames(cross$pheno)[j]))
    cat("\n")
    cat("\n")

    if(is.null(digits)){
        write.table(x=cross$pheno, file=quafile, append=TRUE,
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    } else
        write.table(x=format(x=cross$pheno, digits=digits, trim=TRUE),
                    file=quafile, append=TRUE, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)

    sink()
}

# end of write.cross.mq.R
