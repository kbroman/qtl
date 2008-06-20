######################################################################
#
# countXO.R
#
# copyright (c) 2008, Karl W Broman
# last modified Jun, 2008
# first written Feb, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: countXO
#
######################################################################

######################################################################
#
# countXO: Count number of obligate crossovers for each individual
#          on individual chromosomes or overall           
#
# if bychr=TRUE, return matrix with no. obligate crossovers for each
#                individual on each chromosome
#        =FALSE, return vector with total no. crossovers across the
#                selected chromosomes     
######################################################################

countXO <-
function(cross, chr, bychr=FALSE)
{
  if(!any(class(cross) == "cross")) 
    stop("Input should have class \"cross\".")

  # pull out relevant chromosome
  if(!missing(chr)) cross <- subset(cross,chr=chr)
  chr.name <- names(cross$geno)

  # which type of cross is this?
  type <- class(cross)[1]
  if(type == "f2") {
    if(class(cross$geno[[1]]) == "A") # autosomal
      func <- "R_countXO_f2"
    else func <- "R_countXO_bc"        # X chromsome  
  }
  else if(type == "bc" || type=="riself" || type=="risib" || type=="dh") func <- "R_countXO_bc"
  else if(type == "4way") func <- "R_countXO_4way"
  else 
    stop("ripple not available for cross ", type)

  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  nxo <- matrix(0, ncol=n.chr, nrow=n.ind)
  id <- getid(cross)
  if(is.null(id)) id <- 1:n.ind
  dimnames(nxo) <- list(id, chr.name)

  for(i in 1:n.chr) {
    # data to be input
    genodat <- cross$geno[[i]]$data
    genodat[is.na(genodat)] <- 0
    n.mar <- ncol(genodat)
    
    if(n.mar > 1) {
      z <- .C(func,
              as.integer(n.ind),
              as.integer(n.mar),
              as.integer(genodat),
              oblxo=as.integer(rep(0,n.ind)),
              PACKAGE="qtl")

      nxo[,i] <- z$oblxo
    }
  }
  
  if(!bychr) nxo <- apply(nxo, 1, sum)
  
  nxo
}

# end of countXO.R
