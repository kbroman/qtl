# function under development
# compute standardized information content across the genome rather efficiently

information.content <- function(cross) {
  
  if (!any(class(cross) == 'cross'))
    stop('Must pass in a cross object.')
  
  if (!any(class(cross) == 'f2'))
    stop('InformationContent() currently only implemented for F2 cross type')
  
  
  all.info <- list()
  for (chr.name in names(cross$geno)) {
    
    this.chr <- cross$geno[[chr.name]]
    this.chr.probs <- this.chr[['prob']]
    chr.name <- str_pad(string = chr.name, width = 2, side = 'left', pad = '0')
    
    if (attr(this.chr, 'class') == 'A') {
      
      if (dim(this.chr.probs)[3] != 3)
        stop(paste('Chromosome', chr.name, 'should have 3 genoprobs but doesnt'))
      
      InfoFunc <- function(m) {
        indiv.infos <- apply(X = m, MARGIN = 1, FUN = function(v) {
          return(sum(v * log(v/c(0.25, 0.5, 0.25))))
        })
        return(sum(indiv.infos))
      }
        
      this.chr.info <- data_frame(chr = chr.name,
                                  marker.name = attr(this.chr.probs, 'dimnames')[[2]],
                                  marker.loc = attr(this.chr.probs, 'map'),
                                  info = apply(X = this.chr.probs, MARGIN = c(2), FUN = InfoFunc))
      
    } else if (attr(this.chr, 'class') == 'X') {
      
      if (dim(this.chr.probs)[3] != 2)
        stop(paste('Chromosome', chr.name, 'should have 2 genoprobs but doesnt'))
      
      InfoFunc <- function(m) {
        indiv.infos <- apply(X = m, MARGIN = 1, FUN = function(v) {
          return(sum(v * log(v/c(0.5, 0.5))))
        })
        return(sum(indiv.infos))
      }
      
      this.chr.info <- data_frame(chr = chr.name,
                                  marker.name = attr(this.chr.probs, 'dimnames')[[2]],
                                  marker.loc = attr(this.chr.probs, 'map'),
                                  info = apply(X = this.chr.probs, MARGIN = c(2), FUN = InfoFunc))
      
    } else {
      stop(paste('Chromosome', chr.name, 'has unacceptable "class" attribute. (Must be "A" or "X".)'))
    }
      
    all.info[[chr.name]] <- this.chr.info
  }
  
  return(rbind_all(all.info))
}