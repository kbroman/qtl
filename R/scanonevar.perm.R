#'  @title Conduct Scanonevars on Permuted Genotype Data
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{scanonevar.perm} conducts many \code{scanonevar}s on permuted data
#'    and returns the maximum observed LOD score for each chromosome type in each scan.
#'    The results should be put into \code{convert.scanonevar.to.empirical.ps} with a scan in LODs
#'    to convert that scan to empirical p-values.  It's important that all the parameters used in
#'    \code{scanonevar.perm} are the same as the parameters that were used in the \code{scanonevar}
#'    that they will be used to convert to empirical ps.
#'    
#'  @param n.perms the number of permutations to conduct
#'  @inheritParams scanonevar
#'  
#'  @return Returns a tbl_df of maximum LOD score observed in each genome scan for each chromosome type.
#'    
#'  @seealso  \code{\link{scanonevar}}, \code{\link{convert.scanonevar.to.empirical.ps}}
#'  
#'  @details It is recommended to use approximately 1000 permuted scans to produce highly-replicable,
#'    publication-quality empirial p-values.  For this purpose, users are recommended to dispatch this
#'    function to many computers in parallel, carefully setting the seed on each computer to insure
#'    pseudo-randomness.
#'  
#'  
scanonevar.perm <- function(cross,
                            mean.formula,
                            var.formula,
                            n.perms,
                            chrs = unique(names(cross$geno))) 
{
  
  validated.input <- validate.input.scanonevar(cross = cross, 
                                               mean.formula = mean.formula, 
                                               var.formula = var.formula,
                                               chrs = chrs)
  genoprobs <- validated.input$genoprobs
  mapping.df <- validated.input$mapping.df
  chr.by.marker <- validated.input$chr.by.marker
  pos.by.marker <- validated.input$pos.by.marker
  marker.names <- validated.input$marker.names
  mean.formula <- validated.input$mean.formula
  var.formula <- validated.input$var.formula
  
  all.perms <- NULL
  for (perm.idx in 1:n.perms) {
    perm.scan <- scan.via.dglm(mean.alt.formula = mean.formula,
                               var.alt.formula = var.formula,
                               genoprobs = genoprobs,
                               mapping.df = mapping.df,
                               chr.by.marker = chr.by.marker,
                               pos.by.marker = pos.by.marker,
                               marker.names = marker.names,
                               perm = sample(nrow(genoprobs)))
    
    this.perm <- perm.scan %>%
      group_by(chrtype) %>%
      select(full.lod, mean.lod, var.lod) %>%
      summarise_each(funs(max(., na.rm = TRUE)))
    
    if (is.null(all.perms)) {
      all.perms <- this.perm
    } else {
      all.perms <- rbind(all.perms, this.perm)
    }
  
  }

  return(all.perms)
  
}