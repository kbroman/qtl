filter_.scanonevar <- function(vs, ...) {
  
  out <- vs
  class(out) <- class(out)[-1]
  
  out <- dplyr::filter_(out, ...)
  
  class(out) <- class(vs)
  attr(out, 'method') <- attr(vs, 'method')
  attr(out, 'type') <- attr(vs, 'type')
  attr(out, 'model') <- attr(vs, 'model')
  attr(out, 'dom') <- attr(vs, 'dom')
  attr(out, 'pheno') <- attr(vs, 'pheno')
  attr(out, 'null.effects') <- attr(vs, 'null.effects')
  attr(out, 'units') <- attr(vs, 'units')
  
  return(out)
}
