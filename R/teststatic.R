
# R wrapper for testing static variables within C

teststatic <-
function(start=1, num=3)
{
  z <- .C("R_teststatic",
          as.integer(start),
          as.integer(num),
          PACKAGE="qtl")
  invisible(NULL)
}
