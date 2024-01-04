
#' Take raw Latin Hypercube samples and rescale
#' 
#' @param ranges List of mean, sd for each parameter
#' @param vals Matrix of samples from Latin hypercube uniform distribution.
#'             Dimensions sample by parameter
#' @return Mean, upper and lower limit adjusted
#' 
rescale <- function(ranges, vals) {
  for (i in seq_along(ranges)) {
    rdiff <- ranges[[i]][2] - ranges[[i]][1]
    vals[, i] <- vals[, i]*rdiff + ranges[[i]][1] 
  }
  vals
}


#' Combined named inputs and outputs
#' 
cbind_input_output <- function(param_nm) {
  force(param_nm)
  function(inp, out) {
    cbind(inp, out) |> 
      `colnames<-`(param_nm) |> 
      as.data.frame()
  }
}
