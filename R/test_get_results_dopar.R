
#' Parallelised version
#' 
#' see
#' https://stackoverflow.com/questions/58823105/use-doparallel-with-rcpp-function-on-window-inside-an-r-package
#' 
#' @import parallel foreach Rcpp doParallel
#' @importFrom purrr array_branch
#' @export
#' 
test_get_results_dopar <- function(init_points, indx_in, indx_out) {
  
  # incidence matrix to vector format
  inc_mat_to_vector <- function(x) {
    x_na <- x[, colSums(is.na(x)) == 0]
    # stack years
    c(as.matrix(x_na))
  }
  
  inits_list <- purrr::array_branch(init_points, margin = 1)
  
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  
  on.exit(parallel::stopCluster(cl))
  
  doParallel::registerDoParallel(cl)  # register with foreach package
  
  # # non-parallel
  # init_results <- foreach(i = inits_list) %do%
  #   test_get_results(i, indx_in = indx_in, indx_out = indx_out)
  
  res <- foreach(i = 1:length(inits_list), .combine = 'c',
                 .packages = c("Rcpp", "gonoHistoryMatching"),
                 .export = c('inc_mat_to_vector')
                 # .noexport = c('loadCalibrationParameters',
                 #               'loadCalibratioTargets',
                 #               'loadDemographics',
                 #               'loadParameters',
                 #               'updatedCalibrationParameters')
  ) %dopar% {
    # calibration values not overwritten
    params <- scan(file = "Inputs/Calibration parameters0.txt")
    params[indx_in] <- inits_list[[i]]
    write.table(params, file = "Inputs/Calibration parameters.txt",
                col.names = FALSE, row.names = FALSE)
    runmodel()
    out <- read.delim(file = "Inputs/Calibrated incidence.txt", header = FALSE)
    
    # rearrange to calibration format
    vout <- inc_mat_to_vector(out)
    vout[indx_out]
  }
  
  return(res)
}
