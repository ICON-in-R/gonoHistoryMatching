
#' Parallelised model runner
#' 
#' don't seem to have a read/write race condition problems
#' 
#' see:
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
  
  cl <- parallel::makeCluster(parallel::detectCores() - 1, outfile = "")
  
  on.exit(parallel::stopCluster(cl))
  
  doParallel::registerDoParallel(cl)  # register with foreach package
  
  # # non-parallel
  # init_results <- foreach(i = inits_list) %do%
  #   test_get_results(i, indx_in = indx_in, indx_out = indx_out)
  
  # immutable data
  params0 <- scan(file = "Inputs/Calibration parameters0.txt")
  
  res <- foreach(i = 1:length(inits_list), .combine = 'rbind',
                 .packages = c("Rcpp", "gonoHistoryMatching")#,
                 # .export = c('params0')
  ) %dopar% {
    # calibration values not overwritten
    params <- params0
    params[indx_in] <- inits_list[[i]]
    
    # create separate folder for each run
    dir_name <- paste0("./Inputs/temp/calibration", i, "/")
    dir.create(path = dir_name)
               
    write.table(params, file = paste0(dir_name, "Calibration parameters.txt"),
                col.names = FALSE, row.names = FALSE)
    
    runmodel(dir_name)
    
    out <- read.delim(file = paste0(dir_name, "Calibrated incidence.txt"), header = FALSE)
    
    # rearrange to calibration format
    vout <- inc_mat_to_vector(out)
    vout[indx_out]
    
    # # delete folder
    # unlink(dir_name, recursive = TRUE)
  }
  
  res
}
