
# gonorrhoea transmission dynamic infectious disease model
# history matching calibration using C++ code
# https://danny-sc.github.io/determ_workshop/
Rcpp::compileAttributes()

library(dplyr)
library(purrr)
library(Rcpp)

sourceCpp("src/GonorrheaDTM.cpp", windowsDebugDLL = FALSE)

res <- runmodel()

# original values
# for testing
calib_params <- c(0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0,0,0,0,0,0,0,0,
                  0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0,0,0,0,0,0,0,0,
                  0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0,0,0,0,0,0,0,0)

#
inc_mat_to_vector <- function(x) {
  x_na <- x[ , colSums(is.na(x)) == 0]
  c(t(as.matrix(x_na)))
}

#' DTM helper 
#' @param input transmission rate
#' @value incidence by group
#' @examples
#' out <- get_results(calib_params)
#' 
get_results <- function(input) {
  write.table(input, file = "Inputs/Calibration parameters.txt",
              col.names = FALSE, row.names = FALSE)
  runmodel()
  out <- read.delim(file = "Inputs/Calibrated incidence.txt", header = FALSE)
  
  # rearrange to calibration format
  inc_mat_to_vector(out)
}

#############
# input prep

ethnicity_grps <- c("a","b","c")
sex_grps <- c("male", "female")
sexbeh_grps <- c("heterosexual", "homosexual", "bisexual")
age_grps <- c(0,1,2,3)
inc_years <- 2017:2021

groups_mat <-
  expand.grid(age_grp = age_grps,
              sexbeh = sexbeh_grps,
              sex = sex_grps,
              ethnicity = ethnicity_grps,
              time = inc_years)

groups_out <- do.call(paste, groups_mat)

n_grps_out <- length(groups_out)

groups_in_mat <-
  expand.grid(age_grp = age_grps,
              sexbeh = sexbeh_grps,
              sex = sex_grps,
              ethnicity = ethnicity_grps)

groups_in <- do.call(paste, groups_in_mat)

n_grps_in <- length(groups_in)
ranges_in <- rep(list(c(0, 1)), n_grps_in) |> 
  setNames(groups_in)

# calibration targets

targets_dat <-
  read.delim("Inputs/Calibration targets.txt", sep = "\t", header = FALSE) |> 
  select_if(~ !any(is.na(.)))

# convert to vector
target_val <- data.frame(
  val = inc_mat_to_vector(targets_dat),
  sigma = 0.1)

# convert to list
targets <-
  split(target_val, 1:nrow(target_val)) |> 
  setNames(groups_out) |> 
  map(as.list)

# n_sim <- n_grps_in*10
n_sim <- 2

# latin hypercube design
# cols: parameters
lhs_points <- lhs::maximinLHS(n_sim, n_grps_in)

# rescale
initial_points <- lhs_points
for (i in 1:n_grps_in) {
  rdiff <- ranges_in[[i]][2] - ranges_in[[i]][2]
  initial_points[, i] <- lhs_points[, i]*ranges_in[[i]][2] + rdiff 
}

# test for single input
# initial_results <- get_results(initial_points[i, ])

initial_results <- t(apply(initial_points, 1, get_results))

# all initial values
wave0 <- cbind(initial_points, initial_results) |> 
  `colnames<-`(c(groups_in, groups_out))


###########
# emulator

library(hmer)

##TODO:...

## first wave

ems_wave1 <-
  emulator_from_data(input_data = wave0,
                     output_names = names(targets),
                     ranges = ranges_in, 
                     specified_priors = list(hyper_p = rep(0.55, length(targets))))

emulator_plot(ems_wave1$R200, params = c('beta1', 'gamma'))

plot_actives(ems_wave1)

emulator_plot(ems_wave1$R200, plot_type = 'var', params = c('beta1', 'gamma'))

summary(ems_wave1$R200$model)$adj.r.squared

emulator_plot(ems_wave1$R200, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)

emulator_plot(ems_wave1, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)

space_removed(ems_wave1, targets, ppd=3) + geom_vline(xintercept = 3, lty = 2) + 
  geom_text(aes(x=3, label="x = 3",y=0.33), colour="black", 
            angle=90, vjust = 1.2, text=element_text(size=11))

## second wave

ems_wave1_linear <- emulator_from_data(training, names(targets), 
                                       ranges, quadratic=FALSE,
                                       specified_priors = list(hyper_p = rep(0.55, length(targets))))

R_squared_linear <- list()
for (i in 1:length(ems_wave1_linear)) {
  R_squared_linear[[i]] <- summary(ems_wave1_linear[[i]]$model)$adj.r.squared
}
names(R_squared_linear) <- names(ems_wave1_linear)
unlist(R_squared_linear)

emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

ems_wave1_linear$I200 <-
  ems_wave1_linear$I20$set_hyperparams(
    list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))

emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges), p_size=1)

# Emulator diagnostics

vd <- validation_diagnostics(ems_wave1$R200, validation = validation, targets = targets, plt=TRUE)

sigmadoubled_emulator <- ems_wave1$R200$mult_sigma(2)
vd <- validation_diagnostics(sigmadoubled_emulator, 
                             validation = validation, targets = targets, plt=TRUE)
