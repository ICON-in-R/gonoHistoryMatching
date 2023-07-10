
# gonorrhoea transmission dynamic infectious disease model
# history matching calibration using C++ code
# https://danny-sc.github.io/determ_workshop/
# this is a simplified data testing script

Rcpp::compileAttributes()

library(dplyr)
library(purrr)
library(Rcpp)

sourceCpp("src/GonorrheaDTM.cpp", windowsDebugDLL = FALSE)

res <- runmodel()

# original values for testing
calib_params <-
  c(0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0,0,0,0,0,0,0,0,
    0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0,0,0,0,0,0,0,0,
    0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0,0,0,0,0,0,0,0)

# incidence matrix to vector format
inc_mat_to_vector <- function(x) {
  x_na <- x[ , colSums(is.na(x)) == 0]
  # stack years
  c(as.matrix(x_na))
}

#' DTM helper for testing 
#' 
#' only change a subset of the input calibration parameters
#' and return a subset of the output values
#' for the same groups e.g. age group 1, sexbeh homosexual etc
#' 
#' @param input transmission rates, same length as indx
#' @param indx_in subset input vector index; integer
#' @param indx_out subset output vector index; integer
#' @value incidence by group
#' @examples
#' out <- test_get_results(input = c(0.03, 0.03, 0.03))
#' 
test_get_results <- function(input,
                             indx_in = c(1,2,3,4),
                             indx_out = c(1,2,3,4)) {
  params <- scan(file = "Inputs/Calibration parameters0.txt")
  params[indx_in] <- input
  write.table(params, file = "Inputs/Calibration parameters.txt",
              col.names = FALSE, row.names = FALSE)
  runmodel()
  out <- read.delim(file = "Inputs/Calibrated incidence.txt", header = FALSE)
  
  # rearrange to calibration format
  vout <- inc_mat_to_vector(out)
  vout[indx_out]
}

###################
# input group prep

# create for emulator
#   wave0: data of inputs and outputs from full model
#   ranges_in: upper and lower limits for inputs
#   targets: output calibration central values and sd

# subset vector
indx_in <- c(1,2,3,4)
indx_out <- c(1,2,3,4,      # 2017
              73,74,75,76)  # 2018

ethnicity_grps <- paste0("e", 1:3)
sex_grps <- c("male", "female")
sexbeh_grps <- c("heterosexual", "homosexual", "bisexual")
age_grps <-  paste0("a", 0:3)
inc_years <- paste0("y", 2017:2021)

# all combinations of output covariates
groups_mat <-
  expand.grid(age_grp = age_grps,
              sexbeh = sexbeh_grps,
              sex = sex_grps,
              ethnicity = ethnicity_grps,
              time = inc_years)

# vector of all group names
full_groups_out <- do.call(paste0, groups_mat)

# subset parameters and data for testing
groups_out <- full_groups_out[indx_out]

n_grps_out <- length(groups_out)

# all combinations of input covariates
# without year
groups_in_mat <-
  expand.grid(age_grp = age_grps,
              sexbeh = sexbeh_grps,
              sex = sex_grps,
              ethnicity = ethnicity_grps)

full_groups_in <- do.call(paste0, groups_in_mat)

# subset parameters and data for testing
groups_in <- full_groups_in[indx_in]

n_grps_in <- length(groups_in)

# upper and lower limits for inputs
##TODO: what upper limit?
ranges_in <- 
  rep(list(c(0, 0.1)), n_grps_in) |> 
  setNames(groups_in)

#############################
# create calibration targets

targets_dat <-
  read.delim("Inputs/Calibration targets.txt", sep = "\t", header = FALSE) |> 
  select_if(~ !any(is.na(.)))

# convert to vector
target_val <- data.frame(
  val = inc_mat_to_vector(targets_dat)) |> 
  mutate(sigma = val/500 + 0.1)            ##TODO: this is arbitrary atm

# convert to list
targets <-
  split(target_val, 1:nrow(target_val)) |> 
  setNames(full_groups_out) |> 
  map(as.list)

# subset parameters and data for testing
targets <- targets[indx_out]

# number of full model simulations
n_sim <- n_grps_in*10
n_validation <- 10

# latin hypercube design
# cols: parameters
lhs_points <- lhs::maximinLHS(n_sim, n_grps_in)
lhs_points_validation <- lhs::maximinLHS(n_validation, n_grps_in)

initial_points <- rbind(lhs_points, lhs_points_validation)

# rescale
for (i in 1:n_grps_in) {
  rdiff <- ranges_in[[i]][2] - ranges_in[[i]][2]
  initial_points[, i] <- initial_points[, i]*ranges_in[[i]][2] + rdiff 
}

  
#################
# run full model

# test for single input
# initial_results <- test_get_results(initial_points[i, ], indx_in, indx_out)

initial_results <- t(apply(initial_points, 1,
                           test_get_results, indx_in = indx_in, indx_out = indx_out))

# all initial values
wave0 <-
  cbind(initial_points, initial_results) |> 
  `colnames<-`(c(groups_in, groups_out)) |> 
  as.data.frame()

save(wave0, file = "Outputs/wave0.RData")

# output with known inputs
targets_fake <-
  purrr::map(1:8, ~list(val = 1000,
                        sigma = 100)) |> 
  setNames(groups_out)

# validation set
for (i in 1:n_grps_out) {
  targets_fake[[i]] <- list(val = initial_results[1,i],
                            sigma = 200)
}

###########
# emulator

library(ggplot2)
library(hmer)

##TODO:...

# load("Outputs/wave0.RData")

training <- wave0[1:n_sim, ]
validation <- wave0[(n_sim+1):nrow(wave0), ]

## first wave

ems_wave1 <-
  emulator_from_data(input_data = training,
                     output_names = names(targets),
                     ranges = ranges_in,
                     emulator_type = "deterministic",
                     order = 1,
                     specified_priors = list(hyper_p = rep(0.55, length(targets))))

save(ems_wave1, file = "Outputs/ems_wave1.RData")

emulator_plot(ems_wave1$a0heterosexualmalee1y2017)
emulator_plot(ems_wave1$a1heterosexualmalee1y2017)
emulator_plot(ems_wave1$a2heterosexualmalee1y2017)

# input-output grid
plot_actives(ems_wave1)

# variance
emulator_plot(ems_wave1$a0heterosexualmalee1y2017, plot_type = 'var')

summary(ems_wave1$a0heterosexualmalee1y2017$model)$adj.r.squared

## calibrations plots

emulator_plot(ems_wave1$a0heterosexualmalee1y2017, plot_type = 'imp',
              targets = targets, cb=TRUE)

# implausibility >3 unlikely to give a good fit
emulator_plot(ems_wave1, plot_type = 'imp', targets = targets, cb=TRUE)
emulator_plot(ems_wave1, plot_type = 'nimp', targets = targets, cb=TRUE)  # maximum implausibility

space_removed(ems_wave1, targets, ppd=3) +
  geom_vline(xintercept = 3, lty = 2) + 
  geom_text(aes(x=3, label="x = 3", y=0.33), colour="black", 
            angle=90, vjust = 1.2, text=element_text(size=11))

# fake target data
emulator_plot(ems_wave1, plot_type = 'imp', targets = targets_fake, cb=TRUE)
emulator_plot(ems_wave1, plot_type = 'imp', targets = targets_fake, cb=TRUE,
              params = c('a1heterosexualmalee1', 'a2heterosexualmalee1'))
emulator_plot(ems_wave1, plot_type = 'imp', targets = targets_fake, cb=TRUE,
              params = c('a2heterosexualmalee1', 'a3heterosexualmalee1'))

## second wave

ems_wave1_linear <-
  emulator_from_data(training, names(targets), 
                     ranges, quadratic = FALSE,
                     specified_priors = list(hyper_p = rep(0.55, length(targets))))

R_squared_linear <- list()
for (i in 1:length(ems_wave1_linear)) {
  R_squared_linear[[i]] <- summary(ems_wave1_linear[[i]]$model)$adj.r.squared
}
names(R_squared_linear) <- names(ems_wave1_linear)
unlist(R_squared_linear)

emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c())

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c(), cb = TRUE)

ems_wave1_linear$I200 <-
  ems_wave1_linear$I20$set_hyperparams(
    list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))

emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c())

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c(), cb=TRUE)

wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges), p_size=1)

# Emulator diagnostics

##TODO: errors
vd <- validation_diagnostics(ems_wave1$a3heterosexualmalee1y2017,
                             validation = validation, targets = targets)#, plt=TRUE)

sigmadoubled_emulator <- ems_wave1$a0heterosexualmalee1y2017$mult_sigma(2)
vd <- validation_diagnostics(sigmadoubled_emulator, 
                             validation = validation, targets = targets, plt=TRUE)
