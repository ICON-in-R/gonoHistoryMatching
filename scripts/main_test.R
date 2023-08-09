
# gonorrhoea transmission dynamic infectious disease model
# history matching calibration using C++ code
# input: transmission rate
# output: incidence counts
# https://danny-sc.github.io/determ_workshop/
# this is a simplified data testing script for
# a subset in inputs

Rcpp::compileAttributes()

library(dplyr)
library(purrr)
library(Rcpp)

sourceCpp("src/GonorrheaDTM.cpp", windowsDebugDLL = FALSE)

if (FALSE)
  res <- runmodel()

############
# functions

# incidence matrix to vector format
inc_mat_to_vector <- function(x) {
  x_na <- x[, colSums(is.na(x)) == 0]
  # stack years
  c(as.matrix(x_na))
}

#' DTM helper for testing subset of inputs/outputs
#' 
#' Only change a subset of the input calibration parameters
#' and return a subset of the output values
#' for the same groups e.g. age group 1, sexbeh homosexual etc
#' 
#' @param input transmission rates, same length as indx_in
#' @param indx_in subset of total input vector index; integer vector
#' @param indx_out subset of total output vector index; integer vector
#' 
#' @return incidence by group
#' @examples
#' out <- test_get_results(input = c(0.03, 0.03, 0.03))
#'
test_get_results <- function(input,
                             indx_in = c(1,2,3,4),
                             indx_out = c(1,2,3,4)) {
  # calibration values not overwritten
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
#   - wave0: data of inputs and outputs from full model
#   - ranges_in: upper and lower limits for transmission rates
#   - targets: incidence output calibration central values and sd

# subset vectors
# all ages
# ethnicity group 1
indx_in <- c(1,2,3,4,       # male
             13,14,15,16)   # female

indx_out <- c(1,2,3,4,      # 2017
              13,14,15,16,
              73,74,75,76,  # 2018
              85,86,87,88)

ethnicity_grps <- paste0("e", 1:3)
sex_grps <- c("male", "female")
sexbeh_grps <- c("heterosexual", "homosexual", "bisexual")
age_grps <-  paste0("a", 0:3)
inc_years <- paste0("y", 2017:2021)

# all combinations of _output_ covariates
groups_mat <-
  expand.grid(age_grp = age_grps,
              sexbeh = sexbeh_grps,
              sex = sex_grps,
              ethnicity = ethnicity_grps,
              time = inc_years)

# vector of all output group names
full_groups_out <- do.call(paste0, groups_mat)

# subset parameters and data for testing
groups_out <- full_groups_out[indx_out]

n_grps_out <- length(groups_out)

# all combinations of _input_ covariates
# difference with output is without year
groups_in_mat <-
  expand.grid(age_grp = age_grps,
              sexbeh = sexbeh_grps,
              sex = sex_grps,
              ethnicity = ethnicity_grps)

# vector of all input group names
full_groups_in <- do.call(paste0, groups_in_mat)

# subset parameters and data for testing
groups_in <- full_groups_in[indx_in]

n_grps_in <- length(groups_in)

# upper and lower limits for input transmission rates
##TODO: what upper limit?
# point value in file is 0.022
ranges_in <- 
  rep(list(c(0, 0.05)), n_grps_in) |> 
  setNames(groups_in)

######################################
# create calibration targets (output)

# incidence group by year
targets_dat <-
  read.delim("Inputs/Calibration targets.txt", sep = "\t", header = FALSE) |> 
  select_if(~ !any(is.na(.)))

# convert to vector
# mean and standard deviation
target_val <- data.frame(
  val = inc_mat_to_vector(targets_dat)) |> 
  mutate(sigma = val/200 + 0.1)            ##TODO: this is arbitrary atm

# convert to named list
all_targets <-
  split(target_val, 1:nrow(target_val)) |> 
  setNames(full_groups_out) |> 
  map(as.list)

# subset parameters and data for testing
targets <- all_targets[indx_out]

if (save)
  save(targets, file = "Outputs/targets.RData")

# number of full model simulations
n_sim <- n_grps_in*10
n_validation <- n_grps_in

##################################
# latin hypercube design (inputs)
# cols: parameters
lhs_points <- lhs::maximinLHS(n_sim, n_grps_in)
lhs_points_validation <- lhs::maximinLHS(n_validation, n_grps_in)

# all LHS inputs
init_points <-
  rbind(lhs_points, lhs_points_validation) |> 
  `colnames<-`(groups_in)

# rescale
for (i in 1:n_grps_in) {
  rdiff <- ranges_in[[i]][2] - ranges_in[[i]][2]
  init_points[, i] <- init_points[, i]*ranges_in[[i]][2] + rdiff 
}

if (save)
  save(init_points, file = "Outputs/init_points.RData")

#################
# run full model

# test for single input
if (FALSE) {
  init_results <- test_get_results(init_points[i, ], indx_in, indx_out)
}

# run model for all LHS inputs
init_results <- t(apply(init_points, 1,
                        test_get_results, indx_in = indx_in, indx_out = indx_out))

# all named initial inputs and outputs
wave0 <-
  cbind(init_points, init_results) |> 
  `colnames<-`(c(groups_in, groups_out)) |> 
  as.data.frame()

if (FALSE)
  save(wave0, file = "Outputs/wave0.RData")

# # output incidence with known inputs
# # so can check calibration against
# targets_fake <-
#   purrr::map(1:n_grps_out,
#              ~list(val = 1000,
#                    sigma = 100)) |> 
#   setNames(groups_out)
# 
# # or
# # single validation set list of outputs
# for (i in 1:n_grps_out) {
#   targets_fake[[i]] <- list(val = init_results[1, i],
#                             sigma = 200)
# }

###########
# emulator
###########

library(ggplot2)
library(hmer)

save <- FALSE

if (FALSE)
  load("Outputs/wave0.RData")  # inputs

# split data set
training <- wave0[1:n_sim, ]
validation <- wave0[(n_sim+1):nrow(wave0), ]

#############
# first wave

ems_wave1 <-
  emulator_from_data(input_data = training,            # named inputs and outputs from full model
                     output_names = names(targets),
                     range = ranges_in,               # min, max inputs
                     emulator_type = "deterministic",
                     order = 2) #,                     # of regression 
# specified_priors = list(hyper_p = rep(0.55, length(targets))))

if (save)
  save(ems_wave1, file = "Outputs/ems_wave1.RData")

# contour plot of pair of input parameters
# fills parameter space not sampled in full model
emulator_plot(ems_wave1$a0heterosexualmalee1y2017)
emulator_plot(ems_wave1$a1heterosexualmalee1y2017)
emulator_plot(ems_wave1$a2heterosexualmalee1y2017)
emulator_plot(ems_wave1$a3heterosexualmalee1y2017)

# input-output grid
plot_actives(ems_wave1)

# variance
emulator_plot(ems_wave1$a0heterosexualmalee1y2017, plot_type = 'var')

# R^2 statistic
summary(ems_wave1$a0heterosexualmalee1y2017$model)$adj.r.squared

## calibrations plots

# implausibility
# >3 unlikely to give a good fit
emulator_plot(ems_wave1$a0heterosexualmalee1y2017, plot_type = 'imp',
              targets = targets, cb=TRUE)

emulator_plot(ems_wave1, plot_type = 'imp', targets = targets, cb=TRUE)
emulator_plot(ems_wave1, params = c("a0heterosexualmalee1", "a0heterosexualfemalee1"), plot_type = 'imp', targets = targets, cb=TRUE)
emulator_plot(ems_wave1, params = c("a2heterosexualmalee1", "a2heterosexualfemalee1"), plot_type = 'imp', targets = targets, cb=TRUE)

# nth-maximum implausibility
emulator_plot(ems_wave1, plot_type = 'nimp', targets = targets, cb=TRUE)  # maximum implausibility

space_removed(ems_wave1, targets, ppd=3) +
  geom_vline(xintercept = 3, lty = 2) + 
  geom_text(aes(x=3, label="x = 3", y=0.33),
            colour="black", angle=90, vjust = 1.2, text=element_text(size=11))

# # fake target data
# emulator_plot(ems_wave1, plot_type = 'imp', targets = targets_fake, cb=TRUE)
# emulator_plot(ems_wave1, plot_type = 'imp', targets = targets_fake, cb=TRUE,
#               params = c('a1heterosexualmalee1', 'a2heterosexualmalee1'))
# emulator_plot(ems_wave1, plot_type = 'imp', targets = targets_fake, cb=TRUE,
#               params = c('a2heterosexualmalee1', 'a3heterosexualmalee1'))

# Emulator diagnostics
vd <- validation_diagnostics(ems_wave1,
                             validation = wave0, targets = targets, plt = TRUE)
                             # validation = validation[-8, ], targets = targets, plt = TRUE)

sigmadoubled_emulator <- ems_wave1$a0heterosexualmalee1y2017$mult_sigma(2)

vd <- validation_diagnostics(sigmadoubled_emulator,
                             validation = validation, targets = targets, plt = TRUE)

##TODO: errors
vd <- validation_diagnostics(ems_wave1,
                             validation = validation[-8,], targets = targets, plt=TRUE)

sigmadoubled_emulator <- ems_wave1$a0heterosexualmalee1y2017$mult_sigma(2)
vd <- validation_diagnostics(sigmadoubled_emulator, 
                             validation = validation, targets = targets, plt=TRUE)

##############
# second wave

wave1_points <- generate_new_design(ems_wave1, 60, targets, verbose = TRUE)
# wave1_points <- generate_new_runs(ems_wave1, n_points = 60, targets, verbose = TRUE)

if (save)
  save(wave1_points, file = "Outputs/wave1_points.RData")

# grid of LHS sampled inputs omitting implausible regions
plot_wrap(wave1_points, ranges = ranges_in)

# rerun model for all LHS inputs
wave1_results <- t(apply(wave1_points, 1,
                         test_get_results, indx_in = indx_in, indx_out = indx_out)) #|> 
  # setNames()

wave1 <-
  cbind(wave1_points, wave1_results) |> 
  `colnames<-`(c(groups_in, groups_out)) |> 
  as.data.frame()

if (save)
  save(wave1, file = "Outputs/wave1.RData")

n_sample <- 50
wave1_training <- wave1[1:n_sample, ]
wave1_validation <- wave1[(n_sample+1):nrow(wave1), ]

ems_wave2 <- emulator_from_data(input_data = wave1_training,
                                output_names = names(targets),
                                check.ranges = TRUE,           # update for plausible only
                                ranges = ranges_in,
                                emulator_type = "deterministic",
                                order = 2)
if (save)
  save(ems_wave2, file = "Outputs/ems_wave2.RData")

# contour plot after history matching
emulator_plot(ems_wave2, plot_type = 'imp', targets = targets, cb=TRUE)
emulator_plot(ems_wave2, params = c("a0heterosexualmalee1", "a0heterosexualfemalee1"), plot_type = 'imp', targets = targets, cb=TRUE)
emulator_plot(ems_wave2, params = c("a2heterosexualmalee1", "a2heterosexualfemalee1"), plot_type = 'imp', targets = targets, cb=TRUE)

vd <- validation_diagnostics(ems_wave2, validation = wave1_validation, targets = targets, 
                             plt = TRUE)

# inflate sigma for better fit
ems_wave2_mult_sigma <- ems_wave2$a0heterosexualmalee1y2017$mult_sigma(2)

vd <- validation_diagnostics(ems_wave2_mult_sigma, validation = wave1_validation, targets = targets, 
                             plt = TRUE)

# wave2_points <- generate_new_runs(c(ems_wave2, ems_wave1), 60, targets, verbose = TRUE)
wave2_points <- generate_new_runs(ems_wave2[-c(5,6,13)], 60, targets, verbose = TRUE)

plot_wrap(wave2_points, ranges = ranges_in)

if (save)
  save(wave2_points, file = "Outputs/wave2_points.RData")

##TODO:
## more runs?...

# visualisations of non-implausible space by wave
all_points <- list(as.data.frame(init_points), wave1_points, wave2_points)
wave_points(all_points, input_names = colnames(wave1_points), p_size=1, zero_in = T)
# wave_points(all_points, input_names = names(ranges_in), p_size=1)

all_waves <- list(wave0, wave1) #, wave2)
wave_values(all_waves, targets, l_wid=1, p_size=1)

behaviour_plot(ems_wave1, targets = targets)
behaviour_plot(ems_wave2, targets = targets)

