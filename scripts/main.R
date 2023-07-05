
# gonorrhoea transmission dynamic infectious disease model
# history matching calibration using C++ code
# https://danny-sc.github.io/determ_workshop/

library(dplyr)
library(Rcpp)

Rcpp::compileAttributes()
sourceCpp("src/GonorrheaDTM.cpp", windowsDebugDLL = FALSE)

transmissionRate <- 1
res <- runmodel(transmissionRate)

# inputs
# transmission rates:
#   race 3
#   gender 2
#   sexual behaviour 3

n_grps <- 3*2*3
ranges <- rep(list(c(0, 100)), n_grps)

# calibration targets

targets_dat <-
  read.delim("Inputs/Calibration targets.txt", sep = "\t", header = FALSE) |> 
  select_if(~ !any(is.na(.)))

target_val <- c(t(as.matrix(targets_dat))) 

targets <- list()
for (i in 1:length(target_val)) {
  targets[[i]] <- list(val = target_val[i], sigma = 0.1)
}

# latin hypercube design
lhs_points <- lhs::maximinLHS(n_grps*10, n_grps)

# rescale
initial_points <- lhs_points
for (i in 1:n_grps) {
  initial_points[, i] <- lhs_points[, i]*ranges[[i]][2]
}

initial_results <- runmodel(initial_points)

# initial values
wave0 <- cbind(initial_points, initial_results)


##TODO: write a get_results helper?
# 
# get_results <- function(input, ...) {
#   write.table(input)  
#   runmodel(...)
#   out <- read.delim()
#   # extract, rearrange values for calibration
# }


###########
# emulator

##TODO:...

## first wave

ems_wave1 <- emulator_from_data(wave0, names(targets), ranges, 
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

ems_wave1_linear$I200 <- ems_wave1_linear$I20$set_hyperparams(
  list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))
emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges), p_size=1)