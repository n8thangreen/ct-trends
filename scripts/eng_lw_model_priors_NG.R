### Maximum/Minimum estimates for Prevalence, Incidence and Screening Rate of Chlamydia in England 
### LW model parameterised with LW model prior parameters


dir_joanna <- "C:/Users/ngreen1/Documents/clamydia-Bristol-data/joanna/ct_trends-master/"

################
# read in data #
################

## from 2000 to 2015
# diagnoses
diags_f_max <- read.csv(paste0(dir_joanna, 'data/diagnoses_f_max.csv'))
diags_f_min <- read.csv(paste0(dir_joanna, 'data/diagnoses_f_min.csv'))
diags_m_max <- read.csv(paste0(dir_joanna, 'data/diagnoses_m_max.csv'))
diags_m_min <- read.csv(paste0(dir_joanna, 'data/diagnoses_m_min.csv')) 
# tests
tests_f_max <- read.csv(paste0(dir_joanna, 'data/tests_f_max.csv'))
tests_f_min <- read.csv(paste0(dir_joanna, 'data/tests_f_min.csv'))
tests_m_max <- read.csv(paste0(dir_joanna, 'data/tests_m_max.csv'))
tests_m_min <- read.csv(paste0(dir_joanna, 'data/tests_m_min.csv'))
# population
pop_f <- read.csv(paste0(dir_joanna, 'data/population_f.csv'))
pop_m <- read.csv(paste0(dir_joanna, 'data/population_m.csv'))


########
# prep #
########

### Start Timer
start_time <- Sys.time()

### 10,000 in the LW model paper
n_sample <- 5000

### Organise population data
pop_f_15_19 <- pop_f[,2]
pop_f_20_24 <- pop_f[,3]
pop_m_15_19 <- pop_m[,2]
pop_m_20_24 <- pop_m[,3]

### number of years we have data
years <- length(pop_f_15_19)


### Q: why was this commented out?
### Population stats need to be changed to count only sexually active individuals
# proportion sexually active (from table 2 of 2017 LW model paper)

p_active_f_16_19 <- rbeta(n_sample, shape1 = 532, shape2 = 238)
p_active_f_20_24 <- rbeta(n_sample, shape1 = 569, shape2 = 53.8)
p_active_m_16_19 <- rbeta(n_sample, shape1 = 506, shape2 = 208)
p_active_m_20_24 <- rbeta(n_sample, shape1 = 467, shape2 = 45.8)

# population sexually active in each age-sex group each year
pop_active_f_15_19 <- matrix(0, years, n_sample)
pop_active_f_20_24 <- matrix(0, years, n_sample)
pop_active_m_15_19 <- matrix(0, years, n_sample)
pop_active_m_20_24 <- matrix(0, years, n_sample)

for (i in seq_len(years)){
  pop_active_f_15_19[i,] <- rbinom(n_sample, pop_f_15_19[i], p_active_f_16_19)
  pop_active_f_20_24[i,] <- rbinom(n_sample, pop_f_20_24[i], p_active_f_20_24)
  pop_active_m_15_19[i,] <- rbinom(n_sample, pop_m_15_19[i], p_active_m_16_19)
  pop_active_m_20_24[i,] <- rbinom(n_sample, pop_m_20_24[i], p_active_m_20_24)
}

### Organise diagnoses data
# note: only 15-19 & 20-24 age groups have diagnoses data from 2000 to 2015
diags_f_15_19_max <- diags_f_max[ ,2]
diags_f_15_19_min <- diags_f_min[ ,2]
diags_f_20_24_max <- diags_f_max[ ,3]
diags_f_20_24_min <- diags_f_min[ ,3]
diags_m_15_19_max <- diags_m_max[ ,2]
diags_m_15_19_min <- diags_m_min[ ,2]
diags_m_20_24_max <- diags_m_max[ ,3]
diags_m_20_24_min <- diags_m_min[ ,3]
# create new 2D matrices for the data
diag_rate_f_15_19_max <- matrix(0, years, n_sample)
diag_rate_f_15_19_min <- matrix(0, years, n_sample)
diag_rate_f_20_24_max <- matrix(0, years, n_sample)
diag_rate_f_20_24_min <- matrix(0, years, n_sample)
diag_rate_m_15_19_max <- matrix(0, years, n_sample)
diag_rate_m_15_19_min <- matrix(0, years, n_sample)
diag_rate_m_20_24_max <- matrix(0, years, n_sample)
diag_rate_m_20_24_min <- matrix(0, years, n_sample)

rgamma_diagrate <- function(n_sample, diags, pop)
  rgamma(n_sample, shape = diags*pop/100000)/pop

# note: diagnoses are per 100,000
for (i in seq_len(years)){
  diag_rate_f_15_19_max[i,] <- rgamma_diagrate(n_sample, diags_f_15_19_max[i], pop_f_15_19[i])
  diag_rate_f_15_19_min[i,] <- rgamma(n_sample, diags_f_15_19_min[i]*pop_f_15_19[i]/100000)/pop_f_15_19[i]
  diag_rate_f_20_24_max[i,] <- rgamma(n_sample, diags_f_20_24_max[i]*pop_f_20_24[i]/100000)/pop_f_20_24[i]
  diag_rate_f_20_24_min[i,] <- rgamma(n_sample, diags_f_20_24_min[i]*pop_f_20_24[i]/100000)/pop_f_20_24[i]
  diag_rate_m_15_19_max[i,] <- rgamma(n_sample, diags_m_15_19_max[i]*pop_m_15_19[i]/100000)/pop_m_15_19[i]
  diag_rate_m_15_19_min[i,] <- rgamma(n_sample, diags_m_15_19_min[i]*pop_m_15_19[i]/100000)/pop_m_15_19[i]
  diag_rate_m_20_24_max[i,] <- rgamma(n_sample, diags_m_20_24_max[i]*pop_m_20_24[i]/100000)/pop_m_20_24[i]
  diag_rate_m_20_24_min[i,] <- rgamma(n_sample, diags_m_20_24_min[i]*pop_m_20_24[i]/100000)/pop_m_20_24[i]
}

### Organise testing data 
# note: only 15-19 & 20-24 age groups have diagnoses data from 2000 to 2015
tests_f_15_19_max <- tests_f_max[ ,2]
tests_f_15_19_min <- tests_f_min[ ,2]
tests_f_20_24_max <- tests_f_max[ ,3]
tests_f_20_24_min <- tests_f_min[ ,3]
tests_m_15_19_max <- tests_m_max[ ,2]
tests_m_15_19_min <- tests_m_min[ ,2]
tests_m_20_24_max <- tests_m_max[ ,3]
tests_m_20_24_min <- tests_m_min[ ,3]

# create new 2D matrices for the data
test_rate_f_15_19_max <- matrix(0, years, n_sample)
test_rate_f_15_19_min <- matrix(0, years, n_sample)
test_rate_f_20_24_max <- matrix(0, years, n_sample)
test_rate_f_20_24_min <- matrix(0, years, n_sample)
test_rate_m_15_19_max <- matrix(0, years, n_sample)
test_rate_m_15_19_min <- matrix(0, years, n_sample)
test_rate_m_20_24_max <- matrix(0, years, n_sample)
test_rate_m_20_24_min <- matrix(0, years, n_sample)

# note: tests are % of total population
for (i in seq_len(years)){
  test_rate_f_15_19_max[i,] <- rgamma(n_sample, tests_f_15_19_max[i]*pop_f_15_19[i]/100)/pop_f_15_19[i]
  test_rate_f_15_19_min[i,] <- rgamma(n_sample, tests_f_15_19_min[i]*pop_f_15_19[i]/100)/pop_f_15_19[i]
  test_rate_f_20_24_max[i,] <- rgamma(n_sample, tests_f_20_24_max[i]*pop_f_20_24[i]/100)/pop_f_20_24[i]
  test_rate_f_20_24_min[i,] <- rgamma(n_sample, tests_f_20_24_min[i]*pop_f_20_24[i]/100)/pop_f_20_24[i]
  test_rate_m_15_19_max[i,] <- rgamma(n_sample, tests_m_15_19_max[i]*pop_m_15_19[i]/100)/pop_m_15_19[i]
  test_rate_m_15_19_min[i,] <- rgamma(n_sample, tests_m_15_19_min[i]*pop_m_15_19[i]/100)/pop_m_15_19[i]
  test_rate_m_20_24_max[i,] <- rgamma(n_sample, tests_m_20_24_max[i]*pop_m_20_24[i]/100)/pop_m_20_24[i]
  test_rate_m_20_24_min[i,] <- rgamma(n_sample, tests_m_20_24_min[i]*pop_m_20_24[i]/100)/pop_m_20_24[i]
}

### Parameters for sensitivity (true_pos) and 1-specificity (false_pos)
# beta parameters from 2017 Epidemniology paper table 2:
p_true_pos_f <- rbeta(shape1 = 130, shape2 = 12, n_sample)
p_true_pos_m <- rbeta(shape1 = 33, shape2 = 1, n_sample)
p_false_pos_f <- rbeta(shape1 = 5, shape2 = 2324, n_sample)
p_false_pos_m <- rbeta(shape1 = 3, shape2 = 951, n_sample)

########
# save #
########

write.csv(p_true_pos_f, file = '../../output/eng_lw_model_priors_p_true_pos_f.csv')
write.csv(p_true_pos_m, file = '../../output/eng_lw_model_priors_p_true_pos_m.csv')
write.csv(p_false_pos_f, file = '../../output/eng_lw_model_priors_p_false_pos_f.csv')
write.csv(p_false_pos_m, file = '../../output/eng_lw_model_priors_p_false_pos_m.csv')


### Parameters for rate of spontaneous recovery (spon_rec) per year
# created from an MCMC using STAN software, available from https://github.com/joanna-lewis/ct_surveillance/tree/master/england/stan 
spon_rec_f <- read.csv(paste0(dir_joanna,  file = 'stan/chlamydia_two_exponentials_women.csv'))[1:n_sample, 1]
spon_rec_m <- read.csv(paste0(dir_joanna,  file = 'stan/chlamydia_two_exponentials_men.csv'))[1:n_sample, 1]

# save locally
write.csv(spon_rec_f, file = '../../output/eng_lw_model_priors_spon_rec_f.csv')
write.csv(spon_rec_m, file = '../../output/eng_lw_model_priors_spon_rec_m.csv')


# rewritten start --------------------------------------------------------------


### Parameter for rate of treatment seeking per year (att_symp)
# Mercer Sex. Transm. Infect. (2007) (see table above).

# Find beta distributions corresponding to 95% CIs reported in
library(rriskDistributions)
library(purrr)

partial_beta <- partial(get.beta.par,
                        plot = FALSE, show.output = FALSE)

under1week <- partial_beta(p = c(0.025, 0.975),
                           q = c(0.144, 0.442))
weeks1to2  <- partial_beta(p = c(0.025, 0.975),
                           q = c(0.061, 0.302))
weeks2to4  <- partial_beta(p = c(0.025, 0.975),
                           q = c(0.133, 0.310))
weeks4to6  <- partial_beta(p = c(0.025, 0.975),
                           q = c(0.085, 0.299))
over6weeks <- partial_beta(p = c(0.025, 0.975),
                           q = c(0.055, 0.564))

beta_params <- 
  list(
    "<1 week" = c(
      a = as.numeric(under1week["shape1"]),
      b = as.numeric(under1week["shape2"])),
    "7-13 days" = c(
      a = as.numeric(weeks1to2["shape1"]),
      b = as.numeric(weeks1to2["shape2"])),
    "14-27 days" = c(
      a = as.numeric(weeks2to4["shape1"]),
      b = as.numeric(weeks2to4["shape2"])),
    "28-41 days" = c(
      a = as.numeric(weeks4to6["shape1"]),
      b = as.numeric(weeks4to6["shape2"])),
    "42 days and over" = c(
      a = as.numeric(over6weeks["shape1"]),
      b = as.numeric(over6weeks["shape2"]))
  )
  

## MH for rate of treatment

NEG_NUM <<- -1e10
burnin <- 1000
h_current <- 0.04
sd_h <- 0.05
niter <- n_sample + burnin

# proportion in each time window
tps <- c(0.0, 7.0, 14.0, 28.0, 42.0, Inf)

# partial functions
prop_in_each_period <- prop_periods(tps)


########
# main #
########

res <- runMCMC(ninter,
               h_current,
               sd_h,
               beta_params)

att_symp <- res$h_symp
loglik <- res$loglik

att_symp = att_symp[burnin:niter] # remove burn-in
loglik = loglik[burnin:length(loglik)]
att_symp = att_symp * 365.25 # convert rate from per day to per year

## save
write.csv(att_symp, file = here::here('../../output/eng_lw_model_priors_att_symp.csv'))


### Parameters for proportion of infections asymptomatic (p_asymp)
# data from table 2 of the 2017 LW model paper

p_asymp_f_alpha <-
  get.beta.par(p = c(0.025, 0.975),
               q = c(0.468, 0.752),
               plot = FALSE)[1]
p_asymp_f_beta <-
  get.beta.par(p = c(0.025, 0.975),
               q = c(0.468, 0.752),
               plot = FALSE)[2]
p_asymp_f <-
  rbeta(n_sample, shape1 = p_asymp_f_alpha, shape2 = p_asymp_f_beta)

p_asymp_m_alpha <-
  get.beta.par(p = c(0.025, 0.975),
               q = c(0.264, 0.759),
               plot = FALSE)[1]
p_asymp_m_beta <-
  get.beta.par(p = c(0.025, 0.975),
               q = c(0.264, 0.759),
               plot = FALSE)[2]
p_asymp_m <-
  rbeta(n_sample, shape1 = p_asymp_m_alpha, shape2 = p_asymp_m_beta)

# save
write.csv(p_asymp_f, file = '../../output/eng_lw_model_priors_p_asymp_f.csv')
write.csv(p_asymp_m, file = '../../output/eng_lw_model_priors_p_asymp_m.csv')

# rewritten end ---------------------------------------------------------------


### sol_dyn, model_test_diag, test_diag_fun functions from https://github.com/joostsmid 
library(rriskDistributions)
library(deSolve)
library(rootSolve)
library(ggplot2)
library(grid)

sol_dyn <- list(
  S_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) {
    S <- alpha_AU*alpha_US/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(S)},
  
  U_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) {
    U <- alpha_AU*alpha_SU/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(U)},
  
  A_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) {
    A <- alpha_SU*alpha_UA/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(A)},
  
  prev_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) { 
    # this is just S + A
    prev <- alpha_AU*alpha_US/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA)) +
      alpha_SU*alpha_UA/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(prev)}
)

model_test_diag <- list(
  test_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    tsym <- ( ssym + (1 - A - U)*test_sym ) # number of tests
    return(tsym)},

  diag_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    dsym <- ( A*ssym*true_pos + U*ssym*false_pos + (1 - A - U)*test_sym*true_pos ) # number of diagnoses
    return(dsym)}
)

test_diag_fun <- function(parms, p_symp, self_cure, test_sym, true_pos, false_pos, cov=NULL, adpc=NULL){
  # function that returns number of tests and number of diagnoses, given a certain incidence and screening rate (parms)
  inc <- parms[1] # incidence
  scr <- parms[2] # screening rate
  if (is.null(cov)) cov <- 0
  if (is.null(adpc)) adpc <- 0
  A = sol_dyn$A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, scr*true_pos + test_sym*true_pos)
  U = sol_dyn$U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, scr*true_pos + test_sym*true_pos)
  test <- model_test_diag$test_fun(A, U, scr, test_sym, true_pos, false_pos) - cov
  diag <- model_test_diag$diag_fun(A, U, scr, test_sym, true_pos, false_pos) - adpc
  return(c(test=test, diag=diag))
}


### Create empty matrices for the solutions
# Incidence of chlamydia 
inc_f_15_19_max <- matrix(0, years, n_sample)
inc_f_15_19_min <- matrix(0, years, n_sample)
inc_f_20_24_max <- matrix(0, years, n_sample)
inc_f_20_24_min <- matrix(0, years, n_sample)
inc_m_15_19_max <- matrix(0, years, n_sample)
inc_m_15_19_min <- matrix(0, years, n_sample)
inc_m_20_24_max <- matrix(0, years, n_sample)
inc_m_20_24_min <- matrix(0, years, n_sample)
# force of infection (foi) of chlamydia
foi_f_15_19_max <- matrix(0, years, n_sample)
foi_f_15_19_min <- matrix(0, years, n_sample)
foi_f_20_24_max <- matrix(0, years, n_sample)
foi_f_20_24_min <- matrix(0, years, n_sample)
foi_m_15_19_max <- matrix(0, years, n_sample)
foi_m_15_19_min <- matrix(0, years, n_sample)
foi_m_20_24_max <- matrix(0, years, n_sample)
foi_m_20_24_min <- matrix(0, years, n_sample)
# Screening rate for chlamydia
scr_f_15_19_max <- matrix(0, years, n_sample)
scr_f_15_19_min <- matrix(0, years, n_sample)
scr_f_20_24_max <- matrix(0, years, n_sample)
scr_f_20_24_min <- matrix(0, years, n_sample)
scr_m_15_19_max <- matrix(0, years, n_sample)
scr_m_15_19_min <- matrix(0, years, n_sample)
scr_m_20_24_max <- matrix(0, years, n_sample)
scr_m_20_24_min <- matrix(0, years, n_sample)
# Prevalence of chlamydia
prev_f_15_19_max <- matrix(0, years, n_sample)
prev_f_15_19_min <- matrix(0, years, n_sample)
prev_f_20_24_max <- matrix(0, years, n_sample)
prev_f_20_24_min <- matrix(0, years, n_sample)
prev_m_15_19_max <- matrix(0, years, n_sample)
prev_m_15_19_min <- matrix(0, years, n_sample)
prev_m_20_24_max <- matrix(0, years, n_sample)
prev_m_20_24_min <- matrix(0, years, n_sample)

### nested loop to compute the solutions of the differential equations
for (j in seq_len(years)){
  for (i in seq_len(n_sample)){
    # maximum estimates
    # f_15_19
    solution_f_15_19_max <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_f[i], spon_rec_f[i], att_symp[i], p_true_pos_f[i], p_false_pos_f[i]) - 
        c(test_rate_f_15_19_max[j,i], diag_rate_f_15_19_max[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_f_15_19_max[j,i] <- solution_f_15_19_max[1]
    inc_f_15_19_max[j,i] <- foi_f_15_19_max[j,i]
    scr_f_15_19_max[j,i] <- solution_f_15_19_max[2]
    prev_f_15_19_max[j,i] <- sol_dyn$prev_fun(inc_f_15_19_max[j,i]*p_asymp_f[i], 
                                              spon_rec_f[i] + scr_f_15_19_max[j,i]*p_true_pos_f[i], 
                                              inc_f_15_19_max[j,i]*(1-p_asymp_f[i]), 
                                              scr_f_15_19_max[j,i]*p_true_pos_f[i] + att_symp[i]*p_true_pos_f[i])
    # f_20_24
    solution_f_20_24_max <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_f[i], spon_rec_f[i], att_symp[i], p_true_pos_f[i], p_false_pos_f[i]) - 
        c(test_rate_f_20_24_max[j,i], diag_rate_f_20_24_max[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_f_20_24_max[j,i] <- solution_f_20_24_max[1]
    inc_f_20_24_max[j,i] <- foi_f_20_24_max[j,i]
    scr_f_20_24_max[j,i] <- solution_f_20_24_max[2]
    prev_f_20_24_max[j,i] <- sol_dyn$prev_fun(inc_f_20_24_max[j,i]*p_asymp_f[i], 
                                              spon_rec_f[i] + scr_f_20_24_max[j,i]*p_true_pos_f[i], 
                                              inc_f_20_24_max[j,i]*(1-p_asymp_f[i]), 
                                              scr_f_20_24_max[j,i]*p_true_pos_f[i] + att_symp[i]*p_true_pos_f[i])
    
    # m_15_19
    solution_m_15_19_max <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_m[i], spon_rec_m[i], att_symp[i], p_true_pos_m[i], p_false_pos_m[i]) - 
        c(test_rate_m_15_19_max[j,i], diag_rate_m_15_19_max[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_m_15_19_max[j,i] <- solution_m_15_19_max[1]
    inc_m_15_19_max[j,i] <- foi_m_15_19_max[j,i]
    scr_m_15_19_max[j,i] <- solution_m_15_19_max[2]
    prev_m_15_19_max[j,i] <- sol_dyn$prev_fun(inc_m_15_19_max[j,i]*p_asymp_m[i], 
                                              spon_rec_m[i] + scr_m_15_19_max[j,i]*p_true_pos_m[i], 
                                              inc_m_15_19_max[j,i]*(1-p_asymp_m[i]), 
                                              scr_m_15_19_max[j,i]*p_true_pos_m[i] + att_symp[i]*p_true_pos_m[i])
    # m_20_24
    solution_m_20_24_max <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_m[i], spon_rec_m[i], att_symp[i], p_true_pos_m[i], p_false_pos_m[i]) - 
        c(test_rate_m_20_24_max[j,i], diag_rate_m_20_24_max[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_m_20_24_max[j,i] <- solution_m_20_24_max[1]
    inc_m_20_24_max[j,i] <- foi_m_20_24_max[j,i]
    scr_m_20_24_max[j,i] <- solution_m_20_24_max[2]
    prev_m_20_24_max[j,i] <- sol_dyn$prev_fun(inc_m_20_24_max[j,i]*p_asymp_m[i], 
                                              spon_rec_m[i] + scr_m_20_24_max[j,i]*p_true_pos_m[i], 
                                              inc_m_20_24_max[j,i]*(1-p_asymp_m[i]), 
                                              scr_m_20_24_max[j,i]*p_true_pos_m[i] + att_symp[i]*p_true_pos_m[i])
    # minimum estimates
    # f_15_19
    solution_f_15_19_min <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_f[i], spon_rec_f[i], att_symp[i], p_true_pos_f[i], p_false_pos_f[i]) - 
        c(test_rate_f_15_19_min[j,i], diag_rate_f_15_19_min[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_f_15_19_min[j,i] <- solution_f_15_19_min[1]
    inc_f_15_19_min[j,i] <- foi_f_15_19_min[j,i]
    scr_f_15_19_min[j,i] <- solution_f_15_19_min[2]
    prev_f_15_19_min[j,i] <- sol_dyn$prev_fun(inc_f_15_19_min[j,i]*p_asymp_f[i], 
                                              spon_rec_f[i] + scr_f_15_19_min[j,i]*p_true_pos_f[i], 
                                              inc_f_15_19_min[j,i]*(1-p_asymp_f[i]), 
                                              scr_f_15_19_min[j,i]*p_true_pos_f[i] + att_symp[i]*p_true_pos_f[i])
    # f_20_24
    solution_f_20_24_min <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_f[i], spon_rec_f[i], att_symp[i], p_true_pos_f[i], p_false_pos_f[i]) - 
        c(test_rate_f_20_24_min[j,i], diag_rate_f_20_24_min[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_f_20_24_min[j,i] <- solution_f_20_24_min[1]
    inc_f_20_24_min[j,i] <- foi_f_20_24_min[j,i]
    scr_f_20_24_min[j,i] <- solution_f_20_24_min[2]
    prev_f_20_24_min[j,i] <- sol_dyn$prev_fun(inc_f_20_24_min[j,i]*p_asymp_f[i], 
                                              spon_rec_f[i] + scr_f_20_24_min[j,i]*p_true_pos_f[i], 
                                              inc_f_20_24_min[j,i]*(1-p_asymp_f[i]), 
                                              scr_f_20_24_min[j,i]*p_true_pos_f[i] + att_symp[i]*p_true_pos_f[i])
    # m_15_19
    solution_m_15_19_min <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_m[i], spon_rec_m[i], att_symp[i], p_true_pos_m[i], p_false_pos_m[i]) - 
        c(test_rate_m_15_19_min[j,i], diag_rate_m_15_19_min[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_m_15_19_min[j,i] <- solution_m_15_19_min[1]
    inc_m_15_19_min[j,i] <- foi_m_15_19_min[j,i]
    scr_m_15_19_min[j,i] <- solution_m_15_19_min[2]
    prev_m_15_19_min[j,i] <- sol_dyn$prev_fun(inc_m_15_19_min[j,i]*p_asymp_m[i], 
                                              spon_rec_m[i] + scr_m_15_19_min[j,i]*p_true_pos_m[i], 
                                              inc_m_15_19_min[j,i]*(1-p_asymp_m[i]), 
                                              scr_m_15_19_min[j,i]*p_true_pos_m[i] + att_symp[i]*p_true_pos_m[i])
    # m_20_24
    solution_m_20_24_min <- multiroot(function(parms) 
      test_diag_fun(parms, 1-p_asymp_m[i], spon_rec_m[i], att_symp[i], p_true_pos_m[i], p_false_pos_m[i]) - 
        c(test_rate_m_20_24_min[j,i], diag_rate_m_20_24_min[j,i]), start = c(0.02, 0.44), positive=TRUE)$root
    foi_m_20_24_min[j,i] <- solution_m_20_24_min[1]
    inc_m_20_24_min[j,i] <- foi_m_20_24_min[j,i]
    scr_m_20_24_min[j,i] <- solution_m_20_24_min[2]
    prev_m_20_24_min[j,i] <- sol_dyn$prev_fun(inc_m_20_24_min[j,i]*p_asymp_m[i], 
                                              spon_rec_m[i] + scr_m_20_24_min[j,i]*p_true_pos_m[i], 
                                              inc_m_20_24_min[j,i]*(1-p_asymp_m[i]), 
                                              scr_m_20_24_min[j,i]*p_true_pos_m[i] + att_symp[i]*p_true_pos_m[i])
    # print to console progress of the loop
    print_message <- paste('year',j,'sample',i) 
    print(print_message, quote = FALSE)
  }
}

### Convert matrices into data frame and save to /output
# incidence data frames:
inc_f_15_19_max <- as.data.frame(inc_f_15_19_max*pop_f_15_19)
inc_f_15_19_min <- as.data.frame(inc_f_15_19_min*pop_f_15_19)
inc_f_20_24_max <- as.data.frame(inc_f_20_24_max*pop_f_20_24)
inc_f_20_24_min <- as.data.frame(inc_f_20_24_min*pop_f_20_24)
inc_m_15_19_max <- as.data.frame(inc_m_15_19_max*pop_m_15_19)
inc_m_15_19_min <- as.data.frame(inc_m_15_19_min*pop_m_15_19)
inc_m_20_24_max <- as.data.frame(inc_m_20_24_max*pop_m_20_24)
inc_m_20_24_min <- as.data.frame(inc_m_20_24_min*pop_m_20_24)
rownames(inc_f_15_19_max) <- 2000:2015
rownames(inc_f_15_19_min) <- 2000:2015
rownames(inc_f_20_24_max) <- 2000:2015
rownames(inc_f_20_24_min) <- 2000:2015
rownames(inc_m_15_19_max) <- 2000:2015
rownames(inc_m_15_19_min) <- 2000:2015
rownames(inc_m_20_24_max) <- 2000:2015
rownames(inc_m_20_24_min) <- 2000:2015
  
  
## save
write.csv(inc_f_15_19_max,'output/eng_lw_model_priors_inc_f_15_19_max.csv')
write.csv(inc_f_15_19_min,'output/eng_lw_model_priors_inc_f_15_19_min.csv')
write.csv(inc_f_20_24_max,'output/eng_lw_model_priors_inc_f_20_24_max.csv')
write.csv(inc_f_20_24_min,'output/eng_lw_model_priors_inc_f_20_24_min.csv')
write.csv(inc_m_15_19_max,'output/eng_lw_model_priors_inc_m_15_19_max.csv')
write.csv(inc_m_15_19_min,'output/eng_lw_model_priors_inc_m_15_19_min.csv')
write.csv(inc_m_20_24_max,'output/eng_lw_model_priors_inc_m_20_24_max.csv')
write.csv(inc_m_20_24_min,'output/eng_lw_model_priors_inc_m_20_24_min.csv')

# force of infection data frames
foi_f_15_19_max <- as.data.frame(foi_f_15_19_max)
foi_f_15_19_min <- as.data.frame(foi_f_15_19_min)
foi_f_20_24_max <- as.data.frame(foi_f_20_24_max)
foi_f_20_24_min <- as.data.frame(foi_f_20_24_min)
foi_m_15_19_max <- as.data.frame(foi_m_15_19_max)
foi_m_15_19_min <- as.data.frame(foi_m_15_19_min)
foi_m_20_24_max <- as.data.frame(foi_m_20_24_max)
foi_m_20_24_min <- as.data.frame(foi_m_20_24_min)

rownames(foi_f_15_19_max) <- 2000:2015
rownames(foi_f_15_19_min) <- 2000:2015
rownames(foi_f_20_24_max) <- 2000:2015
rownames(foi_f_20_24_min) <- 2000:2015
rownames(foi_m_15_19_max) <- 2000:2015
rownames(foi_m_15_19_min) <- 2000:2015
rownames(foi_m_20_24_max) <- 2000:2015
rownames(foi_m_20_24_min) <- 2000:2015

## save
write.csv(foi_f_15_19_max,'output/eng_lw_model_priors_foi_f_15_19_max.csv')
write.csv(foi_f_15_19_min,'output/eng_lw_model_priors_foi_f_15_19_min.csv')
write.csv(foi_f_20_24_max,'output/eng_lw_model_priors_foi_f_20_24_max.csv')
write.csv(foi_f_20_24_min,'output/eng_lw_model_priors_foi_f_20_24_min.csv')
write.csv(foi_m_15_19_max,'output/eng_lw_model_priors_foi_m_15_19_max.csv')
write.csv(foi_m_15_19_min,'output/eng_lw_model_priors_foi_m_15_19_min.csv')
write.csv(foi_m_20_24_max,'output/eng_lw_model_priors_foi_m_20_24_max.csv')
write.csv(foi_m_20_24_min,'output/eng_lw_model_priors_foi_m_20_24_min.csv')

# screening Rate data frames:
scr_f_15_19_max <- as.data.frame(scr_f_15_19_max*100)
scr_f_15_19_min <- as.data.frame(scr_f_15_19_min*100)
scr_f_20_24_max <- as.data.frame(scr_f_20_24_max*100)
scr_f_20_24_min <- as.data.frame(scr_f_20_24_min*100)
scr_m_15_19_max <- as.data.frame(scr_m_15_19_max*100)
scr_m_15_19_min <- as.data.frame(scr_m_15_19_min*100)
scr_m_20_24_max <- as.data.frame(scr_m_20_24_max*100)
scr_m_20_24_min <- as.data.frame(scr_m_20_24_min*100)

rownames(scr_f_15_19_max) <- 2000:2015
rownames(scr_f_15_19_min) <- 2000:2015
rownames(scr_f_20_24_max) <- 2000:2015
rownames(scr_f_20_24_min) <- 2000:2015
rownames(scr_m_15_19_max) <- 2000:2015
rownames(scr_m_15_19_min) <- 2000:2015
rownames(scr_m_20_24_max) <- 2000:2015
rownames(scr_m_20_24_min) <- 2000:2015

write.csv(scr_f_15_19_max,'output/eng_lw_model_priors_scr_f_15_19_max.csv')
write.csv(scr_f_15_19_min,'output/eng_lw_model_priors_scr_f_15_19_min.csv')
write.csv(scr_f_20_24_max,'output/eng_lw_model_priors_scr_f_20_24_max.csv')
write.csv(scr_f_20_24_min,'output/eng_lw_model_priors_scr_f_20_24_min.csv')
write.csv(scr_m_15_19_max,'output/eng_lw_model_priors_scr_m_15_19_max.csv')
write.csv(scr_m_15_19_min,'output/eng_lw_model_priors_scr_m_15_19_min.csv')
write.csv(scr_m_20_24_max,'output/eng_lw_model_priors_scr_m_20_24_max.csv')
write.csv(scr_m_20_24_min,'output/eng_lw_model_priors_scr_m_20_24_min.csv')

# prevalence data frames:
prev_f_15_19_max <- as.data.frame(prev_f_15_19_max*100)
prev_f_15_19_min <- as.data.frame(prev_f_15_19_min*100)
prev_f_20_24_max <- as.data.frame(prev_f_20_24_max*100)
prev_f_20_24_min <- as.data.frame(prev_f_20_24_min*100)
prev_m_15_19_max <- as.data.frame(prev_m_15_19_max*100)
prev_m_15_19_min <- as.data.frame(prev_m_15_19_min*100)
prev_m_20_24_max <- as.data.frame(prev_m_20_24_max*100)
prev_m_20_24_min <- as.data.frame(prev_m_20_24_min*100)

rownames(prev_f_15_19_max) <- 2000:2015
rownames(prev_f_15_19_min) <- 2000:2015
rownames(prev_f_20_24_max) <- 2000:2015
rownames(prev_f_20_24_min) <- 2000:2015
rownames(prev_m_15_19_max) <- 2000:2015
rownames(prev_m_15_19_min) <- 2000:2015
rownames(prev_m_20_24_max) <- 2000:2015
rownames(prev_m_20_24_min) <- 2000:2015

write.csv(prev_f_15_19_max,'output/eng_lw_model_priors_prev_f_15_19_max.csv')
write.csv(prev_f_15_19_min,'output/eng_lw_model_priors_prev_f_15_19_min.csv')
write.csv(prev_f_20_24_max,'output/eng_lw_model_priors_prev_f_20_24_max.csv')
write.csv(prev_f_20_24_min,'output/eng_lw_model_priors_prev_f_20_24_min.csv')
write.csv(prev_m_15_19_max,'output/eng_lw_model_priors_prev_m_15_19_max.csv')
write.csv(prev_m_15_19_min,'output/eng_lw_model_priors_prev_m_15_19_min.csv')
write.csv(prev_m_20_24_max,'output/eng_lw_model_priors_prev_m_20_24_max.csv')
write.csv(prev_m_20_24_min,'output/eng_lw_model_priors_prev_m_20_24_min.csv')

### Function to give a summary of the data in the dataframes for easier plotting
summary_stats <- function(dataframe){
  years <- nrow(dataframe)
  summary_dataframe <- matrix(0,years,5)
  for (j in 1:years){
    summary_dataframe[j,1] <- stats::median(as.numeric(dataframe[j,]),na.rm=TRUE)
    quant95 <- stats::quantile(dataframe[j,],c(0.025, 0.975),na.rm=TRUE) # 95% creible interval
    summary_dataframe[j,2] <- quant95[1]
    summary_dataframe[j,3] <- quant95[2]
    quant50 <- stats::quantile(dataframe[j,],c(0.25, 0.75),na.rm=TRUE)
    summary_dataframe[j,4] <- quant50[1]
    summary_dataframe[j,5] <- quant50[2]
  }
  summary_dataframe <- as.data.frame(summary_dataframe)
  rownames(summary_dataframe) <- 2000:2015
  colnames(summary_dataframe) <- c('median','2.5th percentile','97.5th percentile','25th percentile','75th percentile')
  return(summary_dataframe)
}

### Create summary statistics from the data
# incidence
inc_f_15_19_max_stats <- summary_stats(inc_f_15_19_max)
inc_f_15_19_min_stats <- summary_stats(inc_f_15_19_min)
inc_f_20_24_max_stats <- summary_stats(inc_f_20_24_max)
inc_f_20_24_min_stats <- summary_stats(inc_f_20_24_min)
inc_m_15_19_max_stats <- summary_stats(inc_m_15_19_max)
inc_m_15_19_min_stats <- summary_stats(inc_m_15_19_min)
inc_m_20_24_max_stats <- summary_stats(inc_m_20_24_max)
inc_m_20_24_min_stats <- summary_stats(inc_m_20_24_min)

write.csv(inc_f_15_19_max_stats,'output/eng_lw_model_priors_inc_f_15_19_max_stats.csv')
write.csv(inc_f_15_19_min_stats,'output/eng_lw_model_priors_inc_f_15_19_min_stats.csv')
write.csv(inc_f_20_24_max_stats,'output/eng_lw_model_priors_inc_f_20_24_max_stats.csv')
write.csv(inc_f_20_24_min_stats,'output/eng_lw_model_priors_inc_f_20_24_min_stats.csv')
write.csv(inc_m_15_19_max_stats,'output/eng_lw_model_priors_inc_m_15_19_max_stats.csv')
write.csv(inc_m_15_19_min_stats,'output/eng_lw_model_priors_inc_m_15_19_min_stats.csv')
write.csv(inc_m_20_24_max_stats,'output/eng_lw_model_priors_inc_m_20_24_max_stats.csv')
write.csv(inc_m_20_24_min_stats,'output/eng_lw_model_priors_inc_m_20_24_min_stats.csv')

# force of infection
foi_f_15_19_max_stats <- summary_stats(foi_f_15_19_max)
foi_f_15_19_min_stats <- summary_stats(foi_f_15_19_min)
foi_f_20_24_max_stats <- summary_stats(foi_f_20_24_max)
foi_f_20_24_min_stats <- summary_stats(foi_f_20_24_min)
foi_m_15_19_max_stats <- summary_stats(foi_m_15_19_max)
foi_m_15_19_min_stats <- summary_stats(foi_m_15_19_min)
foi_m_20_24_max_stats <- summary_stats(foi_m_20_24_max)
foi_m_20_24_min_stats <- summary_stats(foi_m_20_24_min)

write.csv(foi_f_15_19_max_stats,'output/eng_lw_model_priors_foi_f_15_19_max_stats.csv')
write.csv(foi_f_15_19_min_stats,'output/eng_lw_model_priors_foi_f_15_19_min_stats.csv')
write.csv(foi_f_20_24_max_stats,'output/eng_lw_model_priors_foi_f_20_24_max_stats.csv')
write.csv(foi_f_20_24_min_stats,'output/eng_lw_model_priors_foi_f_20_24_min_stats.csv')
write.csv(foi_m_15_19_max_stats,'output/eng_lw_model_priors_foi_m_15_19_max_stats.csv')
write.csv(foi_m_15_19_min_stats,'output/eng_lw_model_priors_foi_m_15_19_min_stats.csv')
write.csv(foi_m_20_24_max_stats,'output/eng_lw_model_priors_foi_m_20_24_max_stats.csv')
write.csv(foi_m_20_24_min_stats,'output/eng_lw_model_priors_foi_m_20_24_min_stats.csv')

# screening Rate
scr_f_15_19_max_stats <- summary_stats(scr_f_15_19_max)
scr_f_15_19_min_stats <- summary_stats(scr_f_15_19_min)
scr_f_20_24_max_stats <- summary_stats(scr_f_20_24_max)
scr_f_20_24_min_stats <- summary_stats(scr_f_20_24_min)
scr_m_15_19_max_stats <- summary_stats(scr_m_15_19_max)
scr_m_15_19_min_stats <- summary_stats(scr_m_15_19_min)
scr_m_20_24_max_stats <- summary_stats(scr_m_20_24_max)
scr_m_20_24_min_stats <- summary_stats(scr_m_20_24_min)

write.csv(scr_f_15_19_max_stats,'output/eng_lw_model_priors_scr_f_15_19_max_stats.csv')
write.csv(scr_f_15_19_min_stats,'output/eng_lw_model_priors_scr_f_15_19_min_stats.csv')
write.csv(scr_f_20_24_max_stats,'output/eng_lw_model_priors_scr_f_20_24_max_stats.csv')
write.csv(scr_f_20_24_min_stats,'output/eng_lw_model_priors_scr_f_20_24_min_stats.csv')
write.csv(scr_m_15_19_max_stats,'output/eng_lw_model_priors_scr_m_15_19_max_stats.csv')
write.csv(scr_m_15_19_min_stats,'output/eng_lw_model_priors_scr_m_15_19_min_stats.csv')
write.csv(scr_m_20_24_max_stats,'output/eng_lw_model_priors_scr_m_20_24_max_stats.csv')
write.csv(scr_m_20_24_min_stats,'output/eng_lw_model_priors_scr_m_20_24_min_stats.csv')

# prevalence
prev_f_15_19_max_stats <- summary_stats(prev_f_15_19_max)
prev_f_15_19_min_stats <- summary_stats(prev_f_15_19_min)
prev_f_20_24_max_stats <- summary_stats(prev_f_20_24_max)
prev_f_20_24_min_stats <- summary_stats(prev_f_20_24_min)
prev_m_15_19_max_stats <- summary_stats(prev_m_15_19_max)
prev_m_15_19_min_stats <- summary_stats(prev_m_15_19_min)
prev_m_20_24_max_stats <- summary_stats(prev_m_20_24_max)
prev_m_20_24_min_stats <- summary_stats(prev_m_20_24_min)

write.csv(prev_f_15_19_max_stats,'output/eng_lw_model_priors_prev_f_15_19_max_stats.csv')
write.csv(prev_f_15_19_min_stats,'output/eng_lw_model_priors_prev_f_15_19_min_stats.csv')
write.csv(prev_f_20_24_max_stats,'output/eng_lw_model_priors_prev_f_20_24_max_stats.csv')
write.csv(prev_f_20_24_min_stats,'output/eng_lw_model_priors_prev_f_20_24_min_stats.csv')
write.csv(prev_m_15_19_max_stats,'output/eng_lw_model_priors_prev_m_15_19_max_stats.csv')
write.csv(prev_m_15_19_min_stats,'output/eng_lw_model_priors_prev_m_15_19_min_stats.csv')
write.csv(prev_m_20_24_max_stats,'output/eng_lw_model_priors_prev_m_20_24_max_stats.csv')
write.csv(prev_m_20_24_min_stats,'output/eng_lw_model_priors_prev_m_20_24_min_stats.csv')

### Stop timer and print runtime
end_time <- Sys.time()
run_time <- end_time - start_time
print_message <- paste('program runtime:', run_time,'units') 
print(print_message, quote = FALSE)
