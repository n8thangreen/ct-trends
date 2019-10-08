
#' main MCMC MH iterations
#'
#' @param ninter Number of iterations
#' @param h_current Current (constant) hazard estimate
#' @param sd_h Standard deviation of hazard
#' @param beta_params List of beta parameters a, b
#'
#' @return List of all h_symp, loglik, props
#' @export
#'
runMCMC <- function(ninter,
                    h_current,
                    sd_h,
                    beta_params){
  
  ## initialise
  props  <- list()
  h_symp <- rep(NA_real_, niter) # testing rate per person per year
  loglik <- rep(NA_real_, niter)
  
  chain <- list()
  chain$h_symp <- h_current
  chain$props  <- prop_in_each_period(h_current)
  chain$loglik <- sum(beta_loglik(beta_params,
                                  chain$props))
  nperiods <- length(beta_params)
  
  for (i in seq_len(niter)) {
    
    h_star <- proposal(mean = chain$h_symp, sd_h)
    pos_h  <- (h_star > 0)
    
    if (pos_h)
      chain <- MHstep(beta_params, h_star, chain)
    
    h_symp[i]  <- chain$h_symp
    loglik[i]  <- chain$loglik
    props[[i]] <- chain$props
  }
  
  return(
    list(h_symp = h_symp,
         loglik = loglik,
         props  = props))
}


#' Metropolis-Hastings sampler step
#'
#' @param beta_params List of beta parameters a, b
#' @param h_star Proposal rate h
#' @param chain List of current h_symp, loglik, props
#'
#' @return List of current h_symp, loglik, props
#' @export
#'
MHstep <- function(beta_params,
                   h_star,
                   chain){
  
  props_new <- prop_in_each_period(h_star)
  all_pos <- all(props_new > 0)
  
  if (all_pos) {
    
    dbeta_new  <- beta_loglik(beta_params,
                              props_new)
    loglik_new <- sum(dbeta_new)
    log_ratio  <- loglik_new - chain$loglik
    
    accept <- log(runit()) < log_ratio
    
    if (accept) {
      chain$h_symp <- h_star
      chain$loglik <- loglik_new
      chain$props  <- props_new
    }
  }
  
  return(chain)
}


# likelihood function
beta_loglik <- function(beta_params,
                        props){
  
  a <- map(beta_params, "a")
  b <- map(beta_params, "b")
  out <- NULL
  
  for (i in seq_along(a)){
    
    out[i] <- dbeta(props[i],
                    shape1 = as.numeric(a[i]),
                    shape2 = as.numeric(b[i]),
                    log = TRUE)
  }
  return(out)
}


proposal <- function(mean, sd)
  rnorm(1, mean = mean, sd = sd)

runit <- function()
  runif(1, min = 0, max = 1)

prop_periods <- function(tps)
  function(haz) exp(-haz * tps[1:5]) - exp(-haz * tps[2:length(tps)])


