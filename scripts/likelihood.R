###############################################################################################################
# likelihood.R                                                                                                #
#                                                                                                             #
# This script contains the log-likelihood, log-prior, and log-posterior functions.                            #
###############################################################################################################


# log-likelihood function
# inputs:
#   1) pars: vector of parameters
#   2) parnames: vector of parameter names (should be c('a', 'b', 'scale')).
#   3) dat: data frame with columns 'ECS' (response) and 'Psi' (covariate). These
#           should be log-transformed if needed.
#   4) density.name: 'normal,' 'log-normal', or 't' (used to select density function and truncate if necessary)
# returns:
#   1) log-likelihood value
log.lik <- function(pars, parnames, dat, density.name) {
  # extract variables
  a <- pars[match('a', parnames)]
  b <- pars[match('b', parnames)]
  scale <- pars[match('scale', parnames)]
  
  # compute expected value and residuals
  e <- a+b*dat$Psi
  r <- dat$ECS - e
  # return log likelihood of residuals for appropriate density function
  if (density.name == 't') {
    ll <- sum(dt(r, scale, log=T) - pt(-e, scale, log=T, lower.tail=F))
  } else if (density.name == 'normal') {
    ll <- sum(dnorm(r, 0, scale, log=T) - pnorm(-e, 0, scale, log=T, lower.tail=F))
  } else {
    ll <- sum(dnorm(r, 0, scale, log=T))
  }
  ll
}

# log-prior density function
# prior for location parameters ('a': intercept, 'b': linear coefficient) are wide normal priors (sd=10) around
#   the MLE estimate for held-out data.
# prior for scale parameters is gamma(1,2).
# inputs:
#   1) pars: vector of parameters
#   2) parnames: vector of parameter names (should be c('a', 'b', 'scale')
#   3) mle: MLE estimate for held-out data.
# returns:
#   1) log-prior across all parameters
log.pri <- function(pars, parnames, mle) {
  lp <- 0
  for (i in 1:length(parnames)) {
    # scale parameter is given a gamma(1,2) prior
    if (parnames[i] == 'scale') {
      # if proposed value is non-positive, prior density is 0
      if (pars[i] <= 0) {
        lp <- -Inf
      } else {
        lp <- lp + dgamma(pars[i], shape=1, scale=2, log=T)
      }
    } else {
      # otherwise, use a wide normal centered at the mle estimate
      lp <- lp + dnorm(pars[i], mle[i], 10, log=T)
    }
  }
  # return log prior
  lp
}

# log-posterior density function
# inputs:
#   1) pars: vector of parameter values
#   2) parnames: vectof of parameter names (should be c('a', 'b', 'scale'))
#   3) dat: data frame with columns 'ECS' (response) and 'Psi' (covariate). These
#           should be log-transformed if needed.
#   4) mle: MLE estimate for held-out data.
#   5) density.name: 'normal,' 'log-normal', or 't' (used to select density function and truncate if necessary)
log.post <- function(pars, parnames, dat, mle, density.name) {
  # compute log-prior density
  l.pri <- log.pri(pars, parnames, mle)
  # if log-prior density is infinite, don't bother computing likelihood.
  if (is.finite(l.pri)) {
    l.lik <- log.lik(pars, parnames, dat, density.name)
  } else {
    l.lik <- -Inf
  }
  # return non-normalized log posterior density
  l.pri + l.lik
}

# negative log likelihood for DE Optim
neg.log.lik <- function(pars, parnames, dat, density.name) {
  -1*log.lik(pars, parnames, dat, density.name)
}

