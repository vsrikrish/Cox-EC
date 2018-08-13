# define log-likelihood function
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

# define priors
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

# define posterior
log.post <- function(pars, parnames, dat, mle, density.name) {
  l.pri <- log.pri(pars, parnames, mle)
  if (is.finite(l.pri)) {
    l.lik <- log.lik(pars, parnames, dat, density.name)
  } else {
    l.lik <- -Inf
  }
  # return non-normalized log posterior density
  l.pri + l.lik
}

# define negative log likelihood for DE Optim
neg.log.lik <- function(pars, parnames, dat, density.name) {
  -1*log.lik(pars, parnames, dat, density.name)
}

