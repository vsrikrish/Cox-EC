# load libraries
library(DEoptim)
library(adaptMCMC)
library(parallel)

args <- commandArgs(trailingOnly=TRUE)
winlen <- as.numeric(args[1])
aid <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# read data
dat <- read.csv(paste('data/ECS-', winlen, '.csv', sep=''), header=T, row.names=1)

# get number of unique models
nmod <- 16

# odd array ids use log responses
if (aid %% 2 != 0) {
  dat[2:nrow(dat), c('ECS', 'Psi')] <- log(dat[2:nrow(dat), c('ECS', 'Psi')])
  log.name <- 'log-'
} else {
  log.name <- ''
}
# aids 0 and 1 use normal distributions, 2 uses t
if (aid < 2) {
  FUN.density <- dnorm
  density.name <- 'normal'
} else {
  FUN.density <- dt
  density.name <- 't'
}

parnames <- c('a', 'b', 'scale')

source('scripts/likelihood.R')

# compute maximum likelihood estimates over held out data to use for the prior
de_out <- DEoptim(neg.log.lik, lower=c(-5, -20, 0), upper=c(5, 20, 2), control=list(NP=200, itermax=1000, trace=F), parnames=parnames, dat=dat[(nmod+2):nrow(dat),], density.name=density.name)$optim$bestmem
mle_prior <- de_out

# compute maximum likelihood estimates over non-held out data as initial value
de_out <- DEoptim(neg.log.lik, lower=c(-5, -20, 0), upper=c(5, 20, 2), control=list(NP=200, itermax=1000, trace=F), parnames=parnames, dat=dat[2:(nmod+1),], density.name=density.name)$optim$bestmem

init.val <- de_out

# run MCMC chains in parallel
n.chains <- 4
n.iter <- 1000000
n.cpu <- parallel::detectCores()-1
accept_rate_many <- 0.234
accept_rate_few <- 0.44
accept_rate <- accept_rate_many + (accept_rate_few-accept_rate_many)/length(parnames)

amcmc_out <- MCMC.parallel(log.post, n.iter, init.val, n.chain=n.chains, n.cpu=n.cpu, list=T, adapt=T, acc.rate=accept_rate, gamma=0.67, n.start=max(500, round(0.05)*n.iter), parnames=parnames, dat=dat[2:(nmod+1),], density.name=density.name, mle=mle_prior)
amcmc_out[['prior']] <- mle_prior
amcmc_out[['MLE']] <- init.val
amcmc_out[['parnames']] <- parnames
amcmc_out[['dat']] <- dat
amcmc_out[['density.name']] <- density.name
saveRDS(amcmc_out, paste('output/MCMC/mcmc_out-', winlen, '-', log.name, density.name, '.rds', sep=''))s', sep=''))