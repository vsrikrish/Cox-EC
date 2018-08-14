###############################################################################################################
# mcmc_driver.R                                                                                               #
#                                                                                                             #
# This script generates MCMC chains for the Cox et al (Nature, 2018) emergent constraint model.               #
# Accepted command-line arguments are:                                                                        #
#   1) window length (odd integer from 45 to 73 years) for generating the temperature statistic               #
#   2) model type ("normal", "log-normal", "t"). Note that the model type can be controlled using a PBS       #
#       job array, where PBS Array ID 0 = normal, 1 = log-normal, and 2 = t.                                  #
# This script also assumes that the job will be run in parallel, so that multiple chains for the Gelman-Rubin #
#   statistic will be generated. Modifications will be needed to run it in serial.                            #
###############################################################################################################

# load libraries
library(DEoptim)
library(adaptMCMC)
library(parallel)

args <- commandArgs(trailingOnly=TRUE) # read in command line args
winlen <- as.numeric(args[1]) # first argument is the window length (from 45 to 73 years, in odd numbers)

aid <- Sys.getenv('PBS_ARRAYID', unset=NA)  # check if script was called from within a PBS array
# otherwise read in the second command line arg as the name of the density
# needs to be one of "normal", "log-normal", or "t"
if (is.na(aid)) {
  density.name <- args[2]
}

# set script, output, and data paths
R.path <- 'scripts'
out.path <- 'output'
data.path <- 'data'

# read data
dat <- read.csv(file.path(data.path, paste('ECS-', winlen, '.csv', sep='')), header=T, row.names=1)

# get number of unique models
nmod <- 16



source(paste(R.path, 'likelihood.R')) # source likelihood and prior functions

# if aid is not NA, use it to determine the density name and whether there's a log response or not
if (!is.na(aid) {
  aid <- as.numeric(aid)
  # odd array ids use log responses
  if (aid %% 2 != 0) {
    dat[2:nrow(dat), c('ECS', 'Psi')] <- log(dat[2:nrow(dat), c('ECS', 'Psi')])
    log.name <- 'log-'
  } else {
    log.name <- ''
  }
  # aids 0 and 1 use normal distributions, 2 uses t
  if (aid < 2) {
    density.name <- 'normal'
  } else {
    density.name <- 't'
  }
# otherwise parse the density name to determine if a log response is needed
} else {
  if (grepl('log', density.name)) {
    dat[2:nrow(dat), c('ECS', 'Psi')] <- log(dat[2:nrow(dat), c('ECS', 'Psi')])
    log.name <- 'log-'
    density.name <- 'normal'
  } else {
    log.name <- ''
  }
}

parnames <- c('a', 'b', 'scale') # set parameter names


# compute maximum likelihood estimates over held out data to use for the prior
de_out <- DEoptim(neg.log.lik, lower=c(-5, -20, 0), upper=c(5, 20, 2), control=list(NP=200, itermax=1000, trace=F), parnames=parnames, dat=dat[(nmod+2):nrow(dat),], density.name=paste(log.name, density.name))$optim$bestmem
mle_prior <- de_out

# compute maximum likelihood estimates over non-held out data as initial value
de_out <- DEoptim(neg.log.lik, lower=c(-5, -20, 0), upper=c(5, 20, 2), control=list(NP=200, itermax=1000, trace=F), parnames=parnames, dat=dat[2:(nmod+1),], density.name=paste(log.name, density.name))$optim$bestmem

init.val <- de_out

# run MCMC chains in parallel
n.chains <- 4
n.iter <- 1000000
n.cpu <- parallel::detectCores()-1
accept_rate_many <- 0.234
accept_rate_few <- 0.44
accept_rate <- accept_rate_many + (accept_rate_few-accept_rate_many)/length(parnames)

amcmc_out <- MCMC.parallel(log.post, n.iter, init.val, n.chain=n.chains, n.cpu=n.cpu, list=T, adapt=T, acc.rate=accept_rate, gamma=0.67, n.start=max(500, round(0.05)*n.iter), parnames=parnames, dat=dat[2:(nmod+1),], density.name=paste(log.name, density.name), mle=mle_prior)
amcmc_out[['prior']] <- mle_prior
amcmc_out[['MLE']] <- init.val
amcmc_out[['parnames']] <- parnames
amcmc_out[['dat']] <- dat
amcmc_out[['density.name']] <- density.name
saveRDS(amcmc_out, file.path(out.path, 'MCMC', paste('mcmc_out-', winlen, '-', log.name, density.name, '.rds', sep='')))