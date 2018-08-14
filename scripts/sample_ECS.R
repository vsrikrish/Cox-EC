###############################################################################################################
# sample_ECS.R                                                                                                #
#                                                                                                             #
# This script uses MCMC output files to sample from the ECS posterior predictive distribution.                #
# Accepted command-line arguments are:                                                                        #
#   1) window length (odd integer from 45 to 73 years) for generating the temperature statistic               #
# This script then loops through each of the three model types for that window length.                        #
# Saved output:                                                                                               #
#   1) ecs_samp-winlen.rds, storing a named list (by model) of n.samp posterior predictive samples            #
#   2) mcmc_samp-winlen.rds, storing a named list (by model) of n.samp posterior samples.                     #
# We use the Gelman-Rubin statistic to check for MCMC convergence and burn-in chains.                         #
###############################################################################################################


library(coda)
library(reshape2)

# gelman-rubin diagnostic function
gr_diag <- function(chains) {
  if (is.list(chains)) {
    iter_length <- nrow(chains[[1]]$samples)
    test_seq <- seq(1, iter_length, by = 0.05*iter_length)
    sapply(test_seq, function(m) as.numeric(gelman.diag(lapply(chains, function(l) as.mcmc(l$samples[1:m,])))[[2]]))
  } else {
    iter_length <- nrow(chains$samples)
    test_seq <- seq(1, iter_length, by = 0.05*iter_length)
    sapply(test_seq, function(m) as.numeric(gelman.diag(as.mcmc(chains$samples[1:m,]))[[2]]))
  }
}

# read command line argument
args <- commandArgs(trailingOnly=TRUE)
winlen <- as.numeric(args[1])

# set paths
R.path <- 'scripts'
out.path <- 'output'
data.path <- 'data'

# number of unique models
nmod <- 16

n.samp <- 1000000

source(file.path(R.path, 'likelihood.R'))

# load data
dat <- read.csv(file.path(data.path, paste('ECS-', winlen, '.csv', sep='')), header=T, row.names=1)
obs <- dat[1,c('Psi', 'dPsi')]

mcmc_names <- c('normal', 'log-normal', 't')  # set model names
ecs.samp <- vector('list', length(mcmc_names))  # initialize storage for ECS samples
names(ecs.samp) <- mcmc_names

mcmc_samp <- vector('list', length(mcmc_names)) # initialize storage for posterior samples

for (j in 1:length(mcmc_names)) {
  # read MCMC output
  mcmc_out <- readRDS(paste('output/MCMC/mcmc_out-', winlen, '-', mcmc_names[j], '.rds', sep=''))
  gr <- c(gr_diag(mcmc_out[1:4]))   # calculate Gelman-Rubin diagnostic
  bi <- min(which(gr < 1.1))  # find where G-R statistic is below 1.1 (to estimate convergence)
  # compute number of burn-in MCMC samples, prune each chain, and combine chains
  burn.in <- bi*0.05*nrow(mcmc_out[[1]]$samples)
  # sample from combined posterior samples
  samp <- do.call('rbind', lapply(mcmc_out[1:4], function(l) l$samples[(burn.in+1):nrow(l$samples),]))
  idx <- sample(1:nrow(samp), n.samp, replace=T)
  mcmc_samp[[j]] <- samp[idx,]
  
  # identify if log transforms are needed for posterior predictive sampling
  if (grepl('log', mcmc_names[j])) {
    log.response <- 'log'
  } else {
    log.response <- 'identity'
  }
  
  # MC sample prediction distribution
  # sample Psi from observation range
  Psi.samp <- runif(n.samp, as.numeric(obs['Psi']-obs['dPsi']), as.numeric(obs['Psi']+obs['dPsi']))
  
  # if model does not require a log-transformation, compute ECS samples, rejection-sampling values below 0
  if (log.response == 'identity') {
    ev.samp <- mcmc_samp[[j]][, 1] + mcmc_samp[[j]][, 2]*Psi.samp
    n <- 0
    pos <- logical(n.samp)
    r.samp <- numeric(n.samp)
    while (n < n.samp) {
      if (mcmc_out[['density.name']] == 'normal') {
        r.samp[!pos] <- rnorm(n.samp - n, 0, mcmc_samp[[j]][, 3])
      } else {
        r.samp[!pos] <- rt(n.samp - n, mcmc_samp[[j]][, 3])
      }
      pos[!pos] <- (r.samp[!pos] > -ev.samp[!pos])
      n <- sum(pos)
    }
    ecs.samp[[mcmc_names[j]]] <- ev.samp + r.samp
  # otherwise, log-transform covariate and exponentiate response
  } else {
    r.samp <- rnorm(n.samp, 0, mcmc_samp[[j]][, 3])
    ecs.samp[[mcmc_names[j]]] <- exp(mcmc_samp[[j]][, 1] + mcmc_samp[[j]][, 2]*log(Psi.samp) + r.samp)
  }
}


# save output files
saveRDS(mcmc_samp, file.path(out.path, 'samples', paste('mcmc_samp-', winlen, '.rds', sep='')))

saveRDS(ecs.samp, file.path(out.path, 'samples', paste('ecs_samp-', winlen, '.rds', sep='')))
