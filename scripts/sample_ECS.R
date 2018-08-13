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

args <- commandArgs(trailingOnly=TRUE)
winlen <- as.numeric(args[1])

R.path <- 'scripts'
out.path <- 'output'
# number of unique models
nmod <- 16

n.samp <- 1000000

source(file.path(R.path, 'likelihood.R'))

# load data
dat <- read.csv(paste('data/ECS-', winlen, '.csv', sep=''), header=T, row.names=1)
obs <- dat[1,c('Psi', 'dPsi')]

# load MCMC output
mcmc_names <- c('normal', 'log-normal', 't')
ecs.samp <- vector('list', length(mcmc_names))
names(ecs.samp) <- mcmc_names

mcmc_samp <- vector('list', length(mcmc_names))

for (j in 1:length(mcmc_names)) {
  mcmc_out <- readRDS(paste('output/MCMC/mcmc_out-', winlen, '-', mcmc_names[j], '.rds', sep=''))
  gr <- c(gr_diag(mcmc_out[1:4]))
  bi <- min(which(gr < 1.1))
  burn.in <- bi*0.05*nrow(mcmc_out[[1]]$samples)
  samp <- do.call('rbind', lapply(mcmc_out[1:4], function(l) l$samples[(burn.in+1):nrow(l$samples),]))
  idx <- sample(1:nrow(samp), n.samp, replace=T)
  mcmc_samp[[j]] <- samp[idx,]
  if (grepl('log', mcmc_names[j])) {
    log.response <- 'log'
  } else {
    log.response <- 'identity'
  }
  
  # MC sample prediction distribution
  # sample Phi from observation range
  Phi.samp <- runif(n.samp, as.numeric(obs['Psi']-obs['dPsi']), as.numeric(obs['Psi']+obs['dPsi']))
  
  if (log.response == 'identity') {
    ev.samp <- mcmc_samp[[j]][, 1] + mcmc_samp[[j]][, 2]*Phi.samp
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
  } else {
    r.samp <- rnorm(n.samp, 0, mcmc_samp[[j]][, 3])
    ecs.samp[[mcmc_names[j]]] <- exp(mcmc_samp[[j]][, 1] + mcmc_samp[[j]][, 2]*log(Phi.samp) + r.samp)
  }
}


saveRDS(mcmc_samp, file.path(R.path, 'samples', paste('mcmc_samp-', winlen, '.rds', sep='')))

saveRDS(ecs.samp, file.path(out.path, 'samples'paste('ecs_samp-', winlen, '.rds', sep='')))
inlen, '.rds', sep=''))
ep=''))
))
ep=''))
]][2:17, 'ECS'])
  }
  print(paste('WAIC:', waic, sep=' '))
  
  # MC sample prediction distribution
  # sample Phi from observation range
  Phi.samp <- runif(n.samp, as.numeric(obs['Psi']-obs['dPsi']), as.numeric(obs['Psi']+obs['dPsi']))
  
  if (log.response == 'identity') {
    ev.samp <- mcmc_samp[[j]][, 1] + mcmc_samp[[j]][, 2]*Phi.samp
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
  } else {
    r.samp <- rnorm(n.samp, 0, mcmc_samp[[j]][, 3])
    ecs.samp[[mcmc_names[j]]] <- exp(mcmc_samp[[j]][, 1] + mcmc_samp[[j]][, 2]*log(Phi.samp) + r.samp)
  }

  q <- quantile(ecs.samp[[mcmc_names[j]]], c(0, 0.025, 0.05, 0.165, 0.5, 1-0.165, 0.95, 0.975, 1))
  print(q)
  print(paste('Pr(ECS < 1)', sum(ecs.samp[[mcmc_names[j]]] < 1)/length(ecs.samp[[mcmc_names[j]]]), sep=' '))
  print(paste('Pr(ECS < 1.5)', sum(ecs.samp[[mcmc_names[j]]] < 1.5)/length(ecs.samp[[mcmc_names[j]]]), sep=' '))
  print(paste('Pr(ECS > 4.5)', sum(ecs.samp[[mcmc_names[j]]] > 4.5)/length(ecs.samp[[mcmc_names[j]]]), sep=' '))
  print(paste('Pr(ECS > 6)', sum(ecs.samp[[mcmc_names[j]]] > 6)/length(ecs.samp[[mcmc_names[j]]]), sep=' '))


}

stopCluster(cl)

saveRDS(mcmc_samp, paste('output/samples/mcmc_samp-', winlen, '.rds', sep=''))

saveRDS(ecs.samp, paste('output/samples/ecs_samp-', winlen, '.rds', sep=''))
ep=''))