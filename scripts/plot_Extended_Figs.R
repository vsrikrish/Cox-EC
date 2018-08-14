###############################################################################################################
# plot_Extended_Figs.R                                                                                        #
#                                                                                                             #
# This script plots AIC values for all three models ("normal", "log-normal", "heavy-tailed/t"                 #
#   against time window length (Extended Fig 1), as well as plotting tail probabilities for the               #
#    normal and log-normal models against time window length (Extended Fig 2).                                #
#                                                                                                             #
# It assumes that mcmc_driver.R and sample_ECS.R have been run for all window lengths and model types.        #
###############################################################################################################


library(reshape2)
library(ggplot2)
library(RColorBrewer)

R.path <- 'scripts'
out.path <- 'output'
fig.path <- 'figures'

source(file.path(R.path, 'likelihood.R'))

win.lens <- seq(45, 73, by=2)
models <- c('normal', 'log-normal', 't')

model.stats <- data.frame(win.len=rep(win.lens, each=3), model=rep(models, length(win.lens)), AIC=NA, LT=NA, RT=NA)

nmod <- 16 # number of CMIP models used

for (len in win.lens) {
  ecs.samp <- readRDS(file.path(out.path, 'samples', paste('ecs_samp-', len, '.rds', sep='')))

  for (model in models) {
    # get data frame row id to store stats in
    rowid <- which(model.stats$win.len == len & model.stats$model == model)
    # parse model name for a few parameters
    if (grepl('normal', model)) {
      density.name <- 'normal'
    } else {
      density.name <- 't'
    }
    if (grepl('log', model)) {
      log.response <- 'log'
    } else {
      log.response <- 'identity'
    }
    # load MCMC output for AIC calculation
    mcmc_out <- readRDS(file.path(out.path, 'MCMC', paste('mcmc_out-', len, '-', model, '.rds', sep='')))
    # compute AIC
    elpd.aic <- log.lik(mcmc_out[['MLE']], mcmc_out[['parnames']], mcmc_out[['dat']][2:(nmod+1),c('ECS', 'Psi')], mcmc_out[['density.name']]) - length(mcmc_out[['parnames']])
    aic <- -2*elpd.aic
    if (log.response == 'log') {
       aic <- aic + 2*sum(mcmc_out[['dat']][2:(nmod+1),'ECS'])
    }
    model.stats[rowid, 'AIC'] <- aic
    # compute tail probabilities
    model.stats[rowid, 'LT'] <- sum(ecs.samp[[model]] < 1.5)/length(ecs.samp[[model]])
    model.stats[rowid, 'RT'] <- sum(ecs.samp[[model]] > 4.5)/length(ecs.samp[[model]])
  }
}

# melt model
model.melt <- melt(model.stats, id.var=c('win.len', 'model'))

cols <- brewer.pal(4, 'Dark2')

# plot Extended Data Figure 1 (AIC vs. window length)

p <- ggplot(model.melt[model.melt[,'variable'] == 'AIC',]) + geom_line(aes(x=win.len, y=value, color=model), size=1) + geom_point(aes(x=win.len, y=value, color=model), shape=1) +  theme_bw(base_size=7, base_family='sans') + scale_y_continuous('Akaike Information Criteron', expand=c(0, 0), limits=c(28, 45), breaks=seq(25, 45, 5)) + theme(axis.text=element_text(size=6), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position='right', legend.key.width=unit(2, "lines")) + scale_x_continuous('Temperature Series Window Length (yrs)', breaks=seq(45, 75, 5), expand=c(0,0)) + scale_color_manual(name='Distribution', values=c('normal'=cols[1], 'log-normal'=cols[2], 't'=cols[3]), labels=c('Normal', 'Log-normal', 'Heavy-Tailed'), breaks=c('normal', 'log-normal', 't'))

png(file.path(fig.path, 'Extended_Fig_1.png'), height=64, width=89, unit='mm', res=600)
print(p)
dev.off()

pdf(file.path(fig.path, 'Extended_Fig_1.pdf'), height=2.5, width=3.5)
print(p)
dev.off()

# plot Extended Data Figure 2 (LT/RT probabilities vs. window length)

p <- ggplot(model.melt[model.melt[,'variable'] != 'AIC' & model.melt[,'model'] != 't',]) + geom_line(aes(x=win.len, y=value, color=model, linetype=variable), size=1) +  geom_point(aes(x=win.len, y=value, color=model, shape=variable)) + theme_bw(base_size=7, base_family='sans') + scale_y_continuous('Tail Probability', labels=scales::percent, expand=c(0, 0), limits=c(0, .05)) + theme(axis.text=element_text(size=6), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position='right', legend.key.width=unit(2, "lines")) + scale_x_continuous('Temperature Series Window Length (yrs)', breaks=seq(45, 75, 5), expand=c(0,0)) + scale_color_manual(name='Distribution', values=c('normal'=cols[1], 'log-normal'=cols[2]), labels=c('Normal', 'Log-normal'), breaks=c('normal', 'log-normal')) + scale_linetype_manual(name='Tail', values=c('LT'='solid', 'RT'='dashed'), labels=c('ECS < 1.5K', 'ECS > 4.5K'), breaks=c('LT', 'RT')) + scale_shape_manual(name='Tail', values=c('LT'=1, 'RT'=2), labels=c('ECS < 1.5K', 'ECS > 4.5K'), breaks=c('LT', 'RT'))

png(file.path(fig.path, 'Extended_Fig_2.png'), height=64, width=89, unit='mm', res=600)
print(p)
dev.off()

pdf(file.path(fig.path, 'Extended_Fig_2.pdf'), height=2.5, width=3.5)
print(p)
dev.off()

