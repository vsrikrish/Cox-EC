library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

# read data
dat <- read.csv('cmip-dat.csv', header=F)
colnames(dat) <- c('ECS', 'Phi')
# read MCMC output
mcmc_names <- c('normal', 'log-normal', 't')
mcmc_samp <- readRDS('mcmc_samp.rds')

# allocate storage for plots
p <- vector('list', length(mcmc_names))

# generate sequence of Phi values for curve-plotting
Phi <- seq(0.05, 0.30, length=1000)

for (j in 1:length(mcmc_names)) {
  print(mcmc_names[j])
 # mcmc_means <- colMeans(mcmc_samp[[j]])
  if (grepl('log', mcmc_names[j])) {
    log.response <- T
  } else {
    log.response <- F
  }

#  ECS_mean <- mcmc_means[1] + mcmc_means[2]*Phi
#  if (log.response) {
#    ECS_mean <- exp(mcmc_means[1] + mcmc_means[2]*log(Phi))
#  } else {
#    ECS_mean <- mcmc_means[1] + mcmc_means[2]*Phi
#  }
#  mean.df <- data.frame(ECS=ECS_mean, Phi=Phi)

  if (grepl('t', mcmc_names[j])) {
    ECS_resid <- rt(nrow(mcmc_samp[[j]]), df=mcmc_samp[[j]][,3])
  } else {
    ECS_resid <- rnorm(nrow(mcmc_samp[[j]]), mean=0, sd=mcmc_samp[[j]][,3])
  }

  if (log.response) {
    ECS_samp <- exp(sweep(outer(mcmc_samp[[j]][,2], log(Phi), '*'), 1, mcmc_samp[[j]][,1] + ECS_resid, FUN='+'))
  } else {
    ECS_samp <- sweep(outer(mcmc_samp[[j]][,2], Phi, '*'), 1, mcmc_samp[[j]][,1] + ECS_resid, FUN='+')
  }
  
#  ECS_lin_q <- t(apply(ECS_lin, 2, quantile, probs=c(0.025, 0.975)))
#  lin_q.df <- data.frame(Phi=Phi)
#  lin_q.df[,c('ub', 'lb')] <- ECS_lin_q
  
  ECS_samp_q <- t(apply(ECS_samp, 2, quantile, probs=c(0.025, 0.165, 0.50, 0.835,  0.975)))
  samp_q.df <- data.frame(Phi=Phi)
  samp_q.df[,c('lb95', 'lb67', 'median', 'ub67', 'ub95')] <- ECS_samp_q
  
  p[[j]] <- ggplot() + geom_ribbon(data=samp_q.df, aes(x=Phi, ymin=lb95, ymax=ub95, color='95% CI'), alpha=0.2) + geom_ribbon(data=samp_q.df, aes(x=Phi, ymin=lb67, ymax=ub67, color='67% CI'), alpha=0.2) + theme_bw(base_size=7, base_family='sans') + geom_line(data=samp_q.df, aes(x=Phi, y=median, color='Median')) + geom_point(data=dat[1:16,], aes(x=Phi, y=ECS, color='CMIP Data'), shape='x', size=3.5) + coord_cartesian(xlim=c(0.05, 0.30), ylim=c(0.5, 6)) + scale_y_continuous('ECS (K)', breaks=1:6, expand=c(0,0)) + scale_x_continuous(expression(Psi),  breaks=seq(0.05, 0.30, 0.05), expand=c(0, 0)) + theme(panel.grid.major=element_blank(), axis.text.x=element_text(size=6), panel.grid.minor=element_blank(), plot.margin=unit(c(5.5, 7.5, 5.5, 5.5), "points"), legend.position='right', legend.key = element_rect(color='white')) + geom_vline(xintercept=0.13-0.016, color='blue', linetype='dashed') + geom_vline(xintercept=0.13+0.016, color='blue', linetype='dashed') + scale_color_manual('', values=c('CMIP Data'='black', 'Median'='red', '67% CI'='gray74', '95% CI'='gray37'), breaks=c('CMIP Data', 'Median', '67% CI', '95% CI'), guide=guide_legend(override.aes=list(linetype=c('blank', 'solid', 'blank', 'blank'), shape=c('x', NA, NA, NA), fill=c(NA, NA, 'gray37', 'gray74'))))
  
    
  if (j == 1) {
    legend <- gtable_filter(ggplot_gtable(ggplot_build(p[[j]])), "guide-box")
  }
  
  p[[j]] <- p[[j]] + theme(legend.position='none')
  
  
  p[[j]] <- arrangeGrob(p[[j]], top=textGrob(letters[(j+1)], x=unit(0.025, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=8, fontfamily='sans', fontface='bold')))

}

g <- ggplot(dat[1:16,]) + geom_rect(aes(xmin=1.5, xmax=4.5, ymin=0, ymax=Inf), fill='gray', alpha=0.1) + geom_histogram(aes(x=ECS, y=..density..), breaks=seq(2, 5, 0.5), fill='blue', color='white', closed='left', alpha=0.5) + stat_density(aes(x=ECS, y=..density..), geom='line', color='black', bw='nrd0', adjust=.6) + geom_rug(aes(x=ECS, y=ECS), colour='blue', sides='b', position='jitter') + scale_x_continuous('ECS (K)', limits=c(0, 6), breaks=seq(0, 6, 0.5), expand=c(0, 0)) + scale_y_continuous('Density', limits=c(0, 0.8)) + theme_bw(base_size=7, base_family='sans') + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.margin=unit(c(5.5, 7.5, 5.5, 5.5), "points"),  axis.text.x=element_text(size=6), axis.ticks.y=element_blank())

g <- arrangeGrob(g, top=textGrob('a', x=unit(0.025, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=8, fontfamily='sans', fontface='bold')))

print('Plotting...')

pdf('ecs_fit.pdf', height=6, width=7.2)
grid.arrange(arrangeGrob(g,
             p[[1]],
             p[[2]],
             p[[3]],
             ncol=2),
             legend,
             widths=c(1.5, 0.2),
             ncol=2)
dev.off()


png('ecs_fit.png', height=152, width=183, units='mm', res=600)
grid.arrange(arrangeGrob(g,
             p[[1]],
             p[[2]],
             p[[3]],
	     ncol=2),
             legend,
	     widths=c(1.5, 0.2),
	     ncol=2)
dev.off()
.2),
	     ncol=2)
dev.off()
