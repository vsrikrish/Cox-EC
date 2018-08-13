library(reshape2)
library(RColorBrewer)
library(ggplot2)
#library(grid)
#library(gridExtra)
#library(gtable)

out.path <- 'output'
fig.path <- 'figures'

print('Reading in data...')
# read sampled ECS values
win.len <- 55
ecs.samp <- readRDS(file.path(out.path, 'samples', paste('ecs_samp-', win.len, '.rds', sep='')))
samp.melt <- melt(ecs.samp)
colnames(samp.melt) <- c('sample', 'distribution')
samp.melt$distribution <- factor(samp.melt$distribution, levels=c('normal', 'log-normal', 't'))

cols <- brewer.pal(4, 'Dark2')

print('Plotting PDFs...')
# plot entire distribution
p <- ggplot(samp.melt)+ stat_density(aes(x=sample, y=..density.., color=distribution, linetype=distribution), size=1, geom='line', position='identity') +  theme_bw(base_size=7, base_family='sans') + scale_y_continuous('Density', expand=c(0, 0)) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=6), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position='bottom', legend.key.width=unit(2, "lines")) + scale_x_continuous('ECS (K)', breaks=0:8, expand=c(0,0), limits=c(-10, 20)) + annotate('rect', xmin=1.5, xmax=4.5, ymin=0, ymax=0.8, alpha=.1) + scale_color_manual(name='Model', values=c('normal'=cols[1], 'log-normal'=cols[2], 't'=cols[3]), labels=c('Normal', 'Log-normal', 'Heavy-Tailed'), breaks=c('normal', 'log-normal', 't')) + scale_linetype_manual(name='Model', values=c('normal'='solid', 'log-normal'='dashed', 't'='dotted'), labels=c('Normal', 'Log-normal', 'Heavy-Tailed'), breaks=c('normal', 'log-normal', 't'))

p1 <- p + coord_cartesian(xlim=c(0, 8), ylim=c(0, 0.8))

#p1 <- arrangeGrob(p1, top=textGrob('a', x=unit(0.005, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=8, fontfamily='sans', fontface='bold')))

# # plot cdfs
# print('Plotting CDFs...')
# # fit kernel density estimate to each ensemble set
# # we do this to reduce the size of the cdf images as well as to smooth the ecdfs
# dens <- lapply(ecs.samp, density, from=-2, to=10)
# dens.df <- lapply(dens, function(l) data.frame(ECS=l$x, cum.dens=cumsum(l$y)/sum(l$y)))
# dens.melt <- melt(dens.df, id='ECS')

# p.cdf <- ggplot(dens.melt) + geom_line(aes(x=ECS, y=value, color=L1, linetype=L1), size=1) +  theme_bw(base_size=7, base_family='sans') + geom_errorbar(aes(x=1, ymin=0, ymax=0.05, width=0.1, color='IPCC', linetype='IPCC')) + geom_errorbar(aes(x=6, ymin=0.9, ymax=1, width=0.25, color='IPCC', linetype='IPCC'), show.legend=F)  + scale_y_continuous('Cumulative Density', expand=c(0.05, 0.025), limits=c(0,1)) + theme(axis.text.x=element_text(size=6), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position='bottom', legend.key.width=unit(2, "lines")) + scale_x_continuous('ECS (K)', breaks=0:8, expand=c(0,0)) + annotate('rect', xmin=1.5, xmax=4.5, ymin=0, ymax=1.0, alpha=.1) + scale_color_manual(name='Distribution', values=c('normal'=cols[1], 'log-normal'=cols[2], 't'=cols[3], 'IPCC'=cols[4]), labels=c('Normal', 'Log-normal', 't', 'IPCC Range'), breaks=c('normal', 'log-normal', 't', 'IPCC')) + scale_linetype_manual(name='Distribution', values=c('normal'='solid', 'log-normal'='dashed', 't'='dotted', 'IPCC'='solid'), labels=c('Normal', 'Log-normal', 't', 'IPCC Range'), breaks=c('normal', 'log-normal', 't', 'IPCC'))

# e <- 1e-3

# p2 <- p.cdf + coord_cartesian(xlim=c(0, 1.5), ylim=c(0-e, 0.1+e), expand=F)
# p3 <- p.cdf + coord_cartesian(xlim=c(4.5, 8), ylim=c(0.9-e, 1+e), expand=F)

# legend <- gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box")

# p2 <- p2 + theme(legend.position='none')
# p3 <- p3 + theme(legend.position='none')

# p2 <- arrangeGrob(p2, top=textGrob('b', x=unit(0.01, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=8, fontfamily='sans', fontface='bold')))
# p3 <- arrangeGrob(p3, top=textGrob('c', x=unit(0.01, 'npc'), y=unit(-0.25, 'npc'),just=c('left', 'top'), gp=gpar(col='black', fontsize=8, fontfamily='sans', fontface='bold')))


pdf(file.path(fig.path, 'Fig_1.pdf'), height=3, width=5.5)
# grid.arrange(p1,
#           arrangeGrob(p2,
#                       p3,
#                       nrow = 1),
#           legend,
#           heights=c(1.1, 1.1, 0.2),
#           nrow = 3)
print(p1)
dev.off()

png(file.path(fig.path, 'Fig_1.png'), height=76, width=140, unit='mm', res=600)
# grid.arrange(p1,
#           arrangeGrob(p2,
#                       p3,
#                       nrow = 1),
#           legend,
#           heights=c(1.1, 1.1, 0.2),
#           nrow = 3)
print(p1)
dev.off()
