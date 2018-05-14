#This script rarifies networks of cervinus to look at individual-level changes in network topology
if(interactive()==TRUE){
  library(here)
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(grid)
  library(forcats)
  library(reshape2)
  library(lme4)
  library(plyr)
  library(plotrix)
  #library(multcompView)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())

load('data/output_data/hice_stats/rarified_stats.RDS')


#out_36 <- outmat

melted_mat <- melt(outmat[,-c(3,4)], id.vars = c('N_bats','bat ID', 'Network'))


melted_mat$variable <- gsub('^degree$', 'Normalised degree', melted_mat$variable)
melted_mat$variable <- gsub('normalised degree', 'Normalised degree', melted_mat$variable)
melted_mat$variable <- gsub('partner diversity','Partner diversity', melted_mat$variable)
melted_mat$variable <- gsub('proportional similarity','Proportional similarity', melted_mat$variable)



all_summary <- ddply(melted_mat, c("`bat ID`", "N_bats", "variable", 'Network'), summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)),
      count = length(value))

all_summary$Network <- ordered(all_summary$Network, levels=unique(all_summary$Network)[order(as.character(unique(all_summary$Network)))])


grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 3,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

meanplot <- ggplot(all_summary[all_summary$variable=='Proportional similarity',], aes(mean, sd))+ 
  geom_point(aes(colour= N_bats), alpha = 1/2) +
  facet_wrap(~ Network, ncol = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab('Standard deviation')+ xlab('Mean')+
  labs(color='Number of bats')

meanplot


countplot <- ggplot(all_summary[all_summary$variable=='Proportional similarity',], aes(count, sd))+ 
  geom_point(aes(colour= N_bats), alpha = 1/2) + #Don't plot the legend on this one
  facet_wrap(~ Network, ncol = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab('Standard deviation')+ xlab('Number of observations per individual')+
  labs(color='Number of bats')
countplot


grid_arrange_shared_legend(meanplot, countplot)

pdf('plots/Hice/individual_effects_of_rarefaction', onefile = F)
grid_arrange_shared_legend(meanplot, countplot)
dev.off()
