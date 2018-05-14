#This script rarifies networks of cervinus to look at individual-level changes in network topology
if(interactive()==TRUE){
  library(here)
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
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


# IDs <- c()
# n_bats <- c()
# nets <- c()
# mets <- c()
# means <- c()
# ses <- c()
# 
# for(i in 1:length(unique(melted_mat$`bat ID`))){
#   ID <- unique(melted_mat$`bat ID`)[i]
#   temp <- melted_mat[melted_mat$`bat ID`==ID,]
#   for(n in seq(range(temp$N_bats)[1], range(temp$N_bats)[2])){
#     subs <- temp[which(temp$N_bats==n),]
#     for(m in 1:length(unique(subs$variable))){
#       met <- unique(subs$variable)[m]
#       met_set <- subs[subs$variable==met,]
#       IDs <- c(IDs, as.character(ID))
#       n_bats <- c(n_bats, n)
#       nets <- c(nets, subs$Network[1])
#       mets <- c(mets, met)
#       means <- c(means, mean(subs$value))
#       ses <- c(ses, std.error(subs$value))
#     }
#   }
# }

all_summary <- ddply(melted_mat, c("`bat ID`", "N_bats", "variable", 'Network'), summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)),
      count = length(value))

all_summary$Network <- ordered(all_summary$Network, levels=unique(all_summary$Network)[order(as.character(unique(all_summary$Network)))])


meanplot <- ggplot(all_summary[all_summary$variable=='Proportional similarity',], aes(mean, sd))+ 
  geom_point(aes(colour= N_bats), alpha = 1/2) +
  facet_wrap(~ Network, ncol = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab('Standard deviation')+ xlab('Mean')

meanplot


countplot <- ggplot(all_summary[all_summary$variable=='Proportional similarity',], aes(count, sd))+ 
  geom_point(aes(colour= N_bats), alpha = 1/2) +
  facet_wrap(~ Network, ncol = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab('Standard deviation')+ xlab('Number of observations per individual')
countplot

grid.arrange(meanplot, countplot)

# all_melted$habitat_type <- rep(NA, nrow(all_melted))
# all_melted[grep('SAFE', all_melted$SiteAndYear), 'habitat_type'] <- 'Logged'
# all_melted[grep('SBE', all_melted$SiteAndYear), 'habitat_type'] <- 'Logged, replanted'
# all_melted[grep('Danum', all_melted$SiteAndYear), 'habitat_type'] <- 'Primary'
# all_melted[grep('Maliau', all_melted$SiteAndYear), 'habitat_type'] <- 'Primary'
# 
