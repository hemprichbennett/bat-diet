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

n_bats <- c()
it <- c()
n_nets <- c()
n_signif <- c()

#Make an empty list of lists, with the first level being the number of bats, the second being the iteration
rarelist <- list()
for(n in 20:68){
  rarelist[[as.character(n)]] <- list()
  for(i in 1:length(unique(outmat$iteration))){
    rarelist[[as.character(n)]][[i]] <- list()
  }
  }

#Fill the list
for(i in 1: length(unique(outmat$iteration))){
  for(n in 20:68){
    temp <- outmat[which(outmat$N_bats==n & outmat$iteration==i),]
    #cat(dim(temp),'\n')
    for(m in c(1,2,3,4)){
      if(m==1){
        anova_object <- TukeyHSD(aov(temp$`proportional similarity` ~ temp$Network))
        rarelist[[as.character(n)]][[i]][['proportional_similarity']] <- as.matrix(t(anova_object$`temp$Network`))
      }else if(m==2){
        anova_object <- TukeyHSD(aov(temp$`resource range` ~ temp$Network))
        rarelist[[as.character(n)]][[i]][['resource_range']] <- as.matrix(t(anova_object$`temp$Network`))
      }else if(m==3){
        anova_object <- TukeyHSD(aov(temp$`degree` ~ temp$Network))
        rarelist[[as.character(n)]][[i]][['degree']] <- as.matrix(t(anova_object$`temp$Network`))
      }else if(m==4){
        anova_object <- TukeyHSD(aov(temp$`normalised degree` ~ temp$Network))
        rarelist[[as.character(n)]][[i]][['normalised_degree']] <- as.matrix(t(anova_object$`temp$Network`))
      }
    }
    # if(length(which(anova_object$`temp$Network`[,4] <=0.05))>0){
    #   cat('N_bats is', n, 'number of significant differences is', length(which(anova_object$`temp$Network`[,4] <=0.05)), '\n')
    # }
    #n_bats <- c(n_bats, n)
    #it <- c(it, i)
    #n_nets <- c(n_nets, length(unique(temp$Network)))
    #n_signif <- c(n_signif, length(which(anova_object$`temp$Network`[,4] <=0.05)))
  }
}

raremelted <- melt(rarelist)
colnames(raremelted) <- c('statistic', 'Pair of networks', 'value', 'metric', 'iteration', 'n_bats')
raremelted <- raremelted[which(raremelted$statistic=='p adj'),]# Make it so we only have the p-values
raremelted$n_bats <- as.integer(raremelted$n_bats)

pairs <- c()
ns <- c()
percents <- c()
mets <- c()
for(p in 1:length(unique(raremelted$`Pair of networks`))){
  pair <- as.character(unique(raremelted$`Pair of networks`)[p])
  vals <- raremelted[which(raremelted$`Pair of networks`==pair),]
  for(n in seq(range(vals$n_bats)[1], range(vals$n_bats)[2])){# Doing it this way as otherwise the sample unevenness confuses it
    subvals <- vals[which(vals$n_bats==n),]
    for(m in 1:length(unique(raremelted$metric))){
      subvals <- subvals[which(subvals$metric==unique(raremelted$metric)[m]),]
      percent <- length(which(subvals$value <=0.05))/100
      pairs <- c(pairs, pair)
      ns <- c(ns, n)
      mets <- c(mets, unique(raremelted$metric)[m])
      percents <- c(percents, percent)
    }
    
    }
  }

long_df <- data.frame(pairs, ns, mets, percents)
#Getting the factors in the right order takes a little bit of work #
z <- tapply(long_df$ns, long_df$pairs, max)[order(tapply(long_df$ns, long_df$pairs, max))]
long_df$pairs <- factor(long_df$pairs, levels = names(z))

long_df$mets <- gsub('proportional_similarity', 'Proportional similarity', long_df$mets)
long_df$mets <- gsub('resource_range', 'Resource range', long_df$mets)

pairwise_barchart <- ggplot(data=long_df, aes(ns, fct_rev(pairs)))+ geom_tile(colour = 'black', aes(fill=percents))+
  theme_dark()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
   labs(x = 'Number of bats in network', y = "Pair of networks")+
  scale_fill_gradient2(low= 'black', high ='blue',name = 'Percent\nsignificantly\ndifferent')+
  facet_wrap(~ mets)
  
pairwise_barchart
pdf('plots/Hice/rarified_mets.pdf')
pairwise_barchart
dev.off()


