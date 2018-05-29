##################################################
## Project: bat diet (all species)
## Script purpose: checking the beta-diversity of my networks
## Date: 28/05/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

dir <- getwd()
basedir <- strsplit(dir, split ='/')[[1]][2]
print(basedir)
if(grepl('data', basedir)){
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  
}else{
  library('here')
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(betalink)
  library(forcats)
}

setwd(here())
source('scripts/r/r_network_gen.r')


nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T, include_malua = F)
names(nets) <- gsub('.+ ^', '.+^', names(nets))
graphs <- prepare_networks(nets)


beta <- network_betadiversity(graphs)
#Making a new object to allow prettier plotting
temp <- beta
colnames(temp)[c(1,2)] <- c('j', 'i')
for_plot <- rbind(beta, temp)
#Order the factors so that the plot looks nice
for_plot$i <- ordered(for_plot$i, levels = c("DANUM, 2016", "DANUM, 2017", "MALIAU, 2016", 'MALIAU, 2017',
                                             'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017'))
for_plot$j <- ordered(for_plot$j, levels = c("DANUM, 2016", "DANUM, 2017", "MALIAU, 2016", 'MALIAU, 2017',
                                             'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017'))
melted_forplot <- melt(for_plot)


betaplot <- ggplot(melted_forplot, aes(i, fct_rev(j)))+ geom_point(aes(size=value, colour = value))+
  scale_colour_gradient(low = "black",
                      high = "blue")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_wrap(~ variable, labeller = label_bquote(beta [italic(.(variable))]))+
  facet_wrap(~ variable, labeller = label_bquote(italic(beta [.(as.character(variable))])))+
  labs(x= NULL, y = NULL)
betaplot  
pdf('plots/betaplot.pdf')
betaplot
dev.off()
