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
  library(DataExplorer)
}

setwd(here())
source('scripts/r/r_network_gen.r')


nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T, include_malua = F)
names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))

write.csv(unlist(lapply(nets, sum)), 'results/n_links.csv')

graphs <- prepare_networks(nets)


beta <- network_betadiversity(graphs)
#Making a new object to allow prettier plotting
temp <- beta
colnames(temp)[c(1,2)] <- c('j', 'i')
for_plot <- rbind(beta, temp)
#Order the factors so that the plot looks nice
for_plot$i <- ordered(for_plot$i, levels = c("Danum, 2016", "Danum, 2017", "Maliau, 2016", 'Maliau, 2017',
                                             'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017'))
for_plot$j <- ordered(for_plot$j, levels = c("Danum, 2016", "Danum, 2017", "Maliau, 2016", 'Maliau, 2017',
                                             'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017'))
for_plot$STWN <- for_plot$ST/for_plot$WN
melted_forplot <- melt(for_plot)


betaplot <- ggplot(melted_forplot, aes(i, fct_rev(j)))+ geom_point(aes(size=value, colour = value))+
  scale_colour_gradient(low = "white",
                        high = "blue", limits = c(0,1))+
  theme_dark()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')+
  facet_wrap(~ variable, labeller = label_bquote(italic(beta [.(as.character(variable))])),
             nrow = 3)+
  labs(x= NULL, y = NULL)



betaplot  
pdf('plots/betaplot.pdf')
betaplot
dev.off()

