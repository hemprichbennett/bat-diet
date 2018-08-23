##################################################
## Project: All bat diet
## Script purpose: Showing the difference in functional complementarity etc between networks
## Date: 29/05/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################




dir <- getwd()
basedir <- strsplit(dir, split ='/')[[1]][2]
print(basedir)
if(grepl('data', basedir)){
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  
}else{
  library(here)
  library(ggplot2)
  library(tidyverse)
  library(forcats)
  library(reshape2)
  library(LOTUS)
}

setwd(here())
source('scripts/r/r_network_gen.r')
source('scripts/r/hernani_comparisons.R')

signif_check <- function(networks, mars= 'both',mtype="count", index, level = 'higher'){
  if(length(index)>1){stop('this function fucks up if you use more than one index at a time. Sorry')}
  random <- lapply(networks, function(x) replicate(1000, bipartite::networklevel(vegan::permatfull(x, fixedmar= mars, mtype="count",times=1)$perm[[1]],
                                                                        index = index, level = level)))
  quantiles <- lapply(random, function(x) quantile(x, probs=c(0.025, 0.975)))
  real <- lapply(nets, function(x) networklevel(x, index = index, level = level))
  out <- list(random, quantiles, real)
  names(out) <- c('random', 'quantiles', 'real')
  return(out)
}

nets <- r_network_gen(collapse_species = T, desired_species = NULL, filter_species = T,
                      include_malua = F, lulu= T)

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))



## Section: Analysing and plotting differences in f-c from null ####
##################################################
fun_comp <- signif_check(nets, index = c('functional complementarity'), level = 'higher')
melted_fun_comp <- melt.list(fun_comp)
colnames(melted_fun_comp) <- c('value', 'network', 'list item')
melted_fun_comp$network <- as.factor(melted_fun_comp$network)
melted_fun_comp$`list item` <- as.factor(melted_fun_comp$`list item`)


fun_comp_plot <- ggplot(melted_fun_comp[which(melted_fun_comp$`list item`=='random'),], aes(value)) + 
  geom_histogram(binwidth = 0.5) +facet_wrap( ~ network, ncol = 1)+
  geom_vline(data=melted_fun_comp[which(melted_fun_comp$`list item`=='real'),], aes(xintercept = value), colour="red")+
  geom_vline(data=melted_fun_comp[which(melted_fun_comp$`list item`=='quantiles'),], aes(xintercept = value), colour="blue")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x='Functional complementarity', y = 'Count')
fun_comp_plot  
pdf('plots/fun_comp_sig.pdf')
fun_comp_plot  
dev.off()


#### Section: Analysing and plotting differences in niche from null####
##################################################

niche <- signif_check(nets, index = c('niche overlap'), level = 'higher')
melted_niche <- melt.list(niche)
colnames(melted_niche) <- c('value', 'network', 'list item')
melted_niche$network <- as.factor(melted_niche$network)
melted_niche$`list item` <- as.factor(melted_niche$`list item`)


niche_plot <- ggplot(melted_niche[which(melted_niche$`list item`=='random'),], aes(value)) + 
  geom_histogram(binwidth = 0.005) +facet_wrap( ~ network, ncol = 1)+
  geom_vline(data=melted_niche[which(melted_niche$`list item`=='real'),], aes(xintercept = value), colour="red")+
  geom_vline(data=melted_niche[which(melted_niche$`list item`=='quantiles'),], aes(xintercept = value), colour="blue")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x='Niche overlap', y = 'Count')
niche_plot

pdf('plots/niche_sig.pdf')
niche_plot  
dev.off()

## Section: Mostly pointless bit where I facet by metric, but as the bin widths are so different 
## they show little. To be tinkered with later
##################################################
melted_fun_comp$metric <- rep('functional complementarity', nrow(melted_fun_comp))
melted_niche$metric <- rep('Niche overlap', nrow(melted_niche))
sig <- rbind(melted_fun_comp, melted_niche)
sig$metric <- as.factor(sig$metric)


sig_plot <- ggplot(sig[which(sig$`list item`=='random'),], aes(value)) + 
  geom_histogram(binwidth = 0.005) +facet_grid(metric ~ network, scales = 'free')+
  geom_vline(data=sig[which(sig$`list item`=='real'),], aes(xintercept = value), colour="red")+
  geom_vline(data=sig[which(sig$`list item`=='quantiles'),], aes(xintercept = value), colour="blue")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
sig_plot
