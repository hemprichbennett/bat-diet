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

desired_mets <- c('niche overlap', 'functional complementarity')

vals <- lapply(nets, function(x) networklevel(x, index = desired_mets, level = 'higher'))

df <- data.frame(do.call(rbind, vals))


random <- lapply(nets, function(x) replicate(1000, bipartite::networklevel(vegan::permatfull(x, fixedmar= 'both',mtype="count",times=1)$perm[[1]],
                                        index = 'connectance', level = 'higher')))


rbind(sapply(b, function(x) quantile(x, probs=c(0.025, 0.975))),
sapply(nets, function(x) networklevel(x, index = 'connectance')))



fun_comp <- signif_check(nets, index = c('functional complementarity'), level = 'higher')
niche <- signif_check(nets, index = c('niche overlap'), level = 'higher')

melted_fun_comp <- melt.list(fun_comp)
melted_fun_comp$metric <- rep('functional complementarity', nrow(melted_fun_comp))
melted_niche <- melt.list(niche)
melted_niche$metric <- rep('niche overlap', nrow(melted_niche))
sigs <- rbind(melted_fun_comp, melted_niche)
str(sigs)
colnames(sigs) <- c('value', 'network', 'list item', 'metric')
sigs$network <- as.factor(sigs$network)
sigs$`list item` <- as.factor(sigs$`list item`)
sigs$metric <- as.factor(sigs$metric)

ggplot(sigs[which(sigs$`list item`=='random'),], aes(value)) + 
  geom_histogram(binwidth = 0.5) +facet_grid(network ~ metric, scales = 'free')+
  geom_vline(data=sigs[which(sigs$`list item`=='real'),], aes(xintercept = value), colour="red")+
  geom_vline(data=sigs[which(sigs$`list item`=='quantiles'),], aes(xintercept = value), colour="blue")
  