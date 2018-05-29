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

nets <- r_network_gen(collapse_species = T, desired_species = NULL, filter_species = T,
                      include_malua = F, lulu= T)

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))

desired_mets <- c('niche overlap', 'functional complementarity')

vals <- lapply(nets, function(x) networklevel(x, index = desired_mets, level = 'higher'))

df <- data.frame(do.call(rbind, vals))


randomized_ranges(networks = nets, input_format = 'clust_only', indices = 'connectance', network_level = 'higher',
                  sums_to_preserve = 'both')
