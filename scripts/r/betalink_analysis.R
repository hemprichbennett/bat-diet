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
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(reshape2)
  library(corrplot)
  library(iNEXT)
  library(DataExplorer)
  library(betalink)
}

setwd(here())
source('scripts/r/r_network_gen.r')


nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T, include_malua = F)

graphs <- prepare_networks(nets)


network_betadiversity(graphs)

#metaweb(graphs)
