#### Header ####
## Project: Bat-diet
## Script purpose: Calculating how network-level metrics are altered by species removal
## Date: 06/08/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

dir <- getwd()
basedir <- strsplit(dir, split ='/')[[1]][2]
#print(basedir)
if(grepl('data', basedir)){
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  #library(magrittr, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
  library(reshape2, lib.loc = '/data/home/btw863/r_packages/')
  library(netReducer, lib.loc = '/data/home/btw863/r_packages/')
  library(labeling, lib.loc = '/data/home/btw863/r_packages/')
}else{
  library(here)
  #library(magrittr)
  library(ggplot2)
  library(reshape2)
  library(netReducer)
  
}

setwd(here())


source('scripts/r/r_network_gen.r')

nets <- r_network_gen(lulu = T, filter_species = T, include_malua = F)



ind_vec <- c('functional complementarity', 'nestedness', 'mean number of shared partners', 'NODF', 'niche overlap', 'modularity',
             'discrepancy', 'weighted nestedness', 'NODF', 'weighted NODF')

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)

ind <- ind_vec[args]

#for(i in 1:length(ind_vec)){
  starttime <- Sys.time()
  
  cat('Starting ', ind, ' at ', starttime, '\n', sep = '')
  
  #big <- call_reduce(net = nets, ind = index, datatype = 'list')
  big <- removing_and_randomizing(network = nets, index = ind, network_level = 'higher', sums_to_preserve = 'both', datatype='list', nreplicates = 100)
  
  # big$sp %<>%
  #   gsub('Hice', 'Hipposideros cervinus', .)%<>%
  #   gsub('Hidi', 'Hipposideros diadema', .)%<>%
  #   gsub('Hidy', 'Hipposideros dyacorum', .)%<>%
  #   gsub('Hiri', 'Hipposideros ridleyi', .)%<>%
  #   gsub('Keha', 'Kerivoula hardwickii', .)%<>%
  #   gsub('Kein', 'Kerivoula intermedia', .)%<>%
  #   gsub('Kemi', 'Kerivoula minuta', .)%<>%
  #   gsub('Kepa', 'Kerivoula papillosa', .)%<>%
  #   gsub('Rhbo', 'Rhinolophus borneensis', .)%<>%
  #   gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
  #   gsub('Rhtr', 'Rhinolophus trifoliatus', .)
  # 
  endtime <- Sys.time()
  write.csv(big, paste('results/random_removals/', ind, '.csv', sep = ''))
  big$sp <- gsub(' ', '\n', big$sp)
  
  big_actuals <- data.frame(big$sp, big$network, big$actual, big$lower, big$upper)
  colnames(big_actuals) <- c('sp', 'network', 'actual', 'lower', 'upper')
  rand_plot <- ggplot(data=big, aes(rand_vals)) + geom_histogram(binwidth = 3)+ 
    geom_vline(aes(xintercept=actual), data = big_actuals)+
    geom_vline(aes(xintercept=lower), data = big_actuals, linetype = 'dashed')+
    geom_vline(aes(xintercept=upper), data = big_actuals, linetype = 'dashed')+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    facet_grid(sp ~ network)+ ggtitle(ind)
  
  pdf(paste('plots/random_removals/', ind, '.pdf', sep = ''), width = 13)
  print(rand_plot)
  dev.off()

  cat(ind, 'finished, took ', endtime - starttime, '\n')
