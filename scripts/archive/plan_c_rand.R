
##################################################
## Project: bat-diet (for full-species networks)
## Script purpose: Seeing the significance of various metrics across each of my networks SITEWISE
## Date: 12/09/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


if(interactive()==TRUE){
  library(here)
  library(LOTUS)
}else{
  library(reshape2, lib.loc = '/data/home/btw863/r_packages/')
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(LOTUS, lib.loc = '/data/home/btw863/r_packages/')
}



setwd(here())
getwd()

source('scripts/r/r_network_gen.r')
setwd(here())
# getwd()

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)

netlist <- r_network_gen()

ind <- c('functional complementarity',
         'web asymmetry',
         'Alatalo interaction evenness',
         'togetherness',
         'Fisher alpha', 'mean number of shared partners',
         'niche overlap',
         'nestedness',
         'discrepancy',
         'ISA', 'weighted nestedness', 'NODF', 'weighted NODF', 'modularity')

chosen_ind <- ind[args]

print(chosen_ind)



sums_to_preserve <- 'both'
network_level <- 'higher'
n_perm <- 1000


if(chosen_ind != 'modularity'){
  actuals <- lapply(netlist, function(x)
    bipartite::networklevel(x, index = chosen_ind, level = network_level))
  #for(i in 1:length(netlist)){names(actuals)[i] <- names(netlist)[i]}
  names(actuals) <- names(netlist)
  actuals
  
  rand_list <- lapply(netlist, function(x)
    replicate(n_perm, bipartite::networklevel(vegan::permatswap(x, fixedmar=sums_to_preserve,mtype="count",times=1, method="quasiswap")$perm[[1]],
                                              index = chosen_ind, level = network_level)))
  
}else{
  actuals <- lapply(netlist, function(x)
    slot(bipartite::computeModules(web = vegan::permatswap(x, fixedmar=sums_to_preserve,mtype="count",times=1, method="quasiswap")$perm[[1]]), 'likelihood'))
  
  #for(i in 1:length(netlist)){names(actuals)[i] <- names(netlist)[i]}
  names(actuals) <- names(netlist)
  actuals
  
  rand_list <- lapply(netlist, function(x)
    replicate(n_perm, slot(bipartite::computeModules(web = vegan::permatswap(x, fixedmar=sums_to_preserve,mtype="count",times=1, method="quasiswap")$perm[[1]]), 'likelihood')))
  
}
names(rand_list) <- names(netlist)

rand_out <- do.call(cbind, rand_list)

actual_out <- do.call(cbind, actuals)

actual_outname <- paste('results/rand_values_95/actual_', chosen_ind, '.csv', sep = '')
rand_outname <- paste('results/rand_values_95/rand_', chosen_ind, '.csv', sep = '')

write.csv(actual_out, actual_outname, row.names = F)
write.csv(rand_out, rand_outname, row.names = F)
