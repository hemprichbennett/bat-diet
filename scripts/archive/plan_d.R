##################################################
## Project: bat-diet (for full-species networks)
## Script purpose: Seeing the significance of various metrics across each of my networks SITEWISE
## Date: 21/06/18
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

inpath <- 'data/processed_dna_data/lulu/'
filenames <- list.files(pattern = '.csv', path = inpath)
#filenames <- filenames
filenames
filenames <- paste(inpath, filenames, sep = '')
#filenames <- filenames[grep('lulu', filenames)]

rawnets <- lapply(filenames, read.csv, header = F, stringsAsFactors = F, row.names=1)
names(rawnets) <- gsub('.*\\/', '', filenames)
names(rawnets) <- gsub('_.+', '', names(rawnets))
netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))

names(netlists) <- names(rawnets)

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)

ind <- c('functional complementarity',
         'web asymmetry',
         'Alatalo interaction evenness',
         'togetherness',
         'Fisher alpha', 'mean number of shared partners',
         'niche overlap',
         'nestedness',
         'discrepancy',
         'ISA', 'weighted nestedness', 'NODF', 'weighted NODF', 'modularity')

clust <- seq(91, 98)

pairs <- expand.grid(ind, clust)

chosen_pair <- pairs[args,]

chosen_ind <- chosen_pair[,1]

clust <- chosen_pair[,2]

print(chosen_ind)
print(clust)

#small_netlists <- list()
small_netlist <- netlists[[which(names(netlists)==clust)]]

for(i in 1:length(small_netlist)){
  small_netlist[[i]] <- apply(small_netlist[[i]], 2, as.numeric)
  small_netlist[[i]] <- ifelse(small_netlist[[i]]==0,0,1)
}



sums_to_preserve <- 'both'
network_level <- 'higher'
n_perm <- 1000

if(chosen_ind != 'modularity'){
  actuals <- lapply(small_netlist, function(x)
    bipartite::networklevel(x, index = chosen_ind, level = network_level))
  
  names(actuals) <- names(small_netlist)
  actuals
  
  rand_list <- lapply(small_netlist, function(x)
    replicate(n_perm, bipartite::networklevel(vegan::permatswap(x, fixedmar=sums_to_preserve,mtype="count",times=1, method="quasiswap")$perm[[1]],
                                              index = chosen_ind, level = network_level)))
  
}else{
  actuals <- lapply(small_netlist, function(x)
    slot(bipartite::computeModules(web = vegan::permatswap(x, fixedmar=sums_to_preserve,mtype="count",times=1, method="quasiswap")$perm[[1]]), 'likelihood'))
  
  
  names(actuals) <- names(small_netlist)
  actuals
  
  rand_list <- lapply(small_netlist, function(x)
    replicate(n_perm, slot(bipartite::computeModules(web = vegan::permatswap(x, fixedmar=sums_to_preserve,mtype="count",times=1, method="quasiswap")$perm[[1]]), 'likelihood')))
  
}
names(rand_list) <- names(small_netlist)

rand_out <- do.call(cbind, rand_list)

actual_out <- do.call(cbind, actuals)

actual_outname <- paste('results/rand_values/actual_', clust,'_', chosen_ind, '.csv', sep = '')
rand_outname <- paste('results/rand_values/rand_', clust, '_', chosen_ind, '.csv', sep = '')

write.csv(actual_out, actual_outname, row.names = F)
write.csv(rand_out, rand_outname, row.names = F)



