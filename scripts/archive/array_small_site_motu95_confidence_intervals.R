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

chosen_ind <- ind[args]

print(chosen_ind)

small_netlists <- list()
small_netlists[['92']] <- netlists$`92`
small_netlists[['95']] <- netlists$`95`

real_and_errors <- randomized_ranges(small_netlists, indices = chosen_ind, network_level = 'higher', out_format = 'list',  actual_vals = T, summarise = F, n_perm = 1000, quantiles_to_return = NA)

outname <- paste('data/output_data/randomized_ranges/small_', chosen_ind, '.RDA', sep = '')

save.image(outname)
