##################################################
## Project: bat-diet (for full-species networks)
## Script purpose: Seeing the significance of various metrics across each of my networks SITEWISE
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


if(interactive()==TRUE){
  library(here)
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

real_and_errors <- randomized_ranges(netlists, indices = chosen_ind, network_level = 'higher', out_format = 'data.frame', quantiles_to_return = c(0.025, 0.975), actual_vals = T, n_perm = 100)

outname <- paste('data/output_data/randomized_ranges/', chosen_ind, '.csv', sep = '')
write.csv(real_and_errors, outname)

cat(outname, 'written')


#save.image(paste('data/output_data/all_bats/', chosen_ind, '_real_and_error_calcs.RDS', sep = ''))
#plot_str(trial, type = "r")
