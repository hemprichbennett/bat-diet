##################################################
## Project: bat-diet (for full-species networks)
## Script purpose: Seeing the significance of various metrics across each of my networks SITEWISE
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


library(here)
library(reshape2)
library(LOTUS)
library(foreach)
library(doParallel)



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

for(i in 1:length(rawnets)){
  rawnets[[i]][2:nrow(rawnets[[i]]),] <- ifelse(rawnets[[i]][2:nrow(rawnets[[i]]), ] == 0, 0, 1)
}

netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))

names(netlists) <- names(rawnets)

# args = commandArgs(trailingOnly=TRUE)
# 
# cat(args, '\n')
# 
# args <- as.numeric(args)


ind <- c('functional complementarity',
         'weighted NODF', 'modularity')


# I'm doing this locally as the cluster is down, so use some trickery to run
# parallel loops


# Setup parallel loops ----------------------------------------------------
# code modified from https://www.blasbenito.com/post/02_parallelizing_loops_with_r/

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()


start_time <- Sys.time()
print(start_time)

foreach(parallel_i = 1:1000) %dopar%{
  for(i in 1:length(ind)){
    chosen_ind <- ind[i]
    print(chosen_ind)
    
    real_and_errors <- LOTUS::randomized_ranges(netlists, indices = chosen_ind, 
                                         network_level = 'higher', 
                                         out_format = 'list', 
                                         summarise = F, 
                                         actual_vals = F, n_perm = 1,
                                         modularity_nperm = 1000
    )
    
    # reformat the object, as its output as a list which doesn't like
    # being saved as a csv
    temp <- as.data.frame(real_and_errors)
    z <- data.frame(meta = colnames(temp), vals = as.numeric(temp[1,]))
    
    outname <- paste('data/output_data/for_z_scores/', chosen_ind, 
                     '_', parallel_i, '_', i,
                     '.csv', sep = '')
    write.csv(z, outname)
    
    cat(outname, 'written')
}


  
 
}

end_time <- Sys.time()
print(end_time)

print(end_time - start_time)
parallel::stopCluster(cl = my.cluster)