#### Header ####
## Project:
## Script purpose: making sense of the plan_d output
## Date: 12/09/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


library(here)
library(reshape2)
library(ggplot2)
setwd(here())

indir <- 'results/rand_values/'

allfiles <- list.files(path = indir)

actuals <- allfiles[grep('actual', allfiles)]
rands <- allfiles[grep('rand', allfiles)]

# actuals <- read.csv('results/rand_values_95/actual_functional complementarity.csv')
# rand <- read.csv('results/rand_values_95/rand_functional complementarity.csv')

mets <- gsub('rand_', '', rands)
mets <- gsub('\\.csv', '', mets)
#mets <- gsub('.+_', '', mets) #I'm keeping the number in the metrics for the time being, the loop will have more iterations but potentially be less unwieldy



metrics <- c()
clusts <- c()
net1 <- c()
net2 <- c()
net1_actual <- c()
net2_actual <- c()
lower <- c()
upper <- c()
actual <- c()
signif <- c()


for(a in 1: length(mets)){
  met <- mets[a]
  clust <- gsub('_.+', '', met)
  #print(met)
  actual_to_load <- actuals[grep(met, actuals)]
  #print(actual_to_load)
  rand_to_load <- rands[grep(met, rands)]
  
  actual_file <- read.csv(paste(indir, actual_to_load, sep = ''))
  rand_file <- read.csv(paste(indir, rand_to_load, sep = ''))
  
  combns <- t(combn(colnames(rand_file), 2))
  
  for(i in 1:nrow(combns)){
    rand_range <- quantile(rand_file[,combns[i,1]] - rand_file[,combns[i,2]], probs = c(0.025, 0.975))
    #Order it
    rand_range <- rand_range[order(rand_range)]
    actual_diff <- actual_file[,combns[i,1]] - actual_file[,combns[i,2]]
    
    if(actual_diff > rand_range[1] && actual_diff < rand_range[2]){
      
      sig <- F
    }else if(actual_diff == mean(rand_range)){
      sig <- F
    }else{
      
      sig <- T}
    
    metrics <- c(metrics, gsub('.+_', '', met))
    clusts <- c(clusts, clust)
    net1 <- c(net1, combns[i,1])
    net2 <- c(net2, combns[i,2])
    net1_actual <- c(net1_actual, actual_file[,combns[i,1]])
    net2_actual <- c(net2_actual, actual_file[,combns[i,2]])
    lower <- c(lower, rand_range[1])
    upper <- c(upper, rand_range[2])
    actual <- c(actual, actual_diff)
    signif <- c(signif, sig)
    
  }
  
  
}

#This makes a dataframe of the pairwise combinations for each clustering level, and if they're significantly differernt or not
out_df <- data.frame(metrics, clusts, net1, net2, net1_actual, net2_actual, signif, lower, upper, actual)

mets_to_keep <- c('discrepancy', 'mean number of shared partners', "niche overlap", "NODF", 'functional complementarity', 'modularity')

out_df <- out_df[which(out_df$metrics %in% mets_to_keep),]


colnames(out_df) <- c('Metric', 'Clustering', 'Network 1', 'Network 2', 'Network 1 actual value', 'Network 2 actual value', 'Significant',
                      'Lower random difference', 'Higher random difference', 'Measured difference')

out_df[,c(5,6,8,9,10)] <- round(out_df[,c(5,6,8,9,10)], digits = 3)

out_df <- out_df[order(out_df$Metric),]

out_df <- out_df[,c(1,2,3,5,4,6,10,7,8,9)]
sig_df <- out_df[which(out_df$Significant==T),]
sig_df

write.csv(sig_df, 'results/final_random/significant.csv')
write.csv(out_df, 'results/final_random/all.csv')

