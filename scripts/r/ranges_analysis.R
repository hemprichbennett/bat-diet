##################################################
## Project: bat diet
## Script purpose: Analysing how 95% confidence intervals of individual values change over clustering thresholds used
## Date: 04/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: Uses the output of the script 'array_95_confidence_intervals_inter_network.R'
##################################################

library(ggplot2)

files <- list.files(path = 'data/output_data/randomized_ranges/', pattern ='.csv')
files <- paste('data/output_data/randomized_ranges/', files, sep ='')
file_list <- lapply(files, read.csv, row.names= 1)

ranges_df <- do.call(rbind, file_list)
colnames(ranges_df)[c(4,5)] <- c('lower', 'upper')

p_list <- list()
for(i in 1: length(unique(ranges_df$metric))){
  met <- unique(ranges_df$metric)[i]
  pdf(paste('plots/randomized_ranges/', met, '.pdf', sep =''))
  p_list[[met]] <- ggplot(ranges_df[which(ranges_df$metric==met),], aes(clustering, actual, colour = network))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper, colour = network), width=.1)+
    theme_bw()
  print(p_list[[met]])
  dev.off()
}
