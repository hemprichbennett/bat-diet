#### Header ####
## Project: Bat-diet
## Script purpose: Plotting how network-level metrics are altered by species removal, using 
## output from removing_and_randomizing.R
## Date: 25/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(ggplot2)
library(here)
library(magrittr)
library(reshape2)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

allfiles <- list.files(path = 'results/random_removals/', pattern = '.csv')
allfiles <- paste('results/random_removals/', allfiles, sep = '')

setwd(here())

filelist <- lapply(allfiles, read.csv, row.names= 1)

big_df <- do.call(rbind, filelist)

big_df$network %<>%
  gsub('DANUM', 'Danum', .)%<>%
  gsub('MALIAU', 'Maliau', .)

#Now make a df of columns of species counts

sp_info <- read.csv('data/output_data/all_bats/sitewise_all_individual_info.csv', stringsAsFactors = F, row.names = 1)

sp_table <- table(sp_info$Species, sp_info$Site)
sp_df <- melt(sp_table)
colnames(sp_df) <- c('Species', 'Site', 'N_samples')

#Combine the dataframes
sp_df$sp_and_site <- paste(sp_df$Species, sp_df$Site, sep = '_')
big_df$species_and_site <- paste(big_df$sp, big_df$network, sep = '_')

big_df <- merge(x =big_df, y =sp_df, by.x = 'species_and_site', by.y = 'sp_and_site')

big_df$sp %<>%
   gsub('Hice', 'Hipposideros cervinus', .)%<>%
   gsub('Hidi', 'Hipposideros diadema', .)%<>%
   gsub('Hidy', 'Hipposideros dyacorum', .)%<>%
   gsub('Hiri', 'Hipposideros ridleyi', .)%<>%
   gsub('Keha', 'Kerivoula hardwickii', .)%<>%
   gsub('Kein', 'Kerivoula intermedia', .)%<>%
   gsub('Kemi', 'Kerivoula minuta', .)%<>%
   gsub('Kepa', 'Kerivoula papillosa', .)%<>%
   gsub('Rhbo', 'Rhinolophus borneensis', .)%<>%
   gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
   gsub('Rhtr', 'Rhinolophus trifoliatus', .)

big_df$sp <- gsub(' ', '\n', big_df$sp)




#We now need to split our data frame up by metric in a loop, as the appropriate bin width
#varies by metric
plotlist <- list()
for(i in 1: length(unique(big_df$index_used))){
  met <- unique(big_df$index_used)[i]
  bins <- 1
  if(met=='functional complementarity'){
    bins <- 1
  }
  if(met =='niche overlap'){
    bins <- 0.004
  }
  if(met =='mean number of shared partners'){
    bins <- 0.5
  }
  if(met =='modularity'){
    bins <- 0.005
  }
  if(met=='nestedness'){
    bins <- 0.2
  }
  if(met=='weighted nestedness'){
    bins <- 0.01
  }
  if(met=='NODF'){
    bins <- 0.25
  }
  if(met=='weighted NODF'){
    bins <- 0.25
  }
  temp_big <- big_df[which(big_df$index_used==met),]
  big_actuals <- data.frame(temp_big$sp, temp_big$network, temp_big$actual, temp_big$lower, temp_big$upper)
  colnames(big_actuals) <- c('sp', 'network', 'actual', 'lower', 'upper')
  plotlist[[as.character(met)]] <- ggplot(data=temp_big, aes(rand_vals,  fill = N_samples)) + 
    geom_histogram(binwidth = bins)+ 
    scale_fill_gradientn(colours = c('grey', 'red', 'black'),
                        name = "Number of\nsamples")+
    geom_vline(aes(xintercept=actual), colour = 'blue', data = big_actuals)+
    geom_vline(aes(xintercept=lower), colour = 'black', data = big_actuals, linetype = 'dashed')+
    geom_vline(aes(xintercept=upper), colour = 'black', data = big_actuals, linetype = 'dashed')+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       strip.text.y = element_text(size = 6))+
    facet_grid(sp ~ network, scales = 'free_x') + xlab(firstup(as.character(met)))
  
  pdf(paste('plots/random_removals/', met, '.pdf', sep = ''), width = 11)
  print(plotlist[[as.character(met)]])
  dev.off()
  
}

#plotlist$`niche overlap`

ggplot(data=temp_big, aes(rand_vals,  fill = N_samples)) + 
  geom_histogram(binwidth = bins)+ 
  scale_fill_gradientn(colours = c('grey', 'red', 'black'))+
  geom_vline(aes(xintercept=actual), colour = 'red', data = big_actuals)+
  geom_vline(aes(xintercept=lower), colour = 'black', data = big_actuals, linetype = 'dashed')+
  geom_vline(aes(xintercept=upper), colour = 'black', data = big_actuals, linetype = 'dashed')+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(size = 6))+
  facet_grid(sp ~ network, scales = 'free_x') + xlab(firstup(as.character(met)))

