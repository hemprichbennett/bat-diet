#### Header ####
## Project: Bat-diet
## Script purpose: Calculating species-level metrics for each network
## Date: 25/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(here)
library(ggplot2)
library(magrittr)

source(here('scripts', 'r', 'r_network_gen.r'))

nets <- r_network_gen(lulu = T, filter_species = T)

sp_mets <- lapply(nets, function(x) specieslevel(x, level = 'higher'))

sp_df <- do.call(rbind, sp_mets)
sp_df$site <- gsub('\\..+', '', rownames(sp_df))
sp_df$species <- gsub('.+\\.', '', rownames(sp_df))

sp_df$species %<>%
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



write.csv(sp_df, 'results/species_level_data.csv')

sp_df <- read.csv(here('results', 'species_level_data.csv'))

ggplot(sp_df, aes(x=site, y = degree))+ geom_point()+ facet_wrap(~ species)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#Now calculate the degree distribution of the networks
degree_dist <- melt(lapply(nets, function(x) rowSums(x)))
colnames(degree_dist) <- c('Degree', 'Site')

ggplot(degree_dist, aes(x=Degree))+ geom_histogram(binwidth = 1)+ facet_wrap(~Site)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))





  
