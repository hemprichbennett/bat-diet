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
library(reshape2)

source(here('scripts', 'r', 'r_network_gen.r'))

nets <- r_network_gen(lulu = T, filter_species = T)


#####Do the actual calculations ####
sp_mets <- lapply(nets, function(x) specieslevel(x, level = 'higher'))


####Format the data ####
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


####Now calculate the degree distribution of the networks####

#Make individual-based networks

ind_nets <- r_network_gen(collapse_species = F, lulu = T, filter_species = T)
ind_nets <- ind_nets[,-1] #These are the rownames but we dont need them

ind_nets <- ind_nets[,-which(!ind_nets[1,] %in% c('SAFE', 'DANUM', 'MALIAU'))]

ind_netlist <- list()

for(i in 1:length(unique(as.character(ind_nets[1,])))){
  site <- unique(as.character(ind_nets[1,]))[i]
  print(site)
  mat <- ind_nets[3:nrow(ind_nets),which(ind_nets[1,]==site)]
  ind_netlist[[site]] <- apply(mat, 2, as.numeric)
}

degree_dist <- melt(lapply(ind_netlist, function(x) colSums(x)))
colnames(degree_dist) <- c('Degree', 'Site')

ggplot(degree_dist, aes(x=Degree))+ geom_histogram(binwidth = 1)+ facet_wrap(~Site)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

fit <- aov(Degree ~ Site, degree_dist)

plot(fit)

TukeyHSD(fit)
  
