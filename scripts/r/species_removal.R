## @knitr species_removal
#### Header ####
## Project: bat-diet
## Script purpose: Finding the effect of removing different bat species from out networks
## Date: 21/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(ggplot2)
library(magrittr)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots



source('scripts/r/r_network_gen.r')

nets <- r_network_gen(lulu=T, filter_species = T)

ind <- c('functional complementarity',
         'web asymmetry',
         'Alatalo interaction evenness',
         'togetherness',
         'Fisher alpha', 'mean number of shared partners',
         'niche overlap',
         'nestedness',
         'discrepancy',
         'ISA', 'weighted nestedness', 'NODF', 'weighted NODF')

ind <- c('functional complementarity',
         'Fisher alpha', 'mean number of shared partners',
         'nestedness',
         'discrepancy',
         'ISA', 'weighted nestedness', 'NODF', 'weighted NODF')


orig <- lapply(ind, function(i) lapply(nets, function(x) bipartite::networklevel(x, index = i, level = 'higher')))


names(orig) <- ind




mods <- lapply(nets, function(x) slot(bipartite::computeModules(web = x), 'likelihood'))

mods <- melt(mods)
colnames(mods) <- c('value', 'Network')
mods$Metric <- 'modularity'
ind <- c(ind, 'modularity')

orig <- rbind(orig, mods)

orig_melted <- melt(orig)
colnames(orig_melted) <- c('value', 'Network', 'Metric')
orig_melted$minus_species <- "No species removed"


sp_deleter <- function(networks, chosen_index){
  outlist <- list()
  for(i in 1:length(networks)){
    if(chosen_index=='modularity'){
      a <- sapply(seq(1,ncol(networks[[i]])), function(a) slot(bipartite::computeModules(web =networks[[i]][,-a]), 'likelihood'))
    }else{
      a <- sapply(seq(1,ncol(networks[[i]])), function(a) bipartite::networklevel(networks[[i]][,-a], index = chosen_index, level = 'higher'))
    }
    
    a <- data.frame('value'= a, 'minus_species' =colnames(networks[[i]]))
    #print(a)
    outlist[[names(networks)[i]]] <- a
  }
  #cols <- lapply(networks, function(x) colnames(x))
  return(outlist)
}


#slot(bipartite::computeModules(web = nets$SAFE), 'likelihood')

#one_met <- sp_deleter(nets, chosen_index = 'nestedness')
#one_melted <- melt(one_met)

all_mets <- lapply(ind, function(i) sp_deleter(nets, chosen_index = i))
names(all_mets) <- ind
melted_all <- melt(all_mets)
melted_all <- melted_all[,-2]
colnames(melted_all) <- c('minus_species', 'value', 'Network', 'Metric')

#melted_all$Network <- gsub('DANUM', 'Danum', melted_all$Network)
#melted_all$Network <- gsub('MALIAU', 'Maliau', melted_all$Network)

#plot_str(all_mets)


master_df <- rbind(orig_melted, melted_all)
master_df$Network <- gsub('DANUM', 'Danum', master_df$Network)
master_df$Network <- gsub('MALIAU', 'Maliau', master_df$Network)

master_df$minus_species <- as.factor(master_df$minus_species)


master_df$Metric <- firstup(master_df$Metric)

palette <- c("#75aa56",
             "#8259b1",
             "#be7239")

master_df$minus_species %<>% 
  gsub('Hice', 'Hipposideros cervinus', .)%<>%
  gsub('Hidi', 'Hipposideros diadema', .)%<>%
  gsub('Hidy', 'Hipposideros dyacorum', .)%<>%
  gsub('Hiri', 'Hipposideros ridleyi', .)%<>%
  gsub('Keha', 'Kerivoula hardwickii', .)%<>%
  gsub('Kein', 'Kerivoula intermedia', .)%<>%
  gsub('Kepa', 'Kerivoula papillosa', .)%<>%
  gsub('Kepe', 'Kerivoula pellucida', .)%<>%
  gsub('Rhbo', 'Rhinolophus borneensis', .)%<>%
  gsub('Rhse', 'Rhinolophus sedulus', .)%<>%
  gsub('Rhtr', 'Rhinolophus trifoliatus', .)

master_df$minus_species <- as.factor(master_df$minus_species)

master_df$minus_species <- relevel(master_df$minus_species, "No species removed")#Makes 'none' the first factor level

master_df$Metric %<>% 
  gsub('ISA', 'Interaction\nstrength\nasymmetry', .)%<>%
  gsub('Mean number of shared partners', 'Mean number\nof shared\npartners', .)%<>%
  gsub('Weighted nestedness', 'Weighted\nnestedness', .)%<>%
  gsub('Functional complementarity', 'Functional\ncomplementarity', .)#%<>%


#master_df$Metric <- gsub(' ', '\n', master_df$Metric)

sp_plot <- ggplot(master_df, aes(x=minus_species, y = value, colour = Network))+
  geom_point(alpha=0.8)+ scale_color_manual(values=palette, name = 'Site')+
  facet_wrap( ~ Metric, scales = 'free')+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x ='Species removed', y = 'Metric value')
sp_plot

## @knitr species_removal_plot_saving

pdf('plots/species_removal.pdf', height = 11)
sp_plot
dev.off()

