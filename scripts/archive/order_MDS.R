#### Header ####
## Project:
## Script purpose: Using the data from all_species_sitewise_individual_analysis.R to look at species
##                and echolocation guild effects on diet
## Date: 26/07/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(vegan)

load('data/output_data/sitewise_all_individual_info.RDS')

interactions <- all_ecology[,c(23:38, 58)]

interactions[,-ncol(interactions)] <- apply(interactions[,-ncol(interactions)], 2, as.numeric)
#interactions$Species <- all_ecology[,58]

#Sum the interactions by species
interactions <- aggregate(interactions[,seq(1, ncol(interactions)-1)], by=list(Category=interactions$Species), FUN=sum)
rownames(interactions) <- interactions[,1]
interactions <- interactions[,-1]

diet_MDS <- metaMDS(interactions, # Our community-by-species matrix
                                 k=2) # The number of reduced dimensions


plot(diet_MDS)

#make a vector of bat echolocation guilds
echo <- rep('LDC', length(rownames(interactions)))
echo[which(grepl('Hi', rownames(interactions)))] <- 'HDC'
echo[which(grepl('Rh', rownames(interactions)))] <- 'HDC'
#colnames(interactions)[1] <- 'Bat_species'
#rownames(interactions) <- interactions[,17]


plot(diet_MDS, type = 'n')
ordihull(diet_MDS,groups=echo,draw="polygon",col="grey90",label=F)
orditorp(diet_MDS,display="sites",
         air=0.01,cex=1.25)

pdf('plots/order_mds.pdf')
plot(diet_MDS, type = 'n')
ordihull(diet_MDS,groups=echo,draw="polygon",col="grey90",label=F)
orditorp(diet_MDS,display="sites",
         air=0.01,cex=1.25)
title('order')
dev.off()