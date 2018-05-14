#This script rarifies networks of cervinus to look at individual-level changes in network topology
if(interactive()==TRUE){
  library(here)
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(reshape2)
  library(lme4)
  library(plyr)
  #library(multcompView)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

field_data <- read.csv('data/Edited_all_files_with_lat_long_VKedits.csv')
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data_2 <- field_data
#The bit below deals with the inevitable merging issues with some samples being from faeces_no1, some from faeces_no2
field_data$Faeces_no2 <- NULL

field_data_2$Faeces_no1 <- field_data_2$Faeces_no2
field_data_2$Faeces_no2 <- NULL

field_data <- rbind(field_data, field_data_2)
all_interactions <- r_network_gen(collapse_species = F, desired_species = 'Hice', include_malua = T, lulu=T)

#####Make a list with a network for each site####

sites_list <- list()
sites <- unique(colnames(all_interactions))

c_names <- c()


for(i in 1:length(sites)){
  m =  all_interactions[,which(colnames(all_interactions )==sites[i])]
  colnames(m) = m[1,]
  m = m[-c(1,2),]
  m <- m[-which(rowSums(m)==0),]
  print(ncol(m))
  c_names <- c(c_names, colnames(m))
  sites_list[[i]] <- m
}
print(sites)
names(sites_list) <- sites
names(sites_list) <- gsub('DANUM', 'Danum', names(sites_list))
names(sites_list) <- gsub('MALIAU', 'Maliau', names(sites_list))
names(sites_list) <- gsub('MALUA', 'SBE', names(sites_list))



#####Rarify #####

desired_mets <- c('functional complementarity',
                  'web asymmetry',
                  'Alatalo interaction evenness',
                  'togetherness',
                  'Fisher alpha', 'mean number of shared partners',
                  'niche overlap',
                  'nestedness',
                  'discrepancy',
                  'ISA',
                  'robustness')


net_outmat <- matrix(nrow= 0, ncol = 6)



for(i in 1:length(sites_list)){
  site <- names(sites_list)[i]
  for(p in 20: ncol(sites_list[[i]])){
    for(r in 1:100){
      net <- sites_list[[i]][,sample(seq(1, ncol(sites_list[[i]])), replace = F, size = p)]
      netlevel <- networklevel(net, level = 'higher', index = desired_mets)
      n_MOTU <- length(which(rowSums(net)>0))
      out <- cbind(site, p , n_MOTU, r, names(netlevel), netlevel)
      print(out)
      net_outmat <- rbind(net_outmat, out)
    }
    
  }
  
}
colnames(net_outmat) <- c('Network', 'N_bats', 'N_MOTU', 'iteration', 'metric', 'value')

save.image('data/output_data/hice_stats/rarified_network_stats.RDS')
#write.csv(net_outmat, 'shiny/hice_rarifying/hice_network_rarified.csv')

load('data/output_data/hice_stats/rarified_network_stats.RDS')
#As I fucked up and made it a matrix not a data frame, this reformats it
net_outdf <- as.data.frame(net_outmat)
net_outdf$N_bats <- as.numeric(as.character(net_outdf$N_bats))
net_outdf$N_MOTU <- as.numeric(as.character(net_outdf$N_MOTU))
net_outdf$iteration <- as.numeric(as.character(net_outdf$iteration))
net_outdf$value <- as.numeric(as.character(net_outdf$value))

sink('data/output_data/hice_stats/full_network_lme.txt')
for(i in 1:length(unique(net_outdf$metric))){
  met <- as.character(unique(net_outdf$metric)[i])
  print(met)
  subs = net_outdf[which(net_outdf$metric==met),]
  
  mod <- lmer(value ~ N_bats + (1|Network) , data=subs)
  print(summary(mod))
}

sink()

#####Plot with ggplot#####


#net_outdf <- net_outdf[,-c(4,9)]#I don't care about the bat ID or iteration for this bit, it just gets in the way of the melting

# net_outdf$Network <- gsub('SAFE, ', 'SAFE,\n', net_outdf$Network)
# net_outdf$Network <- gsub('Danum, ', 'Danum,\n', net_outdf$Network)
# net_outdf$Network <- gsub('Maliau, ', 'Maliau,\n', net_outdf$Network)
# net_outdf$Network <- gsub('SBE, ', 'SBE,\n', net_outdf$Network)

#net_outdf <- melt(net_outdf, id.vars = c('Network', 'N_bats', 'N_MOTU'))

net_outdf$metric <- gsub('web asymmetry', 'Web\nasymmetry', net_outdf$metric)
net_outdf$metric <- gsub('Alatalo interaction evenness', 'Alatalo\ninteraction\nevenness', net_outdf$metric)
net_outdf$metric <- gsub('discrepancy.HL', 'Discrepancy\nhigher', net_outdf$metric)
net_outdf$metric <- gsub('Fisher alpha', 'Fisher\nalpha', net_outdf$metric)
net_outdf$metric <- gsub('togeherness.HL', 'Togetherness\nhigher', net_outdf$metric)
net_outdf$metric <- gsub('functional.complementarity.higher', 'Functional\ncomplementarity\nhigher', net_outdf$metric)
net_outdf$metric <- gsub('mean.number.of.shared.partners.HL', 'Mean\nnumber\nof shared\npartners', net_outdf$metric)
net_outdf$metric <- gsub('interaction strength asymmetry', 'Interaction\nstrength\nasymmetry', net_outdf$metric)
net_outdf$metric <- gsub('niche.overlap.HL', 'Niche\noverlap', net_outdf$metric)
net_outdf$metric <- gsub('togetherness.HL', 'Togetherness', net_outdf$metric)
net_outdf$metric <- gsub('nestedness', 'Nestedness', net_outdf$metric)
net_outdf$metric <- gsub('functional.complementarity.HL', 'Functional\ncomplementarity', net_outdf$metric)
# net_outdf$variable <- gsub('normalised degree', 'Normalised degree', net_outdf$variable)
# net_outdf$variable <- gsub('^degree$', 'Degree', net_outdf$variable)
# net_outdf$variable <- gsub('resource range', 'Resource range', net_outdf$variable)
# net_outdf$variable <- gsub('proportional similarity', 'Proportional similarity', net_outdf$variable)

#Reorder the factors of the network column to make nicer plots
net_outdf$Network <- ordered(net_outdf$Network, levels = unique(net_outdf$Network)[order(unique(as.character(net_outdf$Network)))])
#temp <- net_outdf[1:400000,]

batfacet <- ggplot(net_outdf, aes(N_bats, value)) + geom_bin2d(bins =70)+
  scale_fill_gradient(low = "blue",
                      high = "red", name = 'Count') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.title.y = element_blank())+ labs(x='Number of bats in network')+
  facet_grid(metric ~ Network, scale = 'free_y', switch = 'y') 
batfacet
pdf('plots/Hice/network_rarifying_bats_facetted.pdf')
batfacet
dev.off()  


motufacet <- ggplot(net_outdf, aes(N_MOTU, value)) + geom_bin2d(bins =70)+
  scale_fill_gradient(low = "blue",
                      high = "red", name = 'Count') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.title.y = element_blank())+ labs(x='Number of MOTU in network')+
  facet_grid(metric ~ Network, scale = 'free_y', switch = 'y') 
motufacet
pdf('plots/Hice/network_rarifying_motu_facetted.pdf')
motufacet 
dev.off()


