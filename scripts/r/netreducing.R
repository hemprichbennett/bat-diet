##################################################
## Project: All bats
## Script purpose: Calculating how network metrics change when we adjust the number of bats in a network
## Date: 15/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: 
##################################################
if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(ggridges)
  library(reshape2)
  library(forcats)
  library(netReducer)
  library(DataExplorer)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(netReducer, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'), stringsAsFactors = F)
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data$Site <- gsub('DVCA', 'Danum', field_data$Site)
field_data$Site <- gsub('DANUM', 'Danum', field_data$Site)
field_data$Site <- gsub('MALIAU', 'Maliau', field_data$Site)

all_interactions <- r_network_gen(collapse_species = F, desired_species = NULL, include_malua = T, filter_species = T, lulu = T)

desired_cols <- c('MOTU', 'DANUM', 'MALIAU', 'SAFE', 'SBE')

all_interactions <- all_interactions[,-which(!all_interactions[1,] %in% desired_cols)]


colnames(all_interactions) <- all_interactions[2,]
all_interactions <- all_interactions[-c(2),]
rownames(all_interactions) <- all_interactions[,1]
all_interactions <- all_interactions[,-1]

locations <- c()

for(i in 1:ncol(all_interactions)){
  locations[i] <- all_interactions[1,i]
  names(locations)[i] <- all_interactions[1,i]
}


row_names <- rownames(all_interactions) #Store these as an object as the apply below kills them
all_interactions <- apply(all_interactions, 2, as.numeric)
rownames(all_interactions) <- row_names

prey_data <- read.csv('data/taxonomy/order.csv')
colnames(prey_data) <- c('MOTU', 'Taxa')
prey_data$MOTU <- as.character(prey_data$MOTU)




###Add the taxonomic information to the interactions for everything
taxa_mat <- matrix(nrow = 0, ncol = ncol(all_interactions))
colnames(taxa_mat) <- colnames(all_interactions)
z <- 1
for(i in 1: nrow(all_interactions)){   ####THIS NEEDS SOME SERIOUS EDITING BEFORE IT'LL WORK
  rowname = rownames(all_interactions)[i]
  if(!rowname %in% prey_data$MOTU){
    next()
  }
  tax = as.character(prey_data[which(prey_data$MOTU == rowname),'Taxa'])
  if(is.null(nrow(taxa_mat))){ #If its the first iteration there won't be any rownames yet, so the next if statement will fail
    taxa_mat <- rbind(taxa_mat, as.numeric(all_interactions[i,]))
    rownames(taxa_mat)[z] <- tax
    z <- z+1
  }
  if(tax %in% rownames(taxa_mat)){
    to_merge = which(rownames(taxa_mat)==tax)
    taxa_mat[to_merge,] <- taxa_mat[to_merge,]+ as.numeric(all_interactions[i,])
  }else{
    taxa_mat <- rbind(taxa_mat, as.numeric(all_interactions[i,]))
    rownames(taxa_mat)[z] <- tax
    z <- z+1
  }
}

#####Make a list with a network for each site####

sites_list <- list()





for(i in 1:length(unique(names(locations)))){
  loc = unique(names(locations))[i]
  sites_list[[i]] <- all_interactions[,which(names(locations)==loc)]
  sites_list[[i]] <- sites_list[[i]][-1,]#Now that we have the site, this row is worthless
  
}
names(sites_list) <- unique(names(locations))

names(sites_list)

sites_list <- sites_list[-which(names(sites_list)=='SBE')] #Theres no point doing this with SBE as it only has 1 species

#Make a vector of sample IDs and species
sampleIDs <- c(field_data$Species, field_data$Species)
names(sampleIDs) <- c(field_data$Faeces_no1, field_data$Faeces_no2)
which(colnames(sites_list$SAFE) %in% names(sampleIDs))

#Rename all the columns to species-site
for(n in 1:length(sites_list)){
  badcols <- c()
  for(i in 1:ncol(sites_list[[n]])){
    pos <- which(names(sampleIDs)==colnames(sites_list[[n]])[i])
    pos <- pos[1]
    #print(colnames(sites_list[[n]])[i])
    if(is.na(pos)){
      badcols <- c(badcols, i)
      next()
    }
    #print(paste(sampleIDs[pos], names(sampleIDs[pos]), sep = '-'))
    colnames(sites_list[[n]])[i] <- paste(sampleIDs[pos], names(sampleIDs[pos]), sep = '-')
  }
  if(length(badcols)>0){
    sites_list[[n]] <- sites_list[[n]][,-badcols]
  }
}



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
         'ISA', 'modularity')

chosen_ind <- ind[args]

#Now do the analysis
out_df <- netreducing(input = sites_list, input_type = 'list', n_iterations = 100, min_nodes = 40, metric_chosen = 'connectance',
                      type_chosen = 'network', level = 'higher')


write.csv(out_df, paste('results/rarifying_networks/reducing_', chosen_ind,'.csv', sep =''))

