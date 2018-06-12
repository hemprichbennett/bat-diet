##################################################
## Project: All bats
## Script purpose: Calculating individual-based metrics and taxonomic composition for individual bats
## Date: 12/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: This is a modified version of the all_species_individual_analysis script, here the analysis looks at the sites with years combined,
## rather than splitting the metaweb by site AND year
##################################################
if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(ggridges)
  library(reshape2)
  library(forcats)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
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
  
}
names(sites_list) <- unique(names(locations))

names(sites_list)

#####Do some ecology ####
all_ecology <- matrix(nrow=0, ncol=1+ncol(specieslevel(matrix(sample(0: 1, size =100, replace = T), nrow = 10, ncol = 10), level = 'higher')))

for(i in 1:length(sites_list)){
  starttime <- Sys.time()
  if(length(which(duplicated(colnames(sites_list[[i]]))))>0){
    sites_list[[i]] <- sites_list[[i]][,-which(duplicated(colnames(sites_list[[i]])))]
  }
  splevel = specieslevel(sites_list[[i]], level = 'higher')
  endtime <- Sys.time()
  cat(names(sites_list)[i],'took', endtime-starttime,'\n')
  all_ecology <- rbind(all_ecology, cbind(rep(names(sites_list)[i], nrow(splevel)), splevel))
}
all_ecology <- cbind(rownames(all_ecology), all_ecology)
#write.csv(all_ecology, 'all_ecology_1.csv')
#all_ecology <- read.csv('all_ecology_1.csv')
#####Add some taxonomic information
colnames(all_ecology)[c(1:2)] <- c('Sample','Site')
t_taxa <- t(taxa_mat)
t_taxa <- cbind(t_taxa, rownames(t_taxa))
t_taxa <- as.data.frame(t_taxa)
colnames(t_taxa)[ncol(t_taxa)] <- 'Sample_no'

all_ecology <- merge(all_ecology, t_taxa, by.x='Sample' , by.y= 'Sample_no')
#write.csv(all_ecology, 'all_ecology_2.csv')
#all_ecology <- read.csv('all_ecology_2.csv')
###Important: this doesn't include EVERY Hice sample, as some either have too many (or more likely too few) nodes for the analysis

all_ecology <- merge(x = all_ecology, y = field_data, by.x = 'Sample', by.y = 'Faeces_no1')

#write.csv(all_ecology, 'all_ecology_3.csv')
#####Work out the nestedness of each bat, then add it to the df#####


all_ecology <- cbind(all_ecology, rep(NA, nrow(all_ecology)))
colnames(all_ecology)[ncol(all_ecology)] <- 'nestedrank'
for(n in 1:length(sites_list)){
  nest <- nestedrank(sites_list[[n]])[['higher level']]
  for(i in 1:length(nest)){
    if(names(nest)[i] %in% all_ecology[,1]){
      val <- nest[i]
      print(val)
      pos <- which(all_ecology[,1]==names(nest)[i])
      print(pos)
      all_ecology[pos, ncol(all_ecology)] <- val
    }
  }
}
write.csv(all_ecology, 'data/output_data/all_bats/sitewise_all_individual_info.csv')
save.image('data/output_data/sitewise_all_individual_info.RDS')
#####Local work ####
load('data/output_data/sitewise_all_individual_info.RDS')


#Look at the occurence of taxa in each species

tax_df <- all_ecology[,c(1, 58, 74, seq(23,39))]
tax_df <- melt(tax_df, id.vars = c('Sample', 'Species', 'SiteAndYear'))
colnames(tax_df)[c(4,5)] <- c('Order', 'Present/absent')
tax_df$`Present/absent` <- as.integer(tax_df$`Present/absent`)
tax_df$`Present/absent` <- ifelse(tax_df$`Present/absent`== 0, 0, 1)

prop_present <- sapply(unique(tax_df[,c('Species', 'SiteAndYear', 'Order')]),  function(x) as.character(x))
prop <- c()
nbats <- c()
#for(i in 1: 1){
for(i in 1: nrow(prop_present)){
  tem <- tax_df[which(tax_df$Species==as.character(prop_present[i,1]) & tax_df$SiteAndYear== as.character(prop_present[i,2])
                      & tax_df$Order==as.character(prop_present[i,3])),]
  prop <- c(prop, sum(tem$`Present/absent`)/nrow(tem)) #The number of bats that consumed the order, divided by total bats
  nbats <- c(nbats, nrow(tem))
}


prop_present <- cbind(prop_present, prop, nbats)
prop_present <- as.data.frame(prop_present)
prop_present$prop <- as.numeric(as.character(prop_present$prop))
prop_present$nbats <- as.integer(as.character(prop_present$nbats))

prop_present$SiteAndYear <- gsub('DANUM', 'Danum', prop_present$SiteAndYear)
prop_present$SiteAndYear <- gsub('DVCA', 'Danum', prop_present$SiteAndYear)
prop_present$SiteAndYear <- gsub('MALIAU', 'Maliau', prop_present$SiteAndYear)
prop_present$Species <- gsub('Hice', 'Hipposideros cervinus', prop_present$Species)
prop_present$Species <- gsub('Hiri', 'Hipposideros ridleyi', prop_present$Species)
prop_present$Species <- gsub('Hidi', 'Hipposideros diadema', prop_present$Species)
prop_present$Species <- gsub('Hidy', 'Hipposideros dyacorum', prop_present$Species)
prop_present$Species <- gsub('Kein', 'Kerivoula intermedia', prop_present$Species)
prop_present$Species <- gsub('Keha', 'Kerivoula hardwickii', prop_present$Species)
prop_present$Species <- gsub('Kepa', 'Kerivoula papillosa', prop_present$Species)
prop_present$Species <- gsub('Rhbo', 'Rhinolophus borneensis', prop_present$Species)
prop_present$Species <- gsub('Rhse', 'Rhinolophus sedulus', prop_present$Species)
prop_present$Species <- gsub('Rhtr', 'Rhinolophus trifoliatus', prop_present$Species)


balloons <- ggplot(data = prop_present[which(prop_present$nbats >9),], aes(y = fct_rev(Order), x =SiteAndYear)) + geom_point(aes(size=prop, colour = prop))+ 
  scale_colour_gradient2(low = 'white',  high = 'steelblue')+
  theme(panel.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill='Proportion of MOTU present',
       x ="Site and year", y = 'Prey taxa')+
  facet_wrap(~Species)
balloons

tiles <- ggplot(data = prop_present[which(prop_present$nbats >5),], aes(y = fct_rev(Order), x =SiteAndYear)) + geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient(low = "white",
                      high = "black", name = 'Proportion') + # #3e0e4c works well
  labs(fill='Proportion of MOTU present',
       x ="Site and year", y = 'Prey taxa')+
  theme(panel.background=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Species)
tiles

pdf('plots/Site comparisons/sitewise_proportion_of_bats_containing.pdf')    
tiles
dev.off()
  

