####This script makes dietary barplots, and runs a chi-squared test to look at the difference between Hice diet
####Between networks

library(stringr)
library('here')
library('bipartite')
library(iNEXT)
library(tidyverse)
library(reshape2)
source(here('scripts/r/iNEXT_prep.R'))
# source('scripts/r/hernani_comparisons.R')

setwd(here())
getwd()
#####Load in the data, reformat it to make it more useful ####

#all_interactions <- read.table(here('data/processed_dna_data/galaxy_r_workflow/95/all_post_QC_otus.txt.table_binary.out'), sep = '\t', header = F, stringsAsFactors = F, row.names = 1)
all_interactions <- read.csv('data/processed_dna_data/lulu/95/lulu_95.csv', header = F, stringsAsFactors = F)#
all_interactions[1,1] <- 'MOTU'
all_interactions[1,] <- gsub('X', '', all_interactions[1,])
all_interactions[2:nrow(all_interactions),2:ncol(all_interactions)] <- ifelse(all_interactions[2:nrow(all_interactions), 2:ncol(all_interactions)]==0,0,1)


#This finds all samples with 'GC' in the name and gives them a useful name
gc <- grep('GC',all_interactions[1,])
for(i in 1:length(gc)){
  #print(str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7])
  temp <- str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7]
  all_interactions[1,gc[i]] <- str_split(temp, pattern='_')[[1]][1]
}

#This finds all samples without 'GC' in the name and gives them a useful name
non_gc <- seq(1, ncol(all_interactions))[-gc]
str_split(all_interactions[1,non_gc[2]], pattern = '\\.')[[1]][1] #Finish this
for(i in 1:length(non_gc)){
  all_interactions[1,non_gc[i]] <- str_split(all_interactions[1,non_gc[i]], pattern = '\\.')[[1]][1]
}

badcols <- c('1774','4437', '2070', '2275')#Sadly these columns match two different samples, so must be removed for now until checked against the field data

all_interactions <- all_interactions[,-which(all_interactions[1,] %in% badcols)]


field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
field_data$Site <- gsub('DVCA', 'Danum', field_data$Site)
field_data$Site <- gsub('DANUM', 'Danum', field_data$Site)
field_data$Site <- gsub('MALIAU', 'Maliau', field_data$Site)
field_data$Site <- gsub('MALUA', 'SBE', field_data$Site)
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)

#####This loop filters for hice and makes a list entry for each siteandyear ####
batnames <- c()
siteslist <- list()
for(i in 1:ncol(all_interactions)){
  batname <- all_interactions[1,i]
  if(!batname %in%field_data$Faeces_no1 && !batname %in% field_data$Faeces_no2){
    next()
  }
  if(batname %in% field_data$Faeces_no1){
    site <- field_data[which(field_data$Faeces_no1==batname), 'SiteAndYear']
    if(length(which(field_data$Faeces_no1==batname))>1){next()}#If theres more than one match
    if(field_data[which(field_data$Faeces_no1==batname), 'Species']!= 'Hice'){#Only allow hice into the upcoming list
      next()
    }
    
  }else if(batname %in% field_data$Faeces_no2){
    site <- field_data[which(field_data$Faeces_no2==batname), 'SiteAndYear']
    if(length(which(field_data$Faeces_no2==batname))>1){next()}
    if(field_data[which(field_data$Faeces_no2==batname), 'Species']!= 'Hice'){#Only allow hice into the upcoming list
      next()
    }
  
    }
  if(site %in% names(siteslist)){
    pos <- which(names(siteslist) == site)
    siteslist[[pos]] <- cbind(siteslist[[pos]], as.numeric(all_interactions[,i]))
  }else{
    siteslist[[site]] <- matrix(nrow = nrow(all_interactions), ncol = 1, as.numeric(all_interactions[,i]))
    rownames(siteslist[[site]]) <- all_interactions[,1]
  }
}


#The bat names are still in row 1
for(i in 1:length(siteslist)){
  colnames(siteslist[[i]]) <- siteslist[[i]][1,]
  siteslist[[i]] <- siteslist[[i]][-1,]
  print(names(siteslist)[i])
  print(ncol(siteslist[[i]]))
}
#####Add the taxonomic information, format the data #####
orderdata <- read.csv('data/taxonomy/order.csv', stringsAsFactors = F)
familydata <- read.csv('data/taxonomy/family.csv', stringsAsFactors = F)

iNEXT_making <- function(prey_data, sites_list){
  colnames(prey_data) <- c('MOTU', 'Taxa')
  taxa_list <- list()
  
  #Put the taxonomic information onto the sites data
  
  for(b in 1: length(siteslist)){
    all_interactions <- siteslist[[b]]
    taxa_mat <- matrix(nrow=0, ncol=ncol(all_interactions))
    
    z <- 1
    for(i in 1: nrow(all_interactions)){ #Make a matrix of simplified interactions
      rowname = rownames(all_interactions)[i]
      if(!rowname %in% prey_data$MOTU){
        next()
      }
      tax = as.character(prey_data[which(prey_data$MOTU == rowname),'Taxa'])
      if(is.null(nrow(taxa_mat))){ #If its the first iteration there won't be any rownames yet, so the next if statement will fail
        taxa_mat <- rbind(taxa_mat, all_interactions[i,])
        rownames(taxa_mat)[z] <- tax
        z <- z+1
      }
      if(tax %in% rownames(taxa_mat)){
        to_merge = which(rownames(taxa_mat)==tax)
        taxa_mat[to_merge,] <- taxa_mat[to_merge,]+ all_interactions[i,]
      }else{
        taxa_mat <- rbind(taxa_mat, all_interactions[i,])
        rownames(taxa_mat)[z] <- tax
        z <- z+1
      }
    }
    taxa_list[[b]] <- taxa_mat
    names(taxa_list)[b] <- names(siteslist)[b]
  }
  return(taxa_list)
 #  out_list <- list() #For now as I'm not using this code for iNEXT, this isn't needed
 #  for(i in 1:length(taxa_list)){ #Simplify the data for iNEXT
 #    nam <- names(taxa_list)[i]
 #    #out_list[[nam]] <- c(ncol(taxa_list[[i]]), rowSums(taxa_list[[i]]))
 #    out_list[[nam]] <-  rowSums(taxa_list[[i]])
 #    #names(out_list[[nam]]) <- NULL
 #    #out_list[[nam]] <- as.incfreq(taxa_list[[i]])
 #  }
 # return(out_list) 
}

orderlist <- iNEXT_making(prey_data = orderdata, sites_list = siteslist)
prop_list <- list()
for(i in 1:length(orderlist)){
  prop_list[[i]] <- data.frame(rowSums(orderlist[[i]])/sum(orderlist[[i]]), rownames(orderlist[[i]]))
}
names(prop_list) <- names(orderlist)



#####Format the data a bit for the chi-squared etc ####
molten_proportions <- melt(prop_list )
molten_proportions <- molten_proportions[,-2] #This column sucks, get rid of it
colnames(molten_proportions) <- c('prey_taxa', 'value', 'SiteAndYear')

#Annoyingly having melted the data we must now cast it to form our contingency table
hice_wide <- acast(molten_proportions, SiteAndYear ~ prey_taxa)
hice_wide[is.na(hice_wide)] <- 0
hice_wide <- t(hice_wide) #Has to be the other way round for the chi-squared test
chisq <- chisq.test(hice_wide)
chisq$expected
#now try it without the rarest 13 taxa
#badtaxa <- names(sort(rowSums(hice_temp)))[1:13]
#chisq.test(hice_wide[-which(rownames(hice_temp) %in% badtaxa),])

#This is the proportion of the OTUs made up by the most common 4 orders
#colSums(hice_wide[-which(rownames(hice_temp) %in% badtaxa),])

#These are the commonest 4 orders
#rownames(hice_wide)[which(!rownames(hice_wide) %in% badtaxa)]



#####Format the data for plotting ####
#Now get the number of samples for use in plotting
s <- c()
counts <- c()
for(i in 1:length(orderlist)){
  s <- c(s, names(orderlist)[i])
  counts <- c(counts, ncol(orderlist[[i]]))
}
count_data <- data.frame(s,counts)

molten_proportions <- merge(molten_proportions, count_data, by.x = 'SiteAndYear', by.y = 's')
molten_proportions[which(molten_proportions$prey_taxa !='Araneae'),'counts'] <- '' #We only want numbers on the first ones


#####Barplot ####
palette <- c("#d05034",
             "#7f60d2",
             "#a6b63d",
             "#c053c0",
             "#5ab150",
             "#c14489",
             "#4cae90",
             "#e24971",
             "#76813a",
             "#667ad2",
             "#d49b3b",
             "#5d9bd3",
             "#9f612d",
             "#835897",
             "#e29473",
             "#d58dc9",
             "#b05560")

molten_proportions <- molten_proportions[-which(molten_proportions$value==0),] #Get rid of the empty rows or they mess up the legend
hice_plot <- ggplot(data = molten_proportions, aes(x = SiteAndYear, y =value, fill = prey_taxa)) + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+ #Get rid of the annoying background formatting, rotate the x labels
  geom_bar(stat = 'identity') + scale_fill_manual(values =palette, name = 'Prey order')+
  xlab("Site and year")+ ylab('Proportion')+ 
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,0.2))+
  geom_text(aes(label=counts), vjust=-0.3, size=3.5, position = position_stack(vjust = 1))

 
hice_plot
pdf('plots/Hice/hice_diet_bars.pdf', width =10, height = 7)
hice_plot
dev.off()
 
png('plots/Hice/hice_diet_bars.png', width = 10, height = 7)
hice_plot
dev.off()


png('~/Desktop/temp.png')
hice_plot
dev.off()

#####Trying out balloon plotting #####
balloon_data <- molten_proportions
colnames(balloon_data)[3] <- 'Proportion of\nMOTU present'
balloon_data$prey_taxa <-  fct_rev(balloon_data$prey_taxa)

balloons <- ggplot(data = balloon_data, aes(x = SiteAndYear, y =prey_taxa)) + geom_point(aes(size=`Proportion of\nMOTU present`))+
  theme(panel.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill='Proportion of MOTU present',
       x ="Site and year", y = 'Prey taxa')
balloons
pdf('plots/Hice/diet_balloons.pdf')
balloons
dev.off()

png('plots/Hice/diet_balloons.png')
balloons
dev.off()

