  if(collapse_species==T && !is.null(desired_species)){
    break('Cannot have false for collapsing species AND have a species desired for selection')
  }
  if(interactive()==TRUE){
    library('bipartite')
    library('stringr')
    library('igraph')
    library('reshape')
  }else{
    require(methods)
    library(network, lib.loc = '/data/home/btw863/r_packages/')
    library(statnet.common, lib.loc = '/data/home/btw863/r_packages/')
    library(sna, lib.loc = '/data/home/btw863/r_packages/')
    library(igraph, lib.loc = '/data/home/btw863/r_packages/')
    library(permute, lib.loc = '/data/home/btw863/r_packages/')
    library(vegan, lib.loc = '/data/home/btw863/r_packages/')
    library(bipartite, lib.loc = '/data/home/btw863/r_packages/')
    library(stringr, lib.loc = '/data/home/btw863/r_packages/')
    library(reshape, lib.loc = '/data/home/btw863/r_packages/')
  }
  
  source('scripts/R/The.matrix.reloader.R')
  source('scripts/R/hernani_comparisons.R')
  
  
  
  #####Reading in data, formatting it for analysis ####
  all_interactions <- read.table('data/processed_dna_data/galaxy_r_workflow/95/all_post_QC_otus.txt.table_binary.out', sep = '\t', header = F, stringsAsFactors = F)#
  

  
  
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
  
  field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'), stringsAsFactors = F)
  
  
  field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
  
  
  samples_1 <- which(field_data$Faeces_no1 %in% all_interactions[1,])
  samples_2 <- which(field_data$Faeces_no2 %in% all_interactions[1,])  
  
  
  all_good_samples <- c(samples_1, samples_2)
  
  s <- data.frame(all_interactions[1,], field_data[all_good_samples, 'SiteAndYear'])
  
  all_sequenced <- list.files('data/processed_dna_data/galaxy_r_workflow/')
  
  #This finds all samples with 'GC' in the name and gives them a useful name
  gc <- grep('GC',all_sequenced)
  for(i in 1:length(gc)){
    #print(str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7])
    temp <- str_split(all_sequenced[gc[i]], pattern = '\\.')[[1]][7]
    all_sequenced[gc[i]] <- str_split(temp, pattern='_')[[1]][1]
  }
  
  #This finds all samples without 'GC' in the name and gives them a useful name
  non_gc <- seq(1, length(all_sequenced))[-gc]
  str_split(all_sequenced[non_gc[2]], pattern = '\\.')[[1]][1] #Finish this
  for(i in 1:length(non_gc)){
    all_sequenced[non_gc[i]] <- str_split(all_sequenced[non_gc[i]], pattern = '\\.')[[1]][1]
  }
  
  sequenced_1 <- which(field_data$Faeces_no1 %in% all_interactions[1,])
  sequenced_2 <- which(field_data$Faeces_no2 %in% all_interactions[1,])  
  
  
  all_good_sequenced <- c(sequenced_1, sequenced_2)
  
  s <- data.frame(all_interactions[1,], field_data[all_good_samples, 'SiteAndYear'])
  unique(field_data[all_good_sequenced, 'SiteAndYear'])
  
  plan <- read.csv('~/Dropbox/Education/PhD/Labwork/DNA Extractions/2017_April_Dave_plates.csv')

  
  plated1 <- which(field_data$Faeces_no1 %in% plan$Sample.br.ID)
  plated2 <- which(field_data$Faeces_no2 %in% plan$Sample.br.ID)  
  plated <- c(plated1, plated2)
  unique(field_data[plated, 'SiteAndYear'])
  