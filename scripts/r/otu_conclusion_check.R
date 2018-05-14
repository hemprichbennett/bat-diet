#setwd()

#Make a nested list, where the main list has a value for each clustering level, and each item within it is a distinct network generated at that level #####

library(here)
library(bipartite)
library(ggplot2)
library(LOTUS)
r_network_gen <- function(input_network, collapse_species = T, desired_species = NULL, filter_species = F, include_malua=F, lulu= F){
  
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
  
  source('scripts/r/The.matrix.reloader.R')
  source('scripts/r/hernani_comparisons.R')
  
  
  
  #####Reading in data, formatting it for analysis ####
  
  all_interactions <- #read.table('data/processed_dna_data/25_april_strict_lengths/for_r/95/all_post_QC_otus.txt.table_binary.out', sep = '\t', header = F, stringsAsFactors = F)# 
    all_interactions <- input_network
  
  desired_colnames <- c("Rhbo","Kein","Hice","Kepa","Rhse","Hiri","Rhtr", "Keha", "Hidi", "Hidy") #If we want to filter out species
  
  
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
  
  badcols <- c('1774','4437', '2070', '2275', '4260', '4531', "1004", "1007", "1107", "1134", "1165", "1180",  "198",  "209",  "210","387",  "426",  "459",  "497",  "536",  "541",  "567",  "591",  "689","796",  "806",  "822",  "841",  "843",  "899",  "910",  "918",  "986","996", "3712", "4341", "4361",'1774','4437', '2070', '2275', '4260', '4531', '841', '843')#Sadly these columns match two different samples, so must be removed for now until checked against the field data
  
  all_interactions <- all_interactions[,-which(all_interactions[1,] %in% badcols)]
  
  field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
  
  field_data$Site <- gsub('DVCA', 'DANUM', field_data$Site)
  field_data$Site <- gsub('MALUA', 'SBE', field_data$Site)
  field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
  field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
  
  field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
  
  #badsites <- c('DV88, 2016', 'MALUA, 2016', 'DV89, 2016') #These ones won't make viable networks for various reasons
  if(include_malua==F){
    badsites <- c('DV88, 2016', 'SBE, 2016', 'DV89, 2016') #These ones won't make viable networks for various reasons
  }else{
    badsites <- c('DV88, 2016', 'DV89, 2016') #These ones won't make viable networks for various reasons
  }
  
  field_data <- field_data[-which(field_data$SiteAndYear %in% badsites),]
  
  if(collapse_species==T){
    
    nets <- lapply(unique(field_data$SiteAndYear), function(i) the.matrix.reloader(master.data = field_data, ID.column.1 = "Faeces_no1", ID.column.2 = 'Faeces_no2', species.column = "Species", split.by.column = "SiteAndYear", split.by.var = i, OTU.matrix = all_interactions))
    
    names(nets) <- unique(field_data$SiteAndYear)
    for(i in 1:length(nets)){
      print(names(nets)[i])
      print(colnames(nets[[i]]))
    }
    
    
    if(filter_species==T){
      for(i in 1: length(nets)){
        if(length(which(!colnames(nets[[i]]) %in% desired_colnames))>0){#If there are species that we don't want
          to_remove <- which(!colnames(nets[[i]]) %in% desired_colnames)
          nets[[i]] <- nets[[i]][,-to_remove]
        }
      }
    }
    
    
    
    return(nets)
  }
  else if(collapse_species==F){
    all_interactions <- rbind(all_interactions[1,], all_interactions)
    if(is.null(desired_species)){
      if(filter_species==F){
        for(i in 1:ncol(all_interactions)){
          if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no1){
            row <- field_data[which(field_data$Faeces_no1==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
          }else if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no2){
            row <- field_data[which(field_data$Faeces_no2==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
          }  
        }
      }else if(filter_species==T){
        badcols <- c()
        for(i in 1:ncol(all_interactions)){
          if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no1){
            row <- field_data[which(field_data$Faeces_no1==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
            if(!row$Species %in% desired_colnames){
              badcols <- c(badcols,i)
            }
          }else if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no2){
            row <- field_data[which(field_data$Faeces_no2==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
            if(!row$Species %in% desired_colnames){
              badcols <- c(badcols,i)
            }
          }  
        }
        all_interactions <- all_interactions[,-badcols]
      }
      
      
      # if(filter_species==T){
      #   
      #   for(i in 1: length(nets)){
      #     if(length(which(!colnames(nets[[i]]) %in% desired_colnames))>0){#If there are species that we don't want
      #       to_remove <- which(!colnames(nets[[i]]) %in% desired_colnames)
      #       nets[[i]] <- nets[[i]][,-to_remove]
      #     }
      #   }
      # }
      return(all_interactions)
    }else{
      all_interactions_with_extra <- rbind(all_interactions[1,], all_interactions)
      all_interactions_with_extra <- the.matrix.reloader(master.data = field_data, ID.column.1 = "Faeces_no1", ID.column.2 = 'Faeces_no2',species.column = "SiteAndYear", split.by.column = "Species", split.by.var = desired_species, OTU.matrix = all_interactions_with_extra, collapse_top_species = F)
      return(all_interactions_with_extra)
    }
    
    
  }
  
}


setwd(here())
getwd()

filenames <- list.files(pattern = 'all_post_QC_otus.txt.table_binary.out', recursive = T)

filenames <- filenames[-grep('galaxy', filenames)]
filenames <- filenames[-grep('lulu', filenames)]

rawnets <- lapply(filenames, read.table,  sep = '\t', header = F, stringsAsFactors = F)
names(rawnets) <- gsub('.*\\/', '', dirname(filenames))
netlists <- lapply(rawnets, r_network_gen, collapse_species = T, filter_species = T)
names(netlists) <- names(rawnets)


#Specify which index(s) to use

ind <- desired_mets <- c('functional complementarity',
                         'web asymmetry',
                         'Alatalo interaction evenness',
                         'togetherness',
                         'Fisher alpha', 'mean number of shared partners',
                         'niche overlap',
                         'nestedness',
                         'discrepancy',
                         'ISA')


m <- metcalcs(lis= netlists, indexes = ind, netlevel = 'higher')

ggplot(m, aes(x = clustering, y = value, color = dataset)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  geom_vline(xintercept = 93)+ #Vertical line
  geom_vline(xintercept = 97)+ #Vertical line
  facet_wrap(~ metric, scales = 'free_y')

#####The le comber function ######
line_plot <- function(df, dataset, clustering, metric, value, plotname = NULL){
  rankings_mat <- matrix(nrow = length(unique(df$dataset)), ncol = length(unique(df$clustering)))
  colnames(rankings_mat) <- unique(df$clustering)
  ms <- as.character(c())
  clusts <- c()
  ns <- c()
  for(i in 1:length(unique(df$metric))){
    #for(i in 1:1){
    
    met <- unique(df$metric)[i]
    metric_subset <- df[which(df$metric==met),]
    #Make the appropriate subset of the data to play with
    for(a in 1:length(unique(metric_subset$clustering))){
      clust <- unique(metric_subset$clustering)[a]
      #print(clust)
      metric_and_cluster_subset <- metric_subset[which(metric_subset$clustering==clust),]
      rankings_mat[,a] <- metric_and_cluster_subset[order(metric_and_cluster_subset$value),'dataset']
    }
    
    for(b in 1:ncol(rankings_mat)){#Starting from 2 as we have to calculate the similarity to the previous col
      cl <- as.numeric(colnames(rankings_mat)[b])
      if(b==1){##We need to do this step as otherwise our counting backwards will crash things: we want to see how similar the values are to the value before them, which is confusing for the first value in the loop
        n <- 7
      }else{
        n <- length(which(rankings_mat[,(b-1)]==rankings_mat[,b]))  
      }
      
      ns <- c(ns,n)
      
      clusts <- c(clusts, cl)
      ms <- c(ms, as.character(met))
    }
  }
  
  out_df <- data.frame(ms, clusts, ns)
  #print(out_df)
  #return(out_df)
  out_df <- out_df[order(out_df$ms, decreasing = TRUE),]
  
  # plot
  # set x limits
  xmin <- 91
  xmax <- 98
  my_rows <- length(unique(out_df$ms))
  my_rows
  #pdf('../Figures/reliable_range_bars.pdf')
  # empty plot
  par(mar=c(5,12,2.5,3))
  plot(1,type="n",xlim=c(xmin,xmax),ylim=c(0,my_rows+1),axes=TRUE,xlab="clustering (%)", yaxt="n", ylab="",frame=FALSE)
  title(main = plotname)
  axis(2, at=1:length(unique(out_df$ms)), labels=unique(out_df$ms), las = 1, cex.axis=0.65)
  #for(a in 1:2){
  for(a in 1:length(unique(out_df$ms))){
    met <- unique(out_df$ms)[a]
    #print(met)
    connectance_df <- out_df[which(out_df$ms == met),]
    connectance_df$ms <- NULL
    # dummy data
    matches <- sample(0:2,100,replace=TRUE)
    thresholds <- seq(90,100,length=100)
    cbind(thresholds,matches)
    
    
    
    # which entries do we wish to count?
    test <- rep(0,nrow(connectance_df))
    #df <- data.frame(thresholds,matches)
    test[which(connectance_df$ns>6)] <- 1
    df <- data.frame(connectance_df,test)
    
    #add dummy entries to ensure we pick up the correct start and end points
    df <- rbind(c(0,0,0),df)
    df <- rbind(df,c(0,0,0))
    
    # look for starts and stops, ensuring we pick up single points
    start_stop <- rep(NA,nrow(df))
    for(i in 2:nrow(df)){
      # start points
      ifelse(df$test[i]==1 & df$test[i-1]==0, start_stop[i] <- "start", start_stop[i] <- "NA")
    }
    # end points
    for(i in 1:(nrow(df)-1))
    {
      if(df$test[i]==1 & df$test[i+1]==0) start_stop[i] <- "stop"
    }
    # single points
    for(i in 2:(nrow(df)))
    {
      if(df$test[i]==1 & df$test[i+1]==0 & df$test[i-1]==0) start_stop[i] <- "single"
    }
    #For the crap bits  
    for(i in 1:(nrow(df)-1))
    {
      if(df$test[i]==0) start_stop[i] <- "crap"
    }  
    
    # add to dataframe
    df <- data.frame(df,start_stop)
    
    # extract start and end points of each sequence
    line_starts <- subset(df,df$start_stop=="start")[,1]
    line_stops <- subset(df,df$start_stop=="stop")[,1]
    bad_points <- subset(df,df$start_stop=="crap")[,1]
    single_points <- subset(df,df$start_stop=='single')[,1]
    line_ends <- cbind(line_starts,line_stops)
    line_ends <- rbind(line_ends, cbind(bad_points, bad_points))
    line_ends <- rbind(line_ends, cbind(single_points, single_points))
    if(0.0 %in% line_ends[,1]){
      line_ends <- line_ends[-which(line_ends[,1]==0.0),] # Remove the zero values which we'd put in earlier when generating df
    }
    
    line_ends
    
    #These rows aren't in order yet, which would allow us to plot but would mess up the colour scheme. The below line sorts that
    if(!is.null(nrow(line_ends))){
      line_ends <- line_ends[order(line_ends[,1]),]
    } #We need the if statement as some lines are nice and have one value throughout, but that means we cannot order their rows as they don't really have any
    
    
    
    #I've now made the empty plot before initiating the loop
    # guideline
    #lines(matrix(c(xmin,xmax,a,a),ncol=2,byrow=FALSE),lwd=0.5,col="red")
    
    # add lines showing matches
    if(is.null(nrow(line_ends))){
      lines(matrix(c(line_ends[1],line_ends[2],a,a),ncol=2,byrow=FALSE),lwd=3,col='black')
    }
    else if(nrow(line_ends>0)){ #Some of the metrics have zero lines, as they're so utterly shit. 
      #The if statement lets us skip them, as otherwise they crash it
      
      for(i in 1:nrow(line_ends))
      {
        lines(matrix(c(line_ends[i,1],line_ends[i,2],a,a),ncol=2,byrow=FALSE),lwd=3,col='black')
        #print(i)
      }
      
    }
    
  }
  abline(v=c(93,97),lty=2,col="gray")
  
}

line_plot(df = m, dataset = 'dataset', clustering = 'clustering', metric = 'metric', value = 'value', plotname = 'Sabah, non-lulu')

#####ANOVA ####

aov(aov(value ~ clustering + dataset + clustering:dataset, data = m[m$metric=='connectance',]))
