
#####beginning properly: laying out data #####

if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(ggridges)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))

field_data$Site <- gsub('DVCA', 'DANUM', field_data$Site)
field_data$Site <- gsub('MALUA', 'SBE', field_data$Site)
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)


#####Small section to find the names of MOTU which were found in more than one network#####
nets <- r_network_gen(collapse_species = T, filter_species = T)

samples <- lapply(nets, rownames)

#overlappingnames <- Reduce(intersect, samples) #This doesn't work, annoyingly and a bit embarresingly

overlappingnames <- c()
for(i in 1: nrow(t(combn(length(samples),2)))){
  net1 <- t(combn(length(samples),2))[i,1]
  net2 <- t(combn(length(samples),2))[i,2]
  overlappingnames <- c(overlappingnames, samples[[net1]][which(samples[[net1]] %in% samples[[net2]])])
}

##### organise our data for the plot ######
all_interactions <- r_network_gen(collapse_species = F, desired_species = NULL, include_malua = F, filter_species = T, lulu = T)

desired_cols <- c('MOTU', 'DANUM, 2016', 'DANUM, 2017', 'MALIAU, 2016', 'MALIAU, 2017', 'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017')
#desired_cols <- c('MOTU', 'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017')
all_interactions <- all_interactions[,-which(!all_interactions[1,] %in% desired_cols)]

#rownames(all_interactions) <- all_interactions[,1]
#all_interactions <- all_interactions[,-1]#we don't need the MOTU names

network_names <- unique(as.character(all_interactions[1,])) #For use later when getting groups
network_names <- network_names[-which(network_names=='MOTU')]
badcols <- c() #Get the species of each column, add it to the matrix
for(i in 2:ncol(all_interactions)){ #Starting at 2 as the first row still contains rownames
  if(all_interactions[2,i] %in% field_data$Faeces_no1){
    all_interactions[2,i] <- as.character(field_data$Species[which(field_data$Faeces_no1==all_interactions[2,i])])
  }else if(all_interactions[2,i] %in% field_data$Faeces_no2){
    all_interactions[2,i] <- as.character(field_data$Species[which(field_data$Faeces_no2==all_interactions[2,i])])
  }else{badcols <- c(badcols,i)}
}
if(length(badcols)>0){all_interactions <- all_interactions[,-badcols]}


colnames(all_interactions) <- paste(all_interactions[1,], all_interactions[2,], sep = '_')
all_interactions <- all_interactions[-c(1,2),]#Remove the ID bits, we don't need them now
rows <- all_interactions[,1]#Converting the matrix to numeric loses these
#all_interactions <- all_interactions[,-1]
all_interactions <- apply(all_interactions, 2, as.numeric) #Make the matrix numeric
rownames(all_interactions) <- rows



#####Merge together bats of the same species and site ####
for_graph <- matrix(ncol= 0, nrow = nrow(all_interactions))
for(i in 2:length(unique(colnames(all_interactions)))){ #start at 2 as we want to avoid the first column, if i remove it it fucks with the column names....
  col <- unique(colnames(all_interactions))[i]
  to_merge <- which(colnames(all_interactions)==col)
  if(length(to_merge)>1){
    for_graph <- cbind(for_graph, apply(all_interactions[,to_merge], 1, sum))
  }else{
    for_graph <- cbind(for_graph, all_interactions[,i])
  }
  
  colnames(for_graph)[i-1] <- col
}
rownames(for_graph) <- rownames(all_interactions)
#grep(colnames(for_graph), network_names[1])


which(colSums(for_graph)==0)

intra_network <- 0
for(i in 1: length(nets)){
  intra_network <- intra_network+ length(which(rowSums(nets[[i]])>1))
}

cat('number of MOTU found in more than one site is', length(unique(overlappingnames)))
cat('number of MOTU found in multiple bats per site is', intra_network)
#all_adjacency <- as.matrix(as.network.matrix(for_graph, matrix.type = 'incidence', bipartite = F))



#####Make the groups for the plotting function
# net_groups <- list()
# for(i in 1:length(network_names)){ 
#   #net_groups[[network_names[i]]] <- colnames(for_graph)[grep(colnames(for_graph), pattern = network_names[i])]
#   net_groups[[network_names[i]]] <- grep(colnames(for_graph), pattern = network_names[i])
# }



#Now do the same for the prey nodes
#WILL THE USE OF I AS AN IDENTIFIER IN THE FINAL LINE WORK, OR SHOULD IT BE ncol(for_graph) + i
# for(i in 1:nrow(for_graph)){
#   if(rownames(for_graph)[i] %in% overlappingnames){
#     next()
#   }
#   #print(i)
#   cols <- which(for_graph[i,]>0)
#   #print(cols)
#   s <- strsplit(names(cols)[1], split = '_')[[1]][1]
#   net_groups[[s]] <- c(net_groups[[s]], ncol(for_graph)+ i)
# }


#for_graph <- for_graph[-which(rowSums(for_graph)==0),]


melted_for_graph <- melt(for_graph)
melted_for_graph <- melted_for_graph[which(melted_for_graph$value>0),c(1,2)] #This returns only the valid edges, with no weighting as igraph doesn't use it in the conversion below
a <- as.matrix(get.adjacency(graph.data.frame(melted_for_graph))) #This gives us a square adjacency matrix instead of the usual rectangular incidence matrix

groupslist <- list()
for(i in 1:nrow(a)){
  if(grepl('denovo', rownames(a)[i])){
    if(rownames(a)[i] %in% overlappingnames){
      #next()
      
      nam <- strsplit(names(which(a[i,]>0)[1]), split = '_')[[1]][1]
      groupslist[[nam]] <- c(groupslist[[nam]], i)
    }else{
      if(length(strsplit(names(which(a[i,]>0)), split = '_'))>1){
        cat(rownames(a)[i], names(which(a[i,]>0)), '\n')
        nam <- strsplit(names(which(a[i,]>0)), split = '_')[[1]][1]
      }
      
      groupslist[[nam]] <- c(groupslist[[nam]], i)
    }
  }else{
    nam <- strsplit(rownames(a)[i], split ='_')[[1]][1]
    groupslist[[nam]] <- c(groupslist[[nam]], i)
  }
}
groupslist <- groupslist[order(names(groupslist))]

#which(rownames(a) %in% overlappingnames)
qgraph(a, groups = groupslist, legend = T, borders = FALSE, layout ='groups', labels = F, directed = F)


pdf('plots/Site comparisons/big_qgraph.pdf')
qgraph(a, groups = groupslist, legend = T, borders = FALSE, layout ='groups', minimum = 0.25, cut = 0.4, vsize = 1.5, labels = F, directed = F)
dev.off()


# 
# pdf('plots/Site comparisons/big_qgraph.pdf')
# qgraph(melted_for_graph, groups = net_groups, legend = T, borders = FALSE, layout ='groups', minimum = 0.25, cut = 0.4, vsize = 1.5, labels = F)
# dev.off()
