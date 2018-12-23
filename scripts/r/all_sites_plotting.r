
library(here)
library(magrittr)  
setwd(here())
getwd()
source('scripts/r/r_network_gen.r')

nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T)

names(nets) %<>%
  gsub('DANUM', 'Danum', .)%<>%
  gsub('MALIAU', 'Maliau', .)
  

for(i in 1:length(nets)){
  colnames(nets[[i]]) %<>%
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
  
}



#####This tiny chunk shows how many MOTU are shared
samples <- lapply(nets, rownames)

Reduce(intersect, samples)

#or as a %
(length(Reduce(intersect, samples))/16148)*100

#####Yay plotting time ####

igraph_list <- list()

for(n in 1:length(nets)){
  net <- nets[[n]]
  molten_net <- melt(net)
  molten_net_2 <- molten_net[(which(molten_net[,3] != 0)),]
  bat_igraph <- graph.data.frame(molten_net_2, directed = F)
  #bat_igraph$layout <- layout_in_circle
  V(bat_igraph)$type <- V(bat_igraph)$name %in% molten_net_2[,1]
  igraph_list[[n]] <- bat_igraph
  
  #Make colours vector
  node_cul <- rep("yellow", length(V(bat_igraph)))
  for(i in 1: length(V(bat_igraph))){
    if(length(grep('denovo',V(bat_igraph)[[i]]$name))>0){ #If the node is an OTU
      node_cul[i] <- "#01c2cd"
    }
    
  }
  
  V(bat_igraph)$name <- gsub('denovo.+', '', V(bat_igraph)$name)
  pdf(paste('plots/Site comparisons/', names(nets)[n], '.pdf', sep =''))
  plot(bat_igraph, vertex.size = 10, vertex.color=node_cul, main = names(nets)[n], width = 0.01, label.dist = 1)
  dev.off()
}

par(mfrow= c(2,3))


pdf('plots/Site comparisons/all_nets.pdf')
par(mfrow= c(3,1))
plot(igraph_list[[1]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[1])
plot(igraph_list[[2]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[2])
plot(igraph_list[[3]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[3])
dev.off()

png('plots/Site comparisons/all_nets.png')
par(mfrow= c(3,1))
plot(igraph_list[[1]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[1])
plot(igraph_list[[2]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[2])
plot(igraph_list[[3]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[3])

dev.off()
