
if(interactive()==TRUE){
  library('here')
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}
setwd(here())
getwd()
source('scripts/r/r_network_gen.r')

nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T)



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
  V(bat_igraph)$type <- V(bat_igraph)$name %in% molten_net_2[,1]
  
  igraph_list[[n]] <- bat_igraph
  
  #Make colours vector
  node_cul <- rep("yellow", length(V(bat_igraph)))
  for(i in 1: length(V(bat_igraph))){
    if(length(grep('denovo',V(bat_igraph)[[i]]$name))>0){ #If the node is an OTU
      node_cul[i] <- "#01c2cd"
    }
  }
  pdf(paste('plots/Site comparisons/', names(nets)[n], '.pdf', sep =''))
  plot(bat_igraph, vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[n])
  dev.off()
}

par(mfrow= c(2,3))


pdf('plots/Site comparisons/all_nets.pdf')
par(mfrow= c(3,3))
plot(igraph_list[[2]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[2])
plot(igraph_list[[6]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[6])
plot(igraph_list[[4]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[4])
plot(igraph_list[[7]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[7])
plot(igraph_list[[1]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[1])
plot(igraph_list[[3]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[3])
plot(igraph_list[[5]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[5])
dev.off()

png('plots/Site comparisons/all_nets.png')
par(mfrow= c(3,3))
plot(igraph_list[[2]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[2])
plot(igraph_list[[6]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[6])
plot(igraph_list[[4]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[4])
plot(igraph_list[[7]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[7])
plot(igraph_list[[1]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[1])
plot(igraph_list[[3]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[3])
plot(igraph_list[[5]], vertex.size = 10, vertex.color=node_cul, vertex.label=NA, main = names(nets)[5])
dev.off()
