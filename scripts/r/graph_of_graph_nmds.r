#Script for plotting a sprung graph

library('here')
library(ggplot2)
library(reshape2)
library(qgraph)
library(igraph)


setwd(here())
source('scripts/r/r_network_gen.r')

nets <- r_network_gen(collapse_species = T, filter_species = T, include_malua = F)

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))

#####Generating all of the actual values for the datasets ####

desired_mets <- c('functional complementarity',
                  'web asymmetry',
                  'Alatalo interaction evenness',
                  'togetherness',
                  'Fisher alpha', 'mean number of shared partners',
                  'niche overlap',
                  'nestedness',
                  'discrepancy',
                  'ISA')



sets <- c()
mets <- c()
vals <- c()
for(n in 1: length(nets)){
  for(m in 1:length(desired_mets)){
    sets <- c(sets, names(nets)[n])
    st <- Sys.time()
    val <- networklevel(nets[[n]], index = desired_mets[m])[1]
    vals <- c(vals, val) #Have to select the first value as some metrics return multiple metrics and break everything
    mets <- c(mets, names(val))
    en <- Sys.time()
    print(val)
    cat(names(nets)[n], desired_mets[m], 'finished, it took ', en-st, '\n')
  }
}

realstats <- data.frame(sets, mets, vals)

realstats$site <- rep(NA, nrow(realstats))
realstats$year <- rep(NA, nrow(realstats))

for(i in 1:nrow(realstats)){
  s <- strsplit(as.character(realstats$sets[i]),', ', fixed = T)
  realstats$site[i] <- s[[1]][1]
  realstats$year[i] <- s[[1]][2]
}


#####plot all metrics in multidimensional space ####
#First normalise them
realstats$scaled <- rep(NA, nrow(realstats))
for(i in 1:length(unique(realstats$mets))){
  # realstats$scaled[realstats$mets==unique(realstats$mets[i])] <- scale(realstats$vals[realstats$mets==unique(realstats$mets[i])], center = F)
  
  vals <- realstats$vals[realstats$mets==unique(realstats$mets[i])]
  realstats$scaled[realstats$mets==unique(realstats$mets[i])] <- (vals - mean(vals)) / sd(vals)+2
  
  
  print(mean(realstats$scaled[realstats$mets==unique(realstats$mets[i])]))
  print(sd(realstats$scaled[realstats$mets==unique(realstats$mets[i])]))
  
} #This gives us a mean and standard deviation of 2. I do this to ensure they're all positive for the NMDS





#make the normalised values into a matrix for the NMDS
scaled_matrix <- dcast(realstats, sets ~ mets)
rownames(scaled_matrix) <- scaled_matrix$sets
scaled_matrix$sets <- NULL

MDS <- metaMDS(scaled_matrix,k=2,trymax=100)


#pdf('plots/Site comparisons/MDS_vegan.pdf')
ordiplot(MDS,type="n")
orditorp(MDS,display="sites",cex=0.7,air=0.01)
ordicluster(MDS,hclust(vegdist(scaled_matrix,"bray")))
#dev.off()
#Or the same sort of thing using ggplot
MDS1 <- MDS$points[,1]
MDS2 <- MDS$points[,2]

NMDS_coords = data.frame(MDS1 = MDS1, MDS2 = MDS2)
NMDS_coords$network <- rownames(NMDS_coords)

g_NMDS <- ggplot(NMDS_coords, aes(x=MDS1, y=MDS2, col=network)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "NMDS Plot")

g_NMDS
pdf('plots/Site comparisons/NMDS_ggplot.pdf')
g_NMDS
dev.off()


dist_list <- melt(as.matrix(dist(NMDS_coords[,-3])))
dist_list <- dist_list[-which(dist_list$X1==dist_list$X2),] #The above line gives us self-matches, which we don't want
dist_list$value <- 1-dist_list$value
#make sure you invert these values somehow later






net_igraph <- graph.data.frame(dist_list, directed = F)
V(net_igraph)$type <- V(net_igraph)$name %in% dist_list[,1]


# 
# #Make colours vector
# node_cul <- rep("yellow", length(V(net_igraph)))
# for(i in 1: length(V(net_igraph))){
#   if(length(grep('denovo',V(net_igraph)[[i]]$name))>0){ #If the node is an OTU
#     node_cul[i] <- "#01c2cd"
#   }
# }

plot(net_igraph, vertex.size = 10, vertex.label=NA)

layout_with_graphopt(net_igraph)

original_dists <- as.matrix(dist(NMDS_coords[,-3]))

for_graph <-  1-original_dists


write.csv( original_dists, 'data/output_data/site_comparisons/NMDS_dists.csv')

a <- qgraph(for_graph, layout = "spring", label.cex = 2)

pdf('plots/Site comparisons/net_of_nets.pdf')
qgraph(for_graph, layout = "spring", label.cex = 2, directed = F)
dev.off()


png('plots/Site comparisons/net_of_nets.png')
qgraph(for_graph, layout = "spring", label.cex = 2, directed = F)
dev.off()


melted_dists <- melt(original_dists)
colnames(melted_dists) <- c('Network 1', 'Network 2', 'value')
tiles <- ggplot(melted_dists, aes(`Network 1`, `Network 2`)) + geom_tile(aes(fill=value), colour = 'white')+
  scale_fill_gradient(low = "white",
                       high = "black", name = 'Distance') + # #3e0e4c works well
  labs(x = '', y = '')+ 
  theme(panel.background = element_blank())
tiles

pdf('plots/Site comparisons/net_tiles.pdf')
tiles
dev.off()

  

