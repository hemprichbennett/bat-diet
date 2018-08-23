####Node sharing calculation

if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(reshape2)
  library(qgraph)
  
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T)

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))

clusts <- c()
inter_network <- c()
n_prey_nodes <- c()

samples <- lapply(nets, rownames)
#####For ggplot2 #####
net1 <- c()
net2 <- c()
overlappingnames <- c()
for(i in 1: nrow(t(combn(length(samples),2)))){
  n1 <- names(nets)[t(combn(length(samples),2))[i,1]]
  n2 <- names(nets)[t(combn(length(samples),2))[i,2]]
  net1 <- c(net1, n1, n2) #Adding them both as otherwise we only get half the values wanted for the plot
  net2 <- c(net2, n2, n1)
  overlap <- length(samples[[n1]][which(samples[[n1]] %in% samples[[n2]])])
  overlappingnames <- c(overlappingnames, overlap, overlap)
}
df <- data.frame(net1,net2, overlappingnames)

tiles <- ggplot(df, aes(net1, net2)) + geom_tile(aes(fill=overlappingnames), colour = 'white')+
  scale_fill_gradient(low = "white",
                      high = "black", name = 'Number of\nshared MOTU', limits=c(0, max(df$overlappingnames))) + # #3e0e4c works well
  labs(x = '', y = '')+ 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
tiles
##I'm not convinced this is the best way to do this: currently it looks at total of shared links, which unfairly
#disadvantages the networks with few nodes. Proportion of shared links could be better
pdf('plots/Site comparisons/node_sharing.pdf')
tiles
dev.off()

png('plots/Site comparisons/node_sharing.png')
tiles
dev.off()

#####For qgraph #####

q_net1 <- c()
q_net2 <- c()
o <- c()
for(i in 1: nrow(t(combn(length(samples),2)))){
  n1 <- names(nets)[t(combn(length(samples),2))[i,1]]
  n2 <- names(nets)[t(combn(length(samples),2))[i,2]]
  q_net1 <- c(q_net1, n1) #Adding them both as otherwise we only get half the values wanted for the plot
  q_net2 <- c(q_net2, n2)
  overlap <- length(samples[[n1]][which(samples[[n1]] %in% samples[[n2]])])
  o <- c(o, overlap)
}
qdf <- data.frame(q_net1,q_net2, o)

##This gives a spring plot where each node is a site, the nodes sharing most nodes are closest together,
#and the node size is the total number of OTUs the site has
qgraph(qdf, layout = "spring", label.cex = 2, directed = F, vsize = (sapply(nets, function(x) nrow(x))/100))


#Now we try to weight the spring graph by the number of potential otus that are shared
totals <- sapply(nets, function(x) nrow(x))
qdf_prop <- qdf

for(i in 1:nrow(qdf_prop)){
  n1 = as.numeric(totals[which(names(totals) == qdf_prop$q_net1[i])])
  n2 = as.numeric(totals[which(names(totals) == qdf_prop$q_net2[i])])
  qdf_prop$o[i] <- 1/(n1+n2)*qdf_prop$o[i]
}

qgraph(qdf_prop, layout = "spring", label.cex = 2, directed = F, vsize = (sapply(nets, function(x) nrow(x))/100))


##Now make a tile plot based on weighted proportions instead

df_prop <- df

for(i in 1:nrow(df_prop)){
  n1 = as.numeric(totals[which(names(totals) == df_prop$net1[i])])
  n2 = as.numeric(totals[which(names(totals) == df_prop$net2[i])])
  df_prop$overlappingnames[i] <- 1/(n1+n2)*df_prop$overlappingnames[i]
}

connectance_tiles <- ggplot(df_prop, aes(net1, net2)) + geom_tile(aes(fill=overlappingnames), colour = 'white')+
  scale_fill_gradient(low = "white",
                      high = "black", name = 'Inter-network\nconnectance', limits=c(0, max(df_prop$overlappingnames))) + # #3e0e4c works well
  labs(x = '', y = '')+ 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
connectance_tiles
