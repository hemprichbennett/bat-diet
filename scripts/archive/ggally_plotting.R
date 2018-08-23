
if(interactive()==TRUE){
  library('here')
  library(GGally)
  library(network)
  library(sna)
  library(ggplot2)
  library(reshape2)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}
setwd(here())
getwd()
source('scripts/r/r_network_gen.r')

bip_init_network <- function (mymat, mode1="P", mode2="A") {
  require(network)
  #require(ggnet)
  if(!is.matrix(mymat)) mymat <- as.matrix(mymat)
  p<- dim(mymat)[1]    # Plants are rows
  a<- dim(mymat)[2]    # Animals are columns
  net<- network::network(mymat,
                         matrix.type = "bipartite",
                         ignore.eval = FALSE,
                         names.eval = "weights")
  net
  network::set.vertex.attribute(net,"mode",c(rep(mode1,p), rep(mode2,a)))
}

nets <- r_network_gen(collapse_species = T)

nets_edgelist <- lapply(nets, melt)
for(i in 1:length(nets_edgelist)){
  emptyrows <- which(nets_edgelist[[i]][,3]==0)
  nets_edgelist[[i]] <- nets_edgelist[[i]][-emptyrows,] #Get rid of all of the empty rows
  nets_edgelist[[i]] <- nets_edgelist[[i]][,-3] #Kill the column of weighting, we don't need it
  nets_edgelist[[i]][,1] <- as.character(nets_edgelist[[i]][,1])
  nets_edgelist[[i]][,2] <- as.character(nets_edgelist[[i]][,2])
}

a <- nets_edgelist[[1]]

#net <- bip_init_network(nets[[1]])


net <- network(a, directed = T)
net %v% 'node_type' = ifelse(grepl('denovo', net %v% 'vertex.names'), 'prey', 'bat')
#pdf('~/Desktop/temp.pdf')

ggnet2(net,  color = 'node_type', palette = 'Set1',
       size = "node_type", size.palette = c("bat" = 1, "prey" = 0.01))

