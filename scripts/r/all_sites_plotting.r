
library(here)
library(magrittr)
library(bipartite)
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

pdf('plots/Site comparisons/all_bipartite.pdf', width = 10)

par(mfrow = c(2,2))
for(i in 1:length(nets)){
  #a bit of trickery to get the loop to do the plots alphabetically
  j <- order(names(nets))[i]
  
  plotweb(nets[[j]], high.lablength = '0', low.lablength = '0')
  title(names(nets)[j])
}
dev.off()