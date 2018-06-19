##################################################
## Project: Bat diet
## Script purpose: Showing the 'full' values for three key network metrics
## Date: 19/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(here)
library(ggplot2)
library(reshape2)

setwd(here())
source('scripts/r/r_network_gen.r')

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots


sites_list <- r_network_gen(filter_species = T)
names(sites_list) <- gsub('DANUM', 'Danum', names(sites_list))
names(sites_list) <- gsub('MALIAU', 'Maliau', names(sites_list))

ind <- c('nestedness', 'niche overlap', 'functional complementarity')

original_mets <- sapply(ind, function(i) sapply(sites_list, function(x) networklevel(x, index = i, level = 'higher')))

melted_original <- melt(original_mets)
colnames(melted_original) <- c('Site', 'Metric', 'value')
melted_original$Site <- gsub('\\..+', '', melted_original$Site)
melted_original$habitat_type <- rep('Primary', nrow(melted_original))
melted_original$habitat_type[grepl('SAFE', melted_original$Site)] <- 'Logged'
melted_original$Metric <- firstup(as.character(melted_original$Metric))


network_vals <- ggplot(melted_original, aes(x=Site, y=value, colour = habitat_type))+ geom_point() + facet_wrap(~ Metric, scales = 'free')+ 
  theme_dark() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c('blue', 'yellow'))
network_vals

pdf('plots/Site comparisons/network_vals.pdf')
network_vals
dev.off()
