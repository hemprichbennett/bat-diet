#### Header ####
## Project: bat-diet
## Script purpose: using the output files from the netreducing array job to make plots FOR MODULARITY ONLY
## Date: 27/07/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(MASS)
library(reshape2)
library(ggplot2)
library(here)

setwd(here())

allfiles <- list.files(path = 'results/rarifying_networks/', pattern = '.csv')
allfiles <- paste('results/rarifying_networks/', allfiles, sep = '')

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots

modfiles <- allfiles[grep('modularity', allfiles)]

#I accidentally created 500 files when I only needed to make 100. Just use the first 100, otherwise reading in 6GB of data will fry the computer
modfiles <- modfiles[c(1:100)]

modlist <- lapply(modfiles, read.csv)

mod <- do.call(rbind, modlist)


mod <- mod[mod$included ==T,]

infilename <- 'modularity'

bigtax <- dcast(mod[,c(2,5,6,7,8,9)], n_used + netnames + metricval ~ Species, fun.aggregate = length)
bigtax$diversity <- sapply(seq(1,nrow(bigtax)), function(x) vegan::diversity(bigtax[x,seq(5, ncol(bigtax)),]))

longtax <- melt(bigtax, id.vars = c('netnames', 'n_used', 'metricval', 'diversity'))
colnames(longtax)[5] <- 'Species'


palette <- c("#75aa56",
             "#8259b1",
             "#be7239")

diversity_scatter <- ggplot(bigtax, aes(x = diversity, y = metricval, colour= netnames))+ 
  geom_point(alpha=0.3)+ scale_color_manual(values=palette, name = 'Site')+
  labs(x='Shannon diversity', y = firstup(infilename))+
  theme_bw() + theme(legend.position="bottom",
                     panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


diversity_scatter  
pdf(paste('plots/netreducing/rarifying_', infilename, '_diversity.pdf', sep = ''), width = 5)
print(diversity_scatter)
dev.off()

sp_scatter <- ggplot(longtax, aes(x = value, y = metricval, colour= netnames))+ 
  geom_point(alpha=0.3)+ scale_color_manual(values=palette, name = 'Site')+
  labs(x='Number of individuals', y = infilename)+
  theme_bw() + theme(legend.position="bottom",
                     panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~ Species, scales = 'free_x')

sp_scatter  
pdf(paste('plots/netreducing/rarifying_', infilename, '_sp.pdf', sep = ''), width = 5)
print(sp_scatter)
dev.off()



individuals_scatter <- ggplot(longtax, aes(x = n_used, y = metricval, colour= netnames))+ 
  geom_point(alpha=0.3)+ scale_color_manual(values=palette, name = 'Site')+
  labs(x='Number of individuals', y = firstup(infilename))+
  theme_bw() + theme(legend.position="bottom",
                     panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
individuals_scatter
pdf(paste('plots/netreducing/rarifying_', infilename, '_individuals.pdf', sep = ''), width = 5)
print(individuals_scatter)
dev.off()

jpeg(paste('plots/netreducing/rarifying_', infilename, '_individuals.jpeg', sep = ''), width = 1800, height =2400, res = 300)
print(individuals_scatter)
dev.off()

tiff(paste('plots/netreducing/rarifying_', infilename, '_individuals.tiff', sep = ''), width = 1800, height =2400, res = 500)
print(individuals_scatter)
dev.off()

tiff(paste('plots/netreducing/rarifying_', infilename, '_diversity.tiff', sep = ''), width = 1800, height =2400, res = 500)
print(diversity_scatter)
dev.off()

jpeg(paste('plots/netreducing/rarifying_', infilename, '_diversity.jpeg', sep = ''), width = 1800, height =2400, res = 300)
print(diversity_scatter)
dev.off()
