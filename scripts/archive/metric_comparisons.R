##################################################
## Project: Bat diet
## Script purpose: Showing the 'full' values for three key network metrics, and randomizing them
## Date: 19/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


dir <- getwd()
basedir <- strsplit(dir, split ='/')[[1]][2]
print(basedir)
if(grepl('data', basedir)){
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  
}else{
  library('here')
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(reshape2)
  library(corrplot)
  library(iNEXT)
  library(DataExplorer)
}

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

start <- Sys.time()

randomized_met_ranges <-lapply(ind, function(i) lapply(sites_list, function(x) replicate(1000, bipartite::networklevel(vegan::permatfull(x, fixedmar='both',mtype="count",times=1)$perm[[1]],
                                                                                                                    index = i, level = 'higher'))))

end <- Sys.time()

cat('randomisations took ', end-start,'\n')

names(randomized_met_ranges) <- ind
random_met_ranges <- melt(randomized_met_ranges)
colnames(random_met_ranges) <- c('Value', 'Site', 'Metric')
random_met_ranges$`2.5%` <- rep(NA, nrow(random_met_ranges))
random_met_ranges$`97.5%` <- rep(NA, nrow(random_met_ranges))
random_met_ranges$original <- rep(NA, nrow(random_met_ranges))
combns <- unique(random_met_ranges[,c('Site', 'Metric')])

for(i in 1:nrow(combns)){
  s <- combns$Site[i]
  m <- combns$Metric[i]
  quants <- quantile(random_met_ranges$Value[which(random_met_ranges$Site==s & random_met_ranges$Metric==m)], probs = c(0.025, 0.975))
  random_met_ranges$`2.5%`[which(random_met_ranges$Site==s & random_met_ranges$Metric==m)] <- quants[1]
  random_met_ranges$`97.5%`[which(random_met_ranges$Site==s & random_met_ranges$Metric==m)] <- quants[2]
  orig <- melted_original$value[which(melted_original$Site==s & melted_original$Metric==m)]
  random_met_ranges$original[which(random_met_ranges$Site==s & random_met_ranges$Metric==m)] <- orig
}


write.csv(random_met_ranges, 'results/randomized_met_ranges.csv')

random_met_ranges <- read.csv('results/randomized_met_ranges.csv', row.names = 1)
colnames(random_met_ranges) <- gsub('X', '', colnames(random_met_ranges))
library(ggplot2)

palette <- c("#75aa56",
             "#8259b1",
             "#be7239")

nest <- ggplot(random_met_ranges[which(random_met_ranges$Metric=='nestedness'),], aes(x=Value, fill=Site, color=Site)) +
    geom_histogram(alpha=0.5, position="identity", binwidth = 0.02)+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x= 'Nestedness', y= 'Count')+
    scale_color_manual(values=palette)+
    scale_fill_manual(values=palette)+
    geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='nestedness'),], aes(xintercept =`2.5.`, color = Site), linetype='dashed')+
    geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='nestedness'),], aes(xintercept =`97.5.`, color = Site), linetype='dashed')+
    geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='nestedness'),], aes(xintercept =original, color = Site), linetype='solid')
nest

pdf('plots/randomized_ranges/hist_nestedness.pdf', width = 11)
nest
dev.off()

func <- ggplot(random_met_ranges[which(random_met_ranges$Metric=='functional complementarity'),], aes(x=Value, fill=Site, color=Site)) +
  geom_histogram(alpha=0.5, position="identity", binwidth = 1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x= 'Functional complementarity', y= 'Count')+
  scale_color_manual(values=palette)+
  scale_fill_manual(values=palette)+
  geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='functional complementarity'),], aes(xintercept =`2.5.`, color = Site), linetype='dashed')+
  geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='functional complementarity'),], aes(xintercept =`97.5.`, color = Site), linetype='dashed')+
  geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='functional complementarity'),], aes(xintercept =original, color = Site), linetype='solid')
func

pdf('plots/randomized_ranges/hist_func.pdf', width = 11)
func
dev.off()

niche <- ggplot(random_met_ranges[which(random_met_ranges$Metric=='niche overlap'),], aes(x=Value, fill=Site, color=Site)) +
  geom_histogram(alpha=0.5, position="identity", binwidth = 0.001)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x= 'Niche overlap', y= 'Count')+
  scale_color_manual(values=palette)+
  scale_fill_manual(values=palette)+
  geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='niche overlap'),], aes(xintercept =`2.5.`, color = Site), linetype='dashed')+
  geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='niche overlap'),], aes(xintercept =`97.5.`, color = Site), linetype='dashed')+
  geom_vline(data = random_met_ranges[which(random_met_ranges$Metric=='niche overlap'),], aes(xintercept =original, color = Site), linetype='solid')
niche

pdf('plots/randomized_ranges/hist_niche.pdf', width = 11)
niche
dev.off()
