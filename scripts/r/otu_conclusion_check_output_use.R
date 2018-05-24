#This script uses the output from the script otu_conclusion_check_generation.R
## @knitr setup
library(here())
library(LOTUS)
library(ggplot2)
library(gridExtra)
#setwd(here())

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots


#Set up the dataframe properly
m <- read.csv('data/output_data/all_bats/otu_conclusions.csv', row.names = 1)
m$site <- unlist(lapply(m$network, function(x) strsplit(as.character(x), split = ', ')[[1]][1]))
m$year <- unlist(lapply(m$network, function(x) strsplit(as.character(x), split = ', ')[[1]][2]))
m$site <- gsub('DANUM', 'Danum', m$site)
m$site <- gsub('MALIAU', 'Maliau', m$site)
m$habitat_type <- as.factor(ifelse(m$site == 'SAFE', 'Logged', 'Primary'))


sitescatter <- ggplot(m , aes(x = clustering, y = value, color = site)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  facet_wrap(~ firstup(gsub('\\.', ' ', metric)), scales = 'free_y')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#print(sitescatter)
#pdf('plots/Site comparisons/MOTU_sitescatter.pdf')
#sitescatter
#dev.off()

netscatter <- ggplot(m , aes(x = clustering, y = value, color = network)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  facet_wrap(~ firstup(gsub('\\.', ' ', metric)), scales = 'free_y')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#print(netscatter)
#pdf('plots/Site comparisons/MOTU_netscatter.pdf')
#netscatter
#dev.off()


habscatter <- ggplot(m , aes(x = clustering, y = value, color = habitat_type)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  facet_wrap(~ firstup(gsub('\\.', ' ', metric)), scales = 'free_y')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#print(habscatter)
#pdf('plots/Site comparisons/MOTU_habscatter.pdf')
#habscatter
#dev.off()
## @knitr all_scatterplots
grid.arrange(netscatter, sitescatter, habscatter, ncol = 2)
## @knitr lineplot
#motuline <- line_plot(input = m, network = 'network', clustering = 'clustering', metric = 'metric', value = 'value', plotname = 'Sabah')
#pdf('plots/Site comparisons/MOTU_line.pdf')
#line_plot(input = m, network = 'network', clustering = 'clustering', metric = 'metric', value = 'value', plotname = 'Sabah')
#dev.off()
line_plot(input = m, network = 'network', clustering = 'clustering', metric = 'metric', value = 'value', plotname = 'Sabah')
#####ANOVA ####
## @knitr anova
m$network <- as.factor(m$network)
aov(aov(value ~ clustering + network + clustering:network, data = m[m$metric=='web asymmetry',]))


## @knitr funct_comp_scatter
ggplot(m[m$metric=='functional.complementarity.HL',] , aes(x = clustering, y = value, color = network)) +
  geom_point()+
  labs(x = 'Clustering', y= 'Fucntional complementarity') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
