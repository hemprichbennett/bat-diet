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
weighted_m <- read.csv('data/output_data/all_bats/weighted_otu_conclusions.csv', row.names = 1)

m$source <- 'frequency'
weighted_m$source <- 'weighted'
metrics <- rbind(m, weighted_m)

metrics$site <- unlist(lapply(metrics$network, function(x) strsplit(as.character(x), split = ', ')[[1]][1]))
metrics$year <- unlist(lapply(metrics$network, function(x) strsplit(as.character(x), split = ', ')[[1]][2]))
metrics$site <- gsub('DANUM', 'Danum', metrics$site)
metrics$site <- gsub('MALIAU', 'Maliau', metrics$site)
metrics$habitat_type <- as.factor(ifelse(metrics$site == 'SAFE', 'Logged', 'Primary'))
metrics$metric <- gsub(' ', '\n', metrics$metric)
metrics$metric <- gsub('\\.', '\n', metrics$metric)


sitescatter <- ggplot(metrics , aes(x = clustering, y = value, color = site)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  facet_grid(firstup(gsub('\\.', ' ', metric)) ~ source, scales = 'free')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(sitescatter)
#pdf('plots/Site comparisons/MOTU_sitescatter.pdf')
#sitescatter
#dev.off()

netscatter <- ggplot(metrics , aes(x = clustering, y = value, color = network)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  facet_grid(firstup(gsub('\\.', ' ', metric)) ~ source, scales = 'free')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(netscatter)
#pdf('plots/Site comparisons/MOTU_netscatter.pdf')
#netscatter
#dev.off()


habscatter <- ggplot(metrics , aes(x = clustering, y = value, color = habitat_type)) +
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
line_plot(input = metrics, network = 'network', clustering = 'clustering', metric = 'metric', value = 'value', plotname = 'Sabah')
#####ANOVA ####
## @knitr anova
metrics$network <- as.factor(metrics$network)
aov(aov(value ~ clustering + network + clustering:network, data = metrics[metrics$metric=='web asymmetry',]))


## @knitr funct_comp_scatter
funct_comp <- ggplot(metrics[metrics$metric=='functional.complementarity.HL',] , aes(x = clustering, y = value, color = network)) +
  geom_point()+
  labs(x = 'Clustering', y= 'Functional complementarity') +
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(funct_comp)

## @knitr mod_scat
mod_scat <- ggplot(metrics[metrics$metric=='modularity',] , aes(x = clustering, y = value, color = habitat_type)) +
  geom_point()+
  labs(x = 'Clustering', y= 'Modularity') +
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(mod_scat)

## @knitr func_scat
func_scat <- ggplot(metrics[metrics$metric=='functional.complementarity.HL',] , aes(x = clustering, y = value, color = habitat_type)) +
  geom_point()+
  labs(x = 'Clustering', y= 'Functional complementarity') +
  scale_x_continuous(breaks = seq(91, 98, 1))+
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(func_scat)


## @knitr GLM writing
##################################################

mod_mod <- (lm(value ~ site * clustering + year , metrics[metrics$metric=='modularity',]))
plot(mod_mod)

fun_mod <- (lm(value ~ site * clustering + year , metrics[metrics$metric=='modularity',]))
