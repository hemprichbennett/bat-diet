##################################################
## Project: bat diet
## Script purpose: Analysing how 95% confidence intervals of individual values change over clustering thresholds used
## Date: 04/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: Uses the output of the script 'array_95_confidence_intervals_inter_network.R'
##################################################

library(ggplot2)
library(gridExtra)

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} #A function to capitalise the metric names when making plots

# cbPalette <- c("#8d5ce4",
#                 "#00a5bc",
#                 "#ef115d")

cbPalette <- c("#5a3fa8",
                "#d4c200",
                "#a4c7ff")

files <- list.files(path = 'data/output_data/randomized_ranges/', pattern ='.csv')
files <- paste('data/output_data/randomized_ranges/', files, sep ='')

files <- files[-grep('smaller', files)]
file_list <- lapply(files, read.csv, row.names= 1)

ranges_df <- do.call(rbind, file_list)
colnames(ranges_df)[c(4,5)] <- c('lower', 'upper')
ranges_df$network <- gsub('DANUM', 'Danum', ranges_df$network)
ranges_df$network <- gsub('MALIAU', 'Maliau', ranges_df$network)

for_facets <- ranges_df
#Calculate if a value is significant or not
for_facets$dummy_for_plotting <- for_facets$actual
for_facets$signif <- 'significant'
for(i in 1:nrow(for_facets)){
  if(for_facets$actual[i] > for_facets$lower[i] & for_facets$actual[i] < for_facets$upper[i]){
    for_facets$signif[i] <- 'non significant'
    for_facets$actual[i] <- NA
  }
}


#for_facets <- rbind(for_facets, for_sig_only_vals)
for_facets$metric <- as.character(for_facets$metric)
for_facets$metric <- sapply(for_facets$metric, firstup)
for_facets$metric <- gsub(' ', '\n', for_facets$metric)




#Remove metrics which are uninteresting #
unique(for_facets$metric)
for_facets <- for_facets[-which(for_facets$metric %in% c('Alatalo\ninteraction\nevenness', 'Web\nasymmetry', "Fisher\nalpha", "Interaction\nstrength\nasymmetry", "Togetherness", 'Nestedness', 'Weighted\nnestedness', 'Weighted\nNODF')),]


#plot
metrics_facet <- ggplot(for_facets, aes(y =clustering, x= dummy_for_plotting, colour = network))+
  geom_errorbarh(aes(xmin=lower, xmax=upper, colour = network), height = 0.4, alpha = 0.8, show.legend = F)+
  geom_point(aes(y = clustering, x = actual))+
  facet_wrap(~metric, scales = 'free', ncol = 3, strip.position="bottom")+
  scale_colour_manual(values=cbPalette, name = 'Observed network\nvalue')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.placement = "outside", 
        strip.background =element_rect(fill="white", colour = 'white'))+
  theme(legend.position="bottom")+
  labs(y = 'Clustering %', x = NULL)
metrics_facet

pdf('plots/randomized_ranges/all_facet.pdf', width = 15)
metrics_facet
dev.off()

  
