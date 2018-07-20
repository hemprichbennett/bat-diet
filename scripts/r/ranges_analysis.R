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
file_list <- lapply(files, read.csv, row.names= 1)

ranges_df <- do.call(rbind, file_list)
colnames(ranges_df)[c(4,5)] <- c('lower', 'upper')
ranges_df$network <- gsub('DANUM', 'Danum', ranges_df$network)
ranges_df$network <- gsub('MALIAU', 'Maliau', ranges_df$network)

#Calculate if a value is significant or not
ranges_df$signif <- rep('significant', nrow(ranges_df))
for(i in 1:nrow(ranges_df)){
  if(ranges_df$actual[i] > ranges_df$lower[i] & ranges_df$actual[i] < ranges_df$upper[i]){
    ranges_df$signif[i] <- 'non significant'
  }
}

##First plot and output all metrics with all values
p_list <- list()
for(i in 1: length(unique(ranges_df$metric))){
  met <- unique(ranges_df$metric)[i]
  pdf(paste('plots/randomized_ranges/', met, '.pdf', sep =''))
  p_list[[met]] <- ggplot(ranges_df[which(ranges_df$metric==met),], aes(clustering, actual, colour = network))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper, colour = network), width=.1, alpha = 0.5)+
    scale_fill_manual(values=cbPalette)+
    scale_colour_manual(values=cbPalette)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(p_list[[met]])
  dev.off()
}


#Now make a set of paired plots, showing significance or lack of it

ranges_df$metric <- gsub(' ', '\n', ranges_df$metric)
# 
# for(i in 1: length(unique(ranges_df$metric))){
# 
#   met <- unique(ranges_df$metric)[i]
#   
#  nolegend <- ggplot(ranges_df[which(ranges_df$metric==met),], aes(clustering, actual, colour = network))+
#     geom_point()+
#     geom_errorbar(aes(ymin=lower, ymax=upper, colour = network), width=.1, alpha = 0.5)+
#     scale_fill_manual(values=cbPalette)+
#     scale_colour_manual(values=cbPalette)+
#     theme_bw()+
#     theme(legend.position="none")+
#     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#    labs(x = 'Clustering', y = firstup(as.character(met)))
#   
# 
#    only_sig <- ggplot(ranges_df[which(ranges_df$metric==met & ranges_df$signif=='significant'),], aes(clustering, actual, colour = network))+
#        geom_point()+
#        scale_fill_manual(values=cbPalette)+
#        scale_colour_manual(values=cbPalette)+
#        theme_bw()+
#        theme(legend.position="none")+
#        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#    labs(x = 'Clustering', y = firstup(as.character(met)))
#  
#   leg <- g_legend(ggplot(ranges_df[which(ranges_df$metric==met & ranges_df$signif=='significant'),], aes(clustering, actual, colour = network))+
#                     geom_point()+
#                     geom_errorbar(aes(ymin=lower, ymax=upper, colour = network), width=.1, alpha = 0.5)+
#                     scale_fill_manual(values=cbPalette)+
#                     scale_colour_manual(values=cbPalette)+
#                     theme_bw()+
#                     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
#                   
#   pdf(paste('plots/randomized_ranges/sig_and_non_', met, '.pdf', sep =''))
#   grid.arrange(nolegend, only_sig, leg, ncol = 3)
#     dev.off()
# }
# 
# dev.off()




all <- ggplot(ranges_df, aes(clustering, actual, colour = network))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper, colour = network), width=.1, alpha = 0.5)+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x = 'Clustering', y = 'Measured value, random values')+
  facet_wrap(~ metric, nrow = 1, scales = 'free_y')+
  theme(legend.position="none")
  

all_only_sig <- ggplot(ranges_df[which(ranges_df$signif=='significant'),], aes(clustering, actual, colour = network))+
  geom_point()+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x = 'Clustering threshold', y = 'Measured value')+
  facet_wrap(~ metric, nrow = 1, scales = 'free_y')+
  theme(legend.position="none")

grid.arrange(all, all_only_sig, leg, nrow = 3)


#Make a df for a facetted version, as thats clearly a better way of doing this

for_sig_only_vals <- ranges_df[which(ranges_df$signif=='significant'),]
for_sig_only_vals$lower <- NA
for_sig_only_vals$upper <- NA
for_sig_only_vals$signif <- 'Only significant values'

for_facets <- ranges_df
for_facets$signif <- 'All values'


for_facets <- rbind(for_facets, for_sig_only_vals)
for_facets$metric <- firstup(for_facets$metric)

#Remove metrics which are uninteresting #

for_facets <- for_facets[-which(for_facets$metric %in% c('Alatalo\ninteraction\nevenness','Togetherness', 'Niche\noverlap', 'Web\nasymmetry')),]


#plot
metrics_facet <- ggplot(for_facets, aes(y =clustering, x=actual, colour = network))+
  geom_point()+
  facet_grid(signif ~ metric, scales = 'free_x')+
  geom_errorbarh(aes(xmin=lower, xmax=upper, colour = network), height = 0.4, alpha = 0.5)+
  scale_colour_manual(values=cbPalette, name = 'Network')+
  scale_linetype_discrete(name = 'Error')+
  theme_dark()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y = 'Clustering %', x = 'Measured value, random values')
metrics_facet

pdf('plots/randomized_ranges/all_facet.pdf', width = 14)
metrics_facet
dev.off()

  