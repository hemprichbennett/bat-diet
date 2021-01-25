##################################################
## Project: bat diet
## Script purpose: Analysing how 95% confidence intervals of individual values change over clustering thresholds used
## Date: 04/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: Uses the output of the script 'array_95_confidence_intervals_inter_network.R'
##################################################

library(tidyverse)
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

files <- list.files(path = 'data/output_data/for_z_scores/', 
                    pattern ='.csv', full.names = T)


file_list <- lapply(files, read.csv, row.names= 1)

ranges_df <- bind_rows(file_list) %>%
  # do some basic formatting as the csvs weren't nicely formatted
  mutate(meta =gsub('functional.complementarity', 'functional complementarity', meta),
         meta =gsub('weighted.NODF', 'weighted NODF', meta),
         meta = gsub('modularity', 'modularity', meta)) %>%
  separate(meta, into =c('metric', 'clustering level', 'dataset'), sep = '\\.') %>%
  mutate(`clustering level` = as.numeric(`clustering level`)) %>%
  # now do the necessary calculations
  group_by(metric, `clustering level`, dataset) %>%
  summarise(lower = quantile(vals, 0.025),
            upper = quantile(vals, 0.975)) 



# Now calculate observed values -------------------------------------------

source('scripts/r/r_network_gen.r')

inpath <- 'data/processed_dna_data/lulu/'
filenames <- list.files(pattern = '.csv', path = inpath)
#filenames <- filenames
filenames
filenames <- paste(inpath, filenames, sep = '')
#filenames <- filenames[grep('lulu', filenames)]

rawnets <- lapply(filenames, read.csv, header = F, stringsAsFactors = F, row.names=1)
names(rawnets) <- gsub('.*\\/', '', filenames)
names(rawnets) <- gsub('_.+', '', names(rawnets))

for(i in 1:length(rawnets)){
  rawnets[[i]][2:nrow(rawnets[[i]]),] <- ifelse(rawnets[[i]][2:nrow(rawnets[[i]]), ] == 0, 0, 1)
}

netlists <- lapply(rawnets, function(x) r_network_gen(input= x,  collapse_species = T, filter_species = T))

names(netlists) <- names(rawnets)



actual_vals_list <- list()
z <- 1
for(ind in c('modularity', 'functional complementarity', 'weighted NODF')){
  for(i in 1:length(netlists)){
    
    if(ind == 'modularity'){
      out_df <- lapply(netlists[[i]], function(x) slot(computeModules(web = x), 
                                                       'likelihood'))
      out_df <- bind_rows(out_df)
    }else{
      out_df <- lapply(netlists[[i]], function(x) networklevel(x, index = ind, level = 'higher'))
      out_df <- do.call(rbind, out_df) %>% 
        as.data.frame() %>% 
        rownames_to_column()
      
      if(ind == 'functional complementarity'){
         out_df <- out_df %>% 
          pivot_wider(names_from = `rowname`, 
                      values_from = functional.complementarity.HL)
      }
      if(ind == 'weighted NODF')
        out_df <- out_df %>% 
        pivot_wider(names_from = `rowname`, 
                    values_from = `weighted NODF`)
    }
    
    
    
    
    out_df$metric <- ind
    out_df$clust <- c(91:98)[i]
    actual_vals_list[[z]] <- out_df
    z <- z + 1
  }
  
}

actual_vals <- bind_rows(actual_vals_list) %>%
  pivot_longer(cols = c('SAFE', 'DANUM', 'MALIAU'), names_to = 'dataset',
               values_to = 'actual') %>%
  rename(`clustering level` = clust)


all_df <- ranges_df %>%
  left_join(actual_vals) %>%
  rename(network = dataset) %>%
  mutate(network = gsub('MALIAU', 'Maliau', network),
         network = gsub('DANUM', 'Danum', network)) %>%
  mutate(metric = firstup(metric))


# 
# colnames(ranges_df)[c(4,5)] <- c('lower', 'upper')
# ranges_df$network <- gsub('DANUM', 'Danum', ranges_df$network)
# ranges_df$network <- gsub('MALIAU', 'Maliau', ranges_df$network)
# 
# for_BCI <- ranges_df[which(ranges_df$clustering==95),]
# for_BCI[,c('clustering', 'lower', 'upper', 'dummy_for_plotting')] <- NULL
# for_BCI$metric <- firstup(as.character(for_BCI$metric))
# 
# write.csv(for_BCI, 'BCI_report/bat_metrics.csv')

for_facets <- all_df
#Calculate if a value is significant or not, remove it if inside the 
# randomised range as we don't want to plot it
for_facets$dummy_for_plotting <- for_facets$actual
for_facets$signif <- 'significant'
for(i in 1:nrow(for_facets)){
  if(for_facets$actual[i] > for_facets$lower[i] & for_facets$actual[i] < for_facets$upper[i]){
    for_facets$signif[i] <- 'non significant'
    for_facets$actual[i] <- NA
  }
}


#for_facets <- rbind(for_facets, for_sig_only_vals)
# for_facets$metric <- as.character(for_facets$metric)
# for_facets$metric <- sapply(for_facets$metric, firstup)
for_facets$metric <- gsub(' ', '\n', for_facets$metric)


# for_facets <- for_facets %>%
#   #Remove metrics which are uninteresting #
#   filter(metric %in% c('Functional\ncomplementarity',
#                       "Modularity", 'Weighted\nNODF')) %>%
#   mutate(metric = gsub('\\n', ' ',metric)) %>%
#   mutate(metric = gsub('Weighted NODF', 'WNODF', metric))


#plot
metrics_facet <- ggplot(for_facets, aes(y =`clustering level`, x= dummy_for_plotting, colour = network))+
  geom_errorbarh(aes(xmin=lower, xmax=upper, colour = network), height = 0.4, alpha = 0.8, show.legend = F)+
  geom_point(aes(y = `clustering level`, x = actual, size = 1.4))+
  facet_wrap(~metric, scales = 'free', ncol = 3, strip.position="bottom")+
  scale_colour_manual(values=cbPalette, name = 'Observed network\nvalue')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.placement = "outside", 
        strip.background =element_rect(fill="white", colour = 'white'),
        text = element_text(size=20))+
  theme(legend.position="bottom")+
  labs(y = 'Clustering %', x = NULL) + 
  # tell ggplot to plot the colour in the legend but not the size
  guides(colour = "legend", size = "none") +
  # make the points in the legend bigger
  guides(colour = guide_legend(override.aes = list(size=4.5)))
metrics_facet

pdf('plots/randomized_ranges/all_facet.pdf', width = 18)
metrics_facet
dev.off()

jpeg("plots/randomized_ranges/all_facet.jpg", units = "in", width = 18, height = 7, res = 300)  
metrics_facet
dev.off()



