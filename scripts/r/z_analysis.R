#### Header ####
## Project: bat diet
## Script purpose: calculating z-scores from the random files created by 
## array_site_motu_95_confidence_intervals
## Date: 2021-01-06
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

#### Setup ####

# Prevent partial-match errors 
options(warnPartialMatchArgs = TRUE)



library(tidyverse)

# Load and format the random data -----------------------------------------

rand_filenames <- list.files('data/output_data/for_z_scores/',
                          full.names = T)

rand_vals <- lapply(rand_filenames, read_csv) %>%
  bind_rows() %>%
  select(-X1) %>%
  mutate(meta = gsub('functional.complementarity', 
                     'functional complementarity', meta), 
         meta = gsub('weighted.NODF', 'weighted NODF', meta),
         meta = gsub('modularity', 'modularity', meta)) %>%
  separate(meta, into = c('metric', 'clustering', 'network'), sep = '\\.') %>%
  mutate(clustering = as.numeric(clustering)) %>%
  rename(random_vals = vals)


# Get the real values -----------------------------------------------

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
    }else{
      out_df <- lapply(netlists[[i]], function(x) networklevel(x, index = ind, level = 'higher'))
    }
    
    
    out_df <- bind_rows(out_df)
    
    out_df$metric <- ind
    out_df$clust <- c(91:98)[i]
    actual_vals_list[[z]] <- out_df
    z <- z + 1
  }
  
}

actual_vals <- bind_rows(actual_vals_list) %>%
  pivot_longer(cols = c('SAFE', 'DANUM', 'MALIAU'), names_to = 'network',
               values_to = 'actual') %>%
  rename(clustering = clust)


# Calculate 'z-scores' ----------------------------------------------------




summary_vals <- rand_vals %>%
  group_by(network, metric, clustering) %>%
  summarise(mean_rand = mean(random_vals),
            sd_rand = sd(random_vals))


z_vals <- summary_vals %>%
  ungroup() %>%
  # combine the datasets
  left_join(actual_vals) %>%
  # calculate the 'z' score
  mutate(z = (actual - mean_rand)/mean_rand) %>%
  # make some of the variables prettier for the outputs
  mutate(network = gsub('DANUM', 'Danum', network),
         network = gsub('MALIAU', 'Maliau', network),
         metric = gsub('functional complementarity', 'Functional complementarity',
                       metric),
         metric = gsub('modularity', 'Modularity', metric),
         metric = gsub('weighted NODF', 'WNODF', metric))

# write a csv of all of the values in their raw form
write_csv(z_vals, 'results/z_scores/all_z_scores.csv')


# then make a csv for each metric, with the values given to 3 decimal places
z_vals %>%
  select(-mean_rand, -sd_rand) %>%
  mutate(actual = round(actual, digits = 2),
         z = round(z, digits = 3),
         Treatment =ifelse(network == 'SAFE', 'Logged', 'Old growth')) %>%
  rename(Network = network, Metric = metric, `Clustering threshold` = clustering, 
         `Observed Value` = actual) %>%
  group_by(Metric) %>%
  do(write_csv(., paste0('results/z_scores/', unique(.$Metric), "_zscores.csv")))


# Finally, one table to rule them all...

met_list <- list()
cols <- c('actual', 'z')
for(met in unique(z_vals$metric)){
  met_list[[met]] <- z_vals %>%
    filter(metric == met) %>%
    rename_at(cols, list( ~paste( ., met, sep = '_') ) ) %>%
    select(-metric, - mean_rand, - sd_rand)
}

mets_df <- left_join(met_list[[1]], met_list[[2]], by = c('network', 'clustering')) %>%
  left_join(met_list[[3]], by = c('network', 'clustering'))

net_list <- list()
cols <- c('actual', 'z')
for(net in unique(z_vals$network)){
  net_list[[net]] <- z_vals %>%
    filter(network == net) %>%
    mutate(actual = round(actual, digits = 3),
    z = round(z, digits = 3)) %>%
    rename_at(cols, list( ~paste( ., net, sep = '_') ) ) %>%
    select(-network, - mean_rand, - sd_rand)
}

nets_df <- left_join(net_list[[1]], net_list[[2]], by = c('metric', 'clustering')) %>%
  left_join(net_list[[3]], by = c('metric', 'clustering'))

write_csv(nets_df, path = 'results/z_scores/grouped_by_metric.csv')

# Plot --------------------------------------------------------------------


ggplot(rand_vals, aes(x = clustering, y = random_vals)) +
  geom_violin() +
  facet_grid(metric ~ network, scales = 'free') +
  theme_classic()

for_plot <- z_vals %>%
  pivot_longer(cols = c('actual', 'z'), names_to = 'var_type', 
               values_to = 'var_value') %>%
  mutate(metric = gsub('Functional complementarity', 
                       'Functional\ncomplementarity', metric))

actual_points <-   ggplot(filter(for_plot, var_type == 'actual'), 
         aes(x = clustering, y = var_value, colour = network)) +
  geom_point() +
  scale_colour_viridis_d()+
  scale_x_continuous(breaks = c(91:98))+
  facet_wrap(. ~ metric, scales = 'free_y', ncol = 1,
             # sort the facet label placement
             strip.position = 'left') +
  theme_bw() +
    theme(legend.position = 'none',
          # sort the facet label placement
          strip.placement = "outside",
          strip.text.y = element_text(size = 8),
          #axis.title.y=element_blank(),
          axis.title.x=element_blank())+
  ylab('Observed value')

actual_points

z_points <-   ggplot(filter(for_plot, var_type == 'z'), 
                          aes(x = clustering, y = var_value, colour = network)) +
  geom_point() +
  scale_colour_viridis_d()+
  scale_x_continuous(breaks = c(91:98))+
  facet_wrap(. ~ metric, ncol = 1) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        # remove facet label
        strip.text = element_blank(),
        #axis.title.y=element_blank()
        ) +
  ylab('Z-value')

z_points

# Function for extracting the legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- ggplot(for_plot, 
                 aes(x = clustering, y = var_value, colour = network)) +
  geom_point() +
  scale_colour_viridis_d(name = 'Network')+
  theme_bw()+
  theme(legend.position = 'bottom')

legend <- g_legend(legend)


pdf('plots/z_plot_unedited.pdf')
grid_plot <- gridExtra::grid.arrange(actual_points, z_points,
                                     grid::textGrob('Clustering threshold'),
                                     legend,
                                     layout_matrix = rbind(c(1,2),
                                                           c(3,3),
                                                           c(4,4)),
                                     heights = unit(c(3.3, 0.3, 0.3), c('in', 'in', 'in')))

dev.off()
