
if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(reshape2)
  
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T)

#####Generating all of the actual values for the datasets ####

desired_mets <- c('functional complementarity',
                  'web asymmetry',
                  'Alatalo interaction evenness',
                  'togetherness',
                  'Fisher alpha', 'mean number of shared partners',
                  'niche overlap',
                  'nestedness',
                  'discrepancy',
                  'ISA')



sets <- c()
mets <- c()
vals <- c()
for(n in 1: length(nets)){
  for(m in 1:length(desired_mets)){
    sets <- c(sets, names(nets)[n])
    st <- Sys.time()
    val <- networklevel(nets[[n]], index = desired_mets[m])[1]
    vals <- c(vals, val) #Have to select the first value as some metrics return multiple metrics and break everything
    mets <- c(mets, names(val))
    en <- Sys.time()
    print(val)
    cat(names(nets)[n], desired_mets[m], 'finished, it took ', en-st, '\n')
  }
}

realstats <- data.frame(sets, mets, vals)

speciesvals <-  lapply(nets, function(x) specieslevel(x, level = 'higher', 
                                                                   index = c('normalised degree', 'proportional similarity', 'resource range')))

for( i in seq_along(speciesvals)){
  
  speciesvals[[i]]$species <- rownames(speciesvals[[i]])
  
}

speciesvals <- melt(speciesvals)

#####making a very time-intensive list of all of the null-model differences ####
temp_desired_mets <- desired_mets[1]



#outlist = mapply(hernani_comparisons, temp_nets, temp_nets, temp_desired_mets, SIMPLIFY = FALSE)

#outlist = mapply(hernani_comparisons, nets, nets, desired_mets, SIMPLIFY  = FALSE) #If simplify is T it makes the output a matrix/vector etc. I want a list

####Generate the combinations of networks we want, so as to make a nice triangle plot later####

sorted_names <- sort(names(nets), decreasing = T)

combinations <- matrix(ncol=2, nrow=0)
for(i in 1:(length(sorted_names)-1)){
  cat('i is ', sorted_names[i],'\n')
  #print(i)
  for(z in seq(i+1,length(sorted_names))){
    print(sorted_names[z])
    combinations <- rbind(combinations,c(sorted_names[i], sorted_names[z]))
  }
}

####Do the null-model stuff####

probs <- matrix(nrow = 0, ncol = 6)
for(i in 1:nrow(combinations)){
  for(m in 1:length(desired_mets)){
    net1 = combinations[i,1]
    net2 = combinations[i,2]
    metric = desired_mets[m]
    p <- hernani_comparisons(net1 = nets[[net1]], net2=nets[[net2]], metric = metric, sums_to_preserve = 'both')
    outbit <- cbind(rownames(p), net1, net2, p)
    probs <- rbind(probs, outbit)
    
  }
}

save.image('data/output_data/full_species_pairwise_diffs.RDS')
load('data/output_data/full_species_pairwise_diffs.RDS') #####DON'T TRUST ANY OF THIS, THERES A BUG IN THE CODE THAT GENERATES IT

##Need to convert the null-model data into a format that we can actually use for plotting etc
significance <- as.data.frame(probs)
colnames(significance)[1] <- 'metric'
significance$`p (rand <= orig)` <- as.numeric(levels(significance$`p (rand <= orig)`))[significance$`p (rand <= orig)`]
significance$`p (rand >= orig)` <- as.numeric(levels(significance$`p (rand >= orig)`))[significance$`p (rand >= orig)`]

significance$is_signif <- rep('Non-significant', nrow(significance))

for(i in 1:nrow(significance)){
  # if(significance[i,4]>=0.95 & significance[i,5]>=0.95){
  #   #significance[i,'is_signif'] <- 'Significantly similar'
  #   cat('row', i, 'True')
  #   next()  #This is unneeded now as Steve LC pointed out that theres no such thing as lower than random, only ns
  # }
  if(significance[i,4]>=0.95){
    significance[i,'is_signif'] <- 'Greater than random'
  }
  if(significance[i,5]>=0.95){
    next()
  }
}





significance$net1 <- factor(significance$net1, levels = levels(significance$net1)[order(as.character(levels(significance$net1)))])
significance$net2 <- factor(significance$net2, levels = levels(significance$net2)[order(as.character(levels(significance$net2)))]) 


#####Now add a column on what it means if theres a sig difference or not ####

significance$verdict <- rep('N.d/N.s.d', nrow(significance))
#significance$for_minuses <- rep('', nrow(significance))

significance <- significance[-which(grepl('LL', significance$metric)),]

for(i in 1:nrow(significance)){
  v <- significance$is_signif[i]
  if(v=="Non-significant"){next()}
  if(v=="Lower than random"){
    significance[i,'verdict'] <- 'Significantly low difference'
    #significance[i, 'for_minuses'] <- '-'
  }
  if(v=="Greater than random"){
    net1 = significance[i,'net1']
    net2 = significance[i,'net2']
    metric = significance[i,'metric']
    #metric <- gsub('\\.', ' ', metric) #The input and output format of bipartite are different and it fills me with rage
    #metric <- gsub(' HL', '', metric)
    #metric <- gsub(' LL', '', metric)
    
    net1val <- realstats[which(as.character(realstats$sets)==net1 & as.character(realstats$mets)==metric),'vals']#Had to use as.character as the dfs are factor-based
    net2val <- realstats[which(as.character(realstats$sets)==net2 & as.character(realstats$mets)==metric),'vals']
    
    if(net1val> net2val){
      biggest <- 'y-axis significantly\n greater bigger than x-axis'
    }else{
      biggest <- 'x-axis significantly\n greater than y-axis'
    }
    significance[i,'verdict'] <- biggest
  }
}

significance$verdict <- as.factor(significance$verdict)
#significance$for_minuses <- as.factor(significance$for_minuses)

significance$metric<- gsub("Alatalo interaction evenness", 'Alatalo\ninteraction\nevenness', significance$metric)  
significance$metric<- gsub("functional.complementarity", paste('Functional\n', 'complementarity', sep = ''), significance$metric)  
significance$metric<- gsub("web.asymmetry", paste('Web\n', 'asymmetry', sep = ''), significance$metric)  
significance$metric<- gsub('mean.number.of.shared.partners', paste('Mean number\n', 'of shared\n', 'partners', sep = ''), significance$metric)  
significance$metric<- gsub('niche.overlap.HL', 'Niche\noverlap HL', significance$metric)  
significance$metric<- gsub('interaction strength asymmetry', 'Interaction\nstrength\nasymmetry', significance$metric)  
significance$metric<- gsub('discrepancy.HL', 'Discrepancy\nHL', significance$metric)  
significance$metric<- gsub('togetherness.HL', 'Togetherness\nHL', significance$metric)  
####We need to make the realstats item plottable too ####
# 
# real_comparisons <- matrix(nrow = 0, ncol = 4)
# for(i in 1:nrow(combinations)){
#   for(m in 1:length(desired_mets)){
#     net1 = combinations[i,1]
#     net2 = combinations[i,2]
#     metric = desired_mets[m]
#     
#     net1val <- realstats[which(realstats$sets==net1 & realstats$mets==metric),'vals']
#     net2val <- realstats[which(realstats$sets==net2 & realstats$mets==metric),'vals']
#     
#     if(net1val> net2val){
#       biggest <- 'y-axis'
#     }else{
#       biggest <- 'x-axis'
#     }
#     
#     if(metric=="Alatalo interaction evenness"){metric <- paste('Alatalo.interaction\n', 'evenness', sep = '')}
#     if(metric=="functional complementarity"){metric <- paste('Functional\n', 'complementarity', sep = '')}  
#     if(metric=="web asymmetry"){metric <- paste('Web\n', 'asymmetry', sep = '')}  
#     if(metric=="mean number of shared partners"){metric <- paste('Mean number\n', 'of shared\n', 'partners', sep = '')}  
#     out <- cbind(metric, net1, net2, biggest)
#     real_comparisons <- rbind(real_comparisons, out)
#     
#   }
# }
# real_comparisons <- as.data.frame(real_comparisons)
# 
# 



####Plotty plotty plotty####



ultimate_plot <- ggplot(significance, aes(x = net2, y = net1)) +
  geom_tile(aes(fill=verdict))+
  facet_grid(metric ~ .)+
  scale_fill_manual(values =c("white", "lightgray","#9f9244","#7f64b9"))+
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.text.x = element_text(size = 6))#Stuff for the facet labels

ultimate_plot

pdf('plots/Site comparisons/diffs_and_significance.pdf')
ultimate_plot
dev.off()

example_tiles <- ggplot(significance[which(significance$metric=="niche.overlap.HL"),], aes(x = net2, y = net1)) +
  geom_tile(aes(fill=verdict))+
  scale_fill_manual(values =c("white", "lightgray","#9f9244","#7f64b9"))+
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
pdf('~/Desktop/temp_tiles.pdf')
example_tiles
dev.off()
# 
# 
# 
# significance_plot <- ggplot(significance, aes(x = net2, y = net1)) +
#   geom_tile(aes(fill=is_signif))+
#   facet_grid(metric ~ .)+
#   guides(fill=guide_legend(title="Difference")) + 
#   scale_fill_manual(values =c("#c36785","#9f9244","#7f64b9"))+
#   theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
# significance_plot
# 
# realdiffs_plot <- ggplot(real_comparisons, aes(x = net2, y = net1)) +
#   geom_tile(aes(fill=biggest))+
#   facet_grid(metric ~ .)+
#   guides(fill=guide_legend(title="Largest value")) + 
#   scale_fill_manual(values =c("#5d0900","#7190e8"))+
#   theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
# realdiffs_plot
# pdf('plots/Site comparisons/crude_tileplot.pdf')
# grid.arrange(realdiffs_plot, significance_plot, ncol = 2)
# dev.off()


# for(i in 1: length(desired_mets)){
#   
#   print(networklevel(nets[[1]], index = desired_mets[i]))
# }

#####Make a long dataset that we can use to visualise all of this #####

monster_stats <- realstats
monster_stats$site <- rep(NA, nrow(monster_stats))
monster_stats$year <- rep(NA, nrow(monster_stats))

for(i in 1:nrow(monster_stats)){
  s <- strsplit(as.character(monster_stats$sets[i]),', ', fixed = T)
  monster_stats$site[i] <- s[[1]][1]
  monster_stats$year[i] <- s[[1]][2]
}

#####Put the original values in a scatterplot ####

for_example <- ggplot(monster_stats[which(monster_stats$mets=='niche.overlap.HL'),], aes(x =year, y = vals, group = site)) + 
  geom_point(aes(color=site))+ 
  geom_line(aes(color=site))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

pdf('~/Desktop/temp_points.pdf')
for_example
dev.off()

point_plot <- ggplot(monster_stats, aes(x =year, y = vals, group = site)) + 
  geom_point(aes(color=site))+ 
  geom_line(aes(color=site))+
  facet_wrap(~ mets, scales = 'free_y', ncol = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

point_plot

pdf('plots/Site comparisons/point_plot.pdf')
point_plot
dev.off()


####Look at some of the specieslevel stuff #####
colnames(speciesvals)[4] <- 'sets'

speciesvals$site <- rep(NA, nrow(speciesvals))
speciesvals$year <- rep(NA, nrow(speciesvals))

for(i in 1:nrow(speciesvals)){
  s <- strsplit(as.character(speciesvals$sets[i]),', ', fixed = T)
  speciesvals$site[i] <- s[[1]][1]
  speciesvals$year[i] <- s[[1]][2]
}

species_plot <- ggplot(speciesvals, aes(x =year, y = value, group = site)) + 
  geom_point(aes(color=site))+ 
  geom_line(aes(color=site))+
  facet_grid(variable ~ species, scales = 'free_y')#+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))

species_plot

pdf('plots/Site comparisons/species_plot.pdf')
species_plot
dev.off()

f()