#This script rarifies networks of cervinus to look at individual-level changes in network topology
if(interactive()==TRUE){
  library(here)
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(reshape2)
  library(lme4)
  library(plyr)
  library(plotrix)
  #library(multcompView)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data_2 <- field_data
#The bit below deals with the inevitable merging issues with some samples being from faeces_no1, some from faeces_no2
field_data$Faeces_no2 <- NULL

field_data_2$Faeces_no1 <- field_data_2$Faeces_no2
field_data_2$Faeces_no2 <- NULL

field_data <- rbind(field_data, field_data_2)
all_interactions <- r_network_gen(collapse_species = F, desired_species = 'Hice', include_malua = T, lulu=T)

#####Make a list with a network for each site####

sites_list <- list()
sites <- unique(colnames(all_interactions))

c_names <- c()


for(i in 1:length(sites)){
  m =  all_interactions[,which(colnames(all_interactions )==sites[i])]
  colnames(m) = m[1,]
  m = m[-c(1,2),]
  m <- m[-which(rowSums(m)==0),]
  print(ncol(m))
  c_names <- c(c_names, colnames(m))
  sites_list[[i]] <- m
}
print(sites)
names(sites_list) <- sites
names(sites_list) <- gsub('DANUM', 'Danum', names(sites_list))
names(sites_list) <- gsub('MALIAU', 'Maliau', names(sites_list))
names(sites_list) <- gsub('MALUA', 'SBE', names(sites_list))



#####Rarify #####

mets <- c('degree', 'normalised degree', 'resource range', 'proportional similarity')

outmat <- matrix(nrow= 0, ncol = 9)



for(i in 1:length(sites_list)){
  site <- names(sites_list)[i]
  for(p in 20: ncol(sites_list[[i]])){
    for(r in 1:1000){
      net <- sites_list[[i]][,sample(seq(1, ncol(sites_list[[i]])), replace = F, size = p)]
      sp <- specieslevel(net, level = 'higher', index = mets)
      n_MOTU <- length(which(rowSums(net)>0))
      out <- cbind(rep(site, nrow(sp)), rep(p, nrow(sp)), rep(n_MOTU, nrow(sp)), rep(r, nrow(sp)), sp, rownames(sp))
      print(out)
      outmat <- rbind(outmat, out)
    }
    
  }
  
}
colnames(outmat) <- c('Network', 'N_bats', 'N_MOTU', 'iteration', 'degree', 'normalised degree', 'resource range', 'proportional similarity', 'bat ID')

save.image('data/output_data/hice_stats/rarified_stats.RDS')
write.csv(outmat, 'shiny/hice_rarifying/hice_rarified.csv')

load('data/output_data/hice_stats/rarified_stats.RDS')


#####Now onto the statistical analysis section #####
rr_mod <- lmer(`resource range` ~ N_bats + (1|Network) + (1|`bat ID`), data=outmat)
 
summary(rr_mod)

nd_mod <- lmer(`normalised degree` ~ N_bats + (1|Network) + (1|`bat ID`), data=outmat)

summary(nd_mod)

ps_mod <- lmer(`proportional similarity` ~ N_bats + (1|Network) + (1|`bat ID`), data=outmat)

summary(ps_mod)
## As ps_mod fails to converge, I then try it without including bat ID as an effect. It then converges but reveals little
ps_mod_2 <- lmer(`proportional similarity` ~ N_bats + (1|Network), data=outmat)
summary(ps_mod_2)

#####Plot with ggplot#####



outdf <- outmat[,-c(4,9)]#I don't care about the bat ID or iteration for this bit, it just gets in the way of the melting

outdf$Network <- gsub('SAFE, ', 'SAFE,\n', outdf$Network)
outdf$Network <- gsub('Danum, ', 'Danum,\n', outdf$Network)
outdf$Network <- gsub('Maliau, ', 'Maliau,\n', outdf$Network)
outdf$Network <- gsub('SBE, ', 'SBE,\n', outdf$Network)

outdf <- melt(outdf, id.vars = c('Network', 'N_bats', 'N_MOTU'))

outdf$variable <- gsub('normalised degree', 'Normalised degree', outdf$variable)
outdf$variable <- gsub('^degree$', 'Degree', outdf$variable)
outdf$variable <- gsub('resource range', 'Resource range', outdf$variable)
outdf$variable <- gsub('proportional similarity', 'Proportional similarity', outdf$variable)

#Reorder the factors of the network column to make nicer plots
outdf$Network <- ordered(outdf$Network, levels = unique(outdf$Network)[order(unique(as.character(outdf$Network)))])



batfacet <- ggplot(outdf[-which(outdf$Network=='SAFE,\n2015'),], aes(N_bats, value)) + geom_bin2d(bins =70)+
  scale_fill_gradient(low = "blue",
                      high = "red", name = 'Count') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.title.y = element_blank())+ labs(x='Number of bats in network')+
  facet_grid(variable ~ Network, scale = 'free_y', switch = 'y') 
batfacet


pdf('plots/Hice/rarifying_bats_facetted.pdf', width = 10)
batfacet
dev.off()  

motufacet <- ggplot(outdf[-which(outdf$Network=='SAFE,\n2015'),], aes(N_MOTU, value)) + geom_bin2d(bins =70)+
  scale_fill_gradient(low = "blue",
                      high = "red", name = 'Count') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.title.y = element_blank())+ labs(x='Number of MOTU in network')+
  facet_grid(variable ~ Network, scale = 'free_y', switch = 'y') 


motufacet
pdf('plots/Hice/rarifying_motu_facetted.pdf', width = 10)
motufacet 
dev.off()

pdf('~/Desktop/temp_ps_example.pdf', width = 10)
ggplot(outdf[which(outdf$variable=='Proportional similarity'),], aes(N_bats, value)) + geom_bin2d(bins =70)+
  scale_fill_gradient(low = "blue",
                      high = "red", name = 'Count') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.title.y = element_blank())+ labs(x='Number of bats in network')+
  facet_grid(variable ~ Network, scale = 'free_y', switch = 'y') 
dev.off()
pdf('~/Desktop/temp_d_example.pdf', width = 10)
ggplot(outdf[which(outdf$variable=='Degree'),], aes(N_bats, value)) + geom_bin2d(bins =70)+
  scale_fill_gradient(low = "blue",
                      high = "red", name = 'Count') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.title.y = element_blank())+ labs(x='Number of bats in network')+
  facet_grid(variable ~ Network, scale = 'free_y', switch = 'y') 
dev.off()



#####Now using the values from outmat we calculate the mean value per bat when we allow 36 bats per network ####
out_36 <- outmat[outmat$N_bats==52,]


out_36_means <- out_36 %>%
  group_by(`bat ID`) %>%
  summarise_at(vars(-c(Network, N_MOTU, iteration, N_bats)), funs(mean(., na.rm=TRUE)))

out_36_means <- merge(x = out_36_means, y = field_data[which(!duplicated(field_data$Faeces_no1)),c('Faeces_no1', 'SiteAndYear')], by.x= "bat ID", by.y = "Faeces_no1", all.y = F)
melted_means <- melt(out_36_means)



out_36_se <- out_36 %>%
  group_by(`bat ID`) %>%
  summarise_at(vars(-c(Network, N_MOTU, iteration, N_bats)), funs(std.error(.)))

out_36_se <- merge(x = out_36_se, y = field_data[which(!duplicated(field_data$Faeces_no1)),c('Faeces_no1', 'SiteAndYear')], by.x= "bat ID", by.y = "Faeces_no1", all.y = F)

melted_se <- melt(out_36_se)
colnames(melted_means)[4] <- 'mean'
colnames(melted_se)[4] <- 'se'

all_melted <- cbind(melted_means, melted_se$se)

colnames(all_melted)[5] <- 'se'

all_melted$variable <- gsub('normalised.degree', 'Normalised degree', all_melted$variable)
all_melted$variable <- gsub('partner.diversity','Partner diversity', all_melted$variable)
all_melted$variable <- gsub('proportional.similarity','Proportional similarity', all_melted$variable)
all_melted$SiteAndYear <- gsub('DVCA', 'Danum', all_melted$SiteAndYear)
all_melted$SiteAndYear <- gsub('DANUM', 'Danum', all_melted$SiteAndYear)
all_melted$SiteAndYear <- gsub('MALIAU', 'Maliau', all_melted$SiteAndYear)
all_melted$SiteAndYear <- gsub('MALUA', 'SBE', all_melted$SiteAndYear)
all_melted$habitat_type <- rep(NA, nrow(all_melted))
all_melted[grep('SAFE', all_melted$SiteAndYear), 'habitat_type'] <- 'Logged'
all_melted[grep('SBE', all_melted$SiteAndYear), 'habitat_type'] <- 'Logged, replanted'
all_melted[grep('Danum', all_melted$SiteAndYear), 'habitat_type'] <- 'Primary'
all_melted[grep('Maliau', all_melted$SiteAndYear), 'habitat_type'] <- 'Primary'


facet_ridge <- ggplot(all_melted, aes (y=SiteAndYear, x =mean, fill = habitat_type)) + 
  geom_density_ridges(scale= 0.5)+ #The scale determines the space between the rows
  theme_ridges()+ #This changes the theme to make it more aesthetically pleasing
  scale_fill_cyclical(values = c("#85d7da","#cfb4de","#d0ca9f"), guide = 'legend', name = 'Habitat type')+
  scale_x_continuous(expand = c(0.01, 0)) + #Make the space between the labels and plot smaller
  scale_y_discrete(expand = c(0.01, 0))+ #Make it so the top series actually fits in the plot
  ylab(NULL)+ xlab(NULL)+
  facet_wrap( ~ variable, ncol=1, scales = 'free_x', strip.position = 'bottom')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"),#strip stuff sorts the facet labels, spacing adjusts the space between facets
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        text = element_text(size=12))+
  theme(legend.text =element_text(size = 10)) #Trying to standardise the sizes across both plots
facet_ridge

anova_list <- list()
for(i in 1:length(unique(all_melted$variable))){ #Could probably just use lapply now that a load of the loop has been killed
  v <- as.character(unique(all_melted$variable)[i])
  s <- all_melted[which(all_melted$variable==v),]
  my_anova <- aov(s$mean ~ s$SiteAndYear)
  summary(my_anova)
  anova_list[[i]] <- TukeyHSD(my_anova)
  names(anova_list)[i] <- v
}

mets <- c()
net1s <- c()
net2s <- c()
signif <- c()
for(l in 1:length(anova_list)){
  lis <- anova_list[[l]]$`s$SiteAndYear`
  for(i in 1:nrow(lis)){
    net1 <- strsplit(rownames(lis),split = '-')[[i]][1]
    net2 <- strsplit(rownames(lis),split = '-')[[i]][2]
    
    if(lis[i,4]<=0.05){
      sig <- 'significant difference'
    }else{sig <- 'non-significant'}
    
    net1s <- c(net1s, net1)
    net2s <- c(net2s, net2)
    mets <- c(mets, names(anova_list)[l])
    signif <- c(signif, sig)
  }
}

sig_df <- data.frame(as.factor(mets), as.factor(net1s), as.factor(net2s), as.factor(signif))
colnames(sig_df) <- c('metrics', 'network_1', 'network_2', 'significance')
# levels(sig_df$network_1) <- unique(all_melted$SiteAndYear)
# levels(sig_df$network_2) <- unique(all_melted$SiteAndYear)
# sig_df$network_1 <- ordered(all_melted$SiteAndYear, levels=unique(all_melted$variable)[order(as.character(unique(all_melted$variable)))])
sig_df$network_2 <- fct_rev(sig_df$network_2) #ggridges plots the factors in an annoying order, this rectifies it

ultimate_plot <- ggplot(sig_df, aes(x = network_1, y = network_2)) +
  geom_tile(aes(fill=significance))+
  scale_fill_manual(values =c("white", "black"), na.value = 'white')+
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill= '#ededed', '#ededed'), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=12))+
  ylab(NULL)+ xlab(NULL)+
  facet_wrap( ~ metrics, ncol=1, scales = 'free_x', strip.position = 'bottom')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(1.7, "lines"),#strip stuff sorts the facet labels, spacing adjusts the space between facets
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
ultimate_plot


grid.arrange(facet_ridge, ultimate_plot, nrow = 1)


trial_subset <- all_melted[all_melted$SiteAndYear=='Danum, 2017',]
trial_subset <- trial_subset[trial_subset$variable=='Proportional similarity',]


ggplot(trial_subset, aes(x=mean, y = mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se))+
  geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"))



trial_subset <- all_melted[all_melted$variable=='Proportional similarity',]


ggplot(trial_subset, aes(x=mean, y = se)) + 
  geom_point(size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill="white"))+
  facet_wrap(~ SiteAndYear)

