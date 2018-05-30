## @knitr inext_setup
##################################################
## Project: bat-diet (all bats)
## Script purpose: using iNEXT to assess the sampling completeness and total diversity of each of my
## sampling events
## Date: 28/05/18
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
all_interactions <- r_network_gen(collapse_species = F, include_malua = F, lulu = T)



bad <- c("MOTU", "3483", "3508", "3510", "3527", "3570", "3575", "3585", "3586", "3623", "3626", "3647", 
 "3673", "3675", "3682", "3684", "3695", "3715", "3747", "3760", "3775", "3781", "3805", "3806", 
 "3823", "3833", "3845", "3847", "3848", "3851", "3857", "3864", "3893", "3898", "3899", "3911", 
 "3925", "3930", "3934", "3939", "3947", "4004", "4115", "4191", "4194", "4289", "4303", "4321", 
 "1028", "3594", "3749", "3759", "3761", "3945", "3960", "2559", "BLANK", "TT105", "TT240", "TT52", 
 "TT63", "2667", "1332", "3685", "3727", "3784")

all_interactions <- all_interactions[,-which(all_interactions[1,] %in% bad)]

#####Make a list with a network for each site####

sites_list <- list()
sites <- unique(unlist(all_interactions[1,]))




for(i in 1:length(sites)){
  m <- all_interactions[2:nrow(all_interactions),which(all_interactions[1,]==sites[i])]
  colnames(m) = m[1,]
  m <- apply(m, 2, as.numeric)
  m = m[-c(1,2),]
  m <- m[-which(rowSums(m)==0),]
  sites_list[[i]] <- m
}
print(sites)
names(sites_list) <- sites


inext_list <- list()

for(i in 1:length(sites_list)){
  nam <- names(sites_list)[i]
  inext_list[[nam]] <- c(ncol(sites_list[[i]]), rowSums(sites_list[[i]]))
  
}

names(inext_list) <- gsub('DANUM', 'Danum', names(inext_list))
names(inext_list) <- gsub('DVCA', 'Danum', names(inext_list))
names(inext_list) <- gsub('MALUA', 'SBE', names(inext_list))
names(inext_list) <- gsub('MALIAU', 'Maliau', names(inext_list))

#For plotting
a <- iNEXT(inext_list, datatype = 'incidence_freq')

inext_list[['all']] <- c(ncol(all_interactions), rowSums(apply(all_interactions[-c(1,2),], 2, as.numeric)))
inc_all <- iNEXT(inext_list, datatype = 'incidence_freq')

#####Here I replicate the gginext command, but tweak it a bit because it doesn't natively allow it ####
## @knitr inext_plot
z <- fortify.iNEXT(a)
z$site_type <- rep(NA, nrow(z))
z[grep('SAFE', z$site),'site_type'] <- 'Logged'
z[grep('Danum', z$site),'site_type'] <- 'Primary'
z[grep('Maliau', z$site),'site_type'] <- 'Primary'
z$site_type <- as.factor(z$site_type)
type = 1
se = TRUE
facet.var = "site"
color.var = "site"
grey = T

if(ncol(z) ==7) {se <- FALSE}
datatype <- unique(z$datatype)
if(color.var=="none"){
  if(levels(factor(z$order))>1 & "site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
    color.var <- "both"
    z$col <- z$shape <- paste(z$site, z$order, sep="-")
    
  }else if("site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
    color.var <- "site"
    z$col <- z$shape <- z$site
  }else if(levels(factor(z$order))>1){
    warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
    color.var <- "order"
    z$col <- z$shape <- factor(z$order)
  }else{
    z$col <- z$shape <- rep(1, nrow(z))
  }
}else if(color.var=="order"){     
  z$col <- z$shape <- factor(z$order)
}else if(color.var=="site"){
  if(!"site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
    z$col <- z$shape <- factor(z$order)
  }
  z$col <- z$shape <- z$site
}else if(color.var=="both"){
  if(!"site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
    z$col <- z$shape <- factor(z$order)
  }
  z$col <- z$shape <- paste(z$site, z$order, sep="-")
}

z$lty <- factor(z$method, c("interpolated", "observed", "extrapolated"), c("interpolation", "observation", "extrapolation"))
z$col <- factor(z$col)
data.sub <- z[which(z$method=="observed"),]

palette <- c("#E69F00", "#56B4E9", "#009E73")

z[grep('observation', z$lty),'lty'] <- NA

g <- ggplot(na.omit(z), aes_string(x="x", y="y", color = 'site_type')) + 
  geom_point(size=3, data=data.sub)+
  ylab('OTU diversity') + xlab('Number of bats sampled')+
  geom_line(aes_string(linetype="lty"), lwd=0.5)+
  scale_color_manual(values=palette)

#g

g <- g +
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #Get rid of a load of the default crap
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=8))


g <- g + facet_wrap( ~ site, nrow=2, strip.position = 'top')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"))#strip stuff sorts the facet labels, spacing adjusts the space between facets
g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr"), alpha=0.2)
## @knitr plotting
g

pdf('plots/all_inext.pdf', width = 10, height = 7)
g
dev.off()
## @knitr inext_nos

####Rearrange our stats ####
asymptote_ests <- inc_all$AsyEst
asymptote_ests <- asymptote_ests[asymptote_ests$Diversity=='Species richness',]
asymptote_ests$Site <- as.character(asymptote_ests$Site)
asymptote_ests <- cbind(asymptote_ests, inc_all$DataInfo$T)
asymptote_ests <- asymptote_ests[order(asymptote_ests$Site),]
colnames(asymptote_ests)[8] <- 'number of samples'
asymptote_ests$percent_completeness <- (asymptote_ests$Observed*100)/asymptote_ests$Estimator
asymptote_ests$N_samples_reqd <- (asymptote_ests$`number of samples`/asymptote_ests$percent_completeness)*100
asymptote_ests <- asymptote_ests[,c(1,3,4,9,8,10,5,6,7)]

## @knitr writing
write.csv(asymptote_ests, 'results/all_bats_inext.csv')
