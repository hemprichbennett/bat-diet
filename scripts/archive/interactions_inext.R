
##################################################
## Project: bat-diet (all bats)
## Script purpose: using iNEXT to assess the sampling completeness and total diversity of each of my
## sampling events FOR INTERACTIONS, not species
## Date: 11/06/18
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



#Make list of networks
nets <- r_network_gen(collapse_species = T, include_malua = F, lulu = T)

names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))

n_nodes <- c(115, 119, 84, 80, 83, 150, 105)
names(n_nodes) <- names(nets)
#Paste together the columns of each network into a vector, so we have a vector of unique interactions
interactionslist <- lapply(nets, function(x) as.numeric(sapply(x, paste)))
#Add the number of samples used
for(i in 1:length(interactionslist)){
  interactionslist[[names(interactionslist)[i]]] <- c(n_nodes[which(names(n_nodes)== names(interactionslist)[i])], interactionslist[[i]])
}





a <- iNEXT(interactionslist, datatype = 'incidence_freq')




#####Here I replicate the gginext command, but tweak it a bit because it doesn't natively allow it ####

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
  ylab('interaction diversity') + xlab('Number of bats sampled')+
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

g
pdf('plots/inext/all_interaction_inext.pdf')
g
dev.off()

asymptote_ests <- a$AsyEst
asymptote_ests <- asymptote_ests[asymptote_ests$Diversity=='Species richness',]
asymptote_ests$Site <- as.character(asymptote_ests$Site)
asymptote_ests <- cbind(asymptote_ests, a$DataInfo$T)
asymptote_ests <- asymptote_ests[order(asymptote_ests$Site),]
colnames(asymptote_ests)[8] <- 'number of samples'
asymptote_ests$percent_completeness <- (asymptote_ests$Observed*100)/asymptote_ests$Estimator
asymptote_ests$N_samples_reqd <- (asymptote_ests$`number of samples`/asymptote_ests$percent_completeness)*100
asymptote_ests <- asymptote_ests[,c(1,3,4,9,8,10,5,6,7)]


write.csv(asymptote_ests, 'results/interactions_inext.csv')

b <- interactionslist$`SAFE, 2015`[2:length(interactionslist$`SAFE, 2015`)]
b <- ifelse(b, 0,1)

