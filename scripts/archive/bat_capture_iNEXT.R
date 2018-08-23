###Plotting the accumulation of bats in Sabah
###I"m pretty sure this is a logical error atm, the x axes of the facets are wrong!
library(stringr)
library('here')
library('bipartite')
library(iNEXT)
library(ggplot2)

setwd(here())

field_data <- read.csv('data/Edited_all_files_with_lat_long_VKedits.csv')
field_data$Site <- gsub('DANUM', 'Danum', field_data$Site)
field_data$Site <- gsub('DVCA', 'Danum', field_data$Site)
field_data$Site <- gsub('MALIAU', 'Maliau', field_data$Site)
field_data$Site <- gsub('MALUA', 'SBE', field_data$Site)
field_data$Species <- gsub('Nyja', 'Nytr', field_data$Species)
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
trap_effort <- read.csv('data/trapping_effort.csv', header = T, stringsAsFactors = T)

effort <- data.frame(trap_campaign=c('SAFE, 2015', 'SAFE, 2017', 'SAFE, 2017', 'Danum, 2016', 'Danum, 2017', 'Maliau, 2016', 'Maliau, 2017', 'SBE, 2016'),
             effort=c(sum(trap_effort$Number.of.traps[which(trap_effort$SAFE_data.Year=='2015')]),
               sum(trap_effort$Number.of.traps[which(trap_effort$SAFE_data.Year=='2016')]),
              sum(trap_effort$Number.of.traps[which(trap_effort$SAFE_data.Year=='2017')]),
              60, 60, 60, 60, 60))

effort <- effort[order(effort$trap_campaign),]

sp_table <- table(field_data$Species, field_data$SiteAndYear)
badsp <- c('0', 'Em--', 'Hi--', 'Hici', 'Ke--', 'Mi--', 'Mu--', 'My--', 'Pi--', 'Pip?', 'Rh--', 'Unknown', 'Scku')
sp_table <- sp_table[-which(rownames(sp_table) %in% badsp),]
badcols <- c('DV88, 2016', 'DV89, 2016')
sp_table <- sp_table[,-which(colnames(sp_table) %in% badcols)]
sp_table

#Add the effort columns as the first row of the table
#sp_table <- rbind(effort$effort, sp_table)

inext_list <- list()
for(i in 1:ncol(sp_table)){
  inext_list[[i]] <- as.numeric(sp_table[which(sp_table[,i]>0),i])
}
names(inext_list) <- colnames(sp_table)

effortvec <- c()

for(i in 1:length(inext_list)){
  effortvec <- c(effortvec, rep(effort$effort[i], length(inext_list[[i]])))
} 

#There is something wrong with how I'm passing sample effort to iNEXT

a <- iNEXT(inext_list, datatype = 'abundance', size = effortvec)


z <- fortify.iNEXT(a)
z$site_type <- rep(NA, nrow(z))
z[grep('SAFE', z$site),'site_type'] <- 'Logged'
z[grep('SBE', z$site),'site_type'] <- 'Logged, replanted'
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

g <- ggplot(z, aes_string(x="x", y="y", color = 'site_type')) + 
  geom_point(size=3, data=data.sub)+
  ylab('Species diversity') + xlab('Number of trap nights')+
  geom_line(aes_string(linetype="lty"), lwd=1.5)+
  scale_color_manual(values=palette)



g <- g +
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #Get rid of a load of the default crap
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=8))


g <- g + facet_wrap( ~ site, nrow=2, scales = 'free_x', strip.position = 'top')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"))#strip stuff sorts the facet labels, spacing adjusts the space between facets
g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr"), alpha=0.2)
g
