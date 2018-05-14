
if(interactive()==TRUE){
  library('here')
  library(dplyr)
  library(ggplot2)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(dplyr, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

nets <- r_network_gen(collapse_species = T, filter_species = T)

####First, work out the sample sizes ####
#This gives a matrix with EVERY interaction in it
all_interactions <- r_network_gen(collapse_species = F, desired_species = NULL, include_malua = T, filter_species = F)
all_interactions <- all_interactions[,-1] #Get rid of the first column as it still has the rownames
#a lot of samples don't match our IDs in the field data (like 120), CHECK THEM
#desired_sites <- c("SAFE, 2015", "SAFE, 2016", 'SAFE, 2017', 'DANUM, 2016', 'DANUM, 2017', 'MALIAU, 2017')
#all_interactions <- all_interactions[,-which(!all_interactions[1,] %in% desired_sites)] #This gets rid of the bad samples, but sort this shit out!

field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'), stringsAsFactors = F)
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)

site <- c()
sp <- c()
for(i in 1:ncol(all_interactions)){
  
  if(all_interactions[2,i] %in% field_data$Faeces_no1){
    r <- which(field_data$Faeces_no1==all_interactions[2,i])
    species <- field_data[r, 'Species']
    if(length(species)>1){next()}#If two rows match this value, ignore it
    #print(species)
    sp <- c(sp, species)
    site <- c(site, all_interactions[1,i])
    }else if(all_interactions[2,1] %in% field_data$Faeces_no2){
      r <- which(field_data$Faeces_no2==all_interactions[2,i])
      species <- field_data[r, 'Species']
      #print(species)
      if(length(species)>1){next()}
      sp <- c(sp, species)
      site <- c(site, all_interactions[1,i])
      }
}

counts <- data.frame(site, sp)
counts <- melt(table(counts))

for(i in 1:length(nets)){
  print(names(nets)[i])
  print(colnames(nets[[i]]))
}



####Now work out the proportions of items in their diet ####
prey_data <- read.csv('data/taxonomy/order.csv')
colnames(prey_data) <- c('MOTU', 'Taxa')
prey_data$MOTU <- as.character(prey_data$MOTU)

taxa_list <- list()

#Put the taxonomic information onto the sites data

for(b in 1: length(nets)){
  all_interactions <- nets[[b]]
  taxa_mat <- matrix(nrow=0, ncol=ncol(all_interactions))
  
  z <- 1
  for(i in 1: nrow(all_interactions)){
    rowname = rownames(all_interactions)[i]
    if(!rowname %in% prey_data$MOTU){
      next()
    }
    tax = as.character(prey_data[which(prey_data$MOTU == rowname),'Taxa'])
    if(is.null(nrow(taxa_mat))){ #If its the first iteration there won't be any rownames yet, so the next if statement will fail
      taxa_mat <- rbind(taxa_mat, all_interactions[i,])
      rownames(taxa_mat)[z] <- tax
      z <- z+1
    }
    if(tax %in% rownames(taxa_mat)){
      to_merge = which(rownames(taxa_mat)==tax)
      taxa_mat[to_merge,] <- taxa_mat[to_merge,]+ all_interactions[i,]
    }else{
      taxa_mat <- rbind(taxa_mat, all_interactions[i,])
      rownames(taxa_mat)[z] <- tax
      z <- z+1
    }
  }
  taxa_list[[b]] <- taxa_mat
  names(taxa_list)[b] <- names(nets)[b]
}









####Convert to proportions #### 
prop_list <- lapply(taxa_list, prop.table, margin = 2) #margin=2 makes it so it generates the column totals
names(prop_list) <- names(taxa_list)


#Format the data a bit for the plot
molten_proportions <- melt(prop_list)
colnames(molten_proportions) <- c('Prey_taxa', 'Bat_sp', 'Proportion', 'Site_and_year')

#####Add the counts data to the molten_proportions data ####
molten_proportions$counts <- rep(NA, nrow(molten_proportions))
for(i in 1:nrow(molten_proportions)){
  bat = as.character(molten_proportions[i, 'Bat_sp'])
  net = as.character(molten_proportions[i, 'Site_and_year'])
  proprow = which(counts$sp==bat & counts$site==net)
  molten_proportions[i, 'counts'] <- counts[proprow,'value']
}

molten_proportions <- molten_proportions[with(molten_proportions, order(Prey_taxa)),]

#This complicated looking command removes all but the first value of each combination of Bat_sp and Site_and_year
first_values <- which(!duplicated(molten_proportions[,c('Bat_sp', 'Site_and_year')]))
molten_proportions[seq(1,nrow(molten_proportions))[-first_values],'counts'] <- ''


unwanted_sp <- c('Hebl', 'Hidy', 'Kepe', 'Muro','Nytr', 'Hidi')

#molten_proportions <- molten_proportions[-which(molten_proportions$Bat_sp %in% unwanted_sp),]




#For some reason the levels of the prey taxa aren't alphabetical. This changes that
molten_proportions$Prey_taxa <- ordered(molten_proportions$Prey_taxa, levels=unique(molten_proportions$Prey_taxa)[order(as.character(unique(molten_proportions$Prey_taxa)))])


#molten_proportions <- molten_proportions[-which(molten_proportions$Proportion==0),] #Doing this makes the table more readable but fucks up the plot

molten_proportions$Bat_sp <- gsub('Hice', 'Hipposideros cervinus', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Hiri', 'Hipposideros ridleyi', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Hidi', 'Hipposideros diadema', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Hidy', 'Hipposideros dyacorum', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Kein', 'Kerivoula intermedia', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Keha', 'Kerivoula hardwickii', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Kepa', 'Kerivoula papillosa', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Rhbo', 'Rhinolophus borneensis', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Rhse', 'Rhinolophus sedulus', molten_proportions$Bat_sp)
molten_proportions$Bat_sp <- gsub('Rhtr', 'Rhinolophus trifoliatus', molten_proportions$Bat_sp)

molten_proportions$Site_and_year <- gsub('DANUM', 'Danum', molten_proportions$Site_and_year)
molten_proportions$Site_and_year <- gsub('MALIAU', 'Maliau', molten_proportions$Site_and_year)


####Get plotting ####
palette <- c("#d05034",
             "#7f60d2",
             "#a6b63d",
             "#c053c0",
             "#5ab150",
             "#c14489",
             "#4cae90",
             "#e24971",
             "#76813a",
             "#667ad2",
             "#d49b3b",
             "#5d9bd3",
             "#9f612d",
             "#835897",
             "#e29473",
             "#d58dc9",
             "#b05560")
#palette <- c("#e2fcff","#f7abb6","#7cdcda","#f6b1cf","#b0e5ba","#e9b9e7","#ffffd5","#e1c4f9","#b2ad80","#a7a8d0","#f4cd9d","#abf3ff","#bda4b8","#fffae7","#92afcb","#ffe4d1","#7eb4c0")
plot <- ggplot(data = molten_proportions, aes(x = Site_and_year, y =Proportion, fill = Prey_taxa)) + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+ #Get rid of the annoying background formatting, rotate the x labels
  geom_bar(stat = 'identity') + facet_wrap( ~ Bat_sp, ncol = 2)+ theme(strip.background = element_rect(fill="white")) + scale_fill_manual(values =palette, name = 'Prey order')+
  xlab("Site and year")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,0.2))+
  geom_text(aes(label=counts), vjust=-0.3, size=3.5, position = position_stack(vjust = 1))
  
plot
pdf('plots/Site comparisons/species_diets.pdf') #The legend isn't alphabetical. This is annoying.
plot
dev.off()

png('plots/Site comparisons/species_diets.png') #The legend isn't alphabetical. This is annoying.
plot
dev.off()
 
hice_only <- molten_proportions[which(molten_proportions$Bat_sp=='Hipposideros cervinus'),]
hice_only <- hice_only[-which(hice_only$Proportion==0),] #Get rid of the empty rows or they mess up the legend
hice_plot <- ggplot(data = hice_only, aes(x = Site_and_year, y =Proportion, fill = Prey_taxa)) + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+ #Get rid of the annoying background formatting, rotate the x labels
  geom_bar(stat = 'identity') + scale_fill_manual(values =palette, name = 'Prey order')+
  xlab("Site and year")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,0.2))+
  geom_text(aes(label=counts), vjust=-0.3, size=3.5, position = position_stack(vjust = 1))

pdf('plots/Hice/hice_diet.pdf')
hice_plot
dev.off()

png('plots/Hice/hice_diet.png')
hice_plot
dev.off()


# 
# ####Now for the stats on the Hice
# 
# 
# colnames(hice_only)[3] <- 'value'
# hice_wide <- cast(hice_only, Site_and_year ~ Prey_taxa)
# hice_wide[is.na(hice_wide)] <- 0
# hice_wide <- t(hice_wide) #Has to be the other way round for the chi-squared test
# chisq.test(hice_wide) 
# 
# #now try it without the rarest 13 taxa
# badtaxa <- names(sort(rowSums(hice_wide)))[1:13]
# chisq.test(hice_wide[-which(rownames(hice_wide) %in% badtaxa)])
