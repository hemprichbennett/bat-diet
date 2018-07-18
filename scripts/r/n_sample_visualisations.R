#### Header ####
## Project: Bat diet
## Script purpose: Calculating the number of individual samples going into each network per species
## Date: 29/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################
library(reshape2)
library(ggplot2)

all_ecology <- read.csv('data/output_data/all_bats/sitewise_all_individual_info.csv')

all_ecology$Site <- gsub('MALUA', 'SBE', all_ecology$Site)

######make a summary table of the number of samples per species and site######

sample_table <- table(all_ecology$Site, all_ecology$Species)
sample_table <- t(sample_table)


write.csv(sample_table, 'results/sample_table.csv')


chisq <- chisq.test(sample_table[,seq(1,3)])
chisq
chisq$expected

percentage_totals <- sample_table
for(i in 1:ncol(sample_table)){
  percentage_totals[,i] <- sample_table[,i] /colSums(sample_table)[i] *100
  print(sum(percentage_totals[,i]))
}

sum(percentage_totals)


percentage_totals

colSums(percentage_totals[grepl('Hi', rownames(percentage_totals)),])
colSums(percentage_totals[grepl('Ke', rownames(percentage_totals)),])
colSums(percentage_totals[grepl('Rh', rownames(percentage_totals)),])

#####Plot it #####

molten <- melt(sample_table)
colnames(molten) <- c('Species', 'Site', 'Count')

cbPalette <- c("#4bae8f",
               "#8a67d2",
               "#70b540",
               "#c05ea4",
               "#6e944a",
               "#cb5570",
               "#c59931",
               "#718cca",
               "#cd5033",
               "#b87c51")

sp_bars <- ggplot() + geom_bar(aes(x = Site, y = Count, fill = Species), data = molten,
                    stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y= 'Number of samples')+
  scale_fill_manual(values = cbPalette)
sp_bars

gen_bars <- ggplot() + geom_bar(aes(x = Site, y = Count, fill = substr(Species, start =1, stop = 2)), data = molten,
                               stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y= 'Number of samples')+
  scale_color_manual(values = cbPalette)
gen_bars

balloons <- ggplot(data = molten[which(molten$Count !=0),], aes(x = Site, y =Species)) + geom_point(aes(size=Count))+
  theme(panel.background=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill='Proportion of MOTU present',
       x ="Site and year", y = 'Prey taxa')
balloons

library(forcats)
tiles <- ggplot(data = molten, aes(y = fct_rev(Species), x =Site)) + geom_tile(aes(fill=Count), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 184/2, name = 'Number of\nsbats') + 
  labs(fill='Number of bats',
       x ="Site", y = 'Bat species')+
  theme(panel.background=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
tiles

