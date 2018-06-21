## @knitr betalink_setup_and_plot
##################################################
## Project: bat diet (all species)
## Script purpose: checking the beta-diversity of my networks
## Date: 28/05/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

dir <- getwd()
basedir <- strsplit(dir, split ='/')[[1]][2]
#print(basedir)
if(grepl('data', basedir)){
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  
}else{
  library('here')
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  library(betalink)
  library(forcats)
  library(DataExplorer)
  library(gridExtra)
}

setwd(here())
source('scripts/r/r_network_gen.r')


nets <- r_network_gen(collapse_species = T, filter_species = T, lulu = T, include_malua = F, split_by = 'site and year')
names(nets) <- gsub('DANUM', 'Danum', names(nets))
names(nets) <- gsub('MALIAU', 'Maliau', names(nets))

write.csv(unlist(lapply(nets, sum)), 'results/n_links.csv')

graphs <- prepare_networks(nets)


beta <- network_betadiversity(graphs)
#Making a new object to allow prettier plotting
temp <- beta
colnames(temp)[c(1,2)] <- c('j', 'i')
for_plot <- rbind(beta, temp)
#Order the factors so that the plot looks nice
for_plot$i <- ordered(for_plot$i, levels = c("Danum, 2016", "Danum, 2017", "Maliau, 2016", 'Maliau, 2017',
                                             'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017'))
for_plot$j <- ordered(for_plot$j, levels = c("Danum, 2016", "Danum, 2017", "Maliau, 2016", 'Maliau, 2017',
                                             'SAFE, 2015', 'SAFE, 2016', 'SAFE, 2017'))
for_plot$STWN <- for_plot$ST/for_plot$WN
melted_forplot <- melt(for_plot)


betaplot <- ggplot(melted_forplot, aes(i, fct_rev(j)))+ geom_point(aes(size=value, colour = value))+
  scale_colour_gradient(low = "white",
                        high = "blue", limits = c(0,1))+
  theme_dark()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')+
  facet_wrap(~ variable, labeller = label_bquote(italic(beta [.(as.character(variable))])),
             nrow = 3)+
  labs(x= NULL, y = NULL)



betaplot  

## @knitr betalink_extras

pdf('plots/beta/betaplot.pdf')
betaplot
dev.off()

jpeg('plots/beta/betaplot.jpg' , units = 'in', width = 9, height = 9, res=300)
betaplot
dev.off()


S_WN <- ggplot(for_plot, aes(S, WN))+ geom_point()+ geom_abline(intercept = 0, slope = 1, linetype="dotted")+ 
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab(expression(italic(beta[S])))+ ylab(expression(italic(beta[WN])))

S_OS <- ggplot(for_plot, aes(S, OS))+ geom_point()+ geom_abline(intercept = 0, slope = 1, linetype="dotted")+ 
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab(expression(italic(beta[S])))+ ylab(expression(italic(beta[OS])))

S_ST <- ggplot(for_plot, aes(S, ST))+ geom_point()+ geom_abline(intercept = 0, slope = 1, linetype="dotted")+ 
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab(expression(italic(beta[S])))+ ylab(expression(italic(beta[ST])))

WN_OS <- ggplot(for_plot, aes(WN, OS))+ geom_point()+ geom_abline(intercept = 0, slope = 1, linetype="dotted")+ 
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab(expression(italic(beta[WN])))+ ylab(expression(italic(beta[OS])))

pdf('plots/beta/betascatter.pdf')
grid.arrange(S_WN, S_OS, S_ST, WN_OS)
dev.off()

pdf('plots/beta/s_wn.pdf')
S_WN
dev.off()

pdf('plots/beta/s_os.pdf')
S_OS
dev.off()

pdf('plots/beta/s_st.pdf')
S_ST
dev.off()

pdf('plots/beta/wn_os.pdf')
WN_OS
dev.off()


sink('results/betadiversity/s_wn.txt')
summary(lm(S ~ WN, data = for_plot))
sink()

sink('results/betadiversity/s_os.txt')
summary(lm(S ~ OS, data = for_plot))
sink()

sink('results/betadiversity/s_st.txt')
summary(lm(S ~ ST, data = for_plot))
sink()

sink('results/betadiversity/wn_os.txt')
summary(lm(WN ~ OS, data = for_plot))
sink()
