#Script by Dave Hemprich-Bennett
#Script gives examples of the number of potential combinations of nodes
library(here)
library(ggplot2)
library(tidyverse)
library(plyr)
library(plotrix)



setwd(here())
node_combs <- c(10,12,14,16,18)
combs_selected <- c(5,6,7,8)

#factorial(node_combs)/(factorial(combs_selected)*(factorial(node_combs-combs_selected)))

possible_combs <- c()
n_available <- c()
n_used <- c()
for(i in 1:length(node_combs)){
  possible_combs <- c(possible_combs, factorial(node_combs[i])/(factorial(combs_selected)*(factorial(node_combs[i]-combs_selected))))
  n_available <- c(n_available, rep(node_combs[i], length(combs_selected)))
  n_used <- c(n_used, combs_selected)
}

df <- data.frame(n_available, n_used, possible_combs)
df$n_used <- as.factor(n_used)

palette <- c("#E69F00", "#56B4E9", "#009E73", 'black')

pdf('plots/Hice/rarefaction_example.pdf')
ggplot(df, aes(n_available, possible_combs, colour=n_used))+ geom_point()+ geom_line()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab('Number of possible combinations')+ xlab('Number of bats prior to rarefaction')+
  scale_color_manual(values=palette, name ='Number of nodes\nrarefied to')
                                                         
dev.off()
