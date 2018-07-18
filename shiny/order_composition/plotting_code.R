library(here)
library(ggridges)
library(reshape2)
library(corrplot)

setwd(here())
taxa_mat<- read.csv('shiny/order_composition/order_mat.csv', row.names = 1) 


n_matches <- 0
for_bigmat <- t(taxa_mat[,which(colSums(taxa_mat)>n_matches)])
if(0 %in% colSums(for_bigmat)){
  for_bigmat <- for_bigmat[,-which(colSums(for_bigmat)==0)]
}

big_cor <- for_bigmat
bigcormat <- round(cor(big_cor),2)
resbig <- cor.mtest(for_bigmat)

####temp #####
n_matches <- 15
for_biggestmat <- t(taxa_mat[,which(colSums(taxa_mat)>n_matches)])
if(0 %in% colSums(for_biggestmat)){
  for_biggestmat <- for_biggestmat[,-which(colSums(for_biggestmat)==0)]
}

print(for_biggestmat)
bigger_cor <- for_biggestmat
biggercormat <- round(cor(bigger_cor),2)
resbiggest <- cor.mtest(for_biggestmat)
corrplot(biggercormat, method = "circle", p.mat = resbiggest$p, sig.level = .05, type = 'upper', order = 'AOE',
         tl.col = "black", tl.srt = 45, insig = 'blank',
         bg = "black", title= n_matches)
####end temp

molten <- melt(for_bigmat)
colnames(molten) <- c('sample', 'order', 'value')
ggplot(molten[which(molten$order %in% c('Lepidoptera', 'Diptera', 'Coleoptera')),], aes(x=sample, y = value, colour=order))+ geom_point()+
 theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
 geom_line()+
  labs(y='N MOTU')

ggplot() + geom_bar(aes(y = value, x = sample, fill = order), data = molten[which(molten$order %in% c('Lepidoptera', 'Diptera', 'Coleoptera')),],
                    stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# ggplot(molten[which(molten$order %in% c('Lepidoptera', 'Diptera', 'Coleoptera')),], aes(x= sample, y=order))+
#   geom_point()

ggplot(molten, aes(molten$value))+ facet_wrap(~ order, ncol = 3)+
  geom_histogram()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# ggplot(as.data.frame(for_bigmat), aes(x= Coleoptera, y=Diptera))+ geom_point()+
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   geom_smooth(method=lm, SE=T)
  

ggplot(molten, aes(x=order, y=value))+ geom_boxplot()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
hist(rowSums(big_cor))
