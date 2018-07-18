#### Header ####
## Project: Bat diet
## Script purpose: plotting the positive and negative correlations of bat prey
## Date: 04/07/2018
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(corrplot)
library(ggplot2)
library(reshape2)

all_ecology <- read.csv('data/output_data/all_bats/sitewise_all_individual_info.csv')
all_ecology <- all_ecology[,c(seq(24,40), 59)]


unique(all_ecology$Species)

neg_sigs <- data.frame('sp'=c(), 'order1'= c(), 'order2' = c(), 'strength' = c())
for(i in 1:length(unique(all_ecology$Species))){
  sp <- unique(all_ecology$Species)[i]
  min <- 5
  
  taxa_mat <- all_ecology[which(all_ecology$Species== sp),]
  taxa_mat$Species <- NULL
  
  for_bigmat <- taxa_mat[,which(colSums(taxa_mat)>min)]
  if(0 %in% colSums(for_bigmat)){
    for_bigmat <- for_bigmat[,-which(colSums(for_bigmat)==0)]
  }
  
  big_cor <- for_bigmat
  bigcormat <- round(cor(big_cor),2)
  resbig <- cor.mtest(for_bigmat) 
  
  # make a melted dataframe of only the significant correlations
  sig_cors <- melt(bigcormat)[which(resbig$p <=0.05),]
  sig_cors <- sig_cors[-which(sig_cors[,1] == sig_cors[,2]),] #Remove the self-matches
  colnames(sig_cors) <- c('order1', 'order2', 'strength')
  if(length(which(sig_cors$strength < 0))){
    sig_cors <- cbind(rep(sp, nrow(sig_cors)), sig_cors)
    colnames(sig_cors)[1] <- 'sp'
    neg_sigs <- rbind(neg_sigs, sig_cors[which(sig_cors$strength < 0),])
  }
}

