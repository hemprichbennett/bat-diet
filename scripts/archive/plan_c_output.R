#### Header ####
## Project:
## Script purpose: making sense of the plan_c output
## Date: 12/09/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################


library(here)
library(reshape2)
library(ggplot2)
setwd(here())


allfiles <- list.files(path = 'results/rand_values_95/')

actuals <- allfiles[grep('actual', allfiles)]
rands <- allfiles[grep('rand', allfiles)]

actuals <- read.csv('results/rand_values_95/actual_functional complementarity.csv')
rand <- read.csv('results/rand_values_95/rand_functional complementarity.csv')

mets <- gsub('rand_', '', rands)
mets <- gsub('\\.csv', '', mets)

colnames(rand)

combns <- t(combn(colnames(rand), 2))

net1 <- c()
net2 <- c()
net1_actual <- c()
net2_actual <- c()
lower <- c()
upper <- c()
actual <- c()
signif <- c()
for(i in 1:nrow(combns)){
  rand_range <- quantile(rand[,combns[i,1]] - rand[,combns[i,2]], probs = c(0.025, 0.975))
  #Order it
  rand_range <- rand_range[order(rand_range)]
  actual_diff <- actuals[,combns[i,1]] - actuals[,combns[i,2]]
  
  if(actual_diff > rand_range[1] && actual_diff < rand_range[2]){
    sig <- F
  }else{sig <- T}
  
  net1 <- c(net1, combns[i,1])
  net2 <- c(net2, combns[i,2])
  net1_actual <- c(net1_actual, actuals[,combns[i,1]])
  net2_actual <- c(net2_actual, actuals[,combns[i,2]])
  lower <- c(lower, rand_range[1])
  upper <- c(upper, rand_range[2])
  actual <- c(actual, actual_diff)
  signif <- c(signif, sig)
  
}


out_df <- data.frame(net1, net2, net1_actual, net2_actual, signif, lower, upper, actual)

out_df


ggplot(out_df, aes(x = net1, y = net2, colour = signif)) + geom_tile() + scale_fill_discrete()
