rm(list=ls())
#To be ran in an array
if(interactive()==TRUE){
  library('here')
  library(here)
  library(bold)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
  library(httpcode, lib.loc = '/data/home/btw863/r_packages/')
  library(urltools, lib.loc = '/data/home/btw863/r_packages/')
  library(bold, lib.loc = '/data/home/btw863/r_packages/')
  library(dplyr, lib.loc = '/data/home/btw863/r_packages/')
  library(ggplot2, lib.loc = '/data/home/btw863/r_packages/')
  library(reshape2, lib.loc = '/data/home/btw863/r_packages/')
  library(data.table, lib.loc = '/data/home/btw863/r_packages/')
}
setwd(here())

args = commandArgs(trailingOnly=TRUE)

cat(args, '\n')

args <- as.numeric(args)

mydata <- read.table("data/processed_dna_data/galaxy_r_workflow/95/reps_95.tsv", header = F, sep = "\t", dec = ".", stringsAsFactors = FALSE)
colnames(mydata) <- c('seqID','seqs')

head(mydata)

mydata <- mydata[seq(nrow(mydata)/23*args, nrow(mydata)/23*(args+1)),] #Select the appropriate subset for the array

#net <- read.table('swarmOTUs/2014_biom_no_singletons.txt', sep = '\t', header = T, row.names = 1)

#make sure headers are not capitalized. We need to use this command to make a named list of the sequences, 
#otherwise our results from BOLD will not give the name of the OTU they correspond to
mydata2 <- as.list(setNames(mydata$seqs, mydata$seqID))  

t <- Sys.time()
output <- bold_identify(sequences = mydata2, db = "COX1_SPECIES", response=FALSE) #This can take several hours to run
t_after <- Sys.time()

cat('generating output took ', t_after - t, '\n')

outtax40 <- lapply(output, head, n=40)
outtax1 <- lapply(output, head, n = 1)
outtax1_df <- do.call('rbind', outtax1)

#a <- bold_identify_parents(output, wide = T)

outtaxframe <- do.call("rbind", lapply(outtax40, data.frame))
write.csv(outtaxframe, paste('data/taxonomic_info', args,'.csv', sep = ''))

# parents_list <- list()
# a <- 1
# for(i in 1:length(output)){
#   print(i)
#   #print(names(output[[1]][i]))
#   if(is.null(output[[i]])){#If a sequence is null, it means that BOLD found no vaguely useful matches for it. I.e. Its a load of crap. Also, if I don't include this if statement the NULLs will crash the code
#     next()
#   }
#   parents_list[[a]] <- bold_identify_parents(outtax40[[i]], wide = T)
#   names(parents_list)[a] <- names(output)[i]
#   a <- a+1
# }
# 
# temp <- parents_list[[1]][[1]][1,]
# temp_rbinded <- do.call("rbind", lapply(parents_list, data.frame))
# 
# tax_df <- matrix(ncol = 2, nrow = 0)
# for(i in 1:length(parents_list)){
#   nam <- str_split(names(parents_list)[[i]], pattern = ' ')[[1]][1]
#   vec <- as.character(parents_list[[i]][[1]][1,"order"])
#   tax_df <- rbind(tax_df, c(nam, vec))
# }
# 
# badrows <- c()
# net_2 <- cbind(rep(NA,nrow(net)), net)
# 
# for(i in 1:nrow(net)){
#   if(rownames(net)[i] %in% tax_df[,1]){
#     r <- which(tax_df[,1]== rownames(net)[i])
#     net_2[i,1] <- tax_df[r,2]
#   }else{
#     badrows <- c(badrows, i)
#   }
# }
# 
# net_2 <- net_2[-badrows,]
# #This summarises the bats by the orders they feed on
# proportions <- net_2 %>% group_by(`rep(NA, nrow(net))`) %>% summarise_all(funs(sum))
# melted_proportions <- melt(proportions)
# colnames(melted_proportions) <- c('Order', 'Bat', 'consumption')
# #melted_proportions$order <- as.factor(melted_proportions$order)
# 
# #barplot(proportions,  ylim = c(0,1),ylab = "Proportion of species in diet", col = setall, las= 3)
# 
# ggplot(data = melted_proportions, 
#        aes(x = Bat, y = consumption, fill = Order))+ geom_bar(stat = 'identity')
# 
