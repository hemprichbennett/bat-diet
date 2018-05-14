
library(lulu)
library(here)
library(seqinr)
setwd(here())

otutab <- read.csv("data/processed_dna_data/lulu/all_post_QC_otus.txt.table.out",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table("data/processed_dna_data/lulu/match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)


####Using the binary matrix resulted in zero sequences being removed, 
####likely as the algorithm relies on the abundances of the OTUs.
####Lets try it with weighted data

curated_result <- lulu(otutab, matchlist)
#Number of OTUs retained out of original
cat('retained',curated_result$curated_count, 'of', nrow(otutab), 'otus')
#The analysis took
curated_result$runtime

save.image('data/output_data/lulu.RDS')
write.csv(curated_result$curated_table,'data/processed_dna_data/lulu/lulu_95.csv')

load('data/output_data/lulu.RDS')
#Lulu kept the following OTUs
kept <- curated_result$curated_otus

fas <- read.fasta('data/processed_dna_data/lulu/95/reps_95.fasta')

to_keep <- fas[which(names(fas) %in% kept)]

write.fasta(to_keep, names = names(to_keep), file.out = 'data/processed_dna_data/lulu/95/lulu_out.fasta')

