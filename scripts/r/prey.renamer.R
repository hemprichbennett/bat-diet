prey.renamer<- function(input.matrix, prey.data, OTU.column, prey.column, collapse.rows=0, bitscore.column = 0, eval.column = 0, align.length.column = 0, min.bitscore = 0, min.eval = 0, min.align.length = 0, percent.ID.column=0, min.percent.ID=0){
  
if((bitscore.column == 0 & min.bitscore != 0)|(bitscore.column != 0 & min.bitscore ==0)){stop("Sorry, bitscore.column and min.bitscore must either both be blank or both have a value")}
if((eval.column == 0 & min.eval != 0)|(eval.column != 0 & min.eval ==0)){stop("Sorry, eval.column and min.eval must either both be blank or both have a value")}
if((align.length.column == 0 & min.align.length != 0)|(align.length.column != 0 & min.align.length ==0)){stop("Sorry, align.length.column and min.align.length must either both be blank or both have a value")}  
if((percent.ID.column == 0 & min.percent.ID != 0)|(percent.ID.column != 0 & min.percent.ID ==0)){stop("Sorry, percent.ID.column and min.percent.ID must either both be blank or both have a value")}    
  badreads <- 0    

#This finds any rows that didn't have a match, and removes them from temp_matrix_1
  temp_matrix_1 <- cbind(column_1 = rownames(input.matrix), input.matrix)
  rownames(temp_matrix_1) <- NULL
  OTUs_present <- c()
  for (i in 1:nrow(temp_matrix_1)){
    {if(temp_matrix_1[i,1] %in% prey.data[,OTU.column]){
      OTUs_present <- append(OTUs_present, i)
    }
      else{ badreads <- badreads+1}  
  }}
  temp_matrix_1 <- temp_matrix_1[OTUs_present,]
  #return(temp_matrix_1)
  
  #The following steps filter by your desired minimum eval
  if (eval.column !=0){
    desired_OTUs <- c()
    for (i in 1:nrow(temp_matrix_1)){
    position <- which(prey.data[,OTU.column] == temp_matrix_1[i,1])
    {if (prey.data[position, eval.column] >= min.eval){
      desired_OTUs <- append(desired_OTUs, i)}
    else{ badreads <- badreads+1}
    }
  }
  temp_matrix_1 <- temp_matrix_1[desired_OTUs,]
  }
  
  #The following steps filter by your desired minimum bitscore
  if (bitscore.column !=0){
    desired_OTUs <- c()
    for (i in 1:nrow(temp_matrix_1)){
      position <- which(prey.data[,OTU.column] == temp_matrix_1[i,1])
      {if (prey.data[position, bitscore.column] >= min.bitscore){
        desired_OTUs <- append(desired_OTUs, i)}
      else{ badreads <- badreads+1}
      }
    }
    temp_matrix_1 <- temp_matrix_1[desired_OTUs,]
  }
  
  #The following steps filter by your desired minimum align.length
  if (align.length.column !=0){
    desired_OTUs <- c()
    for (i in 1:nrow(temp_matrix_1)){
      position <- which(prey.data[,OTU.column] == temp_matrix_1[i,1])
      {if (prey.data[position, align.length.column] >= min.align.length){
        desired_OTUs <- append(desired_OTUs, i)}
      else{ badreads <- badreads+1}
      }
    }
    temp_matrix_1 <- temp_matrix_1[desired_OTUs,]
  }
  
  #The following steps filter by your desired minimum percent.ID
  if (percent.ID.column !=0){
    desired_OTUs <- c()
    for (i in 1:nrow(temp_matrix_1)){
      position <- which(prey.data[,OTU.column] == temp_matrix_1[i,1])
      {if (prey.data[position, percent.ID.column] >= min.percent.ID){
        desired_OTUs <- append(desired_OTUs, i)}
      else{ badreads <- badreads+1}
      }
    }
    temp_matrix_1 <- temp_matrix_1[desired_OTUs,]
  }
  print(paste0("a total of ", badreads, " OTUs were removed from your table"))
    
  #This makes a matrix where all of the OTUS are renamed to have their taxonomic information on
    for (i in 1: nrow(temp_matrix_1)){
    position <- which(prey.data[,OTU.column] == temp_matrix_1[i,1])
    temp_matrix_1[i,1] <- prey.data[position,prey.column]}
  

  
  

  temp_matrix_2 <- matrix(nrow= nrow(temp_matrix_1), ncol = ncol(temp_matrix_1))
  index2 <- 1  
  if(collapse.rows==TRUE){
  for (i in 1:nrow(temp_matrix_1)){
  if(temp_matrix_1[i,1] %in% temp_matrix_2[,1]){
    rowno <- which(temp_matrix_2[,1]==temp_matrix_1[i,1])
    temp_matrix_2[rowno,2:ncol(temp_matrix_2)] <- (as.numeric((temp_matrix_2[rowno,2:ncol(temp_matrix_2)]))) + as.numeric(temp_matrix_1[i,2:ncol(temp_matrix_1)])
  }
  else{temp_matrix_2[index2,] <- temp_matrix_1[i,]}
    index2 <- index2+1
  
  
  }}#This combines the rows by the taxa that they have in column 1, if desired

#return(temp_matrix_2)
  
  #This is all just formatting stuff
  if(collapse.rows ==TRUE){temp_matrix_2 <- temp_matrix_2[complete.cases(temp_matrix_2),]}  
if(collapse.rows==TRUE){colnames(temp_matrix_2) <- colnames(temp_matrix_1)
                        temp_matrix_2 <- temp_matrix_2[1:(length(unique(temp_matrix_2[,1]))),]}
                      
  
if(collapse.rows==FALSE){temp_matrix_2 <- temp_matrix_1}

output_matrix <- matrix(as.numeric(unlist(temp_matrix_2[,2:ncol(temp_matrix_2)])), nrow = nrow(temp_matrix_2))
rownames(output_matrix) <- temp_matrix_2[,1]
colnamevector <- colnames(temp_matrix_2)
colnamevector <- colnamevector[2:length(colnamevector)]
colnames(output_matrix) <- colnamevector

return(output_matrix)
}