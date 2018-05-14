hernani_comparisons <- function(net1, net2, metric, part = NA, nperm =999, sums_to_preserve = 'both'){
  #Code originally by Hernani Oliveira, then altered to make an R function by Dave-Hemprich-Bennett
  
  
  #network_1 is your first network
  #network_2 is your second network
  #metric is the metric you're after
  #part is the network level you're after. Defaults to NA
  #nperm is the number of permutations you want to do
  if(is.na(part)){
    
    participant_metric <-  F
  }else{
    
    participant_metric <-  T
  }
  
  if(metric != 'quanbimo'){
    
    if(participant_metric==FALSE){
      
      if(metric == 'dirt'){
        dirt1 <- computeModules(net1)
        dirt2 <- computeModules(net2)
        orig = abs(slot(dirt1, 'likelihood')- slot(dirt1, 'likelihood'))
      }else{orig=abs(networklevel(net1,index=metric)-networklevel(net2,index=metric))}
      
      
      
    }
    if(participant_metric ==TRUE){
      
      orig=abs(networklevel(net1,index=metric, level = part)-networklevel(net2,index=metric, level = part))
      
    }
    
  }else{
    if(metric == 'quanbimo'){
      quanbimo1 <- computeModules(net1, method='DormannStrauss')
      quanbimo2 <- computeModules(net2, method='DormannStrauss')
      orig = abs(slot(quanbimo1, 'likelihood')- slot(quanbimo2, 'likelihood'))
    }
    
  }
  
  
  #orig
  
  #Calculate the differences in the metric chosen among the randomized networks.
  #Generally gets reliable from 999 permutations.
  #nperm = 999
  i=1
  randomized.patef=matrix(nrow=length(orig),ncol=nperm+1)
  row.names(randomized.patef)=names(orig)
  randomized.patef[,1]=orig
  
  while(i <=nperm){ 
    
    data_aleat=permatfull(net1,fixedmar=sums_to_preserve,mtype="count",times=1)
    data_aleat=data_aleat$perm[[1]]
    net2_aleat=permatfull(net2,fixedmar=sums_to_preserve,mtype="count",times=1)
    net2_aleat=net2_aleat$perm[[1]]
    if(metric != 'quanbimo'){
      if(participant_metric==TRUE){
        linha<-abs(networklevel(data_aleat, index=metric, level = part)-networklevel(net2_aleat, index=metric, level = part))
      }else{linha<-abs(networklevel(data_aleat, index=metric)-networklevel(net2_aleat, index=metric))}
      
    }else{
      if(metric == 'dirt'){
        fake_dirt1 <- computeModules(data_aleat)
        fake_dirt2 <- computeModules(net2_aleat)
        linha = abs(slot(fake_dirt1, 'likelihood')- slot(fake_dirt1, 'likelihood'))
      }
      if(metric == 'quanbimo'){
        fake_quanbimo1 <- computeModules(data_aleat, method='DormannStrauss')
        fake_quanbimo2 <- computeModules(net2_aleat, method='DormannStrauss')
        linha = abs(slot(fake_quanbimo1, 'likelihood')- slot(fake_quanbimo2, 'likelihood'))
      }
      
    }
    randomized.patef[,i+1]=linha
    #print(i)
    i=i+1
    
  } 
  
  #randomized.patef
  
  # Plot and export the comparison between the observed value and the random metric distribution
  niveis<-row.names(randomized.patef)

  
  
  # Calculate the proportion of the randomized differences that was greater than
  #the difference between the original networks
  significance.patef=matrix(nrow=nrow(randomized.patef),ncol=3)
  row.names(significance.patef)=row.names(randomized.patef)
  colnames(significance.patef)=c("p (rand <= orig)", "p (rand >= orig)", "p (rand=orig)")
  
  
  
  signif.sup=function(x) sum(x>=x[1])/length(x) #unicaudal
  signif.inf=function(x) sum(x<=x[1])/length(x) #unicaudal
  signif.two=function(x) ifelse(min(x)*2>1,1,min(x)*2)#bicaudal
  # 
  # signif.sup=function(x) sum(x>=x[1])/length(x)
  # signif.inf=function(x) sum(x<=x[1])/length(x)
  # signif.two=function(x) ifelse(min(x)*2>1,1,min(x)*2)
  
  significance.patef[,1]=apply(randomized.patef,1,signif.inf)
  significance.patef[,2]=apply(randomized.patef,1,signif.sup)
  
  significance.patef2<-data.frame(significance.patef)
  significance.patef2[,3]=apply(significance.patef2[,-3],1,signif.two)
  
  colnames(significance.patef2)=c("p (rand <= orig)", "p (rand >= orig)", "p (rand=orig)")
  
  significance.patef
  
  # #Export Results
  # if( participant_metric == TRUE){
  #   file_name= paste(as.numeric(args),"output_",net1_name, '_', net2_name, '_', clust_level,'_', metric,part,"_nperm.txt", sep="")
  # }
  # if(participant_metric == FALSE){
  #   file_name= paste(as.numeric(args),"output_",net1_name, '_', net2_name, '_', clust_level,'_', metric,"_nperm.txt", sep="")
  # }
  # 
  # 
  # print(file_name)
  return(significance.patef)
  
  # write.table(significance.patef, file=file_name, 
  #             sep=" ",row.names=TRUE,col.names=TRUE)
  # if(participant_metric == TRUE){
  #   write.table(orig, file=paste(as.numeric(args),"orig_",net1_name, '_', net2_name, '_', clust_level,'_', metric,part,"_nperm.txt", sep=""))
  # }
  # if(participant_metric == FALSE){
  #   write.table(orig, file=paste(as.numeric(args),"orig_",net1_name, '_', net2_name, '_', clust_level,'_', metric,"_nperm.txt", sep=""))
  #   
  # }
  
  #So you will run this script comparing your network and mine for linkage density (etc) at each OTU cut off (you'll run this over and over recording the outputs). The idea will be to end up with 4 graphs, one for each metric. In each we would have a series of point corresponding to the metric as calcualted for your network and mine at each OTU level (just as you showed us yesterday). This will show how the metrics change for each OTU level for two datasets. Ideally we will see they are different but consistently so. Then with this script we can actually say what we would have concluded if we were assessing whether the two food webs were ecologically different. We can indicate on the figures whcih cases would have come out as signifiantly different. Hope that makes sense... if you run into problems we can skype. Maybe we can meet on Tuesday to see how it comes out and work out what will go in your talk for Thursday?
  end_time <- Sys.time()
  duration <- end_time - start_time
  cat(metric,' took ', duration, ' to finish \n', sep= '')
  #}
  
  
  
}
