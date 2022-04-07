###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Calculate percentage
###############################################################################################################


IDIRT_Template$set("public","calculatePercentage", function(){
  
  metadata <- self$metadata
  metadata <- metadata[1:38,]
  metadata[!grepl("mix", metadata$Experimrnt_label),"Label"] <- "Heavy"
  metadata[grepl("mix", metadata$Experimrnt_label),"Label"] <- "Mix"
  colnames(metadata)[2] <- "Experiment_label"
  self$metadata <- metadata
  summary <- self$summary
  percentages <- list()
  proteinGroups <- self$proteinGroups
  
  cat("Calculating the percentage")
  
  
  y <- apply(proteinGroups[,paste0("LFQ.intensity.H.",metadata$Expriment[metadata$Label=="Mix"])],1,mean)
  x <- apply(proteinGroups[,paste0("LFQ.intensity.L.",metadata$Expriment[metadata$Label=="Mix"])],1,mean)
  dy <- data.frame(gene =proteinGroups$Majority.protein.IDs , x=x,y=y)
  
  linearMod <- lm(y ~ x, data=dy[!dy$gene %in% c("Q9UN81","O00370"),])
  linearMod
    
  percentages<- proteinGroups[,c("Protein.IDs","Fasta.headers")]
  percentages[["Gene.names"]] <- apply(percentages, 1,function(x) ifelse(!is.na(strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]),
                                                                             paste0(c(strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3],strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]), collapse ="-"),
                                                                             strsplit(strsplit(x["Fasta.headers"],"GN=")[[1]][2]," ")[[1]][1]))
    
  percentages[["uniprotID"]] <- apply(percentages,1, function(x) strsplit(x[["Protein.IDs"]],";")[[1]][1])
    
    
  if(length(metadata$Expriment[metadata$Label == "Heavy"])>0){
    
    per <- apply(proteinGroups[,colnames(proteinGroups) %in% sapply(metadata$Expriment[metadata$Label == "Heavy"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round((y/coef(linearMod)[2])/(1+(y/coef(linearMod)[2])),3),0))
    per[proteinGroups$uniprotID %in% c("Q9UN81","O00370"),] <- apply(proteinGroups[proteinGroups$uniprotID %in% c("Q9UN81","O00370"),colnames(proteinGroups) %in% sapply(metadata$Expriment[metadata$Label %in% c("Heavy")], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(y/(1+y),3),0))
    percentages <- cbind(percentages,per)
      
  }
  
  if(length(metadata$Expriment[metadata$Label == "Mix"])>0){
    per <- apply(proteinGroups[,colnames(proteinGroups) %in% sapply(metadata$Expriment[metadata$Label == "Mix"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round((y/coef(linearMod)[2])/(1+(y/coef(linearMod)[2])),3),0))
    per[proteinGroups$uniprotID %in% c("Q9UN81","O00370"),] <- apply(proteinGroups[proteinGroups$uniprotID %in% c("Q9UN81","O00370"),colnames(proteinGroups) %in% sapply(metadata$Expriment[metadata$Label %in% c("Mix")], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(y/(1+y),3),0))
    percentages <- cbind(percentages,per)
  }
  
  if(length(metadata$Expriment[metadata$Label == "Light"])>0){
    
    per <- apply(proteinGroups[,colnames(proteinGroups) %in% sapply(metadata$Expriment[metadata$Label == "Light"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(1/(1+(y/coef(linearMod)[2])),3),0))
    per[proteinGroups$Majority.protein.IDs %in% c("Q9UN81","O00370"),] <- apply(proteinGroups[proteinGroups$Majority.protein.IDs %in% c("Q9UN81","O00370"),colnames(proteinGroups) %in% sapply(metadata$Expriment[metadata$Label == "Light"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(1/(1+y),3),0))
    percentages <- cbind(percentages,per)
    
  }
    
   colnames(percentages)[grepl('Ratio', colnames(percentages))] <- gsub("Ratio","Percentage",colnames(percentages)[grepl('Ratio', colnames(percentages))])

  self$percentages <- percentages
  
  invisible(self)
})

