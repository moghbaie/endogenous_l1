###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Calculate percentage
###############################################################################################################


Template$set("public","calculatePercentage", function(){
  
  metadata <- self$metadata
  summary <- self$summary
  percentages <- list()
  proteinGroups <- self$proteinGroups
  
  cat("Calculating the percentage")
  
  for(study in names(metadata)){
    
    percentages[[study]] <- proteinGroups[[study]][,c("Protein.IDs","Fasta.headers")]
    percentages[[study]][["Gene.name"]] <- apply(percentages[[study]], 1,function(x) ifelse(!is.na(strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]),
                                                                                            paste0(c(strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3],strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]), collapse ="-"),
                                                                                            strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3]))
    
    if(length(metadata[[study]]$Name[metadata[[study]]$Label == "Heavy"])>0){
      percentages[[study]] <- cbind(percentages[[study]],
                                    apply(proteinGroups[[study]][,colnames(proteinGroups[[study]]) %in% sapply(metadata[[study]]$Name[metadata[[study]]$Label == "Heavy"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(y/(1+y),3),0))
      )
      
    }
    if(length(metadata[[study]]$Name[metadata[[study]]$Label == "Light"])>0){
      percentages[[study]] <- cbind(percentages[[study]],
                                    apply(proteinGroups[[study]][,colnames(proteinGroups[[study]]) %in% sapply(metadata[[study]]$Name[metadata[[study]]$Label == "Light"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(1/(1+y),3),0))
      )
    }
    
    colnames(percentages[[study]])[grepl('Ratio', colnames(percentages[[study]]))] <- gsub("Ratio","Percentage",colnames(percentages[[study]])[grepl('Ratio', colnames(percentages[[study]]))])
  }
  
  self$percentages <- percentages
  
  invisible(self)
})

## Merge mix data

Template$set("public","calculateMixPercentage", function(){
  
  dirs <- list.dirs(file.path(getwd(),self$input.dir), 
                    full.names = FALSE, 
                    recursive = FALSE)
  
  cat("Merging mix files.\n")
  metadata <- self$metadata
  mixdata <- self$mixdata 
  mixpercentages <- list()
  percentages <- self$percentages
  
  for(study in names(mixdata)){
    if(study=="LysoIP1"){
      mixpercentages[[study]] <- mixdata[[study]][,c("Protein.IDs","Fasta.headers")]
      mixpercentages[[study]][["Gene.name"]] <- apply(mixdata[[study]], 1,function(x) ifelse(!is.na(strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]),
                                                                                             paste0(c(strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3],strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]), collapse ="-"),
                                                                                             strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3]))
      
      names <- gsub("Ratio.H.L.normalized.|Ratio.H.L.normalized","",colnames(mixdata[[study]])[grepl("Ratio.H.L.normalized",colnames(mixdata[[study]]))])
      names <- names[names !=""]
      
      n <- dim(self$metadata[[study]])[1]
      metadatamix <- data.frame(`Tube#`= seq(n+1, n+length(names)), Name = names, Label =ifelse(grepl('a',names),"Heavy","Light"))
      colnames(metadatamix)[1] <- "Tube#"
      
      comment <- data.frame(num = c(1,2,3,4,5), treatment = c("Untreated","GA.treated","Torin1.treated","Torin1.PI.treated","PI.treated"))
      metadatamix[["Treatment"]] <- comment$treatment[match(substr(metadatamix[["Name"]], start = 1, stop =1), comment$num)]
      metadatamix[["Experiment_label"]] <- apply(metadatamix,1,function(x) paste0(c(x[["Treatment"]],"mix"), collapse="_"))
      metadatamix[["Base"]] <- NA
      
      if(length(metadatamix$Name[metadatamix$Label == "Heavy"])>0){
        mixpercentages[[study]] <- cbind( mixpercentages[[study]], 
                                          apply(mixdata[[study]][,colnames(mixdata[[study]]) %in% sapply(metadatamix$Name[metadatamix$Label == "Heavy"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(y/(1+y),3),0))
        )
      }
      
      if(length(metadatamix$Name[metadatamix$Label == "Light"])>0){
        mixpercentages[[study]] <- cbind( mixpercentages[[study]],
                                          apply(mixdata[[study]][,colnames(mixdata[[study]]) %in% sapply(metadatamix$Name[metadatamix$Label == "Light"], function(x) paste0("Ratio.H.L.",as.character(x)))],2, function(y) ifelse(!is.na(y), round(1/(1+y),3),0))
        )
      }
      
      metadata[[study]][["Base"]] <- metadatamix$Experiment_label[match(substr(metadata[[study]][["Name"]],1,2),metadatamix$Name)]
      
      metadata[[study]] <- rbind(metadata[[study]],metadatamix)
      
      self$mixpercentages[[study]] <- mixpercentages[[study]]
      
      HUMAN_9606_idmapping <- read.delim("C:/Users/Mehrnoosh/Downloads/HUMAN_9606_idmapping (3).dat", header=FALSE)
      percentages[[study]][["uniprotID"]] <- apply(percentages[[study]],1, function(x) strsplit(x[["Protein.IDs"]],"\\|")[[1]][2])
      mixpercentages[[study]][["uniprotID"]] <- apply(mixpercentages[[study]],1, function(x) strsplit(x[["Protein.IDs"]],"\\|")[[1]][2])
      colnames(mixpercentages[[study]]) <- gsub("Ratio.H.L.","Percentage.H.L.",colnames(mixpercentages[[study]]))
      
      
      percentages[[study]] <- merge(percentages[[study]][,!colnames(percentages[[study]]) %in% c("Protein.IDs","Fasta.headers")],mixpercentages[[study]][,!colnames(mixpercentages[[study]]) %in% c("Protein.IDs","Fasta.headers")],by = c("uniprotID","Gene.name"), all= TRUE) 
      self$percentages[[study]] <- percentages[[study]]
      
      
      self$metadata[[study]] <- metadata[[study]]
      
    }
  }
  
  invisible(self)
  
})