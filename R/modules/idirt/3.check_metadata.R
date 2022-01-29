###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Unifying metadata
###############################################################################################################

Template$set("public","unifyingMetadata", function(){
  cat("Unifying metadatas.\n")
  summary <- self$summary
  metadata <- self$metadata
  
  for(study in names(metadata)){
    if(study =="mTOR"){
      summary[[study]][["Tube#"]] <- apply(summary[[study]], 1,function(x) strsplit(x["Raw.file"],"_")[[1]][2])
      metadata[[study]][["Name"]] <- summary[[study]][["Experiment"]][match( metadata[[study]][["Tube#"]], summary[[study]][["Tube#"]])]
      metadata[[study]][["Label"]] <- "Light"
      metadata[[study]][["Experiment_label"]] <- apply(metadata[[study]],1, function(x) paste0(c(x["Target"],x["Condition"],x["Buffer"]), collapse="_"))
      metadata[[study]][["Base"]] <- NA
      
      }
    
    if(study =="LysoIP1"){
      summary[[study]][["Tube#"]] <- seq(1:dim(summary[[study]])[1])
      metadata[[study]] <-  summary[[study]]%>% select("Tube#","Experiment")
      colnames(metadata[[study]]) <- c("Tube#","Name")
      metadata[[study]][["Label"]] <- ifelse(grepl("a", metadata[[study]][["Name"]]),"Heavy","Light")
      
      comment <- data.frame(num = c(1,2,3,4,5), treatment = c("Untreated","GA.treated","Torin1.treated","Torin1.PI.treated","PI.treated"))
      metadata[[study]][["Treatment"]] <- comment$treatment[match(substr(metadata[[study]][["Name"]], start = 1, stop =1), comment$num)]
      metadata[[study]][["Experiment_label"]] <- metadata[[study]][["Treatment"]]
      metadata[[study]][["Base"]] <- NA
    }
    
    if(study =="LysoIP2"){
      summary[[study]][["Tube#"]] <- apply(summary[[study]], 1,function(x) strsplit(x["Raw.file"],"_")[[1]][3])
      metadata[[study]][["Tube#"]] <- apply(metadata[[study]], 1,function(x) strsplit(x["Name"],"_")[[1]][3])
      metadata[[study]][["Label"]] <- metadata[[study]]$`HA Label`
      metadata[[study]][["Experiment_label"]] <- apply(metadata[[study]],1, function(x) ifelse(is.na(x["Time"])| x["Time"]=="--"|x["Time"]==0,paste0(c(x["Cell Line"],x["Treatment"]), collapse="_"), paste0(c(x["Cell Line"],x["Treatment"],x["Time"]), collapse="_")))
      metadata[[study]][["Experiment_label"]] <- apply(metadata[[study]],1, function(x) ifelse(grepl("DSP|Mix",x["Experiment"] ),paste0(c(x["Experiment_label"],x["Experiment"]), collapse="_"),x["Experiment_label"]))
      metadata[[study]][["Base"]] <- NA

            }
    
    metadata[[study]][["Name"]] <- summary[[study]][["Experiment"]][match(metadata[[study]][["Tube#"]], summary[[study]][["Tube#"]])]
    if(study=="mTOR"){
      metadata[[study]][["Name"]] <- gsub("/",".",metadata[[study]][["Name"]])
    }
  }
  
  self$summary <- summary
  self$metadata <- metadata
  
  invisible(self)
})