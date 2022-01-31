# Mehrnoosh Oghbaie
# 01/09/2021
# Imputingmissing values


LFQ_Template$set("public","imputeMissingValues", function(intensity = "iBAQ"){
  self$significant_tbl <- self$experiments
  metadata <- self$metadata
  comparison <- self$comparison
  for(i in 1:length(self$experiments)){
    print(i)
    name <- names(self$experiments)[i]
    print(name)
    case <- paste0(paste0(intensity,"."), metadata$Experiment[metadata$Experiment_label %in% comparison$Case[comparison$Comparison==name]])
    control <- paste0(paste0(intensity,"."), metadata$Experiment[metadata$Experiment_label %in% comparison$Control[comparison$Comparison==name]])
    
    exp <- self$experiments[[name]]
    
    
    cat("Normalizing each comparison.\n")
    
    exp.copy <- as.matrix(exp%>% dplyr::select(contains(intensity)))
    exp.copy[exp.copy==0] <- NA
    
    exp.norm <- performGlobalRLRNormalization(exp.copy, noLogTransform=TRUE)
    exp.norm <- as.data.frame(exp.norm)

  self$experiments[[name]][,case] <- impute(exp.norm[,case], amm = "2", pmm = "6")
    
  self$experiments[[name]][,control] <- impute(exp.norm[,control], amm = "2", pmm = "6")
 
  }
  
  invisible(self)
})
