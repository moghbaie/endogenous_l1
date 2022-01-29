# Mehrnoosh Oghbaie
# 01/09/2021
# Imputingmissing values


LFQ_Template$set("public","imputeMissingValues", function(intensity = "iBAQ"){
  self$significant_tbl <- self$experiments
  metadata <- self$metadata
  comparison <- self$comparison
  for(i in 1:length(self$experiments)){
    
    name <- names(self$experiments)[i]
    print(name)
    case <- paste0(paste0(intensity,"."), metadata$Experiment[metadata$Experiment_label %in% comparison$Case[comparison$Comparison==name]])
    control <- paste0(paste0(intensity,"."), metadata$Experiment[metadata$Experiment_label %in% comparison$Control[comparison$Comparison==name]])
    
    exp <- self$experiments[[name]]
  

  self$experiments[[name]][,case] <- impute(exp[,case], amm = "2", pmm = "6")
    
  self$experiments[[name]][,control] <- impute(exp[,control], amm = "2", pmm = "6")
 
  }
  
  invisible(self)
})
