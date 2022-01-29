# Mehrnoosh Oghbaie
# 01/09/2021
# Calculate lfc
# Perform ttest
# Make the hit table

LFQ_Template$set("public","fitLinearRegression", function(intensity ="iBAQ"){
  
  cat("Calculate log fold change.\n")
  metadata <- self$metadata
  comparison <- self$comparison
  
  for(comparison in names(self$experiments)){
    
    modeldata <- rbind(
      self$metadata%>% 
        filter(Experiment_label %in%  self$comparison$Case[self$comparison$Comparison == comparison])%>%
        mutate(tag = 1),
      self$metadata%>% 
        filter(Experiment_label %in% self$comparison$Control[self$comparison$Comparison == comparison])%>%
        mutate(tag = -1))
    
    data <- self$experiments[[comparison]]%>%
      select(contains(intensity))
    rownames(data)<- self$experiments[[comparison]]$uniprotID
    
    data1 <- data %>%pivot_longer(cols = everything())%>% filter(value!=0)
    linear.mod <- lm(value~ as.factor(name)-1,data=data1)
    sum <- summary(linear.mod)$coefficients[,"Estimate"]
    names(sum) <- gsub("as.factor\\(name\\)","", names(sum))
    
    cat("Normalizing each comparison.\n")
    data.norm <- t(t(data)-sum)
    
    mod <- model.matrix(~as.factor(modeldata$tag))
    cols <- paste0("LFQ.intensity",".",modeldata$Experiment)
    dy <- data.norm[,cols]
    fit.limma <- limma::lmFit(dy, mod)
    ebayes.limma <- eBayes(fit.limma)
    
    sig_tbl <- topTable(ebayes.limma,number = dim(dy)[1])%>%
      as.data.frame%>%
      select(logFC, adj.P.Val)%>%
      mutate(Significant = ifelse(logFC>1&adj.P.Val<0.05,"Y","N"))%>%
      mutate(uniprotID = rownames(.))
    
    self$significant_tbl[[comparison]] <- merge(self$significant_tbl[[comparison]],sig_tbl, by = "uniprotID")
  }
  
  invisible(self)
})



LFQ_Template$set("public","extractHits", function(){
  
  cat("Extract hits.\n")
  hits <- list()
  
  for(name in names(self$significant_tbl)){
    hits[[name]] <- self$significant_tbl[[name]]%>%
      filter(Significant=="Y" &(Identification.by.MS.MS.case>0 |Identification.by.MS.MS.control>0))%>%
      select(uniprotID, Gene.names)
  }
  
  
  hits_tbl <- foreach(name = names(hits), .combine="rbind")%do%{
    hits[[name]]
  }
  
  hits_tbl <- unique(hits_tbl)
  
  for(name in names(hits)){
    hits_tbl[[name]] <- ifelse(hits_tbl$uniprotID %in% hits[[name]]$uniprotID,1,0)
  }
  
  self$hits_tbl <- hits_tbl
  
  invisible(self)
})

