# Mehrnoosh Oghbaie
# 01/09/2021
# Calculate lfc
# Perform ttest
# Make the hit table

LFQ_Template$set("public","fitLinearRegression", function(intensity ="iBAQ", limma= TRUE){
  
  cat("Calculate log fold change.\n")
  metadata <- self$metadata
  comparison <- self$comparison
  
  for(comparison in names(self$experiments)){
    
    print(comparison)
    modeldata <- rbind(
      self$metadata%>% 
        filter(Experiment_label %in%  self$comparison$Case[self$comparison$Comparison == comparison])%>%
        mutate(tag = 1),
      self$metadata%>% 
        filter(Experiment_label %in% self$comparison$Control[self$comparison$Comparison == comparison])%>%
        mutate(tag = -1))
    
    data <- self$experiments[[comparison]]%>%
      dplyr::select(contains(intensity))
    rownames(data)<- self$experiments[[comparison]]$uniprotID
    
    
#    cat("Normalizing each comparison.\n")
#    
#    data.copy <- as.matrix(data)
#    data.copy[data.copy==0] <- NA
#    data.norm <- performGlobalRLRNormalization(data.copy, noLogTransform=TRUE)
    
    if(limma ){
      
      mod <- model.matrix(~as.factor(modeldata$tag))
      cols <- paste0("LFQ.intensity",".",modeldata$Experiment)
      dy <- data[,cols]
      fit.limma <- limma::lmFit(dy, mod)
      ebayes.limma <- eBayes(fit.limma)
      
      sig_tbl <- topTable(ebayes.limma,number = dim(dy)[1])%>%
        as.data.frame%>%
        select(logFC, adj.P.Val)%>%
        mutate(Significant = ifelse(logFC>1&adj.P.Val<0.05,"Y","N"))%>%
        mutate(uniprotID = rownames(.))
      
    }else{
      
      mod <- model.matrix(~as.factor(modeldata$tag))
      cols <- paste0("LFQ.intensity",".",modeldata$Experiment)
      dy <- data[,cols]
      dy <- dy[apply(dy,1,function(x) sum(x, na.rm=T))>0,]
      fit.lin <- lm.fit(mod,t(dy))
      logFC <- as.data.frame(t(fit.lin$coefficients))[[2]]
      P.Val <- genefilter::rowFtests(as.matrix(dy), as.factor(modeldata$tag))
      P.Val[["adj.P.Val"]] <- p.adjust(P.Val$p.value,method = "BH")
        
      sig_tbl <- data.frame(logFC=logFC, adj.P.Val= P.Val$adj.P.Val, row.names = rownames(dy))%>%
          mutate(Significant = ifelse(logFC>1&adj.P.Val<0.05,"Y","N"))%>%
          mutate(uniprotID = rownames(.))
    }
    
    self$significant_tbl[[comparison]] <- merge(self$significant_tbl[[comparison]],sig_tbl, by = "uniprotID")
  }
  
  invisible(self)
})



LFQ_Template$set("public","extractHits", function(){
  
  cat("Extract hits.\n")
  hits <- list()
  
  for(name in names(self$significant_tbl)){
    hits[[name]] <- self$significant_tbl[[name]]%>%
      filter(Significant=="Y" )%>%
      dplyr::select(uniprotID, Gene.names)
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

