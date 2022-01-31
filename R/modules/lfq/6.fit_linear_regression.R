# Mehrnoosh Oghbaie
# 01/23/2021
# Anova analysis including different buffers


LFQ_Template$set("public","fitGroupLinearRegression", function(intensity = "iBAQ", impute.par = TRUE, limma=TRUE){
   
   self$comparison[["Multi.Comparison"]] <- apply(self$comparison,1, function(x) paste0(c(gsub("_Buffer.+$|_buffer.+","",x["Case"]), gsub("_Buffer.+$|_buffer.+","",x["Control"])), collapse=".vs."))
   
   for(comparison in unique(self$comparison[["Multi.Comparison"]])){
      
      modeldata <- rbind(
         self$metadata%>% 
            filter(Experiment_label %in%  self$comparison$Case[self$comparison$Multi.Comparison == comparison])%>%
            arrange(BufferNum)%>%
            dplyr::select(Experiment, BufferNum)%>%
            mutate(tag = 1),
         self$metadata%>% 
            filter(Experiment_label %in% self$comparison$Control[self$comparison$Multi.Comparison == comparison])%>%
            arrange(BufferNum)%>%
            dplyr::select(Experiment, BufferNum)%>%
            mutate(tag = -1))
      
      case <- modeldata$Experiment[modeldata$tag==1]
      control <- modeldata$Experiment[modeldata$tag==-1]
      
      
      data <- cbind(
         self$proteinGroups[, sapply(case, function(x) paste0(c(intensity, x) ,collapse="."))],
         self$proteinGroups[, sapply(control, function(x) paste0(c(intensity, x) ,collapse="."))]
      )
      
      rownames(data) <- self$proteinGroups$uniprotID
      
      cat("Calculate log2 fold change.\n")
      data <- data[apply(data,1,sum)>0,]
      
      
      data.copy <- as.matrix(data)
      data.copy[data.copy==0] <- NA
      data.norm <- performGlobalRLRNormalization(data.copy)
      
      
      if(impute.par){
         cat("Impute missing values.\n")
         for(buffer in unique(modeldata$BufferNum)){
            
            dat.flag <- data.norm[, sapply(modeldata$Experiment[modeldata$BufferNum==buffer ], function(x) paste0(c(intensity,x), collapse="."))]
            non_zero_records <- rownames(dat.flag)[apply(dat.flag,1,function(x) sum(x, na.rm=T))>0]
            data.backup <-data.norm[non_zero_records, ]
            
            for(tag in unique(modeldata$tag[modeldata$BufferNum==buffer])){
               
               dat <- data.backup[, sapply(modeldata$Experiment[modeldata$BufferNum==buffer & modeldata$tag==tag], function(x) paste0(c(intensity,x), collapse="."))]
               dat.imputed <- impute(dat, amm = "2", pmm = "6")
               data.norm[non_zero_records, sapply(modeldata$Experiment[modeldata$BufferNum==buffer & modeldata$tag==tag], function(x) paste0(c(intensity,x), collapse="."))] <- dat.imputed
            }
         }
      }
      
      if(limma ){
         significant_tbl <- foreach(buffer= unique(modeldata$BufferNum), .combine= rbind)%do%{
            
            print(buffer)
            mod <- model.matrix(~as.factor(modeldata$tag[modeldata$BufferNum==buffer]))
            cols <- paste0("LFQ.intensity",".",modeldata$Experiment[modeldata$BufferNum==buffer])
            dy <- data.norm[,cols]
            fit.limma <- limma::lmFit(dy, mod)
            ebayes.limma <- eBayes(fit.limma)
            
            dt <- topTable(ebayes.limma,number = dim(dy)[1])%>%
               as.data.frame%>%
               select(logFC, adj.P.Val)%>%
               mutate(Significant = ifelse(logFC>1&adj.P.Val<0.05,"Y","N"))%>%
               #rename_with(~paste0(.,".BufferNum.",buffer), everything())%>%
               mutate(uniprotID = rownames(.), Buffer = paste0("BufferNum.", buffer))
         }
      }else{
         
         significant_tbl <- foreach(buffer= unique(modeldata$BufferNum), .combine= rbind)%do%{
            
            print(buffer)
            mod <- model.matrix(~as.factor(modeldata$tag[modeldata$BufferNum==buffer]))
            cols <- paste0("LFQ.intensity",".",modeldata$Experiment[modeldata$BufferNum==buffer])
            dy <- data.norm[,cols]
            dy <- dy[apply(dy,1,function(x) sum(x, na.rm=T))>0,]
            uniprotID <- rownames(dy)
            fit.lin <- lm.fit(mod,t(dy))
            logFC <- as.data.frame(t(fit.lin$coefficients))[[2]]
            P.Val <- genefilter::rowFtests(dy, as.factor(modeldata$tag[modeldata$BufferNum==buffer]))
            P.Val[["adj.P.Val"]] <- p.adjust(P.Val$p.value,method = "BH")
            
            dt <- data.frame(logFC=logFC, adj.P.Val= P.Val$adj.P.Val, row.names = rownames(dy))%>%
               mutate(Significant = ifelse(logFC>1&adj.P.Val<0.05,"Y","N"))%>%
               mutate(uniprotID = rownames(.), Buffer = paste0("BufferNum.", buffer))
         }
      }
      significant_tbl[["Gene.names"]] <- self$proteinGroups$Gene.names[match(significant_tbl$uniprotID, self$proteinGroups$uniprotID)]
      self$significantGroup_tbl[[comparison]] <- significant_tbl
      
   }
   
   invisible(self)
})
