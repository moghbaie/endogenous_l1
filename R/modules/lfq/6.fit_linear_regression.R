# Mehrnoosh Oghbaie
# 01/23/2021
# Group linear regression
# Anova analysis including different buffers


LFQ_Template$set("public","fitGroupLinearRegression", function(intensity = "iBAQ", impute.par = TRUE, limma=FALSE){
   
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
            mod <- model.matrix(~factor(modeldata$tag[modeldata$BufferNum==buffer],levels=c(-1,1)))
            cols <- paste0(intensity,".",modeldata$Experiment[modeldata$BufferNum==buffer])
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
            mod <- model.matrix(~factor(modeldata$tag[modeldata$BufferNum==buffer],levels=c(-1,1)))
            cols <- paste0(intensity,".",modeldata$Experiment[modeldata$BufferNum==buffer])
            dy <- data.norm[,cols]
            dy <- dy[apply(dy,1,function(x) sum(x, na.rm=T))>0,]
            
            case.avg <- apply(dy[,paste0(paste0(intensity,"."),modeldata$Experiment[modeldata$BufferNum==buffer&modeldata$tag==1])],1, function(x) mean(x, na.rm=T))
            control.avg <- apply(dy[, paste0(paste0(intensity,"."),modeldata$Experiment[modeldata$BufferNum==buffer&modeldata$tag==-1])],1, function(x) mean(x, na.rm=T))
            
            uniprotID <- rownames(dy)
            fit.lin <- lm.fit(mod,t(dy))
            logFC <- as.data.frame(t(fit.lin$coefficients))[[2]]
            P.Val <- genefilter::rowFtests(dy, as.factor(modeldata$tag[modeldata$BufferNum==buffer]))
            P.Val[["adj.P.Val"]] <- p.adjust(P.Val$p.value,method = "BH")
            
            dt <- data.frame(case.avg = case.avg, control.avg = control.avg, logFC=logFC, adj.P.Val= P.Val$adj.P.Val, row.names = rownames(dy))%>%
               mutate(Significant = ifelse(logFC>1&adj.P.Val<0.05,"Y","N"))%>%
               mutate(uniprotID = rownames(.), Buffer = paste0("BufferNum.", buffer))
         }
      }
      significant_tbl[["Gene.names"]] <- self$proteinGroups$Gene.names[match(significant_tbl$uniprotID, self$proteinGroups$uniprotID)]
      self$significantGroup_tbl[[comparison]] <- significant_tbl
      
   }
   
   invisible(self)
})



LFQ_Template$set("public","performAnova", function(intensity = "iBAQ", impute.par = TRUE, anova.vars= c("Drug", "BufferNum")){
   
   self$comparison[["Multi.Comparison"]] <- apply(self$comparison,1, function(x) paste0(c(gsub("_Buffer.+$|_buffer.+","",x["Case"]), gsub("_Buffer.+$|_buffer.+","",x["Control"])), collapse=".vs."))
   
   comparison <- self$comparison
   
   case <- unique(comparison$Case[!comparison$Case %in% comparison$Control])
   control <- unique(comparison$Control[comparison$Control %in% comparison$Case])
   
   modeldata <- rbind(self$metadata%>% 
                         select(Experiment, Experiment_label, any_of(anova.vars))%>%
                         filter(Experiment_label %in% case)%>%
                         mutate(tag=1),
                      self$metadata%>% 
                         select(Experiment, Experiment_label, any_of(anova.vars))%>%
                         filter(Experiment_label %in% control)%>%
                         mutate(tag=-1, Drug = "drug12"),
                      self$metadata%>% 
                         select(Experiment, Experiment_label, any_of(anova.vars))%>%
                         filter(Experiment_label %in% control)%>%
                         mutate(tag=-1, Drug = "drug473"))%>%
      mutate(opt = paste0(tag,Drug))
   
   
   data <- self$proteinGroups[, sapply(modeldata$Experiment, function(x) paste0(c(intensity, x) ,collapse="."))]
   
   rownames(data) <- self$proteinGroups$uniprotID
   
   data <- data[apply(data,1,sum)>0,]
   
   data.copy <- as.matrix(data)
   data.copy[data.copy==0] <- NA
   data.norm <- performGlobalRLRNormalization(data.copy)
   
   if(impute.par){
      cat("Impute missing values.\n")
      for(buffer in unique(modeldata$BufferNum)){
         
         print(buffer)
         dat.flag <- data.norm[, sapply(modeldata$Experiment[modeldata$BufferNum==buffer ], function(x) paste0(c(intensity,x), collapse="."))]
         non_zero_records <- rownames(dat.flag)[apply(dat.flag,1,function(x) sum(x, na.rm=T))>0]
         data.backup <-data.norm[non_zero_records, ]
         
         for(opt in unique(modeldata$opt[modeldata$BufferNum==buffer])){
            print(opt)
            dat <- data.backup[, sapply(modeldata$Experiment[modeldata$BufferNum==buffer & modeldata$opt==opt], function(x) paste0(c(intensity,x), collapse="."))]
            dat.imputed <- impute(dat, amm = "2", pmm = "6")
            data.norm[non_zero_records, sapply(modeldata$Experiment[modeldata$BufferNum==buffer & modeldata$opt==opt], function(x) paste0(c(intensity,x), collapse="."))] <- dat.imputed
         }
      }
   }
   
   anova_tbl <- foreach(buffer= unique(modeldata$BufferNum), .combine= rbind)%do%{
      
      print(buffer)
      modeldata_b <- modeldata[modeldata$BufferNum==buffer,]
      mod <- model.matrix(~factor(modeldata_b$tag, levels=c(-1,1)) + as.factor(modeldata_b$Drug))
      
      cols <- paste0(intensity,".",modeldata_b$Experiment)
      dy <- data.norm[,cols]
      dy <- dy[apply(dy,1,function(x) sum(x, na.rm=T))>0,]
      
      case.lfc <- apply(dy[,paste0(paste0(intensity,"."),modeldata_b$Experiment[modeldata_b$Drug=="drug473"&modeldata_b$tag==1])],1, function(x) mean(x, na.rm=T))-
         apply(dy[,paste0(paste0(intensity,"."),modeldata_b$Experiment[modeldata_b$Drug=="drug473"&modeldata_b$tag==-1])],1, function(x) mean(x, na.rm=T))
      
      control.lfc <- apply(dy[, paste0(paste0(intensity,"."),modeldata_b$Experiment[modeldata_b$Drug=="drug12"&modeldata_b$tag==1])],1, function(x) mean(x, na.rm=T))-
         apply(dy[, paste0(paste0(intensity,"."),modeldata_b$Experiment[modeldata_b$Drug=="drug12"&modeldata_b$tag==-1])],1, function(x) mean(x, na.rm=T))
      
      fit.lin <- lm.fit(mod,t(dy))
      logFC <- as.data.frame(t(fit.lin$coefficients))[[2]]
      P.Val <- genefilter::rowFtests(dy, as.factor(modeldata$tag[modeldata$BufferNum==buffer]))
      P.Val[["adj.P.Val"]] <- p.adjust(P.Val$p.value,method = "BH")
      
      dt <- data.frame(case.lfc = case.lfc, control.lfc = control.lfc, logFC=logFC, adj.P.Val= P.Val$adj.P.Val, row.names = rownames(dy))%>%
         mutate(Significant = ifelse(logFC>1 &adj.P.Val<0.05,"Y","N"))%>%
         mutate(uniprotID = rownames(.), Buffer = paste0("BufferNum.", buffer))
      
      
   }
   
   
   
   self$anova_tbl <- anova_tbl
   
   invisible(self)
})
