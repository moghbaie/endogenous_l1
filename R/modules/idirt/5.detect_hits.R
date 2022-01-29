###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## detect hits
###############################################################################################################

## detect the hits by finding the outliers separately (normal distribution outlier / IQR outlier / finding local minimum within specific range)
## detect the hits comparatively by passing the mix (only for experiment that has mix distribution)
## detect the hits by making networks and clustering

Template$set("public","calculatingAvgValues", function(){
  cat("calculating average intensity.\n")
  cat("calculating average percentage.\n")
  
  metadata <- self$metadata
  percentages <- self$percentages
  proteinGroups <- self$proteinGroups
  mixdata <- self$mixdata
  
  for(study in names(metadata)){
    experiment <- unique(metadata[[study]][["Experiment_label"]])
    
    if(study=="LysoIP2"){
      self$percentages[[study]][["uniprotID"]] <- apply(self$percentages[[study]], 1, function(x) strsplit(x[["Fasta.headers"]],"\\|")[[1]][2])
    }
    
    df_avg <- self$percentages[[study]][, c( "Gene.name", "uniprotID")]
    
    for(exp in experiment){
      
      names <- metadata[[study]][["Name"]][metadata[[study]][["Experiment_label"]] == exp]
      labels <- metadata[[study]][["Label"]][metadata[[study]][["Experiment_label"]] == exp]
      dp <- percentages[[study]][,colnames(percentages[[study]]) %in% sapply( names, function(x) paste0("Percentage.H.L.",x))]
      if(grepl("mix|Mix",exp)&study =="LysoIP1"){
        di <- mixdata[[study]][, colnames(mixdata[[study]]) %in% sapply( names, function(x) paste0("iBAQ.",x))]
      } else{
        di <- proteinGroups[[study]][,colnames(proteinGroups[[study]]) %in% sapply( names, function(x) paste0("iBAQ.",x))]
        
      }
      df_avg[[paste0("Avg.Percentage.H.L.",exp)]] <- apply(dp,1,function(x) ifelse(mean(x, na.rm= T)==0|is.na(mean(x, na.rm= T)), 0, mean(x[x!=0&!is.na(x)])))
      df_avg[[paste0("Avg.iBAQ.",exp)]] <- apply(di,1,function(x)ifelse(mean(x)==0,0,log2(mean(x[x!=0]))))[match(df_avg$uniprotID, self$percentages[[study]]$uniprotID)]
    }
    
    self$avg_values[[study]] <- df_avg
  }
  invisible(self)
})


Template$set("public","performNormalityTest", function(){
  cat("performing normality test on average percentage.\n")
  
  normaltest <- data.frame(matrix(NA, ncol=4, nrow=0))
  colnames(normaltest) <- c("study", "experiment","p.value","result")
  metadata <- self$metadata
  
  for(study in names(metadata)){
   
    df_avg <- self$avg_values[[study]] 
    
    for(col in colnames(df_avg)[grepl("Avg.Percentage.H.L.",colnames(df_avg))]){
     
      dx <- df_avg[[col]]
      dx <- dx[dx!=0]
      experiment <- col
      p.value <- shapiro.test(dx)$p.value
      result <- ifelse(shapiro.test(dx)$p.value > 0.05, "Yes","No")
      normaltest <- rbind(normaltest, cbind(study, experiment, p.value,result))
    }
  }
  self$normaltest <- normaltest
  
  invisible(self)
})


Template$set("public","detectOutliers", function(){
  cat("Detecting outliers.\n")
  
  metadata <- self$metadata
  percentages <- self$percentages
  for(study in names(metadata)[names(metadata)!="mTOR"]){
    experiments <- unique(metadata[[study]]$Experiment_label)
    self$outliers[[study]] <- list()
    # experiments that have both light and heavy tag
    
    experiments <- experiments[!grepl("mix|DSP|Mix", experiments)]
    foreach(j=1:length(experiments)) %do%{
      exp <- experiments[j]
      xa <- metadata[[study]]$Name[metadata[[study]]$Experiment_label == exp & metadata[[study]]$Label=="Heavy"]
      xb <-  metadata[[study]]$Name[metadata[[study]]$Experiment_label == exp & metadata[[study]]$Label=="Light"]
      
      dpvalue <- foreach(i=1:dim(percentages[[study]])[1], .combine=rbind)%do%{
        t.test(percentages[[study]][i, paste0("Percentage.H.L.",xa)], percentages[[study]][i, paste0("Percentage.H.L.",xb)])$p.value
      }
      
      outliers <- percentages[[study]]$Gene.name[dpvalue<0.01]
      outliers <- outliers[!is.na(outliers)]
      self$outliers[[study]][[exp]] <- outliers
      
    }
  } 
  invisible(self)
})


Template$set("public","PercentSignificantMatrix", function(){
  
  metadata <- self$metadata
  percentages <- self$percentages
  outliers <- self$outliers
  
  percentageSignificant <- list()
  
  for(study in names(metadata)[names(metadata)!="mTOR"]){
    
    print(study)
    experiments <- unique(metadata[[study]]$Experiment_label)
    experiments <- experiments[!grepl("mix|DSP|Mix", experiments)]
    
    foreach(j = 1:length(experiments))%do%{
      exp <- experiments[j]
      print(exp)
      names <- metadata[[study]]$Name[metadata[[study]]$Experiment_label==exp]
      dp <- percentages[[study]][,paste0("Percentage.H.L.",names)]
      dxp <- percentages[[study]][apply(dp,1,sum)>0 , c("Gene.name",paste0("Percentage.H.L.",names))]
      
      dep <- matrix(nrow= dim(dxp)[1], ncol= dim(dxp)[1])
      colnames(dep) <- rownames(dep) <- dxp$Gene.name
      
      foreach(i = 1: dim(dxp)[1])%do%{
        y <- dxp[i,-1]
        dx <- dxp[-i,-1]
        
        dep[-i,i] <- tryCatch({
          apply(dx,1, function(x,y) ks.test(unname(x)[!is.na(unname(x))], unname(y)[!is.na(unname(y))])$p.value,y = y)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        
      }
      percentageSignificant[[study]][[exp]] <- dep
    }
    
  }
  self$percentageSignificant <- percentageSignificant
  
  invisible(self)
})


Template$set("public","measureRobustness", function(){
  
  metadata <- self$metadata
  percentages <- self$percentages
  
  for(study in names(metadata)[names(metadata)!="mTOR"]){
    
    experiments <- unique(metadata[[study]]$Experiment_label)
    experiments <- experiments[!grepl("mix|DSP|Mix", experiments)]
    
    foreach(j = 1:length(experiments))%do%{
      exp <- experiments[j]
      
      names <- metadata[[study]]$Name[metadata[[study]]$Experiment_label==exp]
      dp <- percentages[[study]][,colnames(percentages[[study]]) %in% paste0("Percentage.H.L.",names)]
      dxp <- percentages[[study]][apply(dp,1,sum)>0& !is.na(apply(dp,1,sum)) ,colnames(percentages[[study]]) %in% c("Gene.name","uniprotID",paste0("Percentage.H.L.",names))]
      
      dk <- self$percentageSignificant[[study]][[exp]]
      diag(dk) <- 0
      
      groups <- cluster::pam(dk,2)$clustering
      groups <- data.frame(names(groups),groups)
      
      dpi <- self$avg_values[[study]][,c("Gene.name","uniprotID",paste0("Avg.Percentage.H.L.",exp), paste0("Avg.iBAQ.",exp))]
      groups[["iBAQ"]] <- dpi$Avg.iBAQ.Untreated[match(groups$names.groups., dpi$uniprotID)]
      groups[["Percentage"]] <- dpi$Avg.Percentage.H.L.Untreated[match(groups$names.groups., dpi$uniprotID)]
      
      dm <- data.frame(Gene.name = dxp$Gene.name,uniprotID = dxp$uniprotID, num.miss = apply(dxp[,!colnames(dxp)%in% c("uniprotID","Gene.name")], 1,function(x) sum(x==0)))
      groups[["num.miss"]] <- dm$num.miss[match(groups$names.groups., dm$uniprotID)]
      groups$Percentage[groups$Percentage==0] <- NA
      groups$iBAQ[groups$iBAQ==0] <- NA
      
      clusters <- groups%>% 
        group_by(groups)%>%
        summarise(robutness = mean(Percentage[Percentage!=0], na.rm=T)/mean(num.miss, na.rm=T))%>%
        arrange(desc(robutness))%>%
        mutate(cluster = paste0("group", 1:n()))
      
      groups[["cluster"]] <- clusters$cluster[match(groups$groups, clusters$groups)]
      
      self$robustness[[study]][[exp]] <- groups
      
    }
    
  }
  
  invisible(self)
})


