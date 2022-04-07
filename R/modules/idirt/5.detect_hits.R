###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## detect hits
###############################################################################################################

## detect the hits by finding the outliers separately (normal distribution outlier / IQR outlier / finding local minimum within specific range)
## detect the hits comparatively by passing the mix (only for experiment that has mix distribution)
## detect the hits by making networks and clustering

IDIRT_Template$set("public","calculatingAvgValues", function(intensity = "iBAQ"){
  cat("calculating average intensity.\n")
  cat("calculating average percentage.\n")
  
  metadata <- self$metadata
  percentages <- self$percentages
  proteinGroups <- self$proteinGroups

  experiments <- unique(metadata[["Experiment_label"]])
  self$experiments <- experiments

  df_avg <- self$percentages[, c( "Gene.names", "uniprotID")]
    
  for(exp in experiments){
      
      names <- metadata$Expriment[metadata[["Experiment_label"]] == exp]
      dp <- percentages[,colnames(percentages) %in% sapply( names, function(x) paste0("Percentage.H.L.",x))]
      di <- proteinGroups[,colnames(proteinGroups) %in% sapply( names, function(x) paste0(c(intensity,x), collapse = "."))]
      df_avg[[paste0("Avg.Percentage.H.L.",exp)]] <- apply(dp,1,function(x) ifelse(mean(x, na.rm= T)==0|is.na(mean(x, na.rm= T)), 0, mean(x[x!=0&!is.na(x)])))
      df_avg[[paste0("Avg.",intensity,".",exp)]] <- apply(di,1,function(x) mean(log2(1+x)))
      
    }
    
    self$avg_values <- df_avg
  
  invisible(self)
})

IDIRT_Template$set("public","detectHits", function(){
  
  metadata <- self$metadata
  percentages <- self$percentages
  df_avg <- self$avg_values
  
  stats <- df_avg%>%
    dplyr::select(contains("mix")&contains("Percentage"))%>%
    filter(.>0)%>%
    summarize(mean = mean(.[[1]]),
              sd = sd(.[[1]]),
              median = median(.[[1]]),
              Q1=quantile(.[[1]], probs = 0.25),
              Q3=quantile(.[[1]], probs = 0.75),
              IQR = Q3-Q1,
              outlier = median+1.5*IQR,
              outlierNormal = mean+3*sd)

  df_significant <- df_avg%>%
    select(uniprotID, Gene.names,!contains("mix")&contains("Percentage"))
  df_significant[,-c(1,2)] <- apply(df_significant[,-c(1,2)],2,function(x) as.numeric(x>stats$outlier))
  
  colnames(df_significant) <- gsub("Avg.Percentage.H.L.", "Significant.",colnames(df_significant))
  df_avg <- cbind(df_avg, df_significant[,-c(1,2)])
  
  self$avg_values <- df_avg
  self$hits <- cbind(df_avg[,c("uniprotID","Gene.names")],df_significant[,-c(1,2)])
  colnames(self$hits) <- gsub("Significant.","",colnames(self$hits))
  
  invisible(self)
})

