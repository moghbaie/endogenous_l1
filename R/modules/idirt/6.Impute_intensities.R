###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Imputing intensities
##############################################################################################################

Template$set("public","compareIntensity2Mix", function(){
  
  metadata <- self$metadata
  proteinGroups <- self$proteinGroups
  mixdata <- self$mixdata
  robustness <- self$robustness
  for(study in names(metadata)[!names(metadata)%in% c("mTOR","LysoIP2")]){
    
    print(study)
    experiments <- unique(metadata[[study]]$Experiment_label[!is.na(metadata[[study]]$Base)])
    
    foreach(j = 1:length(experiments))%do%{
      exp <- experiments[j]
      case <- metadata[[study]]$Name[metadata[[study]]$Experiment_label==exp]
      control_label <- unique(metadata[[study]]$Base[metadata[[study]]$Experiment_label==exp])
      control <- metadata[[study]]$Name[metadata[[study]]$Experiment_label==control_label]
      
      intensities <- proteinGroups[[study]][,c("Protein.IDs","Fasta.headers", paste0("iBAQ.",case))]
      intensities[["Gene.name"]] <- apply(intensities, 1,function(x) ifelse(!is.na(strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]),
                                                                            paste0(c(strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3],strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]), collapse ="-"),
                                                                            strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3]))
      intensities[["uniprotID"]] <- apply(intensities,1, function(x) strsplit(x["Fasta.headers"], "\\|")[[1]][2])
      
      intensities.mix <- mixdata[[study]][,c("Protein.IDs","Fasta.headers", paste0("iBAQ.",control))]
      intensities.mix[["Gene.name"]] <- apply(intensities.mix, 1,function(x) ifelse(!is.na(strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]),
                                                                                    paste0(c(strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3],strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]), collapse ="-"),
                                                                                    strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3]))
      intensities.mix[["uniprotID"]] <- apply(intensities.mix,1, function(x) strsplit(x["Fasta.headers"], "\\|")[[1]][2])
      
      compare_intensity <- merge(intensities[,!colnames(intensities) %in% c("Protein.IDs","Fasta.headers")], intensities.mix[,!colnames(intensities.mix) %in% c("Protein.IDs","Fasta.headers")], by = c("uniprotID","Gene.name"), all= TRUE)
      
      compare_intensity[,!colnames(compare_intensity) %in% c("Gene.name", "uniprotID")] <- apply(compare_intensity[,!colnames(compare_intensity) %in% c("Gene.name","uniprotID")], 2, function(x) ifelse(x==0, NA,log(x)))
      
      index <- apply(compare_intensity[,!colnames(compare_intensity) %in% c("Gene.name", "uniprotID")],1, function(x) sum(x, na.rm=T))>0
      compare_intensity <- compare_intensity[index,]
      
      compare_intensity[,paste0("iBAQ.",case)]  <- impute(compare_intensity[,paste0("iBAQ.",case)], amm = "2", pmm = "6")
      compare_intensity[,paste0("iBAQ.",control)]  <- impute(compare_intensity[,paste0("iBAQ.",control)], amm = "2", pmm = "6")
      
      ## Just added
      #LysoIP$avg_values[[study]]$Gene.name <- ifelse(!is.na(LysoIP$avg_values[[study]]$Gene.name.x), LysoIP$avg_values[[study]]$Gene.name.x,LysoIP$avg_values[[study]]$Gene.name.y)
      
      compare_intensity[["avg.fold.change"]] <- apply(compare_intensity[,paste0("iBAQ.",case)] ,1,mean) - apply(compare_intensity[,paste0("iBAQ.",control)] ,1,mean)
      compare_intensity[[paste0("Avg.Percentage.H.L.",exp)]] <- self$avg_values[[study]][[paste0("Avg.Percentage.H.L.",exp)]][match(compare_intensity$uniprotID, self$avg_values[[study]]$uniprotID)]
      compare_intensity[[paste0("Avg.Percentage.H.L.",exp,"_mix")]] <- self$avg_values[[study]][[paste0("Avg.Percentage.H.L.",exp,"_mix")]][match(compare_intensity$uniprotID, self$avg_values[[study]]$uniprotID)]
      
      
      compare_intensity[["cluster"]] <- robustness[[study]][[exp]]$cluster[match(compare_intensity$uniprotID,robustness[[study]][[exp]]$names.groups. )]
      compare_intensity[["num.miss"]] <- robustness[[study]][[exp]]$num.miss[match(compare_intensity$uniprotID, robustness[[study]][[exp]]$names.groups. )]
      
      self$compare_intensity[[study]][[exp]] <- compare_intensity
      
    }
  }
  
})