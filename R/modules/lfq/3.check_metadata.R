# 01/09/2021
# Mehrnoosh Oghbaie
# Create metadata folder
# download metadata.txt and comparison.txt from dropbox

LFQ_Template$set("public","importMetadata", function(txt_folder, token){

  txt_folder <- self$txt.dir
  
  metadata_folder <- basename(txt_folder)
  print(metadata_folder)
  
  if(!dir.exists(file.path("metadata",metadata_folder))){
    
    cat("Create metadata folder.\n")
    
    dir.create(file.path("metadata",metadata_folder))
  }
  
  
  mq_files <-  tryCatch(
    { 
      drop_dir(txt_folder, dtoken = token, recursive = TRUE)
      
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", txt_folder))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      
    })
  
  
  metadata_files <- mq_files$path_display[grepl("metadata.txt|comparison.txt",mq_files$name)]
  
  print("Downloading metadata.txt and comparison.txt.")
  lapply(metadata_files, function(x) drop_download(x, local_path = file.path("metadata",metadata_folder), overwrite = TRUE))
  
  self$metadata <- read.delim(file.path("metadata",metadata_folder, "metadata.txt"),sep="")
  
  self$comparison <- read.delim(file.path("metadata",metadata_folder, "comparison.txt"),sep="")
  
  invisible(self)
  
})


LFQ_Template$set("public","separateComparisons", function(intensity="iBAQ"){
  
  cat("Separating comparisons intensities.\n")
  cat("Calculate abundance. \n")
  
  comparison <- self$comparison
  metadata <- self$metadata
  proteinGroups <- self$proteinGroups
  experiments <- list()
  
  for(i in 1: dim(comparison)[1]){

    case <- paste0(paste0(intensity,"."),metadata$Experiment[metadata$Experiment_label %in% comparison$Case[i]])
    control <- paste0(paste0(intensity,"."), metadata$Experiment[metadata$Experiment_label %in% comparison$Control[i]])
    
    exp <- proteinGroups[, colnames(proteinGroups) %in%  c("uniprotID", "Gene.names", case, control)]
    
  
    identification.case <- paste0("Identification.type.",metadata$Experiment[metadata$Experiment_label %in% comparison$Case[i]])
    exp[["Identification.by.MS.MS.case"]] <- apply(proteinGroups[,colnames(proteinGroups) %in% identification.case],1,function(x) sum(ifelse(x=="By MS/MS",1,0)))

    identification.control <- paste0("Identification.type.", metadata$Experiment[metadata$Experiment_label %in% comparison$Control[i]])
    exp[["Identification.by.MS.MS.control"]] <- apply(proteinGroups[,colnames(proteinGroups) %in% identification.control],1,function(x) sum(ifelse(x=="By MS/MS",1,0)))
   
    
    exp[, colnames(exp) %in%  c( case, control)] <- apply(exp[, colnames(exp) %in%  c( case, control)],2,function(x) ifelse(x==0,0,log2(x)))
    
    exp <- exp[rownames(exp)[apply(exp[, colnames(exp) %in%  c( case, control)],1,sum)>0],]
    
    
    comparison_label <- comparison$Comparison[i]
    experiments[[comparison_label]] <- exp
    
  }
  
  
  
  self$experiments <- experiments
  
  invisible(self)
})

