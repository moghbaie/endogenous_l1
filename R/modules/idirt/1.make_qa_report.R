# 01/09/2021
# Mehrnoosh Oghbaie
# Download MQ output txt folder in a temp file
# Run PTXQC
# Upload quality report on the original dropbox folder
# Attach report summary table to the object


IDIRT_Template$set("public","runQAReport", function(txt_folder, token, temp_folder){
  
  self$txt.dir <- txt_folder
  
  cat("Checking txt folder.\n")
  
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
  
  txt_files <- mq_files$path_display[grepl(".txt",mq_files$name)]
  
  if(length(txt_files[grepl("Report|report", txt_files)])==0){
    cat("Downloading MQ output files to temp folder.\n")
    lapply(txt_files, function(x) drop_download(x, local_path = temp_folder))
    
    cat("Run quality control report.\n")
    createReport(temp_folder)
    
    
    all_files <- list.files(temp_folder, full.names= TRUE)
    report_files <- all_files[grepl("Report|report", all_files)]
    
    cat("Upload QC reports.\n")
    lapply(report_files, function(x) drop_upload(x, path = txt_folder, dtoken = token))
    
    self$qa_report <-  read.delim(file.path(temp_folder,basename(report_files[grepl("heatmap.txt", report_files)])))
    
  }else{
    drop_download(txt_files[grepl("heatmap.txt", txt_files)], local_path = temp_folder, overwrite= TRUE)
    self$qa_report <-  read.delim(file.path(temp_folder,basename(txt_files[grepl("heatmap.txt", txt_files)])))
  }
  
  invisible(self)
})


