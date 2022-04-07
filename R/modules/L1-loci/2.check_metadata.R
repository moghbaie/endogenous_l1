# 01/09/2021
# Mehrnoosh Oghbaie
# Create metadata folder
# download metadata.txt and comparison.txt from dropbox

PEP_Template$set("public","importMetadata", function(txt_folder, token){
  
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