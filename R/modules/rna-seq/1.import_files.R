# 01/09/2021
# Mehrnoosh Oghbaie
# Import proteinGroups.txt and summary.txt from dropbox folder

TemplateRna$set("public","importCount", function(count_folder, token, temp_folder){


  self$count.dir <- count_folder
  
  rna_files <-  tryCatch(
    { 
      drop_dir(count_folder, dtoken = token, recursive = TRUE)
      
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", count_folder))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      
    })
  
  
  txt_files <- rna_files$path_display[grepl(".txt",rna_files$name)]
  
  # Attaching count.txt to self
  count_file <- txt_files[grepl("count",basename(txt_files))]
  if(!file.exists(file.path(temp_folder,count_file))){
    drop_download(count_file, local_path = temp_folder,overwrite= TRUE)
  }
  
  count <-  read.delim(file.path(temp_folder,basename(count_file)))
  count <- count%>%rename(Ens.ID. = Ens.ID..L1.locus.,
                 Gene.Name=Gene.Name..L1.type.)
  self$count <- count

  
  invisible(self)
})


TemplateRna$set("public","importMetadata", function(count_folder, token, temp_folder){
  
  
  self$count.dir <- count_folder
  
  rna_files <-  tryCatch(
    { 
      drop_dir(count_folder, dtoken = token, recursive = TRUE)
      
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", count_folder))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      
    })
  
  
  txt_files <- rna_files$path_display[grepl(".txt",rna_files$name)]
  
  metadata_folder <- basename(dirname(count_folder))
  print(metadata_folder)
  
  if(!dir.exists(file.path("metadata", metadata_folder))){
    
    cat("Create metadata folder.\n")
    
    dir.create(file.path("metadata",metadata_folder))
  }
  
  metadata_files <- rna_files$path_display[grepl("metadata.txt|comparison.txt",rna_files$name)]
  
  print("Downloading metadata.txt and comparison.txt.")
  lapply(metadata_files, function(x) drop_download(x, local_path = file.path("metadata",metadata_folder), overwrite = TRUE))
  
  # Attaching metadata.txt to self
  metadata_file <- txt_files[grepl("metadata",txt_files)]
  if(!file.exists(file.path(temp_folder,metadata_file))){
    drop_download(metadata_file, local_path = temp_folder,overwrite = TRUE)
  }
  
  metadata <-  read.delim(file.path(temp_folder,basename(metadata_file)),sep="")
  self$metadata <- metadata
  
  comparison_file <- txt_files[grepl("comparison",txt_files)]
  if(!file.exists(file.path(temp_folder,comparison_file))){
    drop_download(comparison_file, local_path = temp_folder, overwrite = TRUE)
  }
  
  comparison <-  read.delim(file.path(temp_folder,basename(comparison_file)),sep="")
  self$comparison <- comparison
  
  invisible(self)
})
