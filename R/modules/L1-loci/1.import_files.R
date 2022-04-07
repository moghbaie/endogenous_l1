# 03/29/2021
# Mehrnoosh Oghbaie
# Import proteinGroups.txt and peptides.txt from dropbox folder

PEP_Template$set("public","importMQOutput", function(txt_folder, token, temp_folder){
  
  cat("Importing proteinGroups.txt.\n")
  cat("Importing peptides.txt.\n")
  cat("Removing the contaminants and reverse sequences.\n")
  cat("Importing summary.txt .\n")
  cat("Importing parameters.txt .\n")
  
  self$txt.dir <- txt_folder
  
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
  
  # Attaching proteinGroups.txt to self
  proteinGroups_file <- txt_files[grepl("proteinGroups",txt_files)]
  if(!file.exists(file.path(temp_folder,"proteinGroups.txt"))){
    drop_download(proteinGroups_file, local_path = temp_folder)
  }
  
  proteinGroups <-  read.delim(file.path(temp_folder,"proteinGroups.txt"))
  proteinGroups <- proteinGroups%>% 
    filter(Potential.contaminant != "+"& Reverse != "+"& Q.value < 0.05 )
  
  proteinGroups[["Gene.names"]] <- apply(proteinGroups, 1,function(x) ifelse(!is.na(strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]),
                                                                             paste0(c(strsplit(x["Fasta.headers"],"_|;|\\|")[[1]][3],strsplit(strsplit(x["Fasta.headers"],";|\\|")[[1]][2],"-")[[1]][2]), collapse ="-"),
                                                                             strsplit(strsplit(x["Fasta.headers"],"GN=")[[1]][2]," ")[[1]][1]))
  
  proteinGroups[["uniprotID"]] <- apply(proteinGroups,1, function(x) strsplit(x[["Protein.IDs"]],";")[[1]][1])
  self$proteinGroups <- proteinGroups
  
  
  # Attach peptides.txt file to self
  peptides_file <- txt_files[grepl("peptides",txt_files)]
  if(!file.exists(file.path(temp_folder,"peptides.txt"))){
    drop_download(peptides_file, local_path = temp_folder)
  }
  
  peptides <- read.delim(file.path(temp_folder,"peptides.txt"))
  peptides <- peptides%>%
    filter(Potential.contaminant != "+"& Reverse != "+")
  
  
  self$peptides <-  peptides
  
  
  # Attach summary.txt file to self
  summary_file <- txt_files[grepl("summary",txt_files)]
  if(!file.exists(file.path(temp_folder,"summary.txt"))){
    drop_download(summary_file, local_path = temp_folder)
  }
  self$summary <-  read.delim(file.path(temp_folder,"summary.txt"))
  
  
  # Attach parameters.txt file to self
  parameter_file <- txt_files[grepl("parameters.txt",txt_files)]
  if(!file.exists(file.path(temp_folder,"parameters.txt"))){
    drop_download(parameter_file, local_path = temp_folder)
  }
  self$parameters <-  read.delim(file.path(temp_folder,"parameters.txt"))
  
  invisible(self)
})