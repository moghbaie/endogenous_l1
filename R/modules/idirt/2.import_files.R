###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Import input
###############################################################################################################

Template$set("public","importInput", function(){
  
  dirs <- list.dirs(file.path(getwd(),self$input.dir), 
                    full.names = FALSE, 
                    recursive = FALSE)
  
  cat("There are", length(dirs), "folders in input folder.\n")
  cat("Importing proteinGroups.txt.\n")
  cat("Removing the contaminants and reverse sequences.\n")
  cat("Importing summary.txt .\n")
  
  pg <- lapply(dirs, function(x) 
    tryCatch({read.delim(file.path(getwd(),self$input.dir, x,"txt/proteinGroups.txt"))%>%
        filter(Potential.contaminant != "+"& Reverse != "+"& Q.value < 0.05 )},
        error = function(e) {cat(paste("There is no proteinGroups.txt in", x, "folder"))}
    ))
  names(pg) <- dirs
  
  sm <- lapply(dirs, function(x) 
    tryCatch({read.delim(file.path(getwd(),self$input.dir, x,"txt/summary.txt"))%>%
        filter(Experiment!="")},
             error = function(e) {cat(paste("There is no summary.txt in", x, "folder"))}
    ))
  names(sm) <- dirs
  
  self$proteinGroups <- pg
  self$summary <- sm
  
  invisible(self)
})

## Importing metadata
Template$set("public","importMetadata", function(){

  dirs <- list.dirs(file.path(getwd(),self$input.dir), 
                    full.names = FALSE, 
                    recursive = FALSE)
  
  cat("Importing metadata files.\n")
  
  metadata <- lapply(dirs, function(x) if(length(list.files(file.path(getwd(),self$input.dir, x,"metadata")))>0)
    read_excel(list.files(file.path(getwd(),self$input.dir, x,"metadata"),full.names = TRUE)))

  names(metadata) <- dirs
  
  self$metadata <- metadata
  
  invisible(self)
})


## Importing mix data
Template$set("public","importMix", function(){
  
  dirs <- list.dirs(file.path(getwd(),self$input.dir), 
                    full.names = FALSE, 
                    recursive = FALSE)
  
  cat("Importing mix files.\n")
  
  mixdata <- lapply(dirs, function(x) if(length(list.files(file.path(getwd(),self$input.dir, x,"txt_mix")))>0)
    read.delim(file.path(getwd(),self$input.dir, x,"txt_mix","proteinGroups.txt"))%>%
      filter(Potential.contaminant != "+"& Reverse != "+"& Q.value < 0.05 ))
  
  names(mixdata) <- dirs
  
  self$mixdata <- mixdata
  
  invisible(self)

})

## Merge mix data
Template$set("public","mergeMix", function(){
  
  dirs <- list.dirs(file.path(getwd(),self$input.dir), 
                    full.names = FALSE, 
                    recursive = FALSE)
  
  cat("Merging mix files.\n")
  
  mixdata <- lapply(dirs, function(x) if(length(list.files(file.path(getwd(),self$input.dir, x,"txt_mix")))>0)
    read.delim(file.path(getwd(),self$input.dir, x,"txt_mix","proteinGroups.txt")))
  
  names(mixdata) <- dirs
  
  self$mixdata <- mixdata
  
  invisible(self)
  
})
