###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 01/09/2021
## Preparing data from MaxQuant or ...
###############################################################################################################
## Defining the class


PEP_Template <- R6::R6Class("PEP_Template",
                            list(
                              txt.dir = NA, # MQ output txt folder location
                              parameters = NA, # MQ parameters
                              proteinGroups = NA, # proteinGroups.txt
                              peptides = NA, # peptides.txt
                              summary = NA, # summary.txt
                              metadata = NA, # metadata file
                              comparison = NA, # comparison table
                              experiments = list(), 
                              plot = list()
                            ))


## Importing input
source(file.path("R/modules/L1_loci/1.import_files.R"), local = TRUE)$value

## Tmport metadata
source(file.path("R/modules/L1_loci/2.check_metadata.R"), local = TRUE)$value

## Match peptides
source(file.path("R/modules/L1_loci/3.match_peptides.R"), local = TRUE)$value

## Visualize
source(file.path("R/modules/L1_loci/4.visualize.R"), local = TRUE)$value