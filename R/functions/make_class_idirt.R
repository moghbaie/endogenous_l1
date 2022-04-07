###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Preparing data from MaxQuant or ...
###############################################################################################################
## Defining the class


IDIRT_Template <- R6::R6Class("IDIRT_Template",
                    list(
                      txt.dir = NA, # MQ output txt folder location
                      qa_report = NA, # QC report summart table
                      parameters = NA, # MQ parameters
                      proteinGroups = NA, # proteinGroups.txt
                      summary = NA, # summary.txt
                      metadata = NA, # metadata file
                      percentages = NA,
                      experiments = list(),
                      avg_values = NA,
                      hits = NA,
                      t.test_result = NA,
                      enrichment = NA
                      ))


## Run QA
source(file.path("R/modules/idirt/1.make_qa_report.R"), local = TRUE)$value

## Importing input
source(file.path("R/modules/idirt/2.import_files.R"), local = TRUE)$value

## unifying and integrating metadata
source(file.path("R/modules/idirt/3.check_metadata.R"), local = TRUE)$value

## Calculate percentage
source(file.path("R/modules/idirt/4.calculate_percentage.R"), local = TRUE)$value

## Detect hits
source(file.path("R/modules/idirt/5.detect_hits.R"), local = TRUE)$value

## Enrichment analysis
source(file.path("R/modules/idirt/6.Enrichment_analysis.R"), local = TRUE)$value

## Make plot
source(file.path("R/modules/idirt/7.make_plot.R"), local = TRUE)$value



