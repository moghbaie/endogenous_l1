###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 01/09/2021
## Preparing data from MaxQuant or ...
###############################################################################################################
## Defining the class


LFQ_Template <- R6::R6Class("LFQ_Template",
                        list(
                          txt.dir = NA, # MQ output txt folder location
                          qa_report = NA, # QC report summart table
                          parameters = NA, # MQ parameters
                          proteinGroups = NA, # proteinGroups.txt
                          summary = NA, # summary.txt
                          metadata = NA, # metadata file
                          comparison = NA, # comparison table
                          experiments = list(), # List of experiments with replicates intensity
                          outliers = list(), # List of outliers in each replicate
                          significant_tbl = list(),# List of significant tables including p.values and lfc
                          significantGroup_tbl = list(), # List of significant group table (including all buffers)
                          hits_tbl = NA, # List of hits in each comparison
                          enrichment = NA, 
                          anova_tbl = NA,
                          plot = list()
                        ))

## Run QA
source(file.path("R/modules/lfq/1.make_qa_report.R"), local = TRUE)$value

## Importing input
source(file.path("R/modules/lfq/2.import_files.R"), local = TRUE)$value

## Tmport metadata
source(file.path("R/modules/lfq/3.check_metadata.R"), local = TRUE)$value

## Impute missing abundance
source(file.path("R/modules/lfq/4.impute_abundance.R"), local = TRUE)$value

## Run ttest
source(file.path("R/modules/lfq/5.make_significant_table.R"), local = TRUE)$value

## Linear regression
source(file.path("R/modules/lfq/6.fit_linear_regression.R"), local = TRUE)$value

## Enrichment analysis
source(file.path("R/modules/lfq/7.Enrichment_analysis.R"), local = TRUE)$value

## Make plots
source(file.path("R/modules/lfq/8.make_plot.R"), local = TRUE)$value
