###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Preparing data from MaxQuant or ...
###############################################################################################################
## Defining the class


Template <- R6::R6Class("Template",
                    list(
                      input.dir = "input",
                      proteinGroups = list(),
                      mixdata = list(),
                      summary = list(),
                      metadata = list(),
                      percentages = list(),
                      mixpercentages = list(),
                      avg_values = list(),
                      normaltest = NA,
                      outliers = list(),
                      percentageSignificant = list(),
                      robustness = list(),
                      compare_intensity = list(),
                      logfoldchange2Base = list(),
                      hits = list(),
                      t.test_result = list(),
                      enrichment = list()
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

## impute intensities
source(file.path("R/modules/idirt/6.Impute_intensities.R"), local = TRUE)$value





