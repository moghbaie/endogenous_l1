## Mehrnoosh Oghbaie
## 03/21/2022
## Differential Expression analysis ...

################################################################################
### RNA -Seq template
################################################################################



TemplateRna <- R6Class("TemplateRna",
                       list(
                         count.dir = NA,
                         count = NA,
                         metadata = NA,
                         comparison = NA,
                         edgeList = list(),
                         edgeListNormalized = list(),
                         dispersion = list(),
                         glmFit = list(),
                         glmLRT = list(),
                         diffExpress = list(),
                         toptable = list(),
                         toptags = list()
))


## Importing input
source(file.path("R/modules/rna-seq/1.import_files.R"), local = TRUE)$value

## Making edgeList
source(file.path("R/modules/rna-seq/2.make_edgelist.R"), local = TRUE)$value

## Calculate differential expression
source(file.path("R/modules/rna-seq/3.diff_expression.R"), local = TRUE)$value
