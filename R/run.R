
set.seed(123)
# Importing the generic functions
source(file.path("R/functions/functions.R"), local = TRUE)
#renv::init(bare = TRUE)
renv::restore()
# Importing libraries
CRAN.packages <- c("readr","readxl","reshape2","dplyr","magrittr","rdrop2","tidyr","tidyverse","RColorBrewer",
                   "ggplot2","R6","ggrepel", "writexl", "foreach", "preprocessCore", "ggridges","ggpubr", "renv")
bioconductor.packages <- c("PTXQC","org.Hs.eg.db","clusterProfiler","ReactomePA","limma")

install.packages.if.necessary(CRAN.packages,bioconductor.packages)
# Connecting to dropbox folder
token <- drop_auth()

# Importing LFQ class generator files
source(file.path("R/functions/make_class_lfq.R"), local = TRUE)
# Importing IDIRT class generator files
source(file.path("R/functions/make_class_idirt.R"), local = TRUE)
# Importing RNA-Seq class generator files
source(file.path("R/functions/make_class_rna.R"), local = TRUE)

# create a temp folder
temp_folder <- tempdir()
renv::snapshot()

##################################################################################################################
txt_folder_CRC_v2 <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/2. MQ output/LFQ_L1_endogenous_CRC/CRC_ver2/"

# Running the project
LFQ_CRC_v2 <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_CRC_v2, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_CRC_v2, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_CRC_v2, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity")$
  extractHits()$
  fitGroupLinearRegression(intensity = "LFQ.intensity", impute.par = TRUE)$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Tissue","Antibody"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Tissue", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Tissue","Antibody"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Tissue","Antibody"))$
  volcanoplot( genes.interest=c("LORF1"), multiwell= TRUE)$
  plotHeatmap(intensity = "LFQ.intensity", n= 50)

save(LFQ_CRC_v2, file= "output/backup/CRC_v2.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

##################################################################################################################
txt_folder_CRC_v1 <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/2. MQ output/LFQ_L1_endogenous_CRC/CRC_ver1/"

# Running the project
LFQ_CRC_v1 <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_CRC_v1, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_CRC_v1, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_CRC_v1, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity")$
  extractHits()$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Tissue","Antibody"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Tissue", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Tissue","Antibody"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Tissue","Antibody"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= FALSE)

save(LFQ_CRC_v1, file= "output/backup/CRC_v1.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

#################################################################################################################
txt_folder_RIP_Seq <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/2. MQ output/LFQ_L1_endogenous_PA1_NTERA_N2102EP/LFQ_RIP_Seq/"

LFQ_RIP_Seq <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_RIP_Seq, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_RIP_Seq, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_RIP_Seq, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity")$
  extractHits()$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Cell_type","Antibody"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Cell_type", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Cell_type","Antibody"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Cell_type","Antibody"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= FALSE)

save(LFQ_RIP_Seq, file= "output/backup/LFQ_RIP_Seq.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

#############################################################################################################################
txt_folder_24well_drug  <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/2. MQ output/LFQ_L1_endogenous_24well_drug/24well_drug"

# Running the project
LFQ_24well_drug <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_24well_drug, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_24well_drug, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_24well_drug, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity")$
  extractHits()$
  fitGroupLinearRegression(intensity = "LFQ.intensity", impute.par = TRUE)$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Cell_type","Drug"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Cell_type", shape = "Drug")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Cell_type","Drug"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Cell_type","Drug"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= TRUE)$
  plotHeatmap(intensity = "LFQ.intensity", n= 50)

save(LFQ_CRC_v2, file= "output/backup/txt_folder_24well_drug.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

############################################################################################################################
txt_folder_N2102EP_24well  <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/2. MQ output/LFQ_L1_endogenous_N2102EP_24well/24well_drug/N210EP_24well"

# Running the project
LFQ_N210EP_24well <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_N2102EP_24well, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_N2102EP_24well, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_N2102EP_24well, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity")$
  extractHits()$
  fitGroupLinearRegression(intensity = "LFQ.intensity", impute.par = TRUE)$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Cell_type","Drug"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Cell_type", shape = "Drug")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Cell_type","Drug"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Cell_type","Drug"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= TRUE)$
  plotHeatmap(intensity = "LFQ.intensity", n= 50)

save(LFQ_CRC_v2, file= "output/backup/txt_folder_N2102EP_24well.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))