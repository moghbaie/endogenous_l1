set.seed(123)
# Importing the generic functions
source(file.path("R/functions/functions.R"), local = TRUE)
#renv::init(bare = TRUE)
#renv::restore()
# Importing libraries
CRAN.packages <- c("readr","readxl","reshape2","dplyr","magrittr","rdrop2","tidyr","tidyverse","RColorBrewer","plotly","gplots",
                   "ggplot2","R6","ggrepel", "writexl", "foreach", "ggridges","ggpubr", "renv", "cowplot","gridExtra", "matrixStats")
bioconductor.packages <- c("PTXQC","org.Hs.eg.db","clusterProfiler","ReactomePA","limma", "preprocessCore", "NormalyzerDE","genefilter","edgeR")

install.packages.if.necessary(CRAN.packages,bioconductor.packages)
# Connecting to dropbox folder
token <- drop_auth(new_user  = TRUE)

# Importing LFQ class generator files
source(file.path("R/functions/make_class_lfq.R"), local = TRUE)
# Importing IDIRT class generator files
source(file.path("R/functions/make_class_idirt.R"), local = TRUE)
# Importing RNA-Seq class generator files
source(file.path("R/functions/make_class_rna.R"), local = TRUE)
# Importing peptides analyzer class generator files
source(file.path("R/functions/make_class_peptide.R"), local = TRUE)

# create a temp folder
temp_folder <- tempdir()
#renv::snapshot()

##################################################################################################################
## LFQ pipeline

txt_folder_CRC_v2 <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Proteomics/2. MQ output/LFQ_L1_endogenous_CRC/CRC_ver2/"

# Running the project
LFQ_CRC_v2 <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_CRC_v2, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_CRC_v2, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_CRC_v2, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity", limma= FALSE)$
  extractHits()$
  fitGroupLinearRegression(intensity = "LFQ.intensity", impute.par = TRUE, limma=FALSE)$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Tissue","Antibody"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Tissue", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Tissue","Antibody"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Tissue","Antibody"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= TRUE)$
  scatterplot()$
  scatterplotGroup()$
  plotHeatmap(intensity = "LFQ.intensity",genes.interest= c("L1RE1"), n= 50, multiwell = TRUE)$
  plotMdsMultiwell(intensity = "LFQ.intensity")

save(LFQ_CRC_v2, file= "output/backup/CRC_v2.RData")
#base::load("output/backup/CRC_v2.RData") 

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

##################################################################################################################
txt_folder_CRC_v1 <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Proteomics/2. MQ output/LFQ_L1_endogenous_CRC/CRC_ver1/"

# Running the project
LFQ_CRC_v1 <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_CRC_v1, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_CRC_v1, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_CRC_v1, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity", limma=FALSE)$
  extractHits()$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Tissue","Antibody"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Tissue", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Tissue","Antibody"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Tissue","Antibody"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= FALSE)$
  scatterplot()$
  plotHeatmap(intensity = "LFQ.intensity",genes.interest="L1RE1", n= 50, multiwell = FALSE)

save(LFQ_CRC_v1, file= "output/backup/CRC_v1.RData")
#base::load("output/backup/CRC_v1.RData") 

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

#################################################################################################################
txt_folder_RIP_Seq <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Proteomics/2. MQ output/LFQ_L1_endogenous_PA1_NTERA_N2102EP/LFQ_RIP_Seq/"

LFQ_RIP_Seq <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_RIP_Seq, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_RIP_Seq, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_RIP_Seq, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity",limma=FALSE)$
  extractHits()$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Cell_type","Antibody"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Cell_type", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Cell_type","Antibody"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Cell_type","Antibody"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= FALSE)$
  scatterplot()$
  plotHeatmap(intensity = "LFQ.intensity",genes.interest="L1RE1", n= 50, multiwell = FALSE)


save(LFQ_RIP_Seq, file= "output/backup/LFQ_RIP_Seq.RData")
#base::load( "output/backup/LFQ_RIP_Seq.RData")
do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

#############################################################################################################################
txt_folder_24well_drug  <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Proteomics/2. MQ output/LFQ_L1_endogenous_24well_drug/24well_drug"

# Running the project
LFQ_24well_drug <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_24well_drug, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_24well_drug, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_24well_drug, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity", limma= FALSE)$
  extractHits()$
  fitGroupLinearRegression(intensity = "LFQ.intensity", impute.par = TRUE, limma=FALSE)$
  performAnova(intensity = "LFQ.intensity", impute.par = TRUE)$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Cell_type","Drug"))$
  pcaplot(intensity = "LFQ.intensity", colour = "Cell_type", shape = "Drug")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Cell_type","Drug"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Cell_type","Drug"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= TRUE)$
  scatterplot()$
  scatterplotGroup()$
  plotHeatmap(intensity = "LFQ.intensity", genes.interest="L1RE1", n= 50, multiwell = TRUE)$
  plotAnova(intensity = "LFQ.intensity", impute.par = TRUE)$
  plotMdsMultiwell(intensity = "LFQ.intensity")

save(LFQ_24well_drug, file= "output/backup/24well_drug.RData")
#base::load("output/backup/24well_drug.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

############################################################################################################################
txt_folder_N2102EP_24well  <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Proteomics/2. MQ output/LFQ_L1_endogenous_N2102EP_24well/N210EP_24well"

# Running the project
LFQ_N210EP_24well <- LFQ_Template$new()$
  runQAReport(txt_folder= txt_folder_N2102EP_24well, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_N2102EP_24well, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_N2102EP_24well, token)$
  separateComparisons(intensity="LFQ.intensity")$
  imputeMissingValues(intensity = "LFQ.intensity")$
  fitLinearRegression(intensity = "LFQ.intensity", limma= FALSE)$
  extractHits()$
  fitGroupLinearRegression(intensity = "LFQ.intensity", impute.par = TRUE, limma= FALSE)$
  enrichmentAnalysis()$
  qaplot(facet.vars= c("Antibody","BufferNum"))$
  pcaplot(intensity = "LFQ.intensity", colour = "BufferNum", shape = "Antibody")$
  boxplot(intensity = "LFQ.intensity",facet.vars= c("Antibody","BufferNum"))$
  imputationplot(intensity = "LFQ.intensity" , facet.vars= c("Antibody","BufferNum"))$
  volcanoplot( genes.interest=c("L1RE1"), multiwell= TRUE)$
  scatterplot()$
  scatterplotGroup()$
  plotHeatmap(intensity = "LFQ.intensity", genes.interest="L1RE1", n= 50)$
  plotMdsMultiwell(intensity = "LFQ.intensity")

save(LFQ_N210EP_24well, file= "output/backup/LFQ_N210EP_24well.RData")
#base::load("output/backup/LFQ_N210EP_24well.RData")

do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

######################################################################################################################
## IDIRT pipeline

txt_folder_IDIRT_ORF2 <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Proteomics/2. MQ output/IDIRT_L1_ectopic/IDIRT_L1"
# Running the project
IDIRT_L1_multiwell <- IDIRT_Template$new()$
  runQAReport(txt_folder= txt_folder_IDIRT_ORF2, token, temp_folder)$
  importMQOutput(txt_folder= txt_folder_IDIRT_ORF2, token, temp_folder)$
  importMetadata(txt_folder =txt_folder_IDIRT_ORF2, token)$
  calculatePercentage()$
  calculatingAvgValues(intensity = "iBAQ")$
  detectHits()$
  CreateScatterPlot(intensity = "iBAQ")$
  CreateDotPlot()$
  CreateTilePlot(top_n = 80)$
  MDSPlot(intensity = "iBAQ")


save(IDIRT_L1_multiwell, file= "output/backup/IDIRT_L1_multiwell.RData")
#base::load("output/backup/IDIRT_L1_multiwell.RData")
do.call(file.remove, list(list.files(temp_folder, full.names = TRUE)))

######################################################################################################################
## Peptide analysis L1 loci



######################################################################################################################
## RNA-Seq pipeline

count_folder <- "CLOUD_store_MS_RAW/L1_endogenouse_ProtoGenomic/Genomics/RNA-Seq/counts"

# Running the project

RIP_Seq <- TemplateRna$new()$
  importCount(count_folder, token, temp_folder)$
  importMetadata(count_folder, token, temp_folder)$
  makeEdgeList(minCPM = 100, minLibrary=2)$
  makeCountBarPlot()$
  makelogCPMBoxPlot()$
  makePCAPlot()$
  diffExpression()$
  plotBCV()$
  plotSmear()$
  plotHeatmap()$
  drawVolcanoPlot()
 
save(RIP_Seq, file= "output/backup/RIP_Seq.RData")
#base::load("output/backup/RIP_Seq.RData")
do.call(file.remove, list(list.files(temp_folder, full.names = TRUE))) 
