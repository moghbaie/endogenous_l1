###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 11/14/2021
## Imputing intensities
##############################################################################################################

IDIRT_Template$set("public","enrichmentAnalysis", function(){
  
  experiments <- self$experiments[!grepl("mix",self$experiments)]
  
  egos_tbl <- foreach(col = experiments,.combine= "rbind")%do%{
    
    uniprotIDs <- self$hits%>%filter_at(vars(col),~ (.x==1))%>%
      .$uniprotID
    
    proteins <- mapIds(org.Hs.eg.db, uniprotIDs, 'ENSEMBL', 'UNIPROT')
    gene.df <- bitr(proteins , fromType = "ENSEMBL",
                    toType = c("ENTREZID","ENSEMBL", "SYMBOL","UNIPROT"),
                    OrgDb = org.Hs.eg.db)
    
    ego_bp <- performEnreachment(gene.df, cat='BP')
    if(dim(ego_bp)[1]>0){
      ego_bp[["category"]] <- "BP"
    }
    
    ego_mf <- performEnreachment(gene.df, cat='MF')
    if(dim(ego_mf)[1]>0){
      ego_mf[["category"]] <- "MF"
    }
    
    ego_cc <- performEnreachment(gene.df, cat='CC')
    if(dim(ego_cc)[1]>0){
      ego_cc[["category"]] <- "CC"    
    }
    
    ego_tbl <- rbind(ego_bp, ego_mf , ego_cc)
    ego_tbl[["cluster"]] <- col
    
    return(ego_tbl)
  }
  
  self[["enrichment"]] <- egos_tbl
  
  invisible(self)
})