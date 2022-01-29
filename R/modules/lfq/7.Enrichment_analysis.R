# Mehrnoosh Oghbaie
# 01/26/2022
# Enrichment analysis with clusterprofiler

LFQ_Template$set("public","enrichmentAnalysis", function(){
  egos_tbl <- foreach(col = names(self$experiments),.combine= "rbind")%do%{
    print(col)
    uniprotIDs <- self$hits_tbl%>%filter_at(vars(col),~ (.x==1))%>%
      .$uniprotID
    
    proteins <- mapIds(org.Hs.eg.db, uniprotIDs, 'ENSEMBL', 'UNIPROT')
    gene.df <- bitr(proteins , fromType = "ENSEMBL",
                    toType = c("ENTREZID","ENSEMBL", "SYMBOL","UNIPROT"),
                    OrgDb = org.Hs.eg.db)
    
    
    ego_bp <- performEnreachment(gene.df, cat='BP')
    ego_bp[["category"]] <- "BP"
    
    ego_mf <- performEnreachment(gene.df, cat='MF')
    ego_mf[["category"]] <- "MF"
    
    ego_cc <- performEnreachment(gene.df, cat='CC')
    ego_cc[["category"]] <- "CC"
    
    ego_tbl <- rbind(ego_bp, ego_mf , ego_cc)
    ego_tbl[["cluster"]] <- col
    
    return(ego_tbl)
  }
  
  self[["enrichment"]] <- egos_tbl
  
invisible(self)
})