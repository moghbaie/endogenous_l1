## Mehrnoosh Oghbaie
## 03/29/2022
###############################################################################
### Calculate differential expression

TemplateRna$set("public","diffExpression", function(){
  print("Calculate differential expression")
  metadata <- self$metadata
  comparison <- self$comparison
  for(name in comparison$Comparison){
    dgList <- self$edgeList[[name]] 
    logCPMcounts <- self$edgeListNormalized[[name]]
    
    cat(paste("Creating design model for",name))
    design <- model.matrix(~ 0 + dgList$samples$Experient_label, data=dgList$samples) # basically making a binary matrix that contains the information which sample is in which experimental group
    
    cat(paste("Estimating dispersion model for",name,"\n"))
    dge.full <- estimateDisp(dgList, design)
    self$dispersion[[name]] <- dge.full
    
    cat(paste("Fitting glm moel on dispersion model for",name,"\n"))
    fit <- glmFit(dge.full, design) # fitting of the glm
    self$glmFit[[name]] <- fit
    
    cat(paste("Setting the model contrast for",name,"\n")) 
    contrast <- ifelse(grepl('IgG|RNA',colnames(design)),-1,1)
    
    cat(paste("Perform ANOVA like test between case and control using glmLRT for",name,"\n"))
    lrt <- glmLRT(fit, contrast = contrast)
    self$glmLRT[[name]] <- lrt
    
    tags <- topTags(lrt, n=dim(dge.full[[1]])[1], adjust.method="BH", sort.by="logFC") #this is the full table for all genes
    self$toptags[[name]] <- tags
    
    cat(paste("Extract significant genes (logFC > 1 and p < 0.05) for",name,"\n"))
    toptable <- tags[[1]][which(abs(tags[[1]]$logFC) > 0.5 & tags[[1]]$FDR < 0.05),] # this is the table with only the significant genes
    self$toptable[[name]] <- toptable
      # print info
    print(summary(decideTests(lrt, p = 0.05, lfc = 1, adjust.method = "BH")))
   
    toptable <- self$toptable[[name]]
   }
  
  invisible(self)
})

TemplateRna$set("public","plotBCV", function(){
  cat("Plot Biological Coefficient of Variation \n")
  comparison <- self$comparison
  for(name in comparison$Comparison){
    
    dge.full <- self$dispersion[[name]]
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"bcvplot"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"bcvplot"))}
    
    png( file= file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"bcvplot",paste0(name,"_bcvplot.png")),  width=550, height=550)
    plotBCV(dge.full) 
    dev.off()
    }
   
  invisible(self)
})

TemplateRna$set("public","plotSmear", function(){
  
  comparison <- self$comparison
  
  for(name in comparison$Comparison){
    
    lrt <- self$glmLRT[[name]]
    toptable <- self$toptable[[name]]
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"smearplot"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"smearplot"))}
    png( file= file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"smearplot",paste0(name,"_smearplot.png")),  width=600, height=600)
   
    cat(paste("Plots log-Fold Change versus log-Concentration for Count Data for",name,"\n"))
    
    plotSmear(lrt, de.tags = toptable$genes)
    abline(h=c(-1, 1), col=2)
    title(paste("DEGS for ", name, sep = ""))

    dev.off()
  }
  
  invisible(self)
})

TemplateRna$set("public","plotHeatmap", function(){
  
  comparison <- self$comparison
  
  for(name in comparison$Comparison){
    
    logCPMcounts <- self$edgeListNormalized[[name]]
    toptable <- self$toptable[[name]]
    # plot DE genes in a heatmap
    n_allgenes <- logCPMcounts[rownames( toptable), ] # here you can also opt to extract the values from the countspermillion matrix, but logCPM is usually nicer
    n_allgenes <- n_allgenes[order(-abs(toptable$logFC)),]
    n_allgenes <- n_allgenes[!grepl("NA.",rownames(n_allgenes)),]
    
    colfunc <- colorRampPalette(c("steelblue", "slategray2", "white", "tan1", "firebrick3")) # these are nice colors
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"heatmap"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"heatmap"))}
    png( file= file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"heatmap",paste0(name,"_heatmap.png")),  width=500, height=1000)
    
    gplots::heatmap.2(n_allgenes, 
              col=colfunc(50), # define number of sections in color scale
              keysize = 1,
              scale="row", # does a row z-score
              trace="none", # is really ugly
              cexRow = 0.7, # font size of the row labels
              cexCol = 0.7, # font size of the col labels
              #Colv = "NA", # do unsupervised, meaning you don't do column clustering
              dendrogram = 'both', # or choose 'row' or 'col or 'none
              key.title = NA,
              main = "Differentially expressed genes")
    
    dev.off()
    
  }
  
  invisible(self)
})

TemplateRna$set("public","drawVolcanoPlot", function(){
  
  comparison <- self$comparison
  
  for(name in comparison$Comparison){
    toptable <- self$toptable[[name]]
    cat(paste("Draw volcano plot for",name,"\n"))

    qp <- ggplot( toptable, aes(x = logFC, y = -log(FDR), label = genes)) +
      geom_point(colour = "grey80", size = 0.5) + 
      geom_hline(yintercept = -log(0.05), colour = "grey80") + # these are the sifnigicane thresholds (p < 0.05)
      geom_vline(xintercept = c(-1,1), colour = "grey80") + # logFC thresholds (logFC > 1 or < -1)
      geom_point(data = subset(toptable,!grepl("L1_", toptable$genes)), aes(x = logFC, y = -log(FDR)), 
                 shape = 21, 
                 colour = "grey50", 
                 fill = "grey", 
                 alpha=0.7) +
      geom_point(data = subset(toptable,grepl("L1_", toptable$genes)), aes(x = logFC, y = -log(FDR)), 
                 shape = 21, 
                 colour = "red2",
                 fill = "red2")+
      geom_text_repel(data = subset(toptable,grepl("L1_", toptable$genes)), color="red2") + # set treshold from when you want to show gene names, or put a # for this line if you dont want to show gene names
      theme_classic() +
      xlim(-max(abs(toptable$logFC)),max(abs(toptable$logFC))) +
      ggtitle(name) 
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"volcanoplot"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"volcanoplot"))}
    ggsave(file=file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"volcanoplot",paste0(name,"_volcanoplot.png")), 
           qp,  
           width=7,
           height=7,
           dpi=200)
    
  }
  
  invisible(self)
})


