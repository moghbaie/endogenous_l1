## Mehrnoosh Oghbaie
## 03/29/2022
###############################################################################
### Make a filtered list EdgeR


obtainFilteredDGE <- function(countMat, targetInfo, genes, minLibrary) {
  # Making sure that a DGEList is made.
  dge <- DGEList(counts=countMat, samples=targetInfo, genes = genes)
  # Determining the expressed genes.
  isExpr <- rowSums(cpm(dge)>=1) >= minLibrary
  # These expressed genes will be saved within the dge dataset.
  dge <- dge[isExpr, ]
  # Calculate the normalization factors.
  dge <- calcNormFactors(dge)
  # The dge list is returned.
  return (dge)
}

TemplateRna$set("public","makeEdgeList", function(minCPM = 100, minLibrary=2){
  print("Making the Edge list and normalize it.")
  metadata <- self$metadata
  count <- self$count
  comparison <- self$comparison
  
  for(name in unique(comparison$Comparison)){
    
    conditions <- comparison%>%
      filter(Comparison %in% name)%>%
      select(Case,Control)%>%
      unname%>%
      unlist
          

    split <- count[,colnames(count) %in% c("Ensembl.ID","Gene.Name",metadata$File_name[metadata$Experient_label %in% conditions])]
    rownames(split) <- make.names(split$Gene.Name, unique=TRUE)

    split[,!colnames(split) %in% c("Ensembl.ID","Gene.Name")] <- apply(split[,!colnames(split) %in% c("Ensembl.ID","Gene.Name")],2, as.numeric)
    split <- split[!is.na(split$Gene.Name),]
    
    dgList <- obtainFilteredDGE( split[,!colnames(split) %in% c("Ensembl.ID","Gene.Name")], 
                                 metadata[metadata$Experient_label %in% conditions,],
                                 rownames(split),
                                 minLibrary)
    
    ###Filtering
    keep <- rowSums(cpm(dgList)>minCPM) >= minLibrary
    dgList <- dgList[keep,]
    
    self$edgeList[[name]] <- dgList
    
    logCPMcounts <- cpm(dgList,log=TRUE)
    
    self$edgeListNormalized[[name]] <- logCPMcounts
    
  }
  
  invisible(self)
})

TemplateRna$set("public","makeCountBarPlot", function(){
  
  print("Making library size and sum count bar plot.")
  comparison <- self$comparison
  metadata <- self$metadata
  
  for(name in unique(comparison$Comparison)){
    
    dgList <- self$edgeList[[name]]
    
    df <- as.data.frame(dgList$samples)
    
    p1 <- ggplot(df, aes(x = as.factor(File_name), y = lib.size)) + # basic of ggplot, define x and y
      geom_bar(stat="identity", fill = "darksalmon") + # add a barplot
      ggtitle("Barplot of library sizes") + # add title
      theme_classic() + # select color theme
      ylab("Library size") + # y axis label
      xlab("Sample") + # x axis label
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) # rotate x labels
    
    
    # barplot of number of expressed genes per sample
    nGene <- data.frame(matrix(NA, nrow = length(colnames(split[,!colnames(split) %in% c("Ensembl.ID","Gene.Name")])), ncol = 1)) # maak leeg dataframe
    rownames(nGene) <- metadata$File_name[metadata$Experient_label %in% conditions] # rownames word sample names
    colnames(nGene) <- "genes"
    for (j in 1:ncol(split[,!colnames(split) %in% c("Ensembl.ID","Gene.Name")])){ # per sample
      nGene[j,1] <- length(which(split[,j] > 0)) # tellen hoeveel genen een expressie hebben van > 0
    }
    
    p2 <- ggplot(nGene, (aes(x = rownames(nGene), y = genes))) +
      geom_bar(stat = "identity", fill = "green4") +
      xlab("sample") +
      ylab("Number of genes expressed") +
      ggtitle("Number of genes with count > 0") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    
    p <-  grid.arrange(arrangeGrob(p1,
                                   p2,
                                   ncol=2, widths=c(1,1)))
    
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir)))))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir)))))}
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"barplot"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"barplot"))}
    ggsave(file=file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"barplot",paste0(name,"_bar.png")), 
           p,  
           width=9,
           height=8,
           dpi=200)
  }
  
  invisible(self)
})

TemplateRna$set("public","makelogCPMBoxPlot", function(){
  print("Drawing logCPM boxplot.")
  comparison <- self$comparison
  
  for(name in unique(comparison$Comparison)){
    logCPMcounts <- self$edgeListNormalized[[name]] 
    df <- melt(logCPMcounts)
    colnames(df) <- c("Gene.name","condition","logCPM")
    df$condition <- as.character(df$condition)
    q <- ggplot(df, aes(x=condition, y=logCPM)) + 
      geom_boxplot()+
      geom_hline(yintercept=median(logCPMcounts),col="blue")+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90))+
      ggtitle("Boxplots of logCPMs")
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"boxplot"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"boxplot"))}
    ggsave(file=file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"boxplot",paste0(name,"_boxplot.png")), 
           q,  
           width=5,
           height=7,
           dpi=200)
  }
    
  
  invisible(self)
})

TemplateRna$set("public","makePCAPlot", function(){
  print("Drawing PCA plot.")
  comparison <- self$comparison
  
  for(name in unique(comparison$Comparison)){
    
    dgList <- self$edgeList[[name]]
    logCPMcounts <- self$edgeListNormalized[[name]] 

    project.pca <- prcomp(t(logCPMcounts)) 
    summary(project.pca)
    
    # plot
    dg <- project.pca$x
    dg <- dg[,1:2]
    dg <- as.data.frame(dg)
    dg <- cbind(dg, dgList$samples)
    
    # PCA plot colored by library size shows that the separation of the samples is not caused by the library size
    # So that is good because that indicates that the samples segregate based on biological variation instead of technical variation
    # so i hope these 3 groups represent your experimental groups :) 
    qq <- ggplot(dg, aes(PC1, PC2, colour = lib.size)) + # you can change lib.size by any column in your target file
      geom_point(size = 3) + 
      scale_shape_manual(values =  c(21, 23)) +
      geom_text_repel(aes(label = File_name), nudge_x = 5) +
      xlab(paste("PC1 (", (round((project.pca$sdev[1]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
      ylab(paste("PC2 (", (round((project.pca$sdev[2]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
      theme_classic() +
      ggtitle("PCA logCPM")
    
    if(!dir.exists(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"pcaplot"))){dir.create(file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"pcaplot"))}
    ggsave(file=file.path(paste0("output/plot/", basename(dirname(self$count.dir))),"pcaplot",paste0(name,"_pcaplot.png")), 
           qq,  
           width=7,
           height=7,
           dpi=200)
  
  }
  
  invisible(self)
})

