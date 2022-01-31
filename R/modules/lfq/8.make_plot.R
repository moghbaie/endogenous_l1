# Mehrnoosh Oghbaie
# 01/25/2021
# Make multiple plots
#1. QC barplot
#2, Pca plot
#3. Boxplot from original log2(intensities)
#4. Plot before and after imputation
#5. Make Volcano plot
#6. Draw heatmap across all buffers in a condition
#7. Draw protein enrichment plot 


#1. QC barplot
LFQ_Template$set("public","qaplot", function(facet.vars= c("Tissue","Antibody") ){
   
   if(!dir.exists(paste0("output/plot/", basename(self$txt.dir)))){dir.create(paste0("output/plot/", basename(self$txt.dir)))}
   
   dx <- self$qa_report
   dx[["Raw.file"]]<- self$summary$Experiment[1:dim(dx)[1]]
   
   #1. QC barplot
   g1 <- dx%>%
      dplyr::select(Raw.file,contains("Contaminants"), contains("Peptide.Intensity"))%>%
      pivot_longer(!`Raw.file`)%>%
      mutate(name = ifelse(grepl("Contaminants",name),"Contaminants","Peptide.Intensity"))%>%
      merge(self$metadata, by.x = "Raw.file",by.y = "Experiment")%>%
      ggplot(aes(y=`Raw.file`,x=value, fill=name))+
      geom_col()+
      facet_grid(formula(paste0(paste0(facet.vars,collapse="+") ,"~name")), space="free", scale="free")+
      theme_bw()+
      theme(axis.text.y = element_text(size=7))
   
   
   if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa"))}
   ggsave(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa", "qaplot.png"),plot=g1,  width=8, height=(2+dim(dx)[1]/10), limitsize = F)
   
   invisible(self)
})

#2, Pca plot
LFQ_Template$set("public","pcaplot", function(intensity = "iBAQ", colour = "Tissue" ,shape = "Antibody" ){
   
   metadata <- self$metadata
   di <- self$proteinGroups %>% 
      dplyr::select(contains(intensity))
   
   di <- di[apply(di,1,sum)>0,]
   di <- log2(di+1)
   df <- data.frame(t(di))
   
   
   pca_res <- prcomp(df, scale. = TRUE)
   
   dx <- data.frame(pca_res$x)%>%
      dplyr::select(PC1,PC2)
   dx[["label"]] <- gsub(paste0(intensity,"."),"",rownames(df))
   dx <- dx%>%
      merge(metadata, by.x ="label",by.y= "Experiment")
   
   g2 <- ggplot(data = dx)+
      geom_point(aes_string(x= "PC1", y= "PC2", colour= colour, shape = shape), stat="identity")+
      theme_bw()
   
   
   if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"pcaplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"pcaplot"))}
   ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"pcaplot", "pcaplot.png"),g2, height=6, width=6, dpi=200)
   
   invisible(self)
})


#3. Boxplot from original log2(intensities)
LFQ_Template$set("public","boxplot", function(intensity = "iBAQ" ,facet.vars= c("Tissue","Antibody")){
   
   palette <- RColorBrewer::brewer.pal(length(unique(self$metadata$BufferNum)), "RdYlGn")
   
   if("BufferNum" %in% colnames(self$metadata)){
      
      dx <- self$proteinGroups%>%
         dplyr::select(contains(intensity))%>%
         pivot_longer(cols = everything())%>%
         filter(value>0)%>%
         mutate(value=log2(value), Experiment = gsub(paste0(intensity,"."),"",name))%>%
         merge(self$metadata, by = "Experiment")%>%
         arrange(BufferNum)%>%
         mutate(BufferNum = as.character(BufferNum))
      
      g2 <- dx%>% ggviolin(x = "name", y = "value",
                           fill = "BufferNum", palette= palette
      )+
         facet_grid(formula(paste("~",paste0(facet.vars , collapse="+"))), scales="free", space="free")+
         theme(axis.text.x = element_blank(),
               axis.title.x= element_blank(),
               axis.ticks.x=element_blank())+
         ylab(paste0("log2(",intensity,")"))+
         xlab("")
      
   } else{
      
      dx <- self$proteinGroups%>%
         dplyr::select(contains(intensity))%>%
         pivot_longer(cols = everything())%>%
         filter(value>0)%>%
         mutate(value=log2(value), Experiment = gsub(paste0(intensity,"."),"",name))%>%
         merge(self$metadata, by = "Experiment")
      
      g2 <- dx%>% ggviolin(x = "name", y = "value",
                           fill = "gray",
                           palette= palette
      )+
         facet_grid(formula(paste("~",paste0(facet.vars , collapse="+"))), scales="free", space="free")+
         theme(axis.text.x = element_blank(),
               axis.title.x= element_blank(),
               axis.ticks.x=element_blank())+
         ylab(paste0("log2(",intensity,")"))+
         xlab("")
   }
   
   
   
   if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"boxplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"boxplot"))}
   ggsave(file = file.path(paste0("output/plot/", basename(self$txt.dir)),"boxplot", "boxplot.png"),g2, height=8, width=(2+dim(dx)[1]/2000), dpi=200, limitsize=F)
   
   invisible(self)
})

#4. Plot before and after imputation

LFQ_Template$set("public","imputationplot", function(intensity = "iBAQ" , facet.vars= c("Tissue","Antibody")){
   
   dx <- self$significant_tbl
   dy <- self$experiments
   
   for(comparison in names(dx)){
      print(comparison)
      dxl <- dx[[comparison]][,grepl(intensity, colnames(dx[[comparison]]))]%>%
         pivot_longer(cols = everything())%>%
         mutate(cat ="pre-Imputation")%>%
         arrange(name)
      
      dyl <- dy[[comparison]][,grepl(intensity, colnames(dy[[comparison]]))]%>%
         pivot_longer(cols = everything())%>%
         mutate(cat ="post-Imputation")%>%
         arrange(name)
      
      df <- rbind(dxl,dyl)%>%
         mutate(cat= factor(cat , levels = c("pre-Imputation","post-Imputation")),
                name = gsub(paste0(intensity,"."),"",name))%>%
         merge(self$metadata, by.x = "name",by.y = "Experiment")
      
      g3 <- ggplot(df, aes(x = value, y = name)) +
         geom_density_ridges(
            scale=0.7,jittered_points = TRUE, position = "raincloud",
            point_size = 1, point_alpha = 1, alpha = 0.7,
         )+
         facet_grid(formula(paste(paste0(facet.vars , collapse="+"), "~ cat ")), scales="free", space="free")+
         theme_bw()+
         ggtitle(paste(c(strsplit(comparison, ".vs.")[[1]][1],strsplit(comparison, ".vs.")[[1]][2]), collapse=" ~ " ))
      
      
      if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"imputation"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"imputation"))}
      ggsave(file.path(paste0("output/plot/", basename(self$txt.dir)),"imputation", paste0(comparison,"_imputation.png")),g3, height=8, width=10, dpi=200)
      
      
   }
   
   invisible(self)
})

#5. Make Volcano plot
LFQ_Template$set("public","volcanoplot", function( genes.interest=c(), interactors=c() , multiwell= TRUE){
   
   if(multiwell){
      for(comp in names(self$significantGroup_tbl)){
         
         for(buffer in unique(self$significantGroup_tbl[[comp]]$Buffer)){
            dx <- self$significantGroup_tbl[[comp]]%>%
               filter(Buffer==buffer)
            shade = data.frame(x1=c( 1), 
                               x2=c( Inf),
                               y1=c( -log10(0.05)), 
                               y2=c( Inf))
            fold_cutoff = 1
            pvalue_cutoff = 0.05
            
            comparison <- paste0(comp, buffer, collapse=".")
            
            g4 <- ggplot(dx) +
               theme_bw()+
               geom_rect(data=shade, 
                         mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
               geom_point(data  = subset(dx, (Significant=="N")), aes(x =logFC, y = -log10(adj.P.Val),colour=Significant), alpha = 0.5,size=1)+
               geom_point(data  = subset(dx, (Significant=="Y")), aes(x =logFC, y = -log10(adj.P.Val),colour=Significant), alpha = 1,size=1)+
               geom_point(data  = subset(dx, (Gene.names %in% genes.interest)), aes(x =logFC, y = -log10(adj.P.Val)),colour="red", alpha = 1,size=3)+
               geom_point(data  = subset(dx, (Gene.names %in% interactors)), aes(x =logFC, y = -log10(adj.P.Val)),colour="green", alpha = 1,size=2)+
               geom_vline(xintercept = fold_cutoff, col = "blue")+
               geom_vline(xintercept = -fold_cutoff, col = "blue")+
               geom_hline(yintercept = -log10(pvalue_cutoff), col = "green2")+
               ggtitle(paste(c(strsplit(comparison, ".vs.")[[1]][1],strsplit(comparison, ".vs.")[[1]][2]), collapse=" ~ " ))+
               scale_colour_manual(values=c("gray","gray30"),
                                   limits = c("N","Y"))+
               geom_text_repel(data  = subset(dx, (Significant=="Y"|Gene.names ==genes.interest)),
                               aes(x=logFC, y=-log10(adj.P.Val),label=Gene.names),
                               segment.alpha =0.35,max.overlaps =25, point.padding = 0.2,
                               size = 2.5, color = ifelse(subset(dx, (Significant=="Y"|Gene.names ==genes.interest))$Gene.names %in% genes.interest, "red","black"))
            
            if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"volcanoplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"volcanoplot"))}
            ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"volcanoplot",paste0(comparison,"_volcanoplot.png")), g4, width=9, height=9, dpi=200)
            
         }
      }
   }else{
      
      for(comparison in names(self$significant_tbl)){
         
         dx <- self$significant_tbl[[comparison]]
         shade = data.frame(x1=c( 1), 
                            x2=c( Inf),
                            y1=c( -log10(0.05)), 
                            y2=c( Inf))
         fold_cutoff = 1
         pvalue_cutoff = 0.05
         
         
         g4 <- ggplot(dx) +
            theme_bw()+
            geom_rect(data=shade, 
                      mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
            geom_point(data  = subset(dx, (Significant=="N")), aes(x =logFC, y = -log10(adj.P.Val),colour=Significant), alpha = 0.5,size=1)+
            geom_point(data  = subset(dx, (Significant=="Y")), aes(x =logFC, y = -log10(adj.P.Val),colour=Significant), alpha = 1,size=1)+
            geom_point(data  = subset(dx, (Gene.names %in% interactors)), aes(x =logFC, y = -log10(adj.P.Val)),colour="green", alpha = 1,size=2)+
            geom_point(data  = subset(dx, (Gene.names %in% genes.interest)), aes(x =logFC, y = -log10(adj.P.Val)),colour="red", alpha = 1,size=3)+
            geom_vline(xintercept = fold_cutoff, col = "blue")+
            geom_vline(xintercept = -fold_cutoff, col = "blue")+
            geom_hline(yintercept = -log10(pvalue_cutoff), col = "green2")+
            ggtitle(paste(c(strsplit(comparison, ".vs.")[[1]][1],strsplit(comparison, ".vs.")[[1]][2]), collapse=" ~ " ))+
            scale_colour_manual(values=c("gray","gray30"),
                                limits = c("N","Y"))+
            geom_text_repel(data  = subset(dx, (Significant=="Y"|Gene.names ==genes.interest)),
                            aes(x=logFC, y=-log10(adj.P.Val),label=Gene.names),
                            segment.alpha =0.35,max.overlaps =25, point.padding = 0.2,
                            size = 2.5, color = ifelse(subset(dx, (Significant=="Y"|Gene.names ==genes.interest))$Gene.names %in% genes.interest, "red","black"))
         
         if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"volcanoplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"volcanoplot"))}
         ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"volcanoplot",paste0(comparison,"_volcanoplot.png")), g4, width=9, height=9, dpi=200)
         
      }
   }
   
   invisible(self)
})

#6. Draw heatmap across all buffers in a condition
LFQ_Template$set("public","plotHeatmap", function(intensity = "iBAQ",genes.interest=c(), n= 50, multiwell = FALSE){
   
   if(multiwell){
      for(comparison in names(self$significantGroup_tbl)){
         print(comparison)
         
         dx1 <- self$significantGroup_tbl[[comparison]]%>%
            mutate(Significant= ifelse(Significant=="Y",1,0))%>%
            dplyr::select(Buffer,Gene.names,uniprotID,Significant,logFC)
         
         
         
         proteins_top <- self$significantGroup_tbl[[comparison]]%>%
            mutate(Gene.names = ifelse(is.na(Gene.names), uniprotID, Gene.names))%>%
            group_by(uniprotID, Gene.names)%>%
            summarize(lfc.total = sum(logFC), Significant.total = sum(ifelse(Significant=="Y",1,0)))%>%
            arrange(desc(lfc.total ))%>%
            filter(Significant.total>0)%>%
            head(n)%>%
            .$uniprotID
         
         
         dx <- dx1%>%
            mutate(Gene.names = ifelse(is.na(Gene.names), uniprotID, Gene.names))%>%
            filter(uniprotID %in% proteins_top)%>%
            mutate(uniprotID = factor(uniprotID, levels = proteins_top))
         
         
         g5 <- ggplot(data= dx)+
            geom_tile(aes(x=Buffer, y = Gene.names, fill= logFC), color="gray30")+
            geom_point(aes(x= Buffer, y = Gene.names, alpha= Significant, size= Significant), 
                       color=ifelse(dx$Gene.names==genes.interest,"green","black"),
                       size=ifelse(dx$Gene.names==genes.interest,2,1))+
            scale_fill_gradient(low="blue",high="orange", na.value = "white")+
            scale_size(range = c(0,1.5))+
            guides(alpha="none")+
            theme_bw()+
            theme(axis.text.y = element_text(size=7),
                  axis.text.x = element_text(angle=90, size=8))+
            ggtitle(paste(paste(c(strsplit(comparison, ".vs.")[[1]][1],strsplit(comparison, ".vs.")[[1]][2]), collapse=" ~ " ),"- Top",n, "hits"))
         
         if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"heatmap"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"heatmap"))}
         ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"heatmap",paste0(comparison,"_heatmap.png")), g5, width=6, height=8, dpi=200)
      }
      
   }else{
      sig_tbl <- data.frame(matrix(NA, nrow=0, ncol= 5))
      colnames(sig_tbl) <- c("Gene.names", "uniprotID", "Significant", "logFC", "comparison")
      
      for(comparison in names(self$significant_tbl)){
         print(comparison)
         
         dx <- self$significant_tbl[[comparison]]%>%
            mutate(Gene.names = ifelse(is.na(Gene.names), uniprotID, Gene.names))%>%
            mutate(Significant= ifelse(Significant=="Y",1,0))%>%
            dplyr::select(Gene.names, uniprotID,Significant,logFC)%>%
            mutate(comparison = comparison)
         
         sig_tbl <- rbind(sig_tbl, dx)
      }
      
      
      proteins_top <- sig_tbl%>%
         group_by(uniprotID, Gene.names)%>%
         summarize(lfc.total = sum(logFC), Significant.total = sum(Significant))%>%
         arrange(desc(lfc.total ))%>%
         filter(Significant.total>0)%>%
         head(n)%>%
         .$uniprotID
      
      sig_tbl <- sig_tbl%>% 
         mutate(Gene.names = ifelse(is.na(Gene.names), uniprotID, Gene.names))%>%
         filter(uniprotID %in% proteins_top)
      
      
      g5 <- ggplot(data= sig_tbl)+
         geom_tile(aes(x=comparison, y = Gene.names, fill= logFC), color="gray30")+
         geom_point(aes(x= comparison, y = Gene.names, alpha= Significant, size= Significant), 
                    color=ifelse(sig_tbl$Gene.names==genes.interest,"green","black"),
                    size=ifelse(sig_tbl$Gene.names==genes.interest,2,1))+
         scale_fill_gradient(low="blue",high="orange", na.value = "white")+
         scale_size(range = c(0,1.5))+
         guides(alpha="none")+
         theme_bw()+
         theme(axis.text.y = element_text(size=7),
               axis.text.x = element_text(angle=90, size=8))+
         ggtitle(paste(paste(c(strsplit(comparison, ".vs.")[[1]][1],strsplit(comparison, ".vs.")[[1]][2]), collapse=" ~ " ),"- Top",n, "hits"))
      
      if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"heatmap"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"heatmap"))}
      ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"heatmap","heatmap.png"), g5, width=6, height=8, dpi=200)
   }
   
   
   invisible(self)
}) 

#7. Draw protein enrichment plot 
