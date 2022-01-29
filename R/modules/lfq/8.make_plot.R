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
      select(`Raw.file`,contains("Contaminants"), contains("Peptide.Intensity"))%>%
      pivot_longer(!`Raw.file`)%>%
      mutate(name = ifelse(grepl("Contaminants",name),"Contaminants","Peptide.Intensity"))%>%
      merge(self$metadata, by.x = "Raw.file",by.y = "Experiment")%>%
      ggplot(aes(y=`Raw.file`,x=value, fill=name))+
      geom_col()+
      facet_grid(formula(paste0(paste0(facet.vars,collapse="+") ,"~name")), space="free", scale="free")+
      theme_bw()+
      theme(axis.text.y = element_text(size=7))
   
   
   if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa"))}
   ggsave(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa", "qa.plot.png"),plot=g1,  width=8, height=(2+dim(dx)[1]/10), limitsize = F)
   
   invisible(self)
})

#2, Pca plot
LFQ_Template$set("public","pcaplot", function(intensity = "iBAQ", colour = "Tissue" ,shape = "Antibody" ){
   
   metadata <- self$metadata
   di <- self$proteinGroups %>% 
      select(contains(intensity))
   
   di <- di[apply(di,1,sum)>0,]
   di <- log2(di+1)
   df <- data.frame(t(di))
   
   
   pca_res <- prcomp(df, scale. = TRUE)
   
   dx <- data.frame(pca_res$x)%>%select(PC1,PC2)
   dx[["label"]] <- gsub(paste0(intensity,"."),"",rownames(df))
   dx <- dx%>%
      merge(metadata, by.x ="label",by.y= "Experiment")
   
   g2 <- ggplot(data = dx)+
      geom_point(aes_string(x= "PC1", y= "PC2", colour= colour, shape = shape))+
      theme_bw()
   
   
   if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"pcaplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"qa"))}
   ggsave(file.path(paste0("output/plot/", basename(self$txt.dir)),"pcaplot", "pcaplot.plot.png"),g2, height=6, width=6, dpi=200)
   
   invisible(self)
})


#3. Boxplot from original log2(intensities)
LFQ_Template$set("public","boxplot", function(intensity = "iBAQ" ,facet.vars= c("Tissue","Antibody")){
   
   palette <- RColorBrewer::brewer.pal(length(unique(self$metadata$BufferNum)), "RdYlGn")
   
   if("BufferNum" %in% colnames(self$metadata)){
      
      g2 <- self$proteinGroups%>%
         select(contains(intensity))%>%
         pivot_longer(cols = everything())%>%
         filter(value>0)%>%
         mutate(value=log2(value), Experiment = gsub(paste0(intensity,"."),"",name))%>%
         merge(self$metadata, by = "Experiment")%>%
         arrange(BufferNum)%>%
         mutate(BufferNum = as.character(BufferNum))%>%
         ggviolin(x = "name", y = "value",
                  fill = "BufferNum", palette= palette
         )+
         facet_grid(formula(paste("~",paste0(facet.vars , collapse="+"))), scales="free", space="free")+
         theme(axis.text.x = element_blank(),
               axis.title.x= element_blank(),
               axis.ticks.x=element_blank())+
         ylab(paste0("log2(",intensity,")"))+
         xlab("")
   } else{
      
      g2 <- self$proteinGroups%>%
         select(contains(intensity))%>%
         pivot_longer(cols = everything())%>%
         filter(value>0)%>%
         mutate(value=log2(value), Experiment = gsub(paste0(intensity,"."),"",name))%>%
         merge(self$metadata, by = "Experiment")%>%
         ggviolin(x = "name", y = "value",
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
   ggsave(file.path(paste0("output/plot/", basename(self$txt.dir)),"boxplot", "boxplot.plot.png"),g2, height=8, width=(2+dim(dx)[1]/10), dpi=200, limitsize=F)
   
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
 
 LFQ_Template$set("public","plotHeatmap", function(intensity = "iBAQ", n= 50){
    
    for(comparison in names(self$significantGroup_tbl)){
       print(comparison)
       
       dx1 <- self$significantGroup_tbl[[comparison]]%>%
          mutate_at(vars(contains("Significant")) ,~ ifelse(.=="Y",1,0))%>%
          filter_at(vars(contains("Significant")),~ sum(.)>0)%>%
          mutate(Gene.names= as.factor(Gene.names))%>%
          select(Gene.names,contains("Significant.Buffer"))%>%
          pivot_longer(!Gene.names)%>%
          mutate(name= gsub("Significant.","", name))%>%
          rename(Significant = value)
       
       dx2 <- self$significantGroup_tbl[[comparison]]%>%
          mutate_at(vars(contains("Significant")) ,~ ifelse(.=="Y",1,0))%>%
          filter_at(vars(contains("Significant")),~ sum(.)>0)%>%
          mutate(Gene.names= as.factor(Gene.names))%>%
          select(Gene.names,contains("lfc.Buffer"))%>%
          pivot_longer(!Gene.names)%>%
          mutate(name= gsub("lfc.","", name))%>%
          rename(lfc = value)
       
       genes_diff <- self$significantGroup_tbl[[comparison]]%>%
          mutate_at(vars(contains("lfc")),~rank(.))%>%
          mutate(totalrank = rowSums(across(contains("lfc"))))%>%
          arrange(desc(totalrank))%>%
          head(n)%>%
          .$Gene.names
       
       genes_top <- self$significantGroup_tbl[[comparison]]%>%
          filter(significant.tag=="Y")%>%
          arrange(desc(lfc.tag))%>%
          head(n)%>%
          .$Gene.names
       
       dx <- cbind(dx1,dx2%>%
                      select(lfc))%>%
          filter(Gene.names %in% genes_top)%>%
          mutate(Gene.names = factor(Gene.names, levels = genes_top))
       
       
       g5 <- ggplot(data= dx)+
          geom_tile(aes(x=name, y = Gene.names, fill= lfc), color="gray")+
          geom_point(aes(x= name, y = Gene.names, alpha= Significant, size= Significant))+
          scale_fill_gradient(low="blue",high="orange")+
          scale_size(range = c(0,1.5))+
          theme_bw()+
          theme(axis.text.y = element_text(size=7),
                axis.text.x = element_text(angle=90, size=8))+
          ggtitle(paste(paste(c(strsplit(comparison, ".vs.")[[1]][1],strsplit(comparison, ".vs.")[[1]][2]), collapse=" ~ " ),"- Top",n, "hits"))
       
       if(!dir.exists("output/plot/self/heatmap")){dir.create("output/plot/self/heatmap")}
       ggsave(file=file.path("output/plot/self/heatmap",paste0(comparison,"_heatmap.png")), g5, width=6, height=12, dpi=200)
       
       
    }
    
    invisible(self)
 })
 
 
 #7. Draw protein enrichment plot 
 