
IDIRT_Template$set("public", "CreateScatterPlot", function(intensity = "iBAQ"){
  metadata <- self$metadata
  df_avg <- self$avg_values
  # get list of conditions from metadata table
  experiments <- self$experiments
  mix_percentage <- df_avg[[paste0("Avg.Percentage.H.L.",unique(metadata$Experiment_label[metadata$Label=="Mix"]))]]
  
  for(exp in experiments[!grepl("mix",experiments)]){
    
    print(exp)
    
    dp_exp <- df_avg[[paste0("Avg.Percentage.H.L.",exp)]]
    di_exp <- df_avg[[paste0(c("Avg",intensity,exp), collapse=".")]]
    
    dx <- data.frame(percentage = dp_exp,
                     log2Intensity = di_exp,
                     Significant = df_avg[[paste0(c("Significant",exp), collapse=".")]],
                     mix_percentage = mix_percentage,
                     uniprotID = df_avg$uniprotID,
                     Gene.names = df_avg$Gene.names
    )%>%
      mutate(in_mix = ifelse(percentage>0&mix_percentage==0, "No","Yes"))
    
    dx$Gene.names[dx$uniprotID=="O00370"] <- "LORF2"
    
    hits <- df_avg$uniprotID[df_avg[[paste0(c("Significant",exp), collapse=".")]]==1]
    n_hits <- length(hits)
    
    stats <- dx%>%
      filter(mix_percentage>0)%>%
      summarize(
        median = median(mix_percentage),
        Q1=quantile(mix_percentage, probs = 0.25),
        Q3=quantile(mix_percentage, probs = 0.75),
        IQR = Q3-Q1,
        outlier = median+1.5*IQR
        
      )
    
    IP_stats <- dx%>%
      filter(percentage>0)%>%
      summarize(
        median = median(percentage),
        Q1=quantile(percentage, probs = 0.25),
        Q3=quantile(percentage, probs = 0.75),
        IQR = Q3-Q1,
        outlier = median+1.5*IQR
        
      )
    
    
    if(n_hits >= 5){ # only plot if there are at least 5 hits
      
      dx0 <- subset(dx,dx$Significant==1)
      
      
      p_scatter <- ggplot(dx, aes(x=percentage,y=log2Intensity)) +
        geom_vline(xintercept= stats$outlier, colour="grey60", lty=2, size=0.8)+
        geom_vline(xintercept= IP_stats$outlier, lty=2, size=0.8)+
        geom_point(colour =ifelse(dx$percentage>stats$outlier,"red","grey80"), size=2) +
        geom_text_repel(data = dx0, aes(label= Gene.names, colour = in_mix),size=3)+
        scale_colour_manual(limits=c("Yes","No"), values = c("black","green2"))+
        xlab("Avg H/L Percentage ")+
        ylab("Avg log LFQ intensity H") +
        xlim(0,1) +
        theme_classic()+
        theme(legend.position = "none")
      
      
      
      
      p_scatter_extend <- ggplot(dx, aes(x=percentage,y=log2Intensity)) +
        geom_vline(xintercept= stats$outlier, colour="grey60", lty=2, size=0.8)+
        geom_vline(xintercept= IP_stats$outlier, lty=2, size=0.8)+
        geom_point(colour =ifelse(dx$percentage>stats$outlier,"red","grey80"), size=2) +
        geom_text_repel(data = dx0, aes(label= Gene.names, colour = in_mix),size=3)+
        scale_colour_manual(limits=c("Yes","No"), values = c("black","green2"))+
        #xlab("Avg H Percentage ")+
        #ylab("Avg log LFQ intensity H") +
        xlim(stats$outlier, 1) +
        theme_classic()+
        theme(legend.position = "none")
      
      mix <- dx$mix_percentage[dx$mix_percentage>0]
      ip <- dx$percentage[dx$percentage>0]
      xmin <- density(mix)$x[density(mix)$y ==max(density(mix)$y)]
      xmax <- density(ip)$x[density(ip)$y ==max(density(ip)$y)]
      ypos <- max(max(density(ip)$y), max(density(mix)$y))
      yseg <- min(max(density(ip)$y), max(density(mix)$y))
      xseg <- ifelse(yseg ==max(density(ip)$y),xmax,xmin)
      
      
      p_density_top <-  ggplot() +
        geom_density(data = subset(dx, percentage>0),aes(x=percentage), alpha=1, size=1)+
        geom_density(data = subset(dx, mix_percentage>0),aes(x=mix_percentage), alpha=1, size=1, colour = "grey")+
        xlab('H/L Percentage')+
        ylab("H/L % density ")+
        theme(legend.position = "none")+
        geom_vline(xintercept= stats$outlierNormal, colour="grey60", lty=2 , size=0.8)+
        geom_vline(xintercept= IP_stats$outlier, colour="black",lty=2, size=0.8)+
        theme_classic()+
        scale_color_grey()+
        geom_segment(aes(x = xmin, y = ypos, xend = xmax, yend = ypos),size=1)+
        geom_segment(aes(x = xseg, y = yseg, xend = xseg, yend = ypos),size=1)+
        annotate("text",x=mean(xmin,xmax),
                 y=(ypos+0.6), 
                 label = paste("p =",signif(t.test(dx$mix_percentage,dx$percentage)$p.value,3)))+
        xlim(c(0,1))+
        ylim(c(0,ypos+1))+
        guides(colour = "none", size = 'none',  fill="none")
      
      
      
      p_density_right = ggplot(subset(dx,log2Intensity>0)) +
        geom_density(aes(x = log2Intensity, color = in_mix)) +
        scale_color_manual(limits=c("Yes","No"), values = c("black","green2"))+
        theme_classic() +
        xlab("") +
        coord_flip() 
      #theme(legend.position = "none")
      
      pempty = ggplot(dx, aes(x=percentage, y=log2Intensity)) +
        geom_blank() +
        theme_void() +
        xlab("")
      
      
      
      if(n_hits<46){
        p = grid.arrange(p_density_top, pempty, p_scatter, p_density_right,
                         ncol = 2, nrow = 2, widths = c(3,1), heights = c(1, 3))
        save_plots=T
        if(save_plots){
          
          if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir))))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir))))}
          
          if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"scatterplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"scatterplot"))}
          ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"scatterplot",paste0(exp,"_scatter.png")), 
                 p,  
                 width=11,
                 height=11,
                 dpi=300)
          
          
        }
        
      }else{
        
        p  <- plot_grid(p_scatter_extend,plot_grid(p_density_top, pempty,p_scatter,p_density_right,  align="hv",
                                                   rel_widths=c(3,1), rel_heights=c(2,3)), rel_widths=c(3,3))
        save_plots=T
        if(save_plots){
          
          if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"scatterplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"scatterplot"))}
          ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"scatterplot",paste0(exp,"_scatter.png")), 
                 p,  
                 width=16,
                 height=11,
                 dpi=300)
          
        }
      }
      
    }
  }
  
  invisible(self)
})

IDIRT_Template$set("public", "CreateDotPlot", function(){
  
  df_avg <- self$avg_values
  metadata <- self$metadata
  experiments <- self$experiments
  
  dt <- df_avg%>%
    filter(uniprotID %in% c("O00370","Q9UN81"))%>%
    select(contains("Avg.Percentage"))%>%
    select(!contains("mix"))
  
  rownames(dt) <- c("LORF2","L1RE1")
  dt <- data.frame(t(dt))
  dt[["p.value"]] <- NA
  for(exp in experiments[!grepl("mix", experiments)]){
    print(exp)
    
    x <- df_avg%>%
      select(contains("Avg.Percentage"))%>%
      select(contains(exp))
    
    y <-  df_avg%>%
      select(contains("Avg.Percentage"))%>%
      select(contains("Mix"))
    
    if(length(x[!is.na(x)])>2){
      dt[rownames(dt)==paste0("Avg.Percentage.H.L.",exp), "p.value"] <- signif(t.test(x[!is.na(x)], y[!is.na(y)])$p.value,3) 
    }
    
  }
  
  
  dt2 <- dt
  dt2[["condition"]] <- paste0(metadata$BufferChar[match(gsub("Avg.Percentage.H.L.","",rownames(dt)), metadata$Experiment_label)]," (",dt2$p.value,")")
  dt2 <- dt2%>%pivot_longer(cols = c("L1RE1","LORF2"),names_to="target", values_to = "avg.percentage")
  
  dotplot <- ggplot(dt2, aes(avg.percentage, condition)) +
    geom_line(aes(group = condition), size=1) +
    geom_point(aes(color = target), size=3)+
    geom_vline(xintercept=mean(dt$L1RE1) , lty=2)+
    scale_y_discrete(position = "right")+
    theme_classic()+ 
    ylab("")+
    xlab("Normalized log LFQ intensity")
  
  if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"dotplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"dotplot"))}
  ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"dotplot","dotplot.png"), 
         dotplot,
         width=11,
         height=7,
         dpi=300)
  
  invisible(self)
})

IDIRT_Template$set("public", "CreateTilePlot", function(top_n = 80){
  
  df_avg <- self$avg_values
  metadata <- self$metadata
  experiments <- self$experiments[!grepl("mix",self$experiments)]
  
  sig <- df_avg%>%
    select(uniprotID, Gene.names, contains("Significant"))
  sig$Gene.names[sig$uniprotID=="O00370"] <- "LORF2"
  gene_list <- sig$Gene.names[order(-apply(sig[,-c(1,2)],1,sum))]
  
  
  
  df_sig <- sig %>% 
    pivot_longer(!c(uniprotID, Gene.names), names_to = "condition", values_to = "Significant")%>%
    filter(Significant>0)%>%
    mutate( Gene.names = factor(Gene.names, levels = unique(gene_list)))%>%
    filter(Gene.names %in% gene_list[1:top_n])
  df_sig$condition <- metadata$BufferChar[match(gsub("Significant.","",df_sig$condition), metadata$Experiment_label)]
  
  tileplot <- ggplot() + 
    geom_tile(data = df_sig, aes(y=Gene.names, x=Significant), fill="red",color="gray85", alpha=1)+
    # scale_fill_manual(values = c("gray","orange","red2"), labels = c("no hit", "hit", "extreme hit"))+
    facet_grid(~condition, scales="free", space="free",shrink=TRUE, switch = "both")+
    theme_bw()+
    theme(
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=8),
      strip.text.x = element_text(
        size = 8, angle=90, face = "bold"
      ),
      strip.text.y.left = element_text(
        size = 8, color = "red", angle= 0, face = "bold"
      )) +xlab("Condition")
  
  
  if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"tileplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"tileplot"))}
  ggsave(file=file.path(paste0("output/plot/", basename(self$txt.dir)),"tileplot","tileplot.png"), 
         tileplot,
         width=8,
         height=12,
         dpi=300)
  
  invisible(self)
})

IDIRT_Template$set("public", "MDSPlot", function(intensity ="iBAQ"){
  
  df_avg <- self$avg_values
  metadata <- self$metadata
  experiments <- self$experiments[grepl("mix",self$experiments)]
  
  sig <- df_avg%>%
    select(uniprotID, Gene.names, contains("Significant"))
  sig$Gene.names[sig$uniprotID=="O00370"] <- "LORF2"
  protein_list <- sig$uniprotID[apply(sig[,-c(1,2)],1,sum)>0]
  
  df_avg[["hit_count"]] <- apply(sig[,-c(1,2)],1,sum)
  
  df.hits <- df_avg%>%
    select(uniprotID, Gene.names,hit_count, contains("Avg.Percentage"), contains(intensity))%>%
    filter(uniprotID %in% protein_list)
  
  # MDS
  dist.hits = dist(df.hits %>% select(-uniprotID, -Gene.names, -hit_count))
  mds.hits = as.data.frame(cmdscale(dist.hits, k=3))
  
  mds.hits <- cbind(df.hits %>% select(uniprotID, Gene.names, hit_count),mds.hits)
  mds.hits$Gene.names[mds.hits$uniprotID=="O00370"] <- "LORF2"
  
  fig.hits.mds = plot_ly(mds.hits, 
                         type = "scatter3d",
                         x = ~V1, 
                         y = ~V2, 
                         z = ~V3,
                         text = ~Gene.names,
                         color = ~hit_count,
                         hoverinfo = 'text'
  )
  fig.hits.mds = fig.hits.mds %>% layout(title = paste0("MDS \n Protein Heavy % + LFQ Intensity in IP Conditions "))
  
  if(!dir.exists(file.path(paste0("output/plot/", basename(self$txt.dir)),"mdsplot"))){dir.create(file.path(paste0("output/plot/", basename(self$txt.dir)),"mdsplot"))}
  
  htmlwidgets::saveWidget(as_widget(fig.hits.mds), file.path(paste0("output/plot/", basename(self$txt.dir)),"mdsplot","mdsplot.html"))
  
  invisible(self)
  
})