###############################################################################################################
## Author: Mehrnoosh Oghbaie
## 01/09/2021
# This file includes most of the general functions that are going to be used in this project

######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p) 
    }
    library(p,character.only=T)
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p, dependencies = TRUE) 	
    }	
    library(p, character.only=T)
  }
}


########################################################################
## Build Complex_protein db
########################################################################
make_Complex_protein <- function(dir = "db/coreComplexes/coreComplexes.txt"){
  coreComplexes <- read.delim(dir)
  complex_protein <- data.frame(matrix(nrow=0, ncol=10))
  colnames(complex_protein) <- c("ComplexName",
                                 "UniprotID",
                                 "Count",
                                 "Identifiers",
                                 "Cell.line",
                                 "Organism",
                                 "GO.ID", 
                                 "GO.description",
                                 "Complex.comment",
                                 "GeneName")
  for(i in 1:dim(coreComplexes)[1]){
    print(i)
    uniprotID <- strsplit(as.character(coreComplexes[i,"subunits.UniProt.IDs."]),";")[[1]]
    GeneName <- strsplit(as.character(coreComplexes[i,"subunits.Gene.name."]),";")[[1]]
    n <- length(uniprotID)
    Identifiers <- paste(uniprotID, collapse="|")
    
    if(length(GeneName)== length(uniprotID)){
      GeneName <- GeneName
    } else{
      GeneName <- NA
    }
    
    complex <- data.frame(ComplexName = rep(coreComplexes[i,"ComplexName"],n),
                          uniprotID = uniprotID,
                          count = rep(n,n),
                          Identifiers = rep(Identifiers,n),
                          Cell.line = rep(coreComplexes[i,"Cell.line"],n),
                          Organism = rep( coreComplexes[i,"Organism"],n),
                          GO.ID = rep( coreComplexes[i,"GO.ID"],n),
                          GO.description = rep( coreComplexes[i,"GO.description"],n),
                          Complex.comment = rep( coreComplexes[i,"Complex.comment"],n),
                          GeneName = GeneName)
    
    complex_protein <- rbind(complex_protein, complex)  
    
  }
  return(complex_protein)
}



################################################################################
### Enrich list of significant proteins by complex protein relationship

enrichComplex <- function(significant_uniprot,complex_protein){
  ### list of annotated proteins
  k <- length(significant_uniprot)
  
  ### list of proteins in human complexes
  N <-length(unique(unlist(complex_protein$uniprotID)))
  
  selected_complex_count <- complex_protein %>%
    dplyr::filter(uniprotID %in% significant_uniprot) %>%
    dplyr::group_by(ComplexName) %>%
    dplyr::summarise(count=as.numeric(as.character(count[1])),
                     selected_count= as.numeric(n()),
                     GeneRatio=selected_count/count,
                     complete_members=Identifiers[1],
                     Identifiers=paste(uniprotID, collapse="|"),
                     Cell.line = Cell.line[1],
                     Organism = Organism[1],
                     GO.ID = GO.ID[1],
                     GO.description = GO.description[1],
                     Complex.comment = Complex.comment[1]) %>%
    dplyr::mutate(pValue = phyper(q=as.numeric(selected_count),m=as.numeric(count), n =N-as.numeric(count), k=k , lower.tail=FALSE))
  
  
  selected_complex_count$p.adj <- p.adjust(selected_complex_count$pValue, "BH")
  selected_complex_count_test <- selected_complex_count %>% filter(p.adj<5e-02&as.numeric(count)>1) %>% arrange(p.adj,desc(selected_count))
  return(selected_complex_count_test)
}


#####################################################################################
### Impute
#####################################################################################
## impute(df, amm = "2", pmm = "6")

count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}

calculate_stats_nonzeros <- function(data){
  ## Predicting all missing replicates values
  colnames <- colnames(data)
  sample <- data
  stats <- foreach(i=1:dim(data)[2], .combine=rbind)%do%{
    cbind(mean(sample[,i], na.rm=T),sd(sample[,i], na.rm=T))
  }
  
  return(stats)
}

impute_all_zeros <- function(data,amm = "2", pmm ="6"){
  ## Predicting all missing replicates values
  colnames <- colnames(data)
  na_zeros <- rowSums(data,na.rm=T)==0
  print(paste("There are", sum(na_zeros), "total zero records in either case or control replicates to be imputed."))
  if(sum(na_zeros)>0){
    min_zero_col <- min_col(apply(data,2,count_zeros))
    if(amm == "2"){
      sample_min_zeros <- unname(unlist(data[,min_zero_col]))
      mu <- mean(sample_min_zeros, na.rm=T)
      sd <- sd(sample_min_zeros, na.rm=T)
      data[na_zeros,colnames] <- matrix(na_zeros_impute(na_zeros* length(colnames) , mu,sd),ncol= length(colnames) , nrow = sum(na_zeros))
    } 
    if(amm == "1"){
      sample <- data
      stats <- calculate_stats_nonzeros(data)
      
      dat <- matrix(NA,nrow= sum(na_zeros),ncol=0)
      for(i in 1:dim(stats)[1]){
        dat <- cbind(dat,na_zeros_impute(na_zeros , stats[1,1],stats[1,2]))
      }
      data[na_zeros,colnames] <- matrix(dat,
                                        ncol= length(colnames),
                                        nrow = sum(na_zeros))
    }
  }
  return(data)
}


impute_partial_zeros <- function(data, count_na, colnames, nb=0){
  if(length(colnames)==2){
    y <- data[count_na>nb,colnames]
  } else{
    y <- data[count_na>nb,colnames]
  }
  
  if(length(colnames)>2){
    for(i in 1:dim(y)[1]){
      #print(i)
      col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
      for(j in col_NAs){
        col_NA <- j
        col_select <- names(y[i,colnames(y)!=col_NA])
        sample <- data[complete.cases(data[,col_select]),col_select]
        delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
        mu <- mean(unlist(unname(delta)), na.rm=T)
        std <- sd(unlist(unname(delta)),na.rm=T)
        cor <- cor(data[count_na==0,colnames])
        mean_cor <- mean(cor[col_NA,col_select])
        deltanew <- rnorm(1,mu, std/(sqrt(2)*mean_cor))
        y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
      }
    }
  }else{
    for(i in 1:dim(y)[1]){
      col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
      for(j in col_NAs){
        col_NA <- j
        col_select <- names(y)
        sample <- data[complete.cases(data[,col_select]),col_select]
        delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
        mu <- mean(unlist(unname(delta)), na.rm=T)
        std <- sd(unlist(unname(delta)),na.rm=T)
        cor <- cor(data[count_na==0,colnames])
        mean_cor <- mean(cor[col_NA,col_select])
        deltanew <- rnorm(1,mu, std/(sqrt(2)*mean_cor))
        y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
      }
    }
  }
  
  if(length(colnames)==2){
    data[count_na>nb,colnames] <- y
  }else{
    data[count_na>nb,colnames] <- y
  }
  
  return(data)
}


impute <- function(data, amm = "2", pmm = "6"){
  data[data==0] <- NA
  stats <- calculate_stats_nonzeros(data)
  data <- impute_all_zeros(data,amm )
  colnames <- colnames(data)
  count_na <- apply(data,1, count_zeros)
  
  if(pmm == "7"){
    data <- impute_partial_zeros(data, count_na, colnames, nb = 0)
  }
  
  if(pmm == "6"){
    data <- impute_partial_zeros(data, count_na, colnames, nb = 1)
    
    for( i in 1:length(colnames)){
      data_missing <- unlist(lapply(data[,colnames[i]],  function(x) is.na(x)))
      data[data_missing,colnames[i]] <- na_zeros_impute(sum(data_missing) , stats[i,1], stats[i,2])
    }
  }
  
  return(data)
}


##########################################################################################
## Enrichment analysis with clusterprofiler
##########################################################################################

performEnreachment <- function(gene.df, cat='BP'){
  
  ego <- enrichGO(gene          = gene.df[["ENTREZID"]],
                  # universe      = gene.df$SYMBOL,
                  OrgDb         = org.Hs.eg.db,
                  ont           = cat,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
  egos <-  clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  
  return(egos@result)
}