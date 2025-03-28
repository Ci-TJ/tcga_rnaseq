setwd("G:/backup/webService/TCGA/tmp")

library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)
library(tidytable)

mat2plot <- function(project=c("TCGA-LUSC"), data_dir="./GDCdata", num_tp=100, num_nt=100,tp_t="TP", tp_n="NT", is_shor=FALSE, save=TRUE, target=c("FAM135B"), candidate="FAM135B"){
  if (file.exists("tmp") == FALSE){
    dir.create("tmp")
  }
  for (p in project){
    if (file.exists(file.path("tmp", p)) == FALSE){
      dir.create(file.path("tmp", p), recursive = TRUE)
      }
        
    query <- GDCquery(project = p,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts")
    samplesDown <- getResults(query,cols=c("cases"))
    #tumor samples
    dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = tp_t)
    #normal samples
    dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = tp_n)
    if (is_short == TRUR){
      if (all(grepl("TCGA",dataSmNT)) & all(grepl("TCGA",dataSmNT))){
        dataSmTP_short <- dataSmTP[1:ifelse(num_tp <= length(dataSmTP), num_tp, length(dataSmTP))]
        dataSmNT_short <- dataSmNT[1:ifelse(num_nt <= length(dataSmNT), num_nt, length(dataSmNT))]
        queryDown <- GDCquery(project = p, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts", 
                            barcode = c(dataSmTP_short, dataSmNT_short))
        } 
      else{
        dataSmTP_short <- dataSmTP[1:length(dataSmTP)]
        dataSmNT_short <- dataSmNT[1:length(dataSmNT)]
        queryDown <- GDCquery(project = p, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts", 
                            barcode = c(dataSmTP_short, dataSmNT_short))
        }} 
    else{
      dataSmTP_short <- dataSmTP
      dataSmNT_short <- dataSmNT
      queryDown <- GDCquery(project = p, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts", 
                            barcode = c(dataSmTP_short, dataSmNT_short))}
    
    dataPrep1 <- GDCprepare(query = queryDown, directory = data_dir, save = save, save.filename = file.path("tmp", p, ".rda"))
      
    #a step to remove sample outliers using pearson correlation
    rownames(dataPrep1) <- rowData(dataPrep1)$gene_name #transfer to gene names
    dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                            cor.cut = 0.6,)
    
    #step with library size and gcContent normalization using EDASeq
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                            geneInfo = geneInfoHT,
                                            method = "gcContent")
    #quantile filtering to remove genes with low count
    dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)

    #voom transformation of the data (log)
    v.dataFilt<-voom(dataFilt)
    #taking log transformed data for exploration of batch effects
    c.dataFilt <- TCGAbatch_Correction(tabDF = v.dataFilt, batch.factor="Plate", adjustment=c("TSS"), is_plot=FALSE)

    if (length(dataSmNT) > 3 ){
      if (candidate %in% rownames(c.dataFilt)){
        DEG <- TCGAanalyze_DEA(
          mat1=c.dataFilt[ifelse(length(target > 0), targe, rownames(c.dataFilt)), dataSmNT_short],
          mat2=c.dataFilt[ifelse(length(target > 0), targe, rownames(c.dataFilt)), dataSmTP_short],
          pipeline="limma",
          Cond1type = "Normal",
          Cond2type = "Tumor",
          method = "glmLRT")}
      else {
        DEG <- TCGAanalyze_DEA(
          mat1=c.dataFilt[ifelse(length(target > 1), targe, rownames(c.dataFilt)), dataSmNT_short],
          mat2=c.dataFilt[ifelse(length(target > 1), targe, rownames(c.dataFilt)), dataSmTP_short],
          pipeline="limma",
          Cond1type = "Normal",
          Cond2type = "Tumor",
          method = "glmLRT")}
      
      fwrite(as_tidytable(DEG, .keep_rownames = "gene_name"), file.path(p, "_deg.csv"))
      }
    else {
      print(paste(p, "It doesn't have enough normal samples!"))
      }
    tmp_mat <- as_tidytable(c.dataFilt, .keep_rownames = "gene_name")
    fwrite(tmp_mat, file.path(p, "_exp.csv"))
    }
  
  return(c.dataFilt)
  }
    

                            
    
                                         
                                         
    
    
  











































      
