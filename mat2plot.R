setwd("G:/backup/webService/TCGA/tmp")

library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)

mat2plot <- function(project=c("TCGA-ACC"), data_dir="./GDCdata", num_tp=100, num_nt=100,tp_t="TP", tp_n="NT", is_shor=FALSE){
  for (p in project){
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
        dataSmTP_short <- dataSmTP[1:ifele(num_tp <= length(dataSmTP), num_tp, ;length(dataSmTP))]
        dataSmNT_short <- dataSmNT[1:ifele(num_tp <= length(dataSmNT), num_tp, ;length(dataSmNT))]
        queryDown <- GDCquery(project = p, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts", 
                            barcode = c(dataSmTP_short, dataSmNT_short))
        } else{
        dataSmTP_short <- dataSmTP[1:length(dataSmTP)]
        dataSmNT_short <- dataSmNT[1:length(dataSmNT)]
        queryDown <- GDCquery(project = p, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts", 
                            barcode = c(dataSmTP_short, dataSmNT_short))
        }} else{
      dataSmTP_short <- dataSmTP
      dataSmNT_short <- dataSmNT
      queryDown <- GDCquery(project = p, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts", 
                            barcode = c(dataSmTP_short, dataSmNT_short)}
    
    queryDown <- GDCquery(project = p, 
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "STAR - Counts", 
                          barcode = c(dataSmTP_short, dataSmNT_short))
    dataPrep1 <- GDCprepare(query = queryDown, directory = data_dir, save = TRUE, save.filename = file.path(p, ".rda"))
      
    #a step to remove sample outliers using pearson correlation
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
    
                                         
                                         
    
    
  











































      
