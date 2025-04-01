#Bug with TCGAbatch_Correction, you need set is_plot for sva to avoid plot.
setwd("G:/backup/webService/TCGA/tmp")

library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)
library(tidytable)
library(ggplot2)
library(dplyr)
library(DGEobj.utils)
library(readxl)

data_pre <- function(df, cut=0.25, is_filt=TRUE){
  dataNorm <- TCGAanalyze_Normalization(tabDF = df,
                                        geneInfo = geneInfoHT,
                                        method = "gcContent")
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
    
  #quantile filtering to remove genes with low count
  if (is_filt){
    dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                          method = "quantile", 
                                          qnt.cut = cut)
  } else {
    dataFilt <- dataNorm
  }
  return(dataFilt)
}


tcga2gtex <- function(project=c("TCGA-LUSC"), data_dir="./GDCdata", num_tp=100, num_nt=100,tp_t="TP", tp_n="NT", 
                     is_short=FALSE, save=TRUE, target=c("FAM135B"), candidate="FAM135B",is_voom=TRUE, is_log=FALSE, norm_method="none",
                     prior.count=0, unit="tpm",cut=0.25, is_filt=TRUE){
  if (file.exists("tmp") == FALSE){
    dir.create("tmp")
  }
  results <- list()

  c2n <- readxl::read_xlsx("tcga2gtex.xlsx") #cancer to normal tissues
  id2tissue <- fread("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
  df_split <- c2n %>% separate_rows(GTEx_SMTSD, sep = ";") 
  if (is_voom){
    gtex_data <- fread("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct", skip = 2, header = TRUE, sep = "\t") #check the path
  } else {
    gtex_data <- fread("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct", skip = 2, header = TRUE, sep = "\t") #check the path
  }
  
  gtex_data <- gtex_data %>% distinct(Description, .keep_all = T) #if by="Description", it will change the colnames to "by"
  gtex_data <- data.frame(gtex_data); rownames(gtex_data) <- gtex_data$Description; gtex_data <- gtex_data[,-c(1,2)]
  colnames(gtex_data) <- gsub("\\.", "-", colnames(gtex_data)) #colnames change after distinct()
  
  #
  for (p in project){
    if (file.exists(file.path("tmp", p)) == FALSE){
      dir.create(file.path("tmp", p), recursive = TRUE)
    }
    #######
    if (file.exists(file.path("tmp", p, paste0(p,"_gtex_exp.csv"))) == FALSE){
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
      if (is_short == TRUE){
        if (all(grepl("TCGA",dataSmNT)) & all(grepl("TCGA",dataSmTP))){
          dataSmTP_short <- dataSmTP[1:ifelse(num_tp <= length(dataSmTP), num_tp, length(dataSmTP))]
          dataSmNT_short <- dataSmNT[1:ifelse(num_nt <= length(dataSmNT), num_nt, length(dataSmNT))]
          queryDown <- GDCquery(project = p,
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification", 
                                workflow.type = "STAR - Counts", 
                                barcode = c(dataSmTP_short, dataSmNT_short))
        } else {
          dataSmTP_short <- dataSmTP[1:length(dataSmTP)]
          dataSmNT_short <- dataSmNT[1:length(dataSmNT)]
          queryDown <- GDCquery(project = p, 
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification", 
                                workflow.type = "STAR - Counts", 
                                barcode = c(dataSmTP_short, dataSmNT_short))
        }
      } else {
        dataSmTP_short <- dataSmTP
        dataSmNT_short <- dataSmNT
        queryDown <- GDCquery(project = p, 
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification", 
                              workflow.type = "STAR - Counts", 
                              barcode = c(dataSmTP_short, dataSmNT_short))
      }
    
      dataPrep1 <- GDCprepare(query = queryDown, directory = data_dir, save = save, save.filename = file.path("tmp", p, paste0(p,"_gtex.rda")))
      
      #a step to remove sample outliers using pearson correlation
      #rownames(dataPrep1) <- rowData(dataPrep1)$gene_name #transfer to gene names
      dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                            cor.cut = 0.6,)
      if (is_voom) {
        #step with library size and gcContent normalization using EDASeq
        dataFilt <- data_pre(df=dataPrep, is_filt=is_filt, cut=cut)
        id2s <- as_tidytable(data.frame(rowData(dataPrep1))) %>% select(gene_id,gene_name) %>% mutate(gene_id = stringr::str_remove(gene_id, "\\..*"))
        id2s <- id2s %>% filter(gene_id %in% rownames(dataFilt)) %>% distinct(gene_name, .keep_all = T)
        dataFilt <- dataFilt[id2s$gene_id,] #filter the genes with redundancy gene names
        #make sure the order
        if (all(rownames(dataFilt) == id2s$gene_id)){
          rownames(dataFilt) <- id2s$gene_name
        } else {
          print("Order is Wrong!!!")
        }

        #voom transformation of the data (log)
        v.dataFilt<-voom(dataFilt)
        #taking log transformed data for exploration of batch effects
        #c.dataFilt <- TCGAbatch_Correction(tabDF = v.dataFilt, batch.factor="Plate", adjustment=c("TSS"), is_plot=FALSE)
        c.dataFilt <- v.dataFilt$E #初始化为voom转换矩阵，确保后续代码可以继续运行
      } else {
        rownames(dataPrep) <- gsub("[.].*", "",rownames(dataPrep))
        ss <- intersect(rownames(dataPrep), rownames(geneInfoHT))
        dataNorm <- dataPrep[ss,]
        len_info <- geneInfoHT[ss,]$geneLength
        v.dataFilt <- convertCounts(
          countsMatrix = dataNorm,
          unit = unit,
          geneLength = len_info,
          log = is_log,
          normalize  = norm_method,
          prior.count = prior.count #note prior.count=1 ！= log2(tpm + prior.count), prior.count=1 but not recemmend for DE analysis!
        )
        
        if (is_log == FALSE){
          v.dataFilt <- log2(v.dataFilt + 1) #note convertCounts is log2(tpm), not log2(tpm+1)
        }
        id2s <- as_tidytable(data.frame(rowData(dataPrep1))) %>% select(gene_id,gene_name) %>% mutate(gene_id = stringr::str_remove(gene_id, "\\..*"))
        id2s <- id2s %>% filter(gene_id %in% rownames(v.dataFilt)) %>% distinct(gene_name, .keep_all = T)
        v.dataFilt <- v.dataFilt[id2s$gene_id,] #filter the genes with redundancy gene names
        #make sure the order
        if (all(rownames(v.dataFilt) == id2s$gene_id)){
          rownames(v.dataFilt) <- id2s$gene_name
        } else {
          print("Order is Wrong!!!")
        }
        #taking log transformed data for exploration of batch effects
        #c.dataFilt <- TCGAbatch_Correction(tabDF = v.dataFilt, batch.factor="Plate", adjustment=c("TSS"), is_plot=FALSE)
        c.dataFilt <- v.dataFilt #初始化为voom转换矩阵，确保后续代码可以继续运行
      }
      
      #############################
      ###tgex
      #############################
      tissue <- df_split %>% filter(project == p) %>% pull(GTEx_SMTSD)
      gtex_id <- id2tissue %>% filter(SMTSD %in% tissue) %>% pull(SAMPID)
      valid_id <- intersect(gtex_id, colnames(gtex_data))  #some samples are not use to rna-seq
      if (length(valid_id) > 3) {
        gtex_normal <- gtex_data[, valid_id]
        if (is_voom) {
          gtex_normal <- data_pre(df=gtex_normal, is_filt=is_filt, cut=cut)
          #voom transformation of the data (log)
          gtex_normal <- voom(gtex_normal)
          #taking log transformed data for exploration of batch effects
          gtex_normal <- gtex_normal$E #初始化为voom转换矩阵，确保后续代码可以继续运行
        } else {
          gtex_normal <- log2(gtex_normal + 1)
        }
        
        print(dim(gtex_normal))
        # retrieve the genes in common between GEO and TCGA-LUAD datasets
        gtex_normal <- gtex_normal[rownames(gtex_normal) %in% intersect(rownames(gtex_normal),rownames(c.dataFilt)),]
        c.dataFilt <- c.dataFilt[rownames(c.dataFilt) %in% intersect(rownames(c.dataFilt),rownames(gtex_normal)),]
        # merge the two counts matrices
        countsTable <- cbind(c.dataFilt,gtex_normal[match(rownames(c.dataFilt), rownames(gtex_normal)),])

        #create dataframe with batch information
        AnnotationCounts <- matrix(0,ncol(countsTable),3)
        colnames(AnnotationCounts) <- c("Samples","Batch","Conditon")
        rownames(AnnotationCounts) <- colnames(countsTable)
        AnnotationCounts <- as.data.frame(AnnotationCounts)
        AnnotationCounts$Samples <- colnames(countsTable)
        AnnotationCounts[colnames(c.dataFilt),"Batch"] <- "TCGA"
        AnnotationCounts[colnames(gtex_normal),"Batch"] <- "GTEX"

        if (length(dataSmNT_short) > 1){
          AnnotationCounts[dataSmNT_short,"Condition"] <- "Normal"
        }
        AnnotationCounts[dataSmTP_short,"Condition"] <- "Tumor"
        AnnotationCounts[colnames(gtex_normal),"Condition"] <- "Normal"

        #AnnotationCounts需要三列，Condition（协变量，可以理解成生物学条件）Batch（批次，这里主要指平台，tcga和gtex）Samples（样本id）
        #############################
        ###tgex
        #############################
        tryCatch({
          c.dataFilt <- TCGAbatch_Correction(tabDF = countsTable,
                                                UnpublishedData = TRUE, 
                                                AnnotationDF = AnnotationCounts)
          # c.dataFilt <- TCGAbatch_Correction(tabDF = v.dataFilt, batch.factor = "Plate", adjustment = "TSS", is_plot = FALSE)
          }, error = function(e) {
          #Catch and skip
          message("TCGAbatch_Correction error", e$message)
          message("Skip!")
        })
        print("Continue!")
        ############################################################################
        dataSmNT_short <- c(colnames(gtex_normal), dataSmNT_short) #Add 20250331
        ############################################################################
      }
      
      ##
      if (length(valid_id) > 3){
        ###
        if (candidate %in% rownames(c.dataFilt)){
          if (length(target) > 0){
            select_row = c(rownames(c.dataFilt)[1:2], target)
          } else {
            select_row = rownames(c.dataFilt)
          }
          DEG <- TCGAanalyze_DEA(
            #add c(rownames(c.dataFilt)[1:2], target) avoid wrong of data format
            mat1=c.dataFilt[select_row, dataSmNT_short], 
            mat2=c.dataFilt[select_row, dataSmTP_short],
            pipeline="limma",
            Cond1type = "Normal",
            Cond2type = "Tumor",
            method = "glmLRT")
          # 绘制分组箱线图并叠加点和标注
          pvalue <- DEG[candidate,]$P.Value
          plegend <- ifelse(pvalue < 0.0001, "****",
                            ifelse(pvalue< 0.001, "***",
                                   ifelse(pvalue < 0.01, "**",
                                          ifelse(pvalue < 0.05, "*", "ns"))))
          
          expression <- c.dataFilt[candidate,c(dataSmTP_short,dataSmNT_short)]
          exp_data <- data.frame(expression); exp_data$sample <- rownames(exp_data)
          exp_data <- exp_data %>%
                        mutate(group = ifelse(sample %in% dataSmTP_short, "tumor",
                                              ifelse(sample %in% dataSmNT_short, "normal", NA)))
    
          pp <- ggplot(exp_data, aes(x = group, y = expression, fill = group)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # 不显示异常值，使叠加更清晰
          geom_jitter(width = 0.2, aes(color = group), size = 2) +  # 叠加点，增加抖动防止重叠
          labs(
            title = paste("Gene Expression of Cancer vs Normal Samples with GTEx in", p, "\n", "\n", "\n", plegend),
            x = "Sample Group",
            y = ifelse(is_voom, "Corrected Voom-transform Value", "log2(TPM + 1)")
          ) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14))  # 隐藏图例（可选）
          ggsave(file.path("tmp", p, paste0(p, "_gtex_Gene_Expression_Boxplot.pdf")), plot = pp, width = 8, height = 6, dpi = 600)
        } else {
          if (length(target) > 1){
              select_row = c(rownames(c.dataFilt)[1:2], target)
          } else {
            select_row = rownames(c.dataFilt)
          }
          DEG <- TCGAanalyze_DEA(
            mat1=c.dataFilt[select_row, dataSmNT_short], 
            mat2=c.dataFilt[select_row, dataSmTP_short],
            pipeline="limma",
            Cond1type = "Normal",
            Cond2type = "Tumor",
            method = "glmLRT")
        }
      
        fwrite(as_tidytable(DEG, .keep_rownames = "gene_name"), file.path("tmp", p, paste0(p,"_gtex_deg.csv")))
        ##
      } else {
        print(paste(p, "It doesn't have enough normal samples!"))
      }
      tmp_mat <- as_tidytable(c.dataFilt, .keep_rownames = "gene_name")
      tmp_mat <- tmp_mat %>% filter(gene_name %in% tmp_mat$gene_name[grep("FAM135", tmp_mat$gene_name)]) #Add to avoid big output
      fwrite(tmp_mat, file.path("tmp", p, paste0(p,"_gtex_exp.csv")))
    
      results[[paste0(p, "_gtex_exp")]] <- c.dataFilt
      results[[paste0(p, "_gtex_deg")]] <- tmp_mat
    
    } else {
      print(paste(p, "is OK!"))
    }
        
  } #for end
  return(results)
}
#test
#linshi <- mat2plot(data_dir = "../GDCdata/", is_short = TRUE)
#project=c("TCGA-LUSC"); data_dir="../GDCdata"; num_tp=100; num_nt=100;tp_t="TP"; tp_n="NT"; 
#is_short=TRUE; save=TRUE; target=c("FAM135B"); candidate="FAM135B";p=project[1]
#voom=TRUE; is_log=FALSE; norm_method="none";prior.count=0; unit="tpm"
