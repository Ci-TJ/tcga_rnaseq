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

mat2plot <- function(project=c("TCGA-LUSC"), data_dir="./GDCdata", num_tp=100, num_nt=100,tp_t="TP", tp_n="NT", 
                     is_short=FALSE, save=TRUE, target=c("FAM135B"), candidate="FAM135B",is_voom=TRUE, is_log=FALSE, norm_method="none",
                     prior.count=0, unit="tpm", cut=0.25, is_filt=TRUE){
  if (file.exists("tmp") == FALSE){
    dir.create("tmp")
  }
  results <- list()
  #
  for (p in project){
    if (file.exists(file.path("tmp", p)) == FALSE){
      dir.create(file.path("tmp", p), recursive = TRUE)
    }
    #######
    if (file.exists(file.path("tmp", p, paste0(p,"_exp.csv"))) == FALSE){
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
      
      dataPrep1 <- GDCprepare(query = queryDown, directory = data_dir, save = save, save.filename = file.path("tmp", p, paste0(p,".rda")))
      
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
        #taking log transformed data for exploration of batch effects
        #c.dataFilt <- TCGAbatch_Correction(tabDF = v.dataFilt, batch.factor="Plate", adjustment=c("TSS"), is_plot=FALSE)
        c.dataFilt <- v.dataFilt #初始化为voom转换矩阵，确保后续代码可以继续运行
      }
      
      tryCatch({
        c.dataFilt <- TCGAbatch_Correction(tabDF = v.dataFilt, batch.factor = "Plate", adjustment = "TSS", is_plot = FALSE)
        }, error = function(e) {
        #Catch and skip
        message("TCGAbatch_Correction error", e$message)
        message("Skip!")
      })
      print("Continue!")
      ##
      if (length(dataSmNT) > 3 ){
        ###
        if (candidate %in% rownames(c.dataFilt)){
          #survival 
          clin <- GDCquery_clinic(p,"clinical")
          flup <- fread(file.path(data_dir,"clinic_data", paste0(p, "_follow_up.tsv"))) #Note the file path!
          texp <- c.dataFilt[candidate,dataSmTP_short]
          cutoff <- median(texp) #use the median, so 50% high + 50% low 
          group <- ifelse(texp > cutoff, "High", "Low")
          group <- as_tidytable(data.frame(group), .keep_rownames = "barcode")
          barcode <- as_tidytable(data.frame(colData(dataPrep1))) %>% select(barcode,patient) %>% filter(barcode %in% dataSmTP_short) %>% inner_join(group,by="barcode")
          colnames(barcode)[2] <- "submitter_id"
          df_clin <- barcode %>% inner_join(clin,by="submitter_id")
          ##############
          for (i in 1:nrow(df_clin)){
            id <- df_clin$submitter_id[i]
            if (id %in% flup$cases.submitter_id){
              #print(paste(id,"in followup!",'\n'))
              fdf <- flup %>% filter(cases.submitter_id == id)
              days_followup <- as.numeric(fdf$follow_ups.days_to_follow_up)
              #print(days_followup)
              days_followup <- days_followup[!is.na(days_followup)]
              #print(days_followup)
              if (length(days_followup) >1){
                last_day_followup <- max(days_followup)
                df_clin$days_to_last_follow_up[i] <- last_day_followup
              }
            }
            else {
              print(paste(id, "is not in follow table!", "\n"))
            }
          }
          ##############
          TCGAanalyze_survival(data.frame(df_clin), clusterCol="group", legend=candidate, main = paste("Kaplan-Meier Overall Survival Curves of", p), 
                               filename = file.path("tmp", p, paste0(p, ".pdf")))
          if (length(target) > 0){
            select_row = c(rownames(c.dataFilt)[1:2], target)
          }
          else {
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
            title = paste("Gene Expression of Cancer vs Normal Samples in", p, "\n", "\n", "\n", plegend),
            x = "Sample Group",
            y = ifelse(is_voom, "Corrected Voom-transform Value", "log2(TPM + 1)")
          ) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14))  # 隐藏图例（可选）
          ggsave(file.path("tmp", p, paste0(p, "_Gene_Expression_Boxplot.pdf")), plot = pp, width = 8, height = 6, dpi = 600)
      
        } ###
        else {
          if (length(target) > 1){
              select_row = c(rownames(c.dataFilt)[1:2], target)}
          else {
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
      
        fwrite(as_tidytable(DEG, .keep_rownames = "gene_name"), file.path("tmp", p, paste0(p,"_deg.csv")))
        ##
      } else {
        print(paste(p, "It doesn't have enough normal samples!"))
      }
      
      tmp_mat <- as_tidytable(c.dataFilt, .keep_rownames = "gene_name")
      fwrite(tmp_mat, file.path("tmp", p, paste0(p,"_exp.csv")))
    
      results[[paste0(p, "_exp")]] <- c.dataFilt
      results[[paste0(p, "_deg")]] <- tmp_mat
    
    } else {
      print(paste(p, "is OK!"))
    }
    
    #for end     
  }
  return(results)
}
#test
#linshi <- mat2plot(data_dir = "../GDCdata/", is_short = TRUE)
#project=c("TCGA-LUSC"); data_dir="../GDCdata"; num_tp=100; num_nt=100;tp_t="TP"; tp_n="NT"; 
#is_short=TRUE; save=TRUE; target=c("FAM135B"); candidate="FAM135B";p=project[1]
#is_voom=TRUE; is_log=FALSE; norm_method="none";prior.count=0; unit="tpm"
    

                            
    
                                         
                                         
    
    
  











































      
