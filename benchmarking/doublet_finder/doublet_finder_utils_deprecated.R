

load_fatemap_into_seurat<-function(data_dir,cell_labels_prefix){
  expression_matrix <- Read10X(data.dir = data_dir)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)

  singlets_file<-glue("{cell_labels_prefix}_singlets.txt")
  multiplets_file<-glue("{cell_labels_prefix}_multiplets.txt")

  singlets<-read.delim(singlets_file,header = F, sep = "\n")[,1]|>
    as.vector()|>
    paste("-1",sep = "")
  multiplets<-read.delim(multiplets_file,header = F, sep = "\n")[,1]|>
    as.vector()|>
    paste("-1",sep = "")
  all_cells<-c(singlets,multiplets)

  seurat_object <- seurat_object[,colnames(seurat_object) %in% all_cells]
  labels<-sapply(seurat_object@meta.data%>%rownames(), function(x){
    if(x %in% singlets){
      cur_label<-"Singlet"
    }else if (x %in% multiplets){
      cur_label<-"Multiplet"
    }else{
      cur_label<-"NA"
    }
  }
  )
  seurat_object<-AddMetaData(seurat_object,labels,col.name="ground_truth")
  # seurat_object <- subset(x = seurat_object, subset = (ground_truth == "Singlet" | ground_truth == "Multiplet"))
  print(table(labels))
  seurat_object
}


run_doublet_finder<-function(doublet.seurat, expected_doublet_rate){
  out<-tryCatch(
  {
    doublet.seurat <- NormalizeData(doublet.seurat)
    doublet.seurat <- ScaleData(doublet.seurat)
    doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
    doublet.seurat <- RunPCA(doublet.seurat)

    sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
    sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
    bcmvn.doublet <- find.pK(sweep.stats.doublet)
    pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
    nExp_poi <- round(expected_doublet_rate*nrow(doublet.seurat@meta.data))
    doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

    attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
    classification_label<-paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')

    score <- doublet.seurat@meta.data[[attribute]]
    fg <- score[doublet.seurat@meta.data$ground_truth=="Multiplet"]
    bg <- score[doublet.seurat@meta.data$ground_truth=="Singlet"]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
    roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc

    singlet_only<-doublet.seurat@meta.data[doublet.seurat@meta.data$ground_truth=="Singlet",]
    true_singlet_detection_rate<-sum(singlet_only[[classification_label]]=="Singlet")/dim(singlet_only)[1]

    multiplet_only<-doublet.seurat@meta.data[doublet.seurat@meta.data$ground_truth=="Multiplet",]
    singlet_in_multiplet_rate<-sum(multiplet_only[[classification_label]]=="Singlet")/dim(multiplet_only)[1]
    output<-list(
      "Expected_Doublet_Rate"=expected_doublet_rate,
      "True_Singlet_Detection_Rate"=true_singlet_detection_rate,
      "Singet_in_Potential_Doublet_Detection_Rate"=singlet_in_multiplet_rate,
      "roc"=roc$auc
    )

    per_cell_info<-doublet.seurat@meta.data|>
      select("ground_truth",dplyr::all_of(attribute), dplyr::all_of(classification_label))
    colnames(per_cell_info)<-c("ground_truth","prediction_score","prediction_label")
    list("predictions"=per_cell_info,
         "stats"=output)
  },
    error=function (cond){
      message(glue("{expected_doublet_rate} as expected doublet rate is not working!"))
      output<-list(
        "Expected_Doublet_Rate"=expected_doublet_rate,
        "True_Singlet_Detection_Rate"=NULL,
        "Singet_in_Potential_Doublet_Detection_Rate"=NULL,
        "roc"=NULL
      )
      per_cell_info<-doublet.seurat@meta.data|>
        select("ground_truth")|>
        data.frame()
      per_cell_info<-cbind(per_cell_info,prediction_score=rep("N/A",nrow(per_cell_info)))
      per_cell_info<-cbind(per_cell_info,prediction_label=rep("N/A",nrow(per_cell_info)))
      return(list("predictions"=per_cell_info,
                  "stats"=output))
    }
  )
  return (out)

}


run_bchmk<-function(seurat_object, estimated_multiplet_rates_to_test, data_sample_ID, output_dir){
  df_ls<-NULL
  dir.create(file.path(output_dir), showWarnings = FALSE, recursive = TRUE)
  for(cur_expected_doublet_rate in estimated_multiplet_rates_to_test){
    cur_res<-run_doublet_finder(seurat_object,cur_expected_doublet_rate)
    df_ls<-rbindlist(list(df_ls,cur_res$stats))
    cur_prediction_results<-cur_res$predictions

    out_file <-glue("{output_dir}/{data_sample_ID}_{cur_expected_doublet_rate}_doublets_predictions.tsv")
    cur_prediction_results <- cbind(rownames(cur_prediction_results), data.frame(cur_prediction_results, row.names=NULL))
    colnames(cur_prediction_results)[1] <- "cellID"
    write.table(cur_prediction_results, file = out_file, quote = FALSE, sep="\t",row.names = FALSE)

    # Not sure if it's necessary to return the stats for auroc curve
    # skip for now
    # stats_key = glue("{data_sample_ID}_{cur_expected_doublet_rate}")


  }
  all_stats<-data.frame(df_ls)
  all_stats<-cbind(ID=rep(data_sample_ID,nrow(all_stats)), all_stats)
  all_stats
}


run_doublet_finder_on_dataset<-function(data_dir, METHOD_ID, out_dir){
  print(glue("Working on dataset: {data_dir}"))
  output_dir <- glue("{out_dir}/doublet_output/{METHOD_ID}")
  dir.create(file.path(output_dir), showWarnings = FALSE, recursive = TRUE)
  dir_10X <- glue("{data_dir}/10X")
  dataset_id <- strsplit(data_dir,"/")[[1]][[length(strsplit(data_dir,"/")[[1]])]]
  sample_dirs<-list.dirs(path = dir_10X, full.names = TRUE, recursive = F)
  stats <- list()
  out_df<-NULL
  labels_prefix <- glue("{data_dir}/fatemapID/{dataset_id}")

  for (cur_sample_dir in sample_dirs){
    cur_seurat_object<-load_fatemap_into_seurat(cur_sample_dir, labels_prefix)

    # if there are less than 10 valid barcodes, there is something wrong with the data, check
    # if 10X barcodes are consistent with fatemapID file's barcodes
    if(ncol(cur_seurat_object)<50){
      next
    }
    cur_sample_id<-strsplit(cur_sample_dir,"/")[[1]][[length(strsplit(cur_sample_dir,"/")[[1]])]]
    cur_data_sample_id<-glue("{dataset_id}_{cur_sample_id}")
    estimated_multiplet_rates_to_test <- seq(0.05, 0.21, 0.01)

    cur_output_dir<-glue("{output_dir}/{cur_data_sample_id}")
    cur_df <- run_bchmk(cur_seurat_object, estimated_multiplet_rates_to_test, cur_data_sample_id, cur_output_dir)
    out_df<-rbindlist(list(out_df, cur_df))|>
      data.frame()
  }
  out_df<-out_df[order(out_df$Expected_Doublet_Rate),]
}
