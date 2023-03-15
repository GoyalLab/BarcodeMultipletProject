load_fatemap_into_seurat<-function(data_dir, cell_labels_prefix, doublets_pct=0.08){
  expression_matrix <- Read10X(data.dir = data_dir)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  singlets_file<-glue("{cell_labels_prefix}_singlets.txt")

  # key file for simulaitng doublets
  singlets_pair_for_doublet_simulation_file<-list.files(data_dir, pattern="corrected_singlet_pairs.csv")

  # extract the singlets
  singlets_pair_for_doublet_simulation_file<-glue("{data_dir}/{singlets_pair_for_doublet_simulation_file}")
  singlets<-read.delim(singlets_file,header = F, sep = "\n")[,1]|>
    as.vector()|>
    paste("-1",sep = "")

  # extract the singlet paris
  singlet_pair_df<-read.csv(singlets_pair_for_doublet_simulation_file, header = F)|>
    dplyr::mutate_each(funs = function(x) str_remove_all(x," "))
  # remove combinations that did not pass QC in seurat
  all_seurat_cells<-seurat_object@meta.data|>row.names()

  # take average of expression to produce doublets
  real_cells_pt1<-singlet_pair_df$V1
  real_cells_pt2<-singlet_pair_df$V2
  assay_for_simulated_doublets<-seurat_object@assays$RNA@counts[,c(real_cells_pt1, real_cells_pt2)]
  avg_assay<-(assay_for_simulated_doublets[, real_cells_pt1] + assay_for_simulated_doublets[, real_cells_pt2])/2


  # create artificial doublet ID
  singlets_in_seurat <- singlets[singlets %in% colnames(seurat_object)]
  doublet_id<-paste(singlet_pair_df$V1,singlet_pair_df$V2, sep = "--")
  colnames(avg_assay)<-doublet_id
  num_doublets<-round(length(singlets_in_seurat)*doublets_pct/(1-doublets_pct) )
  avg_assay<-avg_assay[,1:num_doublets]
  doublet_object<-CreateSeuratObject(counts = avg_assay)

  # tag doublet labels back to the doublet object
  label<-rep("doublet",ncol(doublet_object))|>
    as.data.frame()
  rownames(label)<-colnames(doublet_object)
  doublet_object<-AddMetaData(doublet_object, label, col.name = "label")

  # merge the two objects
  singlet_object<-seurat_object[, singlets_in_seurat]
  label<-rep("singlet",ncol(singlet_object))|>
    as.data.frame()
  rownames(label)<-colnames(singlet_object)
  singlet_object<-AddMetaData(singlet_object, label, col.name = "label")
  doublet.seurat<-merge(singlet_object, y = doublet_object, project = "combined_doublet_singlet")
  doublet.seurat
}

run_doublet_finder<-function(doublet.seurat, exp_dbl){
  doublet.seurat <- NormalizeData(doublet.seurat)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)

  sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE, num.cores = 8)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
  # nExp_poi <- sum(doublet.seurat[["label"]] == "doublet")
  nExp_poi<-round(dim(doublet.seurat)[2]*exp_dbl)
  doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)

  attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
  classification_label<-paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')

  score <- doublet.seurat@meta.data[[attribute]]
  fg <- score[doublet.seurat@meta.data$label=="doublet"]
  bg <- score[doublet.seurat@meta.data$label=="singlet"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_pr<-pr$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  cur_roc<-roc$auc
  output<-list(
    "roc"=cur_roc,
    "pr"=cur_pr
  )
  output
}
run_bchmk<-function(data_dir, prefix, data_sample_ID, num_to_test, exp_dbl, act_dbl){
  df_ls<-NULL
  for(idx in 1:num_to_test){
    seurat_object<-load_fatemap_into_seurat(data_dir, prefix, doublets_pct=act_dbl)
    cur_res<-run_doublet_finder(seurat_object, exp_dbl)
    df_ls<-rbindlist(list(df_ls,cur_res))
  }
  all_stats<-data.frame(df_ls)
  all_stats<-cbind(ID=rep(data_sample_ID,nrow(all_stats)), all_stats)
  all_stats
}

run_doublet_finder_on_dataset<-function(data_dir, METHOD_ID, exp_dbl, actual_dbl, num_to_test=5){
  print(glue("Working on dataset: {data_dir}"))
  dir_10X <- glue("{data_dir}/10X")
  dataset_id <- strsplit(data_dir,"/")[[1]][[length(strsplit(data_dir,"/")[[1]])]]
  sample_dirs<-list.dirs(path = dir_10X, full.names = TRUE, recursive = F)
  stats <- list()
  out_df<-NULL
  labels_prefix <- glue("{data_dir}/fatemapID/{dataset_id}")
  for (cur_sample_dir in sample_dirs){
    # cur_seurat_object<-load_fatemap_into_seurat(cur_sample_dir, labels_prefix, doublets_pct=actual_dbl)
    cur_sample_id<-strsplit(cur_sample_dir,"/")[[1]][[length(strsplit(cur_sample_dir,"/")[[1]])]]
    cur_data_sample_id<-glue("{dataset_id}_{cur_sample_id}")
    cur_df <- run_bchmk(cur_sample_dir, labels_prefix, cur_data_sample_id, num_to_test, exp_dbl = exp_dbl, act_dbl = actual_dbl)
    out_df<-rbindlist(list(out_df, cur_df))|>
      data.frame()
  }
  out_df
}