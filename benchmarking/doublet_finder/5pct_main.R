
require(DoubletFinder)
require(PRROC)
require(Matrix)
require(Seurat)
require(glue)
require(purrr)
require(data.table)
require(dplyr)
require(stringr)
set.seed(2022)
source("doublet_finder_utils.R")

################################
#                              #
#       Input Parameters       #
#                              #
################################


root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
METHOD_ID<-"doublet_finder"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
num_test <- 3
doublet_rates_to_test<-c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25)
ACT_DBL<-0.05
################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores, num_test){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL

  print("INFO: Running on benchmark datasets!")

  all_bchmk_dataset<-list.files(path=root_bchmk_dataset, pattern="*.rds", recursive = T, full.names = T)
  for(cur_bchmk_dataset in all_bchmk_dataset){
    df_ls<-NULL
    for(idx in 1:num_test){
      cur_data<-readRDS(cur_bchmk_dataset)
      cur_data_id<-strsplit(cur_bchmk_dataset, "/")[[1]][[8]]
      count<-cur_data[[1]]
      label<-cur_data[[2]]
      doublet.seurat <- CreateSeuratObject(counts = count, min.cells = 1)
      doublet.seurat <- AddMetaData(doublet.seurat, col.name="label", label)
      cur_dbl_rate <- sum(doublet.seurat[["label"]] == "doublet")/ncol(doublet.seurat)
      cur_res<-run_doublet_finder(doublet.seurat, cur_dbl_rate)
      df_ls<-rbindlist(list(df_ls,cur_res))
    }
    all_stats<-data.frame(df_ls)
    cur_out_df<-cbind(ID=rep(cur_data_id,nrow(all_stats)), all_stats)
    cur_out_df$dbl_exp <- cur_dbl_rate
    cur_out_df$dbl_act <- "NA"
    out_df<-rbindlist(list(out_df, cur_out_df))|>
      data.frame()
    print("INFO: Done with original paper datasets")

  }

  for(expected_doublet_rate in doublet_rates_to_test){
    for(cur_dataset_dir in all_datasets_dirs){
      print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{ACT_DBL}"))
      cur_out_df<-run_doublet_finder_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
                                                exp_dbl = expected_doublet_rate, actual_dbl = ACT_DBL)
      cur_out_df$dbl_exp <- expected_doublet_rate
      cur_out_df$dbl_act <- ACT_DBL
      out_df<-rbindlist(list(out_df, cur_out_df))|>
        data.frame()
      # cur_out_temp_file <- glue("{out_dir}/{METHOD_ID}/tmp/exp_{expected_doublet_rate}_act_{ACT_DBL}.tsv")
      # write.table(out_df, file = cur_out_temp_file, quote = FALSE, sep="\t", row.names = FALSE)
    }
  }


  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates_{ACT_DBL}.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir, n_cores, num_test)
}