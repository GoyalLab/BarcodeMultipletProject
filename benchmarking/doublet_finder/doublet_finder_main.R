#!/usr/bin/Rscript

require(DoubletFinder)
require(PRROC)
require(Matrix)
require(Seurat)
require(glue)
require(purrr)
require(data.table)
require(dplyr)
require(parallel)
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
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/output"
n_cores<-4
num_test <- 2
doublet_rates_to_test<-c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25)
################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores, num_test){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL


  for(actual_doublet_rate in doublet_rates_to_test){
    for(expected_doublet_rate in doublet_rates_to_test){

      # in case the socket creation error occurs
      attempts <- 1
      out_df_ls<-list()
      run_try <- TRUE
      while(run_try && attempts<=3){
        attempts <- attempts + 1
        try({
          out_df_ls<-mclapply(all_datasets_dirs, function(cur_dataset_dir){
            print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{actual_doublet_rate}"))
            cur_out_df<-run_doublet_finder_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
                                                      exp_dbl = expected_doublet_rate, actual_dbl = actual_doublet_rate)
            cur_out_df$dbl_exp <- expected_doublet_rate
            cur_out_df$dbl_act <- actual_doublet_rate
            cur_out_df
          },mc.cores=n_cores)
          saveRDS(out_df_ls, glue("{out_dir}/cur_out_df_ls.rds"))
          print(out_df_ls)
          # check if any is a try-error
          num_failed<-sapply(out_df_ls, function(x){
            inherits(x, "try-error")
          })|>sum()
          if(num_failed > 0){
            run_try <-TRUE
          }else{
            run_try <- FALSE
          }
        })
      }

      cur_mc_out_df<-rbindlist(out_df_ls)|>
        data.frame()
      out_df<-rbindlist(list(out_df, cur_mc_out_df))|>
        data.frame()
      # for(cur_dataset_dir in all_datasets_dirs){
      #   print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{actual_doublet_rate}"))
      #   cur_out_df<-run_doublet_finder_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
      #                                             exp_dbl = expected_doublet_rate, actual_dbl = actual_doublet_rate)
      #   cur_out_df$dbl_exp <- expected_doublet_rate
      #   cur_out_df$dbl_act <- actual_doublet_rate
      #   out_df<-rbindlist(list(out_df, cur_out_df))|>
      #     data.frame()
      #
      # }
    }
  }

  all_bchmk_dataset<-list.files(path=root_bchmk_dataset, pattern="*.rds", recursive = T, full.names = T)
  for(cur_bchmk_dataset in all_bchmk_dataset){
    cur_data<-readRDS(cur_bchmk_dataset)
    cur_data_id<-strsplit(cur_bchmk_dataset, "/")[[1]][[8]]
    count<-cur_data[[1]]
    label<-cur_data[[2]]
    doublet.seurat <- CreateSeuratObject(counts = count, min.cells = 1)
    doublet.seurat <- AddMetaData(doublet.seurat, col.name="label", label)
    cur_dbl_rate <- sum(doublet.seurat[["label"]] == "doublet")/ncol(doublet.seurat)
    print(glue("INFO: Current bechmarking paper dataset doublet rate is {cur_dbl_rate}!"))
    cur_out_df<-run_bchmk(doublet.seurat, cur_data_id, num_to_test = num_test,
                          exp_dbl = cur_dbl_rate)
    cur_out_df$dbl_exp <- "NA"
    cur_out_df$dbl_act <- "NA"
    out_df<-rbindlist(list(out_df, cur_out_df))|>
      data.frame()
  print("INFO: Done with original paper datasets")

  }

  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir, n_cores, num_test)
}