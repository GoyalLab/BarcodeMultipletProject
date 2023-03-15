set.seed(2022)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(PRROC)
library(pbapply)
library(scDblFinder)
require(glue)
require(purrr)
require(data.table)
library(stringr)
require(parallel)

source("doublet_cell_utils.R")
options(error=traceback)
################################
#                              #
#       Input Parameters       #
#                              #
################################


root_dir <-"/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
root_bchmk_dataset<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
METHOD_ID<-"doublet_cell"
out_dir<-"/projects/p31666/zzhang/doublet-bchmk/final"
n_cores<-8
doublet_rates_to_test<-c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25)
num_test <- 3


################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores){

  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL

  all_bchmk_dataset<-list.files(path=root_bchmk_dataset, pattern="*.rds", recursive = T, full.names = T)
  for(cur_bchmk_dataset in all_bchmk_dataset){
    cur_out_df<-run_bchmk_original_paper(cur_bchmk_dataset, num_to_test = num_test)
    cur_out_df$dbl_exp <- "NA"
    cur_out_df$dbl_act <- "NA"
    out_df<-rbindlist(list(out_df, cur_out_df))|>
      data.frame()
  }
  # loop through doublet rates to test for
  for(actual_doublet_rate in doublet_rates_to_test){
    for(expected_doublet_rate in doublet_rates_to_test){
        out_df_ls<-mclapply(all_datasets_dirs, function(cur_dataset_dir){
          print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{actual_doublet_rate}"))
          cur_out_df<-run_doublet_cell_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
                                                  dbl_exp = expected_doublet_rate, dbl_act = actual_doublet_rate)
          cur_out_df$dbl_exp <- expected_doublet_rate
          cur_out_df$dbl_act <- actual_doublet_rate
          cur_out_df
        },mc.cores=n_cores)
        print(out_df_ls)
        cur_mc_out_df<-rbindlist(out_df_ls)|>
          data.frame()
        out_df<-rbindlist(list(out_df, cur_mc_out_df))|>
          data.frame()
      }

      # for(cur_dataset_dir in all_datasets_dirs){
      #   print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{actual_doublet_rate}"))
      #   cur_out_df<-run_doublet_cell_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
      #                                           dbl_exp = expected_doublet_rate, dbl_act = actual_doublet_rate)
      #
      #   cur_out_df$dbl_exp <- expected_doublet_rate
      #   cur_out_df$dbl_act <- actual_doublet_rate
      #   print(cur_out_df)
      #   out_df<-rbindlist(list(out_df, cur_out_df))|>
      #     data.frame()
      # }
    # }
  }


  all_detection_rates_file<-glue("{out_dir}/{METHOD_ID}/stats/all_detection_rates.tsv")
  dir.create(file.path(glue("{out_dir}/{METHOD_ID}/stats")), showWarnings = TRUE, recursive = TRUE)
  write.table(out_df, file = all_detection_rates_file, quote = FALSE, sep="\t",row.names = FALSE)
}

# runs only when script is run by itself
if (sys.nframe() == 0){
  main(root_dir, METHOD_ID, out_dir, n_cores)
}