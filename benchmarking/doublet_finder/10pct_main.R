
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
ACT_DBL<-0.1
################################
#                              #
#             Main             #
#                              #
################################

main<-function (root_dir, METHOD_ID, out_dir, n_cores, num_test){
  all_datasets_dirs<-list.dirs(path = root_dir, full.names = TRUE, recursive = F)
  out_df<-NULL


  for(expected_doublet_rate in doublet_rates_to_test){
      for(cur_dataset_dir in all_datasets_dirs){
        print(glue("INFO: Currently testing expected doublet rate:{expected_doublet_rate}, actual doublet rate:{ACT_DBL}"))
        cur_out_df<-run_doublet_finder_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test = num_test,
                                                  exp_dbl = expected_doublet_rate, actual_dbl = ACT_DBL)
        cur_out_df$dbl_exp <- expected_doublet_rate
        cur_out_df$dbl_act <- ACT_DBL
        out_df<-rbindlist(list(out_df, cur_out_df))|>
          data.frame()

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