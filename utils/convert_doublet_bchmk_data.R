# data_root<-"/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
data_root<-"/projects/b1042/GoyalLab/zzhang/test"
output_root<-"/projects/b1042/GoyalLab/zzhang/test/python_data"
locs<-c(glue("{data_root}/hm-12k.rds"),
        glue("{data_root}/J293t-dm.rds"),
        glue("{data_root}/pbmc-2ctrl-dm.rds"))
for(i in seq(1:length(locs))){
  loc <- locs[i]
  dtID<-gsub(".*/(.+).rds","\\1",loc)
  data<-readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)

  cellIDs<-colnames(count)
  geneIDs<-rownames(count)

  output_prefix<-glue("{output_root}/{dtID}/{dtID}_")
  write_lines(cellIDs,glue("{output_prefix}cell-ID.txt"))
  write_lines(geneIDs,glue("{output_prefix}genes.txt"))
  write_lines(label,glue("{output_prefix}doublets-labels.txt"))
  writeMM(count,glue("{output_prefix}data.mtx"))


}