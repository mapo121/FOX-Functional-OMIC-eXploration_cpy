source('../scenic_process.R')
load('../pbmc_named.RData')#QA/QC corrected version ehe
# fseurat_scenicV2unction(obj_filename, cluster_on, session_name, cores, topWhat) {
load('m.RData')# its the preloaded human database.

scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=18)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

seurat_scenicV3(pbmc, "c_idents", "Aug27.RData", 18, 50)
