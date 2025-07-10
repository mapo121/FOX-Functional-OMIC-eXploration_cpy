library(SCENIC)
load("pt_1_.RData")
exprMat_log <- log2(exprMat+1)
scenicOptions = readRDS("int/scenicOptions.Rds")
scenicOptions@settings$nCores <- 4
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
