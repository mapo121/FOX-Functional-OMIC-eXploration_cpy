library(SCENIC)
source("../scenic_process.R")
import_arb("1.4_GENIE3_linkList-0.tsv")
load("pt_1_.RData")
exprMat_log <- log2(exprMat+1)
scenicOptions = readRDS("int/scenicOptions.Rds")
scenicOptions@settings$nCores <- 4
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50perTarget"))
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)

