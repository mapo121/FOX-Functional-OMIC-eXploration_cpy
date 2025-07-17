library(SCENIC)
#https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

# NOTE (JUL17): 
# see grnboostv2.py : I by default make it look for 1.4_GENIE3_linkList-0.tsv, this filename was unused, just for documentation
import_arb = function(filename) {
	GRNBoost_linkList <- importArboreto("1.4_GENIE3_linkList-0.tsv")
	colnames(GRNBoost_linkList) = c("TF", "Target", "weight")
	GRNBoost_linkList[[4]] = NULL
	saveRDS(GRNBoost_linkList, file = "int/1.4_GENIE3_linkList.Rds")
}

# NOTE (JUL17):
# The parameters 'cluster_on', 'session_name', 'cores', and 'topWhat' are unused in the function logic below.
# They were left intentionally as placeholders to help track metadata during development — such as run timestamps or labeling conventions.
# These parameters are not referenced internally and serve only as run-level context for the developer.

seurat_scenicV3 = function(obj, cluster_on, session_name, cores, topWhat) {
exprMat = as.matrix(obj[['RNA']]$counts)
save(obj, exprMat, file = "pt_1_.RData")
scenicOptions = readRDS("int/scenicOptions.Rds")

### Co-expression network
genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")
exprMat_filtered_log <- log2(exprMat_filtered+1)
save(exprMat_filtered_log, scenicOptions, exprMat_filtered, file = "pt_2.RData")
}
