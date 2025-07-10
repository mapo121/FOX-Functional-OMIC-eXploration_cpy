library(SCENIC)
library(Seurat)
library(SeuratObject)

ResetSCENICAssay = function(obj) {
DefaultAssay(obj) = "RNA"
obj[['SCENIC']] = NULL
return(obj)
}

remove_e = function(X, cutoff) {
	p1 = strsplit(gsub("*g)" ,"", gsub("[(]", "", X)), " ")
	if (!grepl("extended", p1[[1]][1]) && as.numeric(p1[[1]][2]) > cutoff) {
		return(X)
	}
}

removeLowConfidence = function(regulonAUC_names, cutoff) {
	regulonAUC_names = unique(unlist(lapply(regulonAUC_names, remove_e, cutoff)))
	return(regulonAUC_names)
}

AddSCENICAssayV2 = function(obj, path_regulonAUC, cutoff) {
ls = readRDS(path_regulonAUC)
ls = S4ToList(ls)
regulonAUC = ls$assays$data$listData$AUC

saveThese = removeLowConfidence(rownames(regulonAUC),cutoff)
print(saveThese)
regulonAUC = regulonAUC[saveThese, ]
rownames(regulonAUC) = unlist(lapply(rownames(regulonAUC), function(X) strsplit(X, " ")[[1]][1] )) ## remove the space
obj[['SCENIC']] = CreateAssayObject(as.data.frame(regulonAUC))
return(obj)
}

AddSCENICAssay = function(obj, path_regulonAUC){ 
ls = readRDS(path_regulonAUC)
ls = S4ToList(ls)
regulonAUC = ls$assays$data$listData$AUC
rownames(regulonAUC) = unlist(lapply(rownames(regulonAUC), function(X) strsplit(X, " ")[[1]][1] ))
regulonAUC = as.data.frame(regulonAUC)
regulonAUC$Names = rownames(regulonAUC)
regulonAUC = regulonAUC[!grepl("extended", regulonAUC$Names), ]
regulonAUC$Names = NULL
obj[['SCENIC']] = CreateAssayObject(regulonAUC)
return(obj)
save(obj, file = 'Seurat_SCENIC_assayObject.RData')
}

calcFuncMod = function(obj, dir)  {
a = readRDS(dir)
# Add module scores to Seurat object
obj <- AddModuleScore(obj, features = a, name = "regulon")
# Calculate average regulon module scores per cell type
average_scores <- obj@meta.data %>%
    group_by(Condition) %>%
    summarise(across(starts_with("regulon"), mean))
write.csv(as.data.frame(average_scores), file = "module_scores.csv")
write.csv(names(a),file = 'regulon_N.csv')
return(obj)
}


CalculateUMAP = function(obj, dims, res, assay="SCENIC") {
DefaultAssay(obj) = assay
obj <- NormalizeData(object = obj)
obj <- FindVariableFeatures(object = obj)
obj <- ScaleData(object = obj)
obj <- RunPCA(object = obj)
obj <- FindNeighbors(object = obj, dims = 1:dims)
obj <- FindClusters(object = obj, resolution=res)
obj <- RunUMAP(object = obj, dims = 1:dims)
save(obj, file = paste0("Seurat_", assay ,"_assayObject_Updated", "res: ", res, "dims: ", dims, ".RData"))
return(obj)
}

CalculateOverlayClusters  = function(obj, SCENIC_clusters, RNA_clusters) {
library(dplyr)
s_t = obj@meta.data
by_regulonAUC_res0.1 = s_t %>% group_by_at(c(SCENIC_clusters, RNA_clusters))
return(by_regulonAUC_res0.1 %>% tally())
}


