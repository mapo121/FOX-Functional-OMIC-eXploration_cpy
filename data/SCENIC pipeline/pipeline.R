library(dplyr)
library(Seurat)
library(patchwork)


# taken from the seurat website

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

png("image_umnapverifcation.png", width = 500, height = 500)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
dev.off()

save(pbmc, file = "pbmc_new.RData")

## Note Jul 10 note - earlier pipeline was ran 2/3/2025, we ran a R_session that did this later

### labeling done in post R-session (around 2/3/2025), note is here just for reference
## same biomarkers found in PBMC3k tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial
## new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
##    "NK", "DC", "Platelet")
## names(new.cluster.ids) <- levels(pbmc)
## pbmc <- RenameIdents(pbmc, new.cluster.ids)
## DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
## subsetted object was done with pbmc_noPlatelet = subset(pbmc, idents=c("Platelet"), invert=TRUE)
