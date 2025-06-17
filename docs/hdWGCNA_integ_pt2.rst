AUCell & SCENIC Workflow
========================

This section demonstrates how to perform gene set analysis using R, binarize clusters, and compute Jensen-Shannon Divergence (JSD) for comparing modules and clusters in gene expression data.

hub_df was generated from following the tutorial here: `hdWGCNA <https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html>`_. Scroll down to section called "Getting the module assignment table".

You need to download AUCell, dplyr, SCENIC, and BiocParallel.

You must install these R packages:

.. code-block:: r

   library(AUCell)
   library(SCENIC)
   library(dplyr)

Prepare expression matrix from your Seurat object:

.. code-block:: r

   exprMatrix <- Seurat_obj[["RNA"]]$counts

Previous Steps in the Workflow
------------------------------

**Step 1: Build cell rankings**

.. code-block:: r

   cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)

**Step 2: Calculate AUC**

Use hdWGCNA output or define manually.

Calculate the hub_df from the vignette and plug in your module eigengenes here!

.. code-block:: r

   #  assuming you calculated the gene sets and have a copy of them... see GetHubGenes(seurat_obj, n_hubs = 10) from hdWGCNA
   #  shape of hub_df

   #   gene_name  module       kME
   # 1      Gene1 black    0.3711414
   # 2      Gene2 blue     0.3694937
   # 3      Gene3 blue     0.3318094
   # 4      Gene4 black    0.3304103
   # 5      Gene5 black    0.3312313

   # df_grouped = hub_df %>% group_by(module) %>% summarize(genes_in_module = list(gene_name), .groups = "drop")

    # Create gene sets for each module
    # geneSets <- lapply(1:nrow(df_grouped), function(i) {
    #    GeneSet(df_grouped$genes_in_module[[i]], setName = as.character(df_grouped$module[i]))
    # })

   ## this is how you upload them to AUCell... you can plug in each your gene sets here
   geneSets <- list(
     black = c("Gene1", "Gene5", "Gene4"),
     blue  = c("Gene3", "Gene2")
   )

   cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

**Step 3: Assign cells**

.. code-block:: r

   # Plot histograms and obtain thresholds
   set.seed(123)
   thresholds <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE)

   # Assign cells based on thresholds
   cells_assignment <- AUCell_exploreThresholds(
     cells_AUC,
     plotHist=FALSE,
     assignCells=TRUE
   )

   # Extract threshold values for each regulon
   Thresholds_forAUCell <- getThresholdSelected(cells_assignment)

Export Thresholds to .tsv
-------------------------

.. code-block:: r

   regulon_df <- data.frame(
     regulon = names(Thresholds_forAUCell),
     threshold = as.numeric(Thresholds_forAUCell),
     stringsAsFactors = FALSE
   )

   # Must include "3.5_" in the filename for SCENIC compatibility
   write.table(
     regulon_df,
     file = "3.5_regulon_scores_thresholds.tsv",
     sep = "\t",
     row.names = FALSE,
     quote = FALSE
   )

Get AUC and Generate RSS
------------------------

.. code-block:: r

   cells_test_RAS <- getAUC(cells_AUC)

   # Take labels from Seurat object
   cellInfo <- data.frame(seuratCluster = Idents(Seurat_obj))

   # Optional: Remove low-confidence regulons
   cells_AUC <- cells_AUC[!grepl("extended", rownames(cells_AUC)), ]

   # Calculate RSS
   rss <- calcRSS(
     AUC = getAUC(cells_AUC),
     cellAnnotation = cellInfo[colnames(cells_AUC), "seuratCluster"]
   )

   write.csv(rss, file = "rss_values_.csv")

   # Merge RAS with metadata
   pbmc_cpy@meta.data <- cbind(pbmc_cpy@meta.data, t(cells_test_RAS))
   write.csv(pbmc_cpy@meta.data, file = "RAS_values_dataset.csv")


Usage Example
-------------
To run FOX, you'll need to prepare your data (such as RSS matrices and metadata) and pass it to the class. Here's an example of how to initialize and use FOX:

.. code-block:: python

   from FOXREG import ComparisonTree
   import pandas as pd
   import warnings
   warnings.filterwarnings("ignore")

    # Read in the data
    data = pd.read_csv("rss_values_.csv")  # RSS values
    df_RAS = pd.read_csv("RAS_values_dataset.csv")  # AUC metadata
    
    # Define labels for your comparison
    other_clusters_to_compare = data.columns[1:].tolist()

    # Initialize the ComparisonTree with your data
    comparison = ComparisonTree(
        "<baseline cluster>", df_RAS, "newLabels", data, other_clusters_to_compare, "Unnamed: 0", 
        "3.5_regulon_scores_thresholds.tsv"
    )