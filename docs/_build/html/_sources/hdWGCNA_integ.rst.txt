R Code Example: Gene Set Analysis and Binarization (hdWGCNA integration)
========================================================================

This section demonstrates how to perform gene set analysis using R, binarize clusters, and compute Jensen-Shannon Divergence (JSD) for comparing modules and clusters in gene expression data.

hub_df was generated from following the tutorial here: `hdWGCNA <https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html>`_. Scroll down to section called "Getting the module assignment table".

You need to download AUCell, dplyr, and BiocParallel.

Step 1: Define `calcRSS` Function
---------------------------------

First, we define the `calcRSS` function to calculate the Regulon Specificity Score (RSS) for each cell type. This was taken straight from the code from SCENIC documentation.

.. code-block:: r

    calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
    {
      if(any(is.na(cellAnnotation))) stop("NAs in annotation")
      if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
      normAUC <- AUC/rowSums(AUC)
      if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
      # 
      ctapply <- lapply
      if(require('BiocParallel')) ctapply <- bplapply

      rss <- ctapply(cellTypes, function(thisType)
        sapply(rownames(normAUC), function(thisRegulon)
        {
          pRegulon <- normAUC[thisRegulon,]
          pCellType <- as.numeric(cellAnnotation==thisType)
          pCellType <- pCellType/sum(pCellType)
          .calcRSS.oneRegulon(pRegulon, pCellType)
        })
      )
      rss <- do.call(cbind, rss)
      colnames(rss) <- cellTypes
      return(rss)
    }

Step 2: Internal Functions for JSD and RSS Calculation
------------------------------------------------------

We define the internal functions for calculating Jensen-Shannon Divergence (JSD) and the Regulon Specificity Score (RSS).

.. code-block:: r

    .H <- function(pVect){
      pVect <- pVect[pVect>0] # /sum(pVect) ??
      - sum(pVect * log2(pVect))
    }

    # Jensen-Shannon Divergence (JSD)
    calcJSD <- function(pRegulon, pCellType)
    {
      (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
    }

    # Regulon specificity score (RSS)
    .calcRSS.oneRegulon <- function(pRegulon, pCellType)
    {
      jsd <- calcJSD(pRegulon, pCellType)
      1 - sqrt(jsd)
    }

Step 3: Create Gene Sets for Each Module
----------------------------------------

We now create gene sets for each module (regulon) from the hub data frame.

.. code-block:: r

    # module is your regulon (preprocessing!)
    df_grouped = hub_df %>% group_by(module) %>% summarize(genes_in_module = list(gene_name), .groups = "drop")

    # Create gene sets for each module
    all.sets <- lapply(1:nrow(df_grouped), function(i) {
      GeneSet(df_grouped$genes_in_module[[i]], setName = as.character(df_grouped$module[i]))
    })

    # View the result
    print(all.sets)
    all.sets <- GeneSetCollection(all.sets)
    exprMatrix = seurat_obj[["RNA"]]$counts

Step 4: Calculate AUC and Assign Thresholds
-------------------------------------------

We calculate the Area Under the Curve (AUC) for the gene sets and explore the thresholds using the `AUCell` package.

.. code-block:: r

    cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)

    cells_AUC <- AUCell_calcAUC(all.sets, cells_rankings)
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
    cellInfo = data.frame(seuratCluster=Idents(seurat_obj))
    rss = calcRSS(AUC=getAUC(cells_AUC), cellAnnotation=cellInfo[colnames(cells_AUC), "seuratCluster"])


Step 5: Plug into FOX.
------------------------------------------

Take the cells_AUC and rss and export them to .csv files. Then plug them into FOX with your desired cluster names.


.. code-block:: r

    InputForFOX = function(obj , cells_AUC, rss, cells_assignment) {
      a  = S4ToList(cells_AUC)
      a = a$assays$data$listData$AUC
      a = a[!grepl("extended", rownames(a)), ]
      a = t(a)
      obj@meta.data =cbind(obj@meta.data, a)
      write.csv(obj@meta.data, "RAS_matrix_foxInput.csv")

      rss = rss[!grepl("extended", rownames(rss)), ]
      write.csv(rss, file = "rss_FOX-input.csv")

      for (i in as.character(rownames(cells_AUC))) {
      	threshold = as.double(cells_assignment[[i]]$aucThr$selected)
  	nCells = cells_assignment[[i]]$aucThr$thresholds["minimumDens", "nCells"]
  	threshold_table <- rbind(threshold_table, data.frame(regulon = i, threshold = threshold, 	nCellsAssigned =nCells ))
      }

      # Export the threshold_table to a .tsv file
      write.table(threshold_table, file = "3.5_AUCell-thresholds.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


    }

    InputForFOX(seurat_obj, cells_AUC, rss, cells_assignment)



Conclusion
------------------------------------------

This analysis demonstrates the process of binarizing clusters and regulons, followed by computing Jensen-Shannon Divergence to compare gene sets. These steps are useful for understanding the relationships between gene modules and cell clusters, especially in single-cell RNA-seq data. This can be applied to a variety of different gene network libraries.


