R Code Example: Gene Set Analysis and Binarization (hdWGCNA integration)
=================================================

This section demonstrates how to perform gene set analysis using R, binarize clusters, and compute Jensen-Shannon Divergence (JSD) for comparing modules and clusters in gene expression data.

hub_df was generated from following the tutorial here: `hdWGCNA <https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html>`_. Scroll down to section called "Getting the module assignment table".


Step 1: Create Gene Sets for Each Module

First, we group the data by modules and create gene sets for each module.

.. code-block:: r

    # module is your regulon
    df_grouped = hub_df %>% group_by(module) %>% summarize(genes_in_module = list(gene_name), .groups = "drop")

    # Create gene sets for each module
    all.sets <- lapply(1:nrow(df_grouped), function(i) {
      GeneSet(df_grouped$genes_in_module[[i]], setName = as.character(df_grouped$module[i]))
    })

    # View the result
    print(all.sets)
    all.sets <- GeneSetCollection(all.sets)

Step 2: Prepare Expression Matrix
Next, we extract the expression matrix from a Seurat object, which contains gene expression data.

.. code-block:: r

    exprMatrx = seurat_obj[["RNA"]]$counts



Step 3: Calculate AUC and Assign Thresholds
Using the `AUCell` package, we calculate the Area Under the Curve (AUC) for the gene sets and explore the thresholds.

.. code-block:: r

    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)

Step 4: Binarize Clusters Based on Conditions
We binarize clusters in the Seurat object based on cell type or condition.

.. code-block:: r

    # Loop through each cell type to binarize clusters (remember, pick your column)
    binarize_clusters = lapply(as.character(unique(seurat_obj@meta.data$cell_type)), 
        function(i) {
            cluster = ifelse(seurat_obj@meta.data$cell_type == i, 1, 0)
            return(cluster)
    })
    names(binarize_clusters) = unique(seurat_obj@meta.data$cell_type)

Step 5: Binarize Regulons Based on AUC Thresholds
Now, we binarize the regulons (gene modules) based on their AUC values and the selected thresholds.

.. code-block:: r

    binarize_regulons = lapply(as.character(rownames(cells_AUC)), 
        function(i) {
            threshold = as.double(cells_assignment[[i]]$aucThr$selected)
            b_values <- ifelse(cells_AUC[i, ] > threshold, 1, 0)
            return(b_values)
    })
    names(binarize_regulons) = rownames(cells_AUC)

Step 6: Compute Jensen-Shannon Divergence (JSD)
Finally, we compute the Jensen-Shannon Divergence between the binarized clusters and regulons using the `philentropy <https://cran.r-project.org/web/packages/philentropy/index.html>`_ package.

.. code-block:: r

    library(philentropy)

    # Initialize a matrix to store the JSD values
    jsd_matrix <- matrix(0, nrow = length(binarize_clusters), ncol = length(binarize_regulons))

    # Set row and column names of the matrix
    rownames(jsd_matrix) <- names(binarize_clusters)
    colnames(jsd_matrix) <- names(binarize_regulons)

    # Loop through each combination of binarize_clusters and binarize_regulons
    for (i in names(binarize_clusters)) {
      for (i2 in names(binarize_regulons)) {
        t <- length(binarize_clusters[[i]])  # Length of the current vector
        x_count <- rbind(binarize_clusters[[i]] / t, as.double(binarize_regulons[[i2]]) / t)  # Normalize the vectors
        jsd_matrix[i, i2] <- 1 - sqrt(JSD(x_count, unit = "log2"))
      }
    }
   write.csv(jsd_matrix, file = 'RSS_scores.csv')

Step 7: Plug into FOX.

Take the cells_AUC and jsd_matrix and export them to .csv files. Then plug them into FOX with your desired cluster names.


 Conclusion
This analysis demonstrates the process of binarizing clusters and regulons, followed by computing Jensen-Shannon Divergence to compare gene sets. These steps are useful for understanding the relationships between gene modules and cell clusters, especially in single-cell RNA-seq data. This can be applied to a variety of different gene network libraries.


