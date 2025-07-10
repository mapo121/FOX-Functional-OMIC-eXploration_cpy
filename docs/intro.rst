Python Example: Running the Hypergeometric Test
================================================

FOX requires two inputs: The AUC matrix provided by AUCell and the RSS matrix calculated by SCENIC.  
The functions `returnThreshold_dictionary` and `rundifferential_test_AUC` explain how we test for the expression of regulons.

Hypergeometric testing by FOX (internals)
-----------------
First, we define the essential helper functions needed to conduct differential testing.

.. code-block:: python

    def returnThreshold_dictionary(threshold_file, regulonC = None):
        """
        Returns a dictionary of thresholds based on the input threshold file.  
        The file is processed depending on whether it's a 3.5 format or not.
        
        :param threshold_file: Path to the file containing threshold data.
        :param regulonC: The column/row of the regulon to use.
        :return: A dictionary mapping regulon names to threshold values.
        """
        thresholds, regulonNames, regulonThresholds = None, [], []
        
        # Handle case where threshold_file is in 3.5 format
        if "3.5_" in threshold_file:
            thresholds = pd.read_csv(threshold_file, sep="\t").T
            regulonNames = list(map(lambda x: x.split(" ")[0], thresholds.loc["regulon"].tolist()))
            regulonThresholds = thresholds.loc[regulonC]  # Row-based threshold extraction
        else:
            thresholds = pd.read_csv(threshold_file)
            regulonNames = list(map(lambda x: x.split(" ")[0], thresholds[regulonC]))
            regulonThreshold = thresholds[regulonC]  # Column-based threshold extraction
        
        threshold_dict = {}
        
        for i in range(len(regulonNames)):
            threshold_dict[regulonNames[i]] = regulonThresholds[i]
            i += 1
        
        return threshold_dict

    def rundifferential_test_AUC(df, metaClusterLabel, control, treatment, regulon_name, threshold_file, alpha = .05):
        """
        This function performs a hypergeometric test on the differential expression of regulons 
        between control and treatment groups.
        
        :param df: DataFrame containing the AUC values for all cells.
        :param metaClusterLabel: Label for the meta-cluster (e.g., cell type or sample group).
        :param control: Label for the control group.
        :param treatment: Label for the treatment group.
        :param regulon_name: The name of the regulon being tested.
        :param threshold_file: Path to the file containing thresholds for regulons.
        :param alpha: The significance level for testing.
        
        :return: Tuple containing the regulon name and p-values for both control and treatment groups.
        """
        df.columns = list(map(lambda x: x.split(" ")[0], df.columns))  # Clean up column names
        
        if threshold_file is None:
            calculate_mean = df[regulon_name].mean()
            df['is_successful_{}'.format(regulon_name)] = df[regulon_name] > df[regulon_name].mean()
        else:
            threshold = returnThreshold_dictionary(threshold_file, 'threshold')
            df['is_successful_{}'.format(regulon_name)] = df[regulon_name] > threshold[regulon_name]

        # Count successful control cells
        control_successful_cells = df[(df[metaClusterLabel] == control) & (df['is_successful_{}'.format(regulon_name)])].shape[0]
        
        # Count successful treatment cells
        treatment_successful_cells = df[(df[metaClusterLabel] == treatment) & (df['is_successful_{}'.format(regulon_name)])].shape[0]
        
        total_successful_cells = df[df['is_successful_{}'.format(regulon_name)]].shape[0]

        # Total population size
        total_population = df.shape[0]
        n_1 = df[df[metaClusterLabel] == control].shape[0]
        n_2 = df[df[metaClusterLabel] == treatment].shape[0]

        # Calculate p-values using the hypergeometric distribution
        p_value_1 = hypergeom.sf(control_successful_cells - 1, total_population, total_successful_cells, n_1)
        p_value_2 = hypergeom.sf(treatment_successful_cells - 1, total_population, total_successful_cells, n_2)

        return (regulon_name, round(p_value_1, 2), round(p_value_2, 2))

Preparation of data for FOX
---------------------------
To run FOX, we're assuming that you are using Seurat. We are using SCENIC, after you run it you should see a file named 
3.4_regulonAUC.Rds with a file named 3.5_AUCellThresholds_Info.tsv. You want to save both of those!

.. code-block:: r
    
    library(SCENIC)
    library(AUCell)

    regulonAUC = readRDS("3.4_regulonAUC.Rds")
    
    ## optional step... remove the extended "low-confidence"
    cells_AUC = regulonAUC[!grepl("extended", rownames(regulonAUC)), ]
    
    cells_AUC = getAUC(cells_AUC)
    seurat_obj@meta.data <- cbind(seurat_obj@meta.data, t(cells_AUC))

    write.csv(seurat_obj@meta.data, file = "RAS_AUC_values.csv")
    ## calculate then the RSS file! 
    
    cellInfo <- data.frame(seuratCluster=Idents(seurat_obj))
    rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])

    # if you removed the extended earlier, remove them down here
    rss = rss[!grepl("extended", rownames(rss)), ]

    write.csv(rss, "rss_values_.csv")



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
    df_RAS = pd.read_csv("RAS_AUC_values.csv")  # AUC metadata
    
    # Define labels for your comparison
    labels = data.columns[1:].tolist()

    # Initialize the ComparisonTree with your data
    comparison = ComparisonTree(
        "<Baseline >", df_RAS, "newLabels", data, labels, "Unnamed: 0", 
        "3.5_AUCellThresholds_Info.tsv"
    )

    # Build the tree and run analyses
    comparison.construct_tree() 
    p_vals = comparison.plotRSS_NMF("B", drawQuadrants=True, include_pvals=True)
    comparison.plot_3dEmbedding(rawRSS=False)
    comparison.analyze_factors("B", percentages=True)
    comparison.compareLayers("B", "Naive CD4 T", 0.055)
    tr = comparison.create_global_tree()
    tree, dict_save = tr

Explanation of Functions
------------------------
- **returnThreshold_dictionary**: This function parses the threshold file and returns a dictionary that maps regulon names to threshold values. It handles different file formats.
  
- **rundifferential_test_AUC**: This function runs a differential test on regulon activity between control and treatment groups, using the hypergeometric distribution to compute the statistical significance of observed differences.

---


