
---
![alt text](docs/image_cover.jpg)

Let me know if you need further customization or changes! Email: mapostol@unmc.edu or mapo3126@gmail.com

# F🦊X: **F**unctional OMIC  e**X**ploration of Gene Regulatory Networks

FOX is a highly **modular** and **flexible** methodology for analyzing and comparing gene-regulatory networks, especially in single-cell gene expression data. It integrates several advanced tools, including **SCENIC**, **NMF**, and **Kendall's Tau**, to provide deep insights into gene regulation. FOX can be used to visualize, compare, and analyze the structure and activity of gene regulatory networks under different conditions.

[📘 Website Documentation](https://howard-fox-lab.github.io/FOX-Functional-OMIC-eXploration/) • [🦊 Quick examples (From paper, PBMC3K)](https://howard-fox-lab.github.io/FOX-Functional-OMIC-eXploration/Test_group_final.html) 

# Data Citations:

- Quick examples data is taken from the Seurat website (PBMC3K):  
  [PBMC3K Dataset from 10x Genomics](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)

- Supplemental information data is from:  
  [Human and mouse single-nucleus transcriptomics reveal TREM2-dependent and TREM2-independent cellular responses in Alzheimer’s disease](https://www.nature.com/articles/s41591-019-0695-9)


To see examples of formatted data, see data/

<p align="center">
  <img src="https://em-content.zobj.net/source/twitter/376/fox_1f98a.png" width="100" alt="FOXREG logo"/>
</p>

## Installation (test)

	pip install FOXREG


## Usage

To run FOX, you'll need to prepare your data (such as RSS matrices and metadata) and pass it to the class. Here's an example of how to initialize and use FOX:

```python
	from FOXREG import ComparisonTree
	import pandas as pd
	import warnings
	warnings.filterwarnings("ignore")

        data = pd.read_csv("QA_QC_PBMC_rss_values_Feb3.csv") ## this would be one comparison (RSS)
        df_RAS = pd.read_csv("obj_AUC_metadata2_PBMC.csv") ## grab this from your SCENIC stuff, include ALL METADATA AUC AND cellLabels

	labels = [
    		"B",
    		"CD14+ Mono",
    		"NK",
    		"CD8 T",
    		"FCGR3A+ Mono",
    		"DC",
    		"Memory CD4 T"
	]

        # your new labels here is your "tissue" or "cell" column
        comparison = ComparisonTree("Naive CD4 T", df_RAS, "newLabels", data, labels, "Unnamed: 0", "3.5_AUCellThresholds_Info_PVMC_QA_QC.tsv")


        comparison.construct_tree() 
        p_vals = comparison.plotRSS_NMF("B", drawQuadrants=True, include_pvals=True)
        comparison.plot_3dEmbedding(rawRSS=False)
        comparison.analyze_factors("B", percentages=True)
        comparison.compareLayers("B", "Naive CD4 T", 0.055)
        tr = comparison.create_global_tree()
```


### Example Workflow:
1. **Prepare your single-cell gene expression data** (e.g., CSV format).
2. **Initialize FOX** with the necessary data, including control and treatment conditions.
3. **Compare gene-regulatory layers** across conditions using the `compareLayers` function.
4. **Visualize the network structure** using 2D and 3D plots.
5. **Assess the reproducibility** of the regulatory network using the global tree structure.
6. **Analyze factors and clusters** with advanced statistical methods and visualize the results.

## Contributions

Contributions are welcome! Feel free to fork the repository and submit pull requests for bug fixes, new features, or improvements. Help us improve FOX!

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
Emoji artwork © Twitter, used under CC-BY 4.0 via [Twemoji](https://twemoji.twitter.com/)

---

