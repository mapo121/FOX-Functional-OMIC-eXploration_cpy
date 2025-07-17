import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

# NOTE (JUL17)
# Although the output file is named GENIE3_linkList, we used GRNBoost2 (from Arboreto) to generate the network,
# following SCENIC’s standard input file naming convention.


n = 1
tf_path = "int/1.1_inputTFs.txt"
ex_path = "int/1.1_exprMatrix_filtered_t.txt"
    # ex_matrix is a DataFrame with gene names as column names
if __name__ == '__main__':
    ex_matrix = pd.read_csv(ex_path, sep='\t')
    # tf_names is read using a utility function included in Arboreto
    tf_names = load_tf_names(tf_path)
    for i in range(0, n):
        network = grnboost2(expression_data=ex_matrix, tf_names=tf_names)
        network.to_csv('1.4_GENIE3_linkList-' + str(i) + ".tsv", sep='\t', index=False, header=True)
