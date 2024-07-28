'''

Process data from the PanglaoDB

Input:
 - data_file: data from the PanglaoDB in RData format

Output:
 - processed_data: processed data in seurat object format

'''

import anndata
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


pandas2ri.activate()

ro.r.source('PanglaoDB_proc_R.R')
process_panglaodb = ro.r['process_panglaodb']

def process_PanglaoDB(data_path):
    object = process_panglaodb(data_path)
    
    counts_matrix = np.array(object.rx2('counts_matrix'))
    col_data_df = pandas2ri.rpy2py(object.rx2('col_data_df'))
    row_data_df = pandas2ri.rpy2py(object.rx2('row_data_df'))
    
    # processed_data definition
    
    processed_data = anndata.AnnData(X=counts_matrix)
    processed_data.obs = col_data_df
    processed_data.var = row_data_df
    
    return processed_data