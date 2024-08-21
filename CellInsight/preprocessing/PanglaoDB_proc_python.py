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
from rpy2.robjects import pandas2ri, conversion, default_converter


def process_PanglaoDB(data_path, additional_data):
    '''
    data_path : EX) cellinsight-bucket/PanglaoDB/SRA203368/PanglaoDB/SRA203368_SRS866906.sparse.RData
    additional_data : list consisting of SRA, SRS, SRR, Species, Tumor, Protocol, Instrument, Full-length mRNA-seq, Number of cells,
                      Number of exp. genes, Number of clusters, Tissue, Cell line (Y/N), Primary adult tissue (Y/N), and Target cell population
    
    '''
    with conversion.localconverter(default_converter):
        pandas2ri.activate()

        ro.r.source('./PanglaoDB_proc_R.R')
        process_panglaodb = ro.r['process_panglaodb']

        object = process_panglaodb(data_path)
    
    counts_matrix = np.array(object.rx2('counts_matrix'))
    col_data_df = pandas2ri.rpy2py(object.rx2('col_data_df'))
    row_data_df = pandas2ri.rpy2py(object.rx2('row_data_df'))
    
    # processed_data definition
    
    processed_data = anndata.AnnData(X=counts_matrix)
    processed_data.obs['cell_ids'] = col_data_df
    processed_data.var['gene_ids'] = row_data_df
    
    processed_data.obs['species__ontology_label'] = additional_data['Species']
    processed_data.obs['library_preparation_protocol__ontology_label'] = additional_data['Protocol']
    processed_data.obs['organ__ontology_label'] = additional_data['Tissue']
    processed_data.obs['tumor'] = additional_data['Tumor']
    processed_data.obs['cell_line'] = additional_data['Cell line']
    processed_data.obs['primary_adult_tissue'] = additional_data['Primary adult tissue']
    
    processed_data.uns['SRA'] = additional_data['SRA']
    processed_data.uns['SRS'] = additional_data['SRS']
    processed_data.uns['SRR'] = additional_data['SRR']
    processed_data.uns['instrument'] = additional_data['Instrument']
    
    return processed_data