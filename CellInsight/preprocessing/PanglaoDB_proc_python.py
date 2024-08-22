'''

Process data from the PanglaoDB

Input:
 - data_file: data from the PanglaoDB in RData format

Output:
 - processed_data: processed data in seurat object format

'''
import io
import anndata
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, conversion, default_converter

from migrations import models

def import_additional_data(data_path):
    '''
    data_path : EX) PanglaoDB/SRA621638/SRA621638_SRS2610285.sparse.RData
    '''
    
    parts = data_path.split('/')
    csv_path = f'PanglaoDB/{parts[2]}/{parts[2]}_data.csv'
    unprocessed_data = models.get_s3_objects('cellinsight-bucket', csv_path)
    csv = pd.read_csv(io.StringIO(unprocessed_data.decode('utf-8')))

    #EX) SRA621638_SRS2610285.sparse.RData -> SRA621638_SRS2610285
    
    data_name = (parts[2].split('.'))[0] 
    
    if len(data_name) != 9: #EX) data_name : SRA621638_SRS2610285
        SSRA, SRS = data_name.split('_')

        return csv[(csv['SRA'] == SRA) & (csv['SRS'] == SRS)]

    else: #EX) data_name : SRA621638
        SRA = data_name

        return csv[csv['SRA'] == SRA]

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

