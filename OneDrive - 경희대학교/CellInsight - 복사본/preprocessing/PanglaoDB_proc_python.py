'''

Process data from the PanglaoDB

Input:
 - data_file: data from the PanglaoDB in RData format

Output:
 - processed_data: processed data in seurat object format

'''
import io
import os
import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
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
        SRA, SRS = data_name.split('_')

        return csv[(csv['SRA'] == SRA) & (csv['SRS'] == SRS)]

    else: #EX) data_name : SRA621638
        SRA = data_name

        return csv[csv['SRA'] == SRA]

def process_dgCMatrix(object):
    i = np.array(object.rx2('i'))
    p = np.array(object.rx2('p'))
    x = np.array(object.rx2('x'))
    dims = np.array(object.rx2('Dim'))

    dimnames = object.rx2('Dimnames')
    rownames = [str(name) for name in dimnames[0]]
    colnames = [str(name) for name in dimnames[1]]

    sparse_matrix = sp.csc_matrix((x, i, p), shape=(dims[0], dims[1]))

    processed_data = anndata.AnnData(X=sparse_matrix)
    processed_data.obs_names = rownames
    processed_data.var_names = colnames

    return processed_data

def process_SingleCellExperiment(object):
    counts_matrix = np.array(object.rx2('counts_matrix'))
    col_data_df = pandas2ri.rpy2py(object.rx2('col_data_df'))
    row_data_df = pandas2ri.rpy2py(object.rx2('row_data_df'))
    
    processed_data = anndata.AnnData(X=counts_matrix)

    processed_data.obs['cell_ids'] = col_data_df
    processed_data.var['gene_ids'] = row_data_df

    return processed_data

def add_information(processed_data, additional_data):
    processed_data.uns['species__ontology_label'] = additional_data['Species']
    processed_data.uns['library_preparation_protocol__ontology_label'] = additional_data['Protocol']
    processed_data.uns['organ__ontology_label'] = additional_data['Tissue']
    processed_data.uns['tumor'] = additional_data['Tumor']
    processed_data.uns['cell_line'] = additional_data['Cell line (Y/N)']
    processed_data.uns['primary_adult_tissue'] = additional_data['Primary adult tissue (Y/N)']
    
    processed_data.uns['SRA'] = additional_data['SRA']
    processed_data.uns['SRS'] = additional_data['SRS']
    processed_data.uns['SRR'] = additional_data['SRR']
    processed_data.uns['instrument'] = additional_data['Instrument']
    
    return processed_data

def process_PanglaoDB(data_path, additional_data):
    '''
    data_path : EX) cellinsight-bucket/PanglaoDB/SRA203368/PanglaoDB/SRA203368_SRS866906.sparse.RData
    additional_data : list consisting of SRA, SRS, SRR, Species, Tumor, Protocol, Instrument, Full-length mRNA-seq, Number of cells,
                      Number of exp. genes, Number of clusters, Tissue, Cell line (Y/N), Primary adult tissue (Y/N), and Target cell population
    
    '''

    with conversion.localconverter(default_converter):
        pandas2ri.activate()

        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script_path = os.path.join(script_dir, 'PanglaoDB_proc_R.R')

        ro.r.source(r_script_path)
        process_panglaodb = ro.r['process_panglaodb']

        object = process_panglaodb(data_path)
    
    if 'dgCMatrix' in object.rclass:
        processed_data = process_dgCMatrix(object)
    else:
        processed_data = process_SingleCellExperiment(object)
    
    processed_data = add_information(processed_data, additional_data)
    
    return processed_data

def make_anndata(file_path):
    additional_data = import_additional_data(file_path)
    preprocessed_data = process_PanglaoDB(file_path, additional_data)

    return preprocessed_data

def save_preprocessed_data(preprocessed_data):
    temp_path = os.path.join('temp', 'preprocessed_data.h5ad')
    preprocessed_data.write(temp_path)

    return temp_path
