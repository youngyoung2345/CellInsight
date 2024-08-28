import io
import os

import anndata
import numpy as np
import pandas as pd

import scipy.sparse as sp
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, conversion, default_converter

from migrations import models

pandas2ri.activate()

script_dir = os.path.dirname(os.path.abspath(__file__))
r_script_path = os.path.join(script_dir, 'PanglaoDB_proc_R.R')

uns_keys = {
    'species__ontology_label': 'Species',
    'library_preparation_protocol__ontology_label': 'Protocol',
    'organ__ontology_label': 'Tissue',
    'tumor': 'Tumor',
    'cell_line': 'Cell line',
    'primary_adult_tissue': 'Primary adult tissue',
    'Number of clusters': 'Number of clusters',
    'SRA': 'SRA',
    'SRS': 'SRS',
    'SRR': 'SRR',
    'instrument': 'Instrument'
}

def get_first_value(series):
        return series.values[0] if not pd.isna(series).all() else pd.NA 

def set_uns_value(processed_data, key, data, column_name):
    processed_data.uns[key] = get_first_value(data.get(column_name, pd.Series([pd.NA])))

def import_additional_data(data_path):
    '''
    parameter

    data_path : EX) PanglaoDB/SRA621638/SRA621638_SRS2610285.sparse.RData

    - load additional data of a file of PanglaoDB
    
    '''
    
    # Extract study name
    # EX : PanglaoDB/SRA621638/SRA621638_SRS2610285.sparse.RData -> SRA621638

    parts = data_path.split('/')
    sra_id = parts[1]
    csv_path = f'PanglaoDB/{sra_id}/{sra_id}_data.csv'  

    unprocessed_data = models.get_s3_objects(models.bucket_name, csv_path)
    
    if unprocessed_data is None:
        return None
    
    csv = pd.read_csv(io.StringIO(unprocessed_data.decode('utf-8')))

    # Extract file name 
    # EX : SRA621638_SRS2610285.sparse.RData -> SRA621638_SRS2610285

    SRA, SRS = None, None
    file_name = (parts[2].split('.'))[0] 
    
    if '_' in file_name: 
        # EX : SRA621638_SRS2610285
        SRA, SRS = file_name.split('_', 1)      
        return csv[(csv['SRA'] == SRA) & (csv['SRS'] == SRS)]

    else:
        # EX : SRA621638
        SRA = file_name
        return csv[csv['SRA'] == SRA]
    
def process_PanglaoDB(data_path, additional_data):
    '''
    parameter

    data_path : EX) cellinsight-bucket/PanglaoDB/SRA203368/PanglaoDB/SRA203368_SRS866906.sparse.RData
    additional_data : list consisting of SRA, SRS, SRR, Species, Tumor, Protocol, Instrument, Full-length mRNA-seq, Number of cells,
                      Number of exp. genes, Number of clusters, Tissue, Cell line (Y/N), Primary adult tissue (Y/N), and Target cell population
    
    '''

    with conversion.localconverter(default_converter):
        pandas2ri.activate()

        ro.r.source(r_script_path)
        process_panglaodb = ro.r['process_panglaodb']
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"RData file not found: {data_path}")

        object = process_panglaodb(data_path)

        if 'dgCMatrix' in object.rclass:
            # dgCMatrix 데이터 처리
            i = np.array(object.rx2('i'))
            p = np.array(object.rx2('p'))
            x = np.array(object.rx2('x'))
            dims = np.array(object.rx2('Dim'))
            dimnames = object.rx2('Dimnames')
            rownames = [str(name) for name in dimnames[0]]
            colnames = [str(name) for name in dimnames[1]]

            # 희소 행렬 복원
            sparse_matrix = sp.csc_matrix((x, i, p), shape=(dims[0], dims[1]))

            # AnnData 객체 생성
            processed_data = anndata.AnnData(X=sparse_matrix)
            processed_data.obs_names = rownames
            processed_data.var_names = colnames
        else:
    
            counts_matrix = np.array(object.rx2('counts_matrix'))
            col_data_df = pandas2ri.rpy2py(object.rx2('col_data_df'))
            row_data_df = pandas2ri.rpy2py(object.rx2('row_data_df'))
        
            # processed_data definition
            processed_data = anndata.AnnData(X=counts_matrix)
            processed_data.obs['cell_ids'] = col_data_df
            processed_data.var['gene_ids'] = row_data_df
    
    # add additional data to processed data

    for key, value in uns_keys.items():
        set_uns_value(processed_data, key, additional_data, value)

    return processed_data








