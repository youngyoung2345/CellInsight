'''

Process data from the Single Cell Portal

Input:
 - data_file: data from the Single Cell Portal in h5ad format

Output:
 - processed_data: processed data in anndata format

'''

import scanpy

def process_Single_Cell_Portal(data_file):
    processed_data = scanpy.read_10x_h5(data_file)
    
    return processed_data
def process_Single_Cell_Portal_csv(data_file):
    processed_data = scanpy.read_csv(data_file, delimiter='\t')
    return processed_data

def process_Single_Cell_Portal_CellRanger(data_file):
    processed_data = scanpy.read_10x_mtx(data_file)
    return processed_data