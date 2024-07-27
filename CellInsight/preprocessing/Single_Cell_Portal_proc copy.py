'''

Process data from the Single Cell Portal.

Input:
 - data_file: data from the Single Cell Portal in h5ad format.

Output:
 - processed_data: processed data in anndata format.

'''

import scanpy

def process_Single_Cell_Portal(data_file):
    processed_data = scanpy.read_h5ad(data_file)
    
    return processed_data