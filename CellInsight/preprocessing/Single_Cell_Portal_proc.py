'''

Process data from the Single Cell Portal

Input:
 - data_file: data from the Single Cell Portal in h5ad format

Output:
 - processed_data: processed data in anndata format

'''

import scanpy
from migrations import models

def process_Single_Cell_Portal(file_path, data_type):
    object = models.get_s3_objects('cellinsight-bucket', file_path)

    match data_type:
        case '10x_h5':
            processed_data = scanpy.read_10x_h5(object)
        case 'csv':
            processed_data = scanpy.read_csv(object)
        case '10x_mtx':
            processed_data = scanpy.read_10x_mtx(object)
        case _:
            return False
    
    return processed_data