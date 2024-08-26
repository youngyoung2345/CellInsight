
'''
Process data from the Single Cell Portal

Input:
 - data_file: data from the Single Cell Portal in h5ad format

Output:
 - processed_data: processed data in anndata format

'''

import scanpy as sc
from migrations import models

def process_Single_Cell_Portal(file_path, data_type):
    # study 하나 통째로 받아오기 -> expression, cluster 따로따로
    file_list = models.read_singlecellportal_from_s3('cellinsight-bucket', file_path)

    # expression으로 anndata읽기
    adata = file_list[0]['data']
    
    # cluster 정보 추가하기
    clusters = file_list[1]['data']
    
    # 클러스터링 결과가 세포와 일치하는지 확인 (index가 동일한지 확인)
    if not all(clusters.index == adata.obs.index):
        raise ValueError("The index of cluster file and the index of adata.obs do not match.")

    # 클러스터링 결과를 adata의 obs에 추가합니다.
    adata.obs['cluster'] = clusters['cluster']

    return adata