
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

    adata_list = []
    
    for item in file_list:
        # expression 파일 존재 확인
        if 'expression' in item['type']:
            adata = item['data']
            adata_list.append(adata)
            
            for item in file_list:
                # cluster 파일 존재 확인
                if 'cluster' in item['type']:
                    clusters = item['data']
                    
                    # 클러스터링 결과가 세포와 일치하는지 확인 (index가 동일한지 확인)
                    if clusters.index == adata.obs.index:
                        
                        # 클러스터링 결과 추가
                        if 'cluster' in clusters.columns:
                            adata.obs['cluster'] = clusters['cluster']
                            
                        # UMAP 또는 TSNE 좌표 추가
                        if 'X' in clusters.columns and 'Y' in clusters.columns:
                            adata.obsm['X_umap'] = clusters[['X', 'Y']].values
                            
                        return adata_list  # expression 파일로 만들고 cluster 정보 추가된 adata들 return
                else:
                    return adata_list  # expression 파일로 만든 adata들만 return
        else:
            raise ValueError("Expression file not found in the provided study.")

