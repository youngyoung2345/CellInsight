import numpy as np

bucket_name = ''
file_key = ''

# 특정 유전자 발현 여부를 확인하는 함수
def check_gene_presence(adata, gene_name):
    for gene in adata.var_names:
        if gene_name in gene:
            return True
        else:
            return False

# 특정 유전자가 발현된 세포 클러스터의 개수를 얻는 함수
def count_clusters_with_gene(adata, gene_name, cluster_key='leiden', threshold=0):
    count = 0
    for gene in adata.var_names:
        if gene_name in gene:
            # 유전자 발현 데이터 추출
            gene_expression = adata[:, gene].X
            
            # 필요한 경우 2차원 배열을 1차원으로 변환
            if gene_expression.ndim > 1:
                gene_expression = gene_expression.toarray().flatten()
        
            # 클러스터 정보 추출
            clusters = adata.obs[cluster_key]
            
            # threshold를 초과하는 발현 값을 가진 세포들 선택
            cells_with_gene = gene_expression > threshold

            # 해당 세포들이 속한 클러스터만 선택
            clusters_with_gene = clusters[cells_with_gene]

            # 고유 클러스터 수를 계산
            unique_clusters = np.unique(clusters_with_gene)
            
            count += len(unique_clusters)
            
    return count