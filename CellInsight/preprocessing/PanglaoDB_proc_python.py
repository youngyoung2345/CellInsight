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
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd
import os
import scipy.sparse as sp
#from search import study_search
from migrations import models
from rpy2.robjects import pandas2ri, conversion, default_converter
pandas2ri.activate()

# 현재 스크립트의 디렉토리 경로 가져오기
script_dir = os.path.dirname(os.path.abspath(__file__))

# 절대 경로를 이용해 R 스크립트 경로 설정
r_script_path = os.path.join(script_dir, 'PanglaoDB_proc_R.R')
print(f"R script path: {r_script_path}")

# R 스크립트 로드
ro.r.source(r_script_path)
process_panglaodb = ro.r['process_panglaodb']

def import_additional_data(data_path):
    '''
    data_path : EX) PanglaoDB/SRA621638/SRA621638_SRS2610285.sparse.RData
    '''
    
     # 경로에서 필요한 부분만 추출
    parts = data_path.split('/')
    sra_id = parts[1]  # 예: SRA621638
    csv_path = f'PanglaoDB/{sra_id}/{sra_id}_data.csv'  # 올바른 경로로 변경
    
    print(f"Attempting to fetch data from S3 with path: {csv_path}")
    try:
        # S3에서 CSV 파일 데이터를 가져옴
        unprocessed_data = models.get_s3_objects('cellinsight-bucket', csv_path)
        if unprocessed_data is None:
            raise ValueError(f"Failed to fetch CSV data from S3: {csv_path}")
        
        # 데이터를 판다스 데이터프레임으로 변환
        csv = pd.read_csv(io.StringIO(unprocessed_data.decode('utf-8')))
        
        print(f"Data fetched successfully. Size: {len(unprocessed_data)} bytes")
        
    except Exception as e:
        print(f"Error fetching CSV data from S3: {str(e)}")
        return None  # 오류 발생 시 None 반환
    #EX) SRA621638_SRS2610285.sparse.RData -> SRA621638_SRS2610285
     # SRA 변수 초기화
    SRA = None
    SRS = None

    data_name = (parts[2].split('.'))[0] 
    
    if '_' in data_name:  # 언더스코어가 있는지 확인
        try:
            SRA, SRS = data_name.split('_', 1)  # 최대 두 부분으로 나누기
            return csv[(csv['SRA'] == SRA) & (csv['SRS'] == SRS)]
        except ValueError:
            print(f"Error splitting data_name: {data_name}")
            return None
    else:
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
    
    
        processed_data.obs['species__ontology_label'] = additional_data['Species']
        processed_data.obs['library_preparation_protocol__ontology_label'] = additional_data['Protocol']
        processed_data.obs['organ__ontology_label'] = additional_data['Tissue']
        processed_data.obs['tumor'] = additional_data.get('Tumor', pd.NA)
        processed_data.obs['cell_line'] = additional_data.get('Cell line', pd.NA)
        processed_data.obs['primary_adult_tissue'] = additional_data.get('Primary adult tissue', pd.NA)
        
        processed_data.uns['SRA'] = additional_data.get('SRA', pd.NA)
        processed_data.uns['SRS'] = additional_data.get('SRS', pd.NA)
        processed_data.uns['SRR'] = additional_data.get('SRR', pd.NA)
        processed_data.uns['instrument'] = additional_data.get('Instrument', pd.NA)
        
    return processed_data









