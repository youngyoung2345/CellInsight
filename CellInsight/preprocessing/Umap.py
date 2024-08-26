import boto3
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from botocore.client import Config
import io
import os
import numpy as np
import preprocessing.PanglaoDB_proc_python as PanglaoDB_proc_python
import migrations.models as models
def fetch_s3_folder_list():
    bucket_name = 'cellinsight-bucket'
    s3_client = boto3.client(
        's3',
        aws_access_key_id='',  
        aws_secret_access_key='',
        endpoint_url='https://kr.object.ncloudstorage.com',
        region_name='kr-standard',
        config=Config(signature_version='s3v4')
    )

    response = s3_client.get_object(Bucket=bucket_name, Key='singlecellportal/fourth_name.csv')
    content = response['Body'].read()
    excel_data = pd.read_csv(io.BytesIO(content))
    
    # B열의 폴더 이름과 C열의 표시 이름 매핑
    name_mapping = dict(zip(excel_data['B'], excel_data['C']))


    # S3의 모든 폴더 목록 가져오기
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix='singlecellportal/', Delimiter='/')
    
    # 폴더 목록 추출
    folders = [prefix['Prefix'] for prefix in response.get('CommonPrefixes', [])]
    
    display_folders = []
    for folder in folders:
        folder_key = folder.split('/')[-2]  # 'singlecellportal/scp10/'에서 'scp10' 추출
        if folder_key in name_mapping:
            display_folders.append((folder, name_mapping[folder_key]))
        else:
            display_folders.append((folder, folder_key))  # 매핑이 없으면 폴더 이름 사용

    return display_folders

def fetch_cluster_files(folder_name):
    bucket_name = 'cellinsight-bucket'
    cluster_prefix = folder_name.rstrip('/') + '/cluster/'  # 클러스터 파일 경로

    s3_client = boto3.client(
        's3',
        aws_access_key_id='',  
        aws_secret_access_key='',  
        endpoint_url='https://kr.object.ncloudstorage.com',
        region_name='kr-standard',
        config=Config(signature_version='s3v4')
    )

    # 선택된 폴더의 cluster 디렉토리 내의 파일 목록 가져오기
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=cluster_prefix)
    cluster_files = [item['Key'] for item in response.get('Contents', []) if item['Key'].endswith(('.csv', '.tsv', '.txt'))]

    # cluster_prefix 콘솔에 출력
    print(f"Cluster Prefix: {cluster_prefix}")
 # cluster_prefix 및 cluster_files 콘솔에 출력
    print(f"Cluster Files: {cluster_files}")
    # 파일이 없을 경우 경고 출력
    if not cluster_files:
        print(f"No cluster files found in {cluster_prefix}")
    cluster_files2 = cluster_files[0]
    print(f"Cluster Files2: {cluster_files2}")
    return cluster_files2


def fetch_and_process_file(cluster_files2, delimiter):
    bucket_name = 'cellinsight-bucket'

    s3_client = boto3.client(
        's3',
        aws_access_key_id='',  
        aws_secret_access_key='',  
        endpoint_url='https://kr.object.ncloudstorage.com',
        region_name='kr-standard',
        config=Config(signature_version='s3v4')
    )
    print(f"Cluster Files2: {cluster_files2}")
    try:
        # S3에서 클러스터 파일 가져오기
        response = s3_client.get_object(Bucket=bucket_name, Key=cluster_files2)
        content = response['Body'].read()
    except s3_client.exceptions.NoSuchKey:
        print(f"Error: The specified key {cluster_files2} does not exist.")
        return None
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

    # 파일 확장자에 따라 구분자를 설정
    _, file_extension = os.path.splitext(cluster_files2)
    if file_extension == '.csv':
        delimiter = ','
    elif file_extension in ['.txt', '.tsv']:
        delimiter = '\t'
    else:
        raise ValueError("Unsupported file format.")

    # CSV/TSV/TXT 파일 읽기
    exp = pd.read_csv(io.BytesIO(content), delimiter=delimiter)
    exp = exp.iloc[1:].copy()  # 첫 번째 행 제거

    # 'X'와 'Y' 열을 float으로 변환
    exp['X'] = pd.to_numeric(exp['X'], errors='coerce')
    exp['Y'] = pd.to_numeric(exp['Y'], errors='coerce')

    # Scanpy의 AnnData 객체 생성
    adata = sc.AnnData(X=exp[['X', 'Y']].values)

    # 이웃 계산 (필수 단계)
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X')

    # Leiden 클러스터링 실행
    sc.tl.leiden(adata, resolution=0.5)

    # UMAP 실행
    sc.tl.umap(adata)

    # UMAP 좌표를 데이터프레임에 추가
    adata.obsm['X_umap'] = exp[['X', 'Y']].values

    figures_path = 'media/figures'
    if not os.path.exists(figures_path):
        os.makedirs(figures_path)
    
    
    umap_plot_path = os.path.join(figures_path, 'cluster_files2.png')

    # UMAP 시각화 및 저장
    sc.pl.umap(adata, color='leiden', show=False)
    plt.savefig(umap_plot_path)
    print(f"UMAP plot saved to {umap_plot_path}")
    return umap_plot_path



def fetch_s3_folder_list_Panglao():
    bucket_name='cellinsight-bucket'
    prefix = 'PanglaoDB/'
    response = models.list_s3_objects(bucket_name, prefix, delimiter=True)
    folders = [prefix['Prefix'] for prefix in response.get('CommonPrefixes', [])]
    return folders

def fetch_first_file_in_folder(folder_name):
    bucket_name = 'cellinsight-bucket'
    

    # 선택된 폴더 내 첫 번째 파일 가져오기
    response =  models.list_s3_objects(bucket_name, folder_name)
    files_in_folder = response.get('Contents', [])

    if not files_in_folder:
        print(f"No files found in {folder_name}")
        return None
    
    first_file_key = files_in_folder[0]['Key']
    print('first_file_key',first_file_key)
    return first_file_key


def fetch_and_process_file_Panglao(folder_name):
    # 선택된 폴더에서 첫 번째 파일 가져오기
    first_file_key = fetch_first_file_in_folder(folder_name)
    if not first_file_key:
        return None, None, None
    print('first_file_key2',first_file_key)
    bucket_name = 'cellinsight-bucket'
    figures_path = 'media/figures'
    local_file_path = os.path.join(figures_path, os.path.basename(first_file_key)).replace("\\", "/")
    print(local_file_path)
    # 파일 다운로드
    models.download_s3_objects(bucket_name, first_file_key, local_file_path)
    print(f"File downloaded to {local_file_path}")
    
    # 추가 데이터 가져오기
    try:
        additional_data = PanglaoDB_proc_python.import_additional_data(first_file_key)
    except ValueError as e:
        print(f"Error importing additional data: {str(e)}")
        return None, None, None

    # RData 파일을 읽고 UMAP 생성
    try:
        processed_data = PanglaoDB_proc_python.process_PanglaoDB(local_file_path, additional_data)
    except Exception as e:
        print(f"Error processing RData file: {str(e)}")
        return None, None, None
    # Debugging information
    print("Processed data observation data:")
    print(processed_data.obs.head())
    
    print("Processed data uns keys:")
    print(processed_data.uns.keys())
    # UMAP 생성
    print(additional_data)

    print("Processed data observation data:")

    # 제외할 열과 uns 키 목록
    columns_to_exclude = [
        'cell_ids', 
        'gene_ids',
        'species__ontology_label', 
        'library_preparation_protocol__ontology_label',
        'organ__ontology_label', 
        'tumor', 
        'cell_line', 
        'primary_adult_tissue'
    ]

    uns_keys_to_exclude = ['SRA', 'SRS', 'SRR', 'instrument']

    # 제외할 열을 제거
# 제외할 열을 제거
# 제외할 열을 제거 (존재할 경우에만 제거)
 

    # 확인을 위한 출력
    print("Remaining uns keys:", processed_data.uns.keys())
    print("Remaining obs columns:", processed_data.obs.columns)
    print("Current obs keys:")
    print(processed_data.obs.keys())
    pd.set_option('display.max_rows', 10)
    pd.set_option('display.max_columns', 10)

    # Debugging information
    print("Processed data observation data (first 100 rows):")
    print(processed_data.obs.head(10))
    
    print("Processed data variable data (first 100 rows):")
    print(processed_data.var.head(10))
    processed_data.layers['counts'] = processed_data.X.copy()
    sc.pp.normalize_total(processed_data, target_sum=1e4)
    sc.pp.log1p(processed_data)
    sc.tl.pca(processed_data, n_comps=50)
    processed_data.raw = processed_data
   

    sc.pp.neighbors(processed_data, n_pcs=50, n_neighbors=15, use_rep='X')

    sc.tl.leiden(processed_data, resolution=0.5)
    
    
    
    sc.tl.umap(processed_data)

   
    figures_path = 'media/figures'
    if not os.path.exists(figures_path):
        os.makedirs(figures_path)


    umap_plot_path = os.path.join(figures_path, 'PanglaoDB_Umap.png')
    sc.pl.umap(processed_data, color='leiden', show=False)
    plt.savefig(umap_plot_path)
    print(f"UMAP plot saved to {umap_plot_path}")
    
    # obs와 uns 데이터를 반환하여 프론트에서 설명으로 사용
    return umap_plot_path, processed_data.obs, processed_data.uns