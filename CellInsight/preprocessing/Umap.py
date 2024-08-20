import io
import os

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

    # S3의 모든 폴더 목록 가져오기
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix='singlecellportal/', Delimiter='/')

    folders = [prefix['Prefix'] for prefix in response.get('CommonPrefixes', [])]
    
    return folders

def fetch_cluster_files(folder_name):
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