
import pandas as pd
import boto3
from botocore.client import Config
import io
import json

def fetch_markers():
    bucket_name = 'cellinsight-bucket'
    file_key = 'PanglaoDB/markers/PanglaoDB_markers_27_Mar_2020.tsv'

    s3_client = boto3.client(
        's3',
        aws_access_key_id='',  # AWS Access Key
        aws_secret_access_key='',  # AWS Secret Access Key
        endpoint_url='https://kr.object.ncloudstorage.com',
        region_name='kr-standard',
        config=Config(signature_version='s3v4')
    )

    # S3에서 파일 가져오기
    response = s3_client.get_object(Bucket=bucket_name, Key=file_key)
    content = response['Body'].read()

    # 데이터를 DataFrame으로 읽기
    markergene = pd.read_csv(io.BytesIO(content), delimiter='\t')

    return markergene
