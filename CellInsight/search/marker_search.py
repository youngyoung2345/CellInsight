
import pandas as pd
import boto3
from botocore.client import Config
import io
import json

from migrations import models

def fetch_markers():
    bucket_name = 'cellinsight-bucket'
    file_key = 'PanglaoDB/markers/PanglaoDB_markers_27_Mar_2020.tsv'

    # S3에서 파일 가져오기
    response = models.s3_client.get_object(Bucket=bucket_name, Key=file_key)
    content = response['Body'].read()

    # 데이터를 DataFrame으로 읽기
    markergene = pd.read_csv(io.BytesIO(content), delimiter='\t')

    return markergene