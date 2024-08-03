import boto3
from botocore.client import Config
import pandas as pd
import io

import models

bucket_name = ''
file_key = ''

local_file_path = 'PanglaoDB_markers_27_Mar_2020.tsv'

# upload a local file to database
s3_client.upload_file(local_file_path, bucket_name, file_key)

# get a file from database
response = s3_client.get_object(Bucket=bucket_name, Key=file_key)
content = response['Body'].read()

# read file contents in memory as dataframe
markergene = pd.read_csv(io.BytesIO(content), delimiter='\t')

