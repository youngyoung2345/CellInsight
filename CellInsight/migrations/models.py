import boto3
from botocore.client import Config

import io
import anndata

from rpy2 import robjects
from rpy2.robjects.packages import importr

bucket_name='cellinsight-bucket'

# boto3.resource provides low-level interface 

s3_resource = boto3.resource(
    's3',
    aws_access_key_id='',  
    aws_secret_access_key='',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)

# boto3.client provide high-levee interface(object-oriented interface)

s3_client = boto3.client(
    's3',
    aws_access_key_id='',  
    aws_secret_access_key='',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)

def list_s3_objects(bucket_name, prefix, delimiter=False):
    try:
        params = {
            'Bucket': bucket_name,
            'Prefix': prefix,
        }
        
        if delimiter:
            params['Delimiter'] = '/'

        response = s3_client.list_objects_v2(**params)
        return response
    except:
        return None

def get_s3_objects(bucket_name, file_key):
    try:
        response = s3_client.get_object(Bucket=bucket_name, Key=file_key)
        unprocessed_data = response['Body'].read()
        return unprocessed_data
    except:
        return None
    
def download_s3_objects(bucket_name, file_key, file_path):
    try:
        s3_client.download_file(bucket_name, file_key, file_path)
        return 
    except:
        return False

def read_RData(file_data):
    base = importr('base')
    
    with robjects.conversion.localconverter(robjects.default_converter + robjects.numpy2ri.converter):
        r_object = base.readRDS(io.BytesIO(file_data))
        return robjects.conversion.rpy2py(r_object)

def read_h5ad(file_data):
    with io.BytesIO(file_data) as f:
        data = anndata.read_h5ad(f)
        return data
    
def read_PanglaoDB_from_s3(bucket_name, prefix):
    response = list_s3_objects(bucket_name, prefix)
    
    if response == None:
        return None
    
    file_list = [] 
    
    for file in response['Contents']:
        file_key = file['Key']
        unprocessed_data = get_s3_objects(bucket_name, file_key)
        
        if file_key.endswith('.RData'):
            processed_data = read_RData(unprocessed_data)

            parts = file_key.split('/')
            information = {'study' : parts[1], 'file_name' : parts[2], 'data': processed_data}
                
            file_list.append(information)
    
    return file_list

def read_singlecellportal_from_s3(bucket_name, prefix):
    
    return file_list
