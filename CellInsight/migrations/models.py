import boto3
from botocore.client import Config

import io
import anndata

from rpy2 import robjects
from rpy2.robjects.packages import importr

# boto3.resource provides low-level interface 

s3_resource = boto3.resource(
    's3',
    aws_access_key_id='your-access-key',  
    aws_secret_access_key='your-secret-key',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)

# boto3.client provide high-levee interface(object-oriented interface)

s3_client = boto3.client(
    's3',
    aws_access_key_id='your-access-key',  
    aws_secret_access_key='your-secret-key',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)

def read_RData(file_data):
    base = importr('base')
    
    with robjects.conversion.localconverter(robjects.default_converter + robjects.numpy2ri.converter):
        r_object = base.readRDS(io.BytesIO(file_data))
        return robjects.conversion.rpy2py(r_object)

def read_h5ad(file_data):
    with io.BytesIO(file_data) as f:
        data = anndata.read_h5ad(f)
        return data
    
def list_s3_objects(bucket_name, prefix):
    try:
        response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
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

def read_file_from_s3(bucket_name, prefix):
    response = list_s3_objects(bucket_name, prefix)
    
    if response == None:
        return None
        
    file_list = [] 
    
    for file in response['Contents']:
        file_key = file['Key']
        unprocessed_data = get_s3_objects(bucket_name, file_key)
        
        if file_key.endswith('.RData'):
            data_type = 'RData'
            processed_data = read_RData(unprocessed_data)
        elif file_key.endswith('.h5ad'):
            data_type = 'h5ad'
            processed_data = read_h5ad(unprocessed_data)
        else:
            # print for debugging
            # !!! Delete this print after completing this project !!!
            print(f"Unsupported file type for key: {file_key}") 
            
            continue
        
        file_list.append({'key' : file_key, 'data_type' : data_type,'data': processed_data})
    
    return file_list