import boto3
from botocore.client import Config

# boto3.resource provides low-level interface 
# boto3.client provide high-levee interface(object-oriented interface)

s3_resource = boto3.resource(
    's3',
    aws_access_key_id='your-access-key',  
    aws_secret_access_key='your-secret-key',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)



s3_client = boto3.client(
    's3',
    aws_access_key_id='your-access-key',  
    aws_secret_access_key='your-secret-key',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)


def read_file_from_s3(bucket_name, database_name, case_name):
    
    
    
    
    return file


