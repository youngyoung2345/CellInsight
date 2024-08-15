import migrations.models as models

def search_study(bucket_name, database_name, study_name):
    try:
        prefix = f"{database_name}/"
        
        response = models.s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
        
        if 'Contents' in response:
            for obj in response['Contents']:
                if obj['Key'] == f"{prefix}{study_name}/":
                    return bucket_name, obj['Key']
            return False
    except:
        return False
    
    
    
    