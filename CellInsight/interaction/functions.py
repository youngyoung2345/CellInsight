import io
import pandas as pd 

from migrations import models

def load_singlecellportal_folder_list():
    raw_information = models.get_s3_object(models.bucket_name, 'singlecellportal/information.csv')
    information = pd.read_csv(io.BytesIO(raw_information))

    raw_prefix = models.list_s3_objects(models.bucket_name, Prefix='singlecellportal/', Delimiter='/')
    prefix_list = [prefix['Prefix'] for prefix in raw_prefix.get('CommonPrefixes', [])]
    
    name_mapping = dict(zip(information['B'], information['C']))

    mapping = []

    for prefix in prefix_list:
        study_number = prefix.split('/')[1]

        if study_number in name_mapping:
            mapping.append((prefix, name_mapping[prefix]))
        else:
            mapping.append((prefix, study_number))

    return mapping

def load_cluster_file_list(folder_name):
    cluster_prefix = folder_name.rstrip('/') + '/cluster/'

    raw_prefix = models.list_s3_objects(models.bucket_name, cluster_prefix)
    prefix = [item['Key'] for item in raw_prefix.get('Contents', []) if item['Key'].endswith(('.csv', '.tsv', '.txt'))]

    return prefix
