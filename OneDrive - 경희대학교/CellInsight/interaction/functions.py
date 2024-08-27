import io
import pandas as pd 

from migrations import models
from preprocessing import PanglaoDB_proc_python, Single_Cell_Portal_proc 

def load_prefix_list(prefix, delimiter):
    raw_prefix = models.list_s3_objects(models.bucket_name, prefix = prefix, delimiter = delimiter)
    print(raw_prefix)
    prefix_list = [prefix['Prefix'] for prefix in raw_prefix.get('CommonPrefixes', [])]

    return prefix_list

def load_singlecellportal_folder_list():
    raw_information = models.get_s3_object(models.bucket_name, 'singlecellportal/fourth_name.csv')
    information = pd.read_csv(io.BytesIO(raw_information))

    prefix_list = load_prefix_list('singlecellportal/', '/')
    
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
    prefix_list = [item['Key'] for item in raw_prefix.get('Contents', []) if item['Key'].endswith(('.csv', '.tsv', '.txt'))]

    return prefix_list
