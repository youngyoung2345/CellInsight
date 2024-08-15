import migrations.models as models
import pandas as pd
import io

def search_study(bucket_name, database_name, study_name):
    try:
        prefix = f"{database_name}/"

        response = models.s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
        
        if 'Contents' in response:
            for obj in response['Contents']:
                if obj['Key'] == f"{prefix}{study_name}/":
                    return bucket_name, obj['Key']
            return None
    except:
        return False
    
def search_study_case(bucket_name, filtered_data):

    before_SRA = 0
    key_list = []

    for index, row in filtered_data.iterrows():
        SRA, SRS, SRR = row['SRA'], row['SRS'], row['SRR']

        if before_SRA == SRA:
            for obj in before_response['Content']:
                    if f'SRS{SRS}' in obj['Key']:
                        key_list.append(obj['Key'])
                        break
        else:
            response = models.list_s3_objects(bucket_name, f'PanglaoDB/SRA{SRA}/PanglaoDB/')
            
            if SRR != 'notused': 
                key_list.append(response['Content']['Key'])
            else:
                for obj in response['Content']:
                    if f'SRS{SRS}' in obj['Key']:
                        key_list.append(obj['Key'])
                        break
            
        before_SRA = SRA
        before_response = response

    return key_list

def search_PanglaoDB_study_attribute(bucket_name, attribute, input):
    try:
        response = models.get_s3_objects(bucket_name, 'PanglaoDB/data.csv')
        file_content = response['Body'].read()

        data = pd.read_csv(io.BytesIO(file_content))
        filtered_data = data[data[attribute].str.contains(input, na=False, case=False)]

        if filtered_data.empty:
            return None
        
        key_list = search_study_case(bucket_name, filtered_data[['SRA', 'SRS', 'SRR']])
        return key_list
    except:
        return False

def search_singlecellportal_study_attribute(bucket_name, attribute, input):

    return