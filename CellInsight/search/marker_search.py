import pandas as pd
import boto3
from botocore.client import Config
import io

def fetch_markers(selected_organ='All', selected_cell_type='All'):
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

    response = s3_client.get_object(Bucket=bucket_name, Key=file_key)
    content = response['Body'].read()
    dataset = pd.read_csv(io.BytesIO(content), delimiter='\t')

    # Organ별로 cell type 목록 생성
    organ_cell_types = dataset.groupby('organ')['cell type'].apply(list).to_dict()

    # Cell type 이름과 개수 집계
    for organ in organ_cell_types:
        cell_type_counts = pd.Series(organ_cell_types[organ]).value_counts()
        organ_cell_types[organ] = [f"{cell_type} ({count})" for cell_type, count in cell_type_counts.items()]

    # Official gene symbol에 해당하는 cell type 종류 계산
    official_gene_symbol_counts = dataset.groupby('official gene symbol')['cell type'].nunique()

    if selected_organ != 'All':
        dataset = dataset[dataset['organ'] == selected_organ]
    if selected_cell_type != 'All':
        dataset = dataset[dataset['cell type'].str.contains(selected_cell_type.split(' (')[0])]

    # HTML 코드 생성
    html_code = """
    <table id="data-table">
        <thead>
            <tr>
                <th>Species</th>
                <th>Official gene symbol</th>
                <th>UI</th>
                <th>Sensitivity (human)</th>
                <th>Sensitivity (mouse)</th>
                <th>Specificity (human)</th>
                <th>Specificity (mouse)</th>
                <th>Marker count</th>
                <th>Cell type</th>
                <th>Germ layer</th>
                <th>Organ</th>
                <th>Aliases</th>
                <th>Product description</th>
                <th>Action</th>
            </tr>
        </thead>
        <tbody>
    """

    for _, row in dataset.iterrows():
        html_code += "<tr>"
        html_code += f"<td>{row['species']}</td>"
        html_code += f"<td>{row['official gene symbol']}</td>"
        html_code += f"<td>{round(row['ubiquitousness index'], 3)}</td>"
        html_code += f"<td>{round(row['sensitivity_human'], 3)}</td>"
        html_code += f"<td>{round(row['sensitivity_mouse'], 3)}</td>"
        html_code += f"<td>{round(row['specificity_human'], 3)}</td>"
        html_code += f"<td>{round(row['specificity_mouse'], 3)}</td>"
        marker_count = official_gene_symbol_counts[row['official gene symbol']]
        html_code += f"<td>{marker_count}</td>"
        html_code += f"<td>{row['cell type']}</td>"
        html_code += f"<td>{row['germ layer']}</td>"
        html_code += f"<td>{row['organ']}</td>"
        html_code += f"<td>{row['nicknames']}</td>"
        html_code += f"<td>{row['product description']}</td>"
        html_code += "<td><button onclick='alert(\"Action\")'>flag</button></td>"
        html_code += "</tr>"

    html_code += """
        </tbody>
    </table>
    """


    return organ_cell_types, html_code
