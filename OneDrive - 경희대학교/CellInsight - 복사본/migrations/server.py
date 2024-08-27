import paramiko
from preprocessing import PanglaoDB_proc_python

information = {
    'IP': '211.188.54.91', 
    'port': 22,
    'username': 'root',
    'key_path': 'Cellinsight-key.pem'
}

def connect_server():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(information['IP'], port=information['port'], username=information['username'], key_filename=information['key_path'])

    return ssh

def upload_anndata_to_server(ssh, temp_path, server_path):
    sftp = ssh.open_sftp()
    sftp.put(temp_path, server_path)
    sftp.close()

    return

def compute_on_server(ssh, option, script_path, preprocessed_data_path, server_path):
    command = f"python3 {script_path} --option{option} --preprocessed_data{preprocessed_data_path} --server_path{server_path}"
    
    stdin, stdout, stderr = ssh.exec_command(command)

    return stdout

def download_file_from_server(ssh, server_path, local_file_path):
    sftp = ssh.open_sftp()
    sftp.get(server_path, local_file_path)
    sftp.close()

preprocessed_data = PanglaoDB_proc_python.make_anndata('PanglaoDB/SRA322626/SRA322626.sparse.RData')   
temp_path = PanglaoDB_proc_python.save_preprocessed_data(preprocessed_data)

ssh = connect_server()
upload_anndata_to_server(ssh, temp_path, 'home/root/temp/preprocessed_data.h5ad')
compute_on_server(ssh, 'violin_plot', 'processing/proc.py', temp_path, 'home/root/temp/preprocessed_data.h5ad')
download_file_from_server(ssh, 'home/root/temp/preprocessed_data.h5ad', 'temp/violin_plot.png')
