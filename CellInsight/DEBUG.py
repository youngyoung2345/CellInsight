from preprocessing import PanglaoDB_proc_python
from migrations import server

preprocessed_data = PanglaoDB_proc_python.make_anndata('PanglaoDB/SRA322626/SRA322626.sparse.RData')  
temp_path = PanglaoDB_proc_python.save_preprocessed_data(preprocessed_data)

ssh = server.connect_server()
server.upload_anndata_to_server(ssh, temp_path, 'home/root/temp/preprocessed_data.h5ad')
path = server.compute_on_server(ssh, 'umap_plot', 'processing/proc.py', temp_path, 'home/root/temp')
server.download_file_from_server(ssh, path, 'temp/umap_plot.png')
