import os
import matplotlib
matplotlib.use('Agg')
import scanpy as sc
import anndata as ad
import scrublet as scr
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



def preprocess_and_visualize(file_path, file_format, min_counts, min_genes, max_genes, pct_counts_mt):
    # 파일 형식에 따라 데이터를 로드
    if file_format == 'panglaodb':
        import PanglaoDB_proc_python as pdl
        adata = pdl.process_PanglaoDB(file_path)
    elif file_format == 'h5ad':
        from . import Single_Cell_Portal_proc_copy as scp  # 상대 경로로 수정
        adata = scp.process_Single_Cell_Portal(file_path)
    # 시각화 및 저장 경로 설정
    figures_path = 'media/figures'
    if not os.path.exists(figures_path):
        os.makedirs(figures_path)
    
    violin_plot_path = os.path.join(figures_path, 'violin_plot.png')
    umap_plot_path = os.path.join(figures_path, 'umap_plot.png')
    
    
    # QC 및 전처리
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scr.Scrublet(adata.X).scrub_doublets(verbose=False)
    adata.obs['predicted_doublets'] = adata.obs['predicted_doublets'].astype(str)
    # 사용자로부터 입력받은 QC 값을 사용
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs['pct_counts_mt'] < pct_counts_mt, :]
    adata = adata[adata.obs['predicted_doublets'] == 'False', :]
 # 시각화 파일 저장
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    plt.savefig(violin_plot_path)
    plt.close()
    
    # Normalization
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata)
    sc.pp.scale(adata, max_value=10)

    #Umap

    #cell_cycle_genes = [x.strip() for x in open('/content/drive/MyDrive/2023_BIML/data/regev_lab_cell_cycle_genes.txt')]
    #s_genes = cell_cycle_genes[:43]
    #g2m_genes = cell_cycle_genes[43:]
    #sc.tl.score_genes_cell_cycle(adata_1, s_genes=s_genes, g2m_genes=g2m_genes)
    #cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_1.var_names]
    #adata_1.var.loc[cell_cycle_genes, 'highly_variable'] = False

    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=10)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

     
   
    
    sc.pl.umap(adata, color=['leiden'])
    plt.savefig(umap_plot_path)
    plt.close()

    return violin_plot_path, umap_plot_path