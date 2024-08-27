import os
import scanpy as sc
import matplotlib.pyplot as plt

def draw_and_save_violin_plot(preprocessed_data, figures_path='media/figures'):
    if not os.path.exists(figures_path):
        os.makedirs(figures_path)

    violin_plot_path = os.path.join(figures_path, 'violin_plot.png')

    sc.pl.violin(preprocessed_data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    
    plt.savefig(violin_plot_path)
    plt.close()

    return violin_plot_path

def draw_and_save_umap_plot(preprocessed_data, figures_path='media/figures'):
    if not os.path.exists(figures_path):
        os.makedirs(figures_path)
        
    umap_plot_path = os.path.join(figures_path, 'umap_plot.png')

    sc.tl.pca(preprocessed_data)
    sc.pp.neighbors(preprocessed_data, n_pcs=10)
    sc.tl.umap(preprocessed_data)
    sc.tl.leiden(preprocessed_data, resolution=0.5)

    sc.pl.umap(preprocessed_data, color=['leiden'])

    plt.savefig(umap_plot_path)
    plt.close()

    return umap_plot_path
