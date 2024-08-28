import os
import numpy as np
import scanpy as sc
import anndata as ad

import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt

def draw_and_save_violin_plot(preprocessed_data, figures_path='media/figures'):
    violin_plot_path = os.path.join(figures_path, 'violin_plot.png')

    sc.pl.violin(preprocessed_data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    
    plt.savefig('media/figures/violin_plot.png')
    plt.close()

    return violin_plot_path

def draw_and_save_umap_plot(preprocessed_data, figures_path='media/figures'):
    umap_plot_path = os.path.join(figures_path, 'umap_plot.png')

    sc.tl.pca(preprocessed_data)
    sc.pp.neighbors(preprocessed_data, n_pcs=10)
    sc.tl.umap(preprocessed_data)
    sc.tl.leiden(preprocessed_data, resolution=0.5)

    sc.pl.umap(preprocessed_data, color=['leiden'])

    plt.savefig('media/figures/umap_plot.png')
    plt.close()

    return umap_plot_path

def draw_and_save_violin_plot_on_server(preprocessed_data_path, server_path):
    preprocessed_data = ad.read_h5ad(preprocessed_data_path)
    
    sc.pl.violin(preprocessed_data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

    if not os.path.exists(server_path):
        os.makedirs(server_path)

    violin_plot_path = os.path.join(server_path, 'violin_plot.png')

    plt.savefig(violin_plot_path)
    plt.close()

    return violin_plot_path


def draw_and_save_umap_plot_on_server(preprocessed_data_path, server_path):
    preprocessed_data = ad.read_h5ad(preprocessed_data_path)

    sc.tl.pca(preprocessed_data)
    sc.pp.neighbors(preprocessed_data, n_pcs=10)
    sc.tl.umap(preprocessed_data)
    sc.tl.leiden(preprocessed_data, resolution=0.5)

    sc.pl.umap(preprocessed_data, color=['leiden'])

    if not os.path.exists(server_path):
        os.makedirs(server_path)

    umap_plot_path = os.path.join(server_path, 'umap_plot.png')

    plt.savefig(umap_plot_path)
    plt.close()

    return umap_plot_path

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--option', required=True)
    parser.add_argument('--preprocessed_data', required=True)
    parser.add_argument('--server_path', required=True)

    args = parser.parse_args()

    match args.option:
        case 'umap_plot':
            draw_and_save_umap_plot(args.preprocessed_data, args.server_path)
        case 'violin_plot':
            draw_and_save_violin_plot(args.preprocessed_data, args.server_path)
