import os
import matplotlib
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt

from . import PanglaoDB_proc_python as pdl
from . import Single_Cell_Portal_proc as scp

matplotlib.use('Agg')

def load_data(file_path, file_format):
    match file_format:
        case 'RData':
            preprocessed_data = pdl.make_anndata(file_path)
        case _:
            preprocessed_data = scp.process_Single_Cell_Portal(file_path, file_format)

    return preprocessed_data

def preprocess_data(preprocessed_data, min_counts, min_genes, max_genes, pct_counts_mt):
    preprocessed_data.var_names_make_unique()
    preprocessed_data.var['mt'] = preprocessed_data.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(preprocessed_data, qc_vars=['mt'], percent_top=None, inplace=True)
    preprocessed_data.obs['doublet_scores'], preprocessed_data.obs['predicted_doublets'] = scr.Scrublet(preprocessed_data.X).scrub_doublets(verbose=False)
    preprocessed_data.obs['predicted_doublets'] = preprocessed_data.obs['predicted_doublets'].astype(str)

    sc.pp.filter_cells(preprocessed_data, min_counts=min_counts)
    sc.pp.filter_cells(preprocessed_data, min_genes=min_genes)
    sc.pp.filter_cells(preprocessed_data, max_genes=max_genes)
    
    preprocessed_data = preprocessed_data[preprocessed_data.obs['pct_counts_mt'] < pct_counts_mt, :]
    preprocessed_data = preprocessed_data[preprocessed_data.obs['predicted_doublets'] == 'False', :]

    return preprocessed_data

def normalize_data(preprocessed_data):
    preprocessed_data.layers['counts'] = preprocessed_data.X.copy()
    sc.pp.normalize_total(preprocessed_data, target_sum=1e4)
    sc.pp.log1p(preprocessed_data)
    preprocessed_data.raw = preprocessed_data
    sc.pp.highly_variable_genes(preprocessed_data)
    sc.pp.scale(preprocessed_data, max_value=10)

    return preprocessed_data