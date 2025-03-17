import sklearn
import sklearn.decomposition
from sklearn.decomposition import non_negative_factorization
import numpy as np
import scanpy as sc
import csv
import scipy
import pandas as pd
import h5py
import anndata as ad
import os
import tempfile
import gzip
import shutil

def compress_file(input_path, output_path):
    with open(input_path, 'rb') as f_in, gzip.open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def process_annotation_h5ad_data(data_file_path):
    print("Reading h5ad data...")
    X = sc.read_h5ad(data_file_path)
    
    X.obs_names_make_unique()

    X2 = X.X.toarray()

    X3 = pd.DataFrame(data=X2, columns=X.var_names, index=X.obs.index)

    H = np.load("cnmf_run.spectra.k_18.dt_0_015.consensus.df.npz", allow_pickle=True)

    H2 = pd.DataFrame(H['data'], columns=H['columns'], index=H['index'])

    H3 = H2.filter(items=X3.columns)

    X4 = X3.filter(items=H3.columns)
    
    X5 = X4.values

    H4 = H3.to_numpy()

    X6 = X5.astype(np.float64)

    test = non_negative_factorization(X6, W=None, H=H4, n_components=18, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, random_state=None, verbose=0, shuffle=False)

    test2 = list(test)

    print("Calculation Completed...")

    processed = pd.DataFrame(test2[0], columns=H['index'], index=X4.index)

    row_sums = processed.sum(axis=1)

    processed_data = (processed.div(row_sums, axis=0) * 100)
    
    new_column_names = ['Tcells', 'AC', 'NPC1_OPC', 'Microglia', 'MES2', 'Vascular_MES1', 'Oligodendrocytes', 'MES1', 'CD14_Mono', 'cDC', 'Neutrophils', 'NPC2', 'Giant_Cell_GBM', 'Cycling', 'Pericytes', 'Plasma', 'Endothelial', 'Mast']
    processed_data.columns = new_column_names

    return processed_data

def process_annotation_csv_data(data_file_path):
    print("Reading csv data...")

    X = pd.read_table(data_file_path, index_col=0, sep=',')

    X2 = X.T

    H = np.load("cnmf_run.spectra.k_18.dt_0_015.consensus.df.npz", allow_pickle=True)

    H2 = pd.DataFrame(H['data'], columns=H['columns'], index=H['index'])

    H3 = H2.filter(items=X2.columns)
    
    X4 = X2.filter(items=H3.columns)

    H4 = H3.to_numpy()

    X5 = X4.values
    
    X6 = X5.astype(np.float64)

    test = non_negative_factorization(X6, W=None, H=H4, n_components=18, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=1000, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, random_state=None, verbose=0, shuffle=False)

    test2 = list(test)

    print("Calculation Completed...")

    processed = pd.DataFrame(test2[0], columns=H['index'], index=X4.index)

    row_sums = processed.sum(axis=1)

    processed_data = (processed.div(row_sums, axis=0) * 100)
    
    new_column_names = ['Tcells', 'AC', 'NPC1_OPC', 'Microglia', 'MES2', 'Vascular_MES1', 'Oligodendrocytes', 'MES1', 'CD14_Mono', 'cDC', 'Neutrophils', 'NPC2', 'Giant_Cell_GBM', 'Cycling', 'Pericytes', 'Plasma', 'Endothelial', 'Mast']
    processed_data.columns = new_column_names

    return processed_data

def compress_file(input_path, output_path):
    with open(input_path, 'rb') as f_in, gzip.open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def process_annotation_mtx_data(mtx_file_path, barcodes_file_path, features_file_path):
    print("MTX file path received:", mtx_file_path)
    print("Barcodes file path received:", barcodes_file_path)
    print("Features file path received:", features_file_path)

    with tempfile.TemporaryDirectory() as temp_dir:
        mtx_dest_path = os.path.join(temp_dir, 'matrix.mtx.gz')
        barcodes_dest_path = os.path.join(temp_dir, 'barcodes.tsv.gz')
        features_dest_path = os.path.join(temp_dir, 'features.tsv.gz')

        # Compress files if they are not already compressed
        if mtx_file_path.endswith('.gz'):
            shutil.copy(mtx_file_path, mtx_dest_path)
        else:
            compress_file(mtx_file_path, mtx_dest_path)

        if barcodes_file_path.endswith('.gz'):
            shutil.copy(barcodes_file_path, barcodes_dest_path)
        else:
            compress_file(barcodes_file_path, barcodes_dest_path)

        if features_file_path.endswith('.gz'):
            shutil.copy(features_file_path, features_dest_path)
        else:
            compress_file(features_file_path, features_dest_path)

        print("Reading features file to check format...")
        features = pd.read_csv(features_dest_path, sep='\t', header=None)
        print(features.head())

        if features.shape[1] < 2:
            print("Features file does not have enough columns. Expected at least 2 columns.")
            raise ValueError("Invalid features file format")

        print("Reading mtx data from temporary directory...")
        X = sc.read_10x_mtx(temp_dir, var_names='gene_symbols', cache=True)
      
    X.var_names_make_unique()
    
    print("MTX data loaded successfully")

    X2 = X.X.toarray()

    X3 = pd.DataFrame(data=X2, columns=X.var_names, index=X.obs.index)

    H = np.load("cnmf_run.spectra.k_18.dt_0_015.consensus.df.npz", allow_pickle=True)

    H2 = pd.DataFrame(H['data'], columns=H['columns'], index=H['index'])

    H3 = H2.filter(items=X3.columns)
    
    X4 = X3.filter(items=H3.columns)

    H4 = H3.to_numpy()

    X5 = X4.values
    
    X6 = X5.astype(np.float64)

    test = non_negative_factorization(X6, W=None, H=H4, n_components=18, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=1000, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, random_state=None, verbose=0, shuffle=False)

    test2 = list(test)

    print("Calculation Completed...")

    processed = pd.DataFrame(test2[0], columns=H['index'], index=X4.index)

    row_sums = processed.sum(axis=1)

    processed_data = (processed.div(row_sums, axis=0) * 100)
    
    new_column_names = ['Tcells', 'AC', 'NPC1_OPC', 'Microglia', 'MES2', 'Vascular_MES1', 'Oligodendrocytes', 'MES1', 'CD14_Mono', 'cDC', 'Neutrophils', 'NPC2', 'Giant_Cell_GBM', 'Cycling', 'Pericytes', 'Plasma', 'Endothelial', 'Mast']
    processed_data.columns = new_column_names

    return processed_data

def process_myeloid_program_h5ad_data(data_file_path):
    print("Reading h5ad data...")
    X = sc.read_h5ad(data_file_path)
    
    X.obs_names_make_unique()

    X2 = X.X.toarray()

    X3 = pd.DataFrame(data=X2, columns=X.var_names, index=X.obs.index)

    H = pd.read_table("Myeloid_NMF_Average_Gene_Spectra.txt", sep="\t", index_col=0)

    H2 = H.T

    H3 = H2.filter(items=X3.columns)

    X4 = X3.filter(items=H3.columns)
    
    X5 = X4.values
    
    H4 = H3.to_numpy()

    X6 = X5.astype(np.float64)

    test = non_negative_factorization(X6, W=None, H=H4, n_components=14, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=1000, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, random_state=None, verbose=0, shuffle=False)

    test2 = list(test)

    print("Calculation Completed...")

    processed = pd.DataFrame(test2[0], columns=H.columns, index=X.obs.index)

    row_sums = processed.sum(axis=1)

    processed_data = (processed.div(row_sums, axis=0) * 100)

    return processed_data

def process_myeloid_program_csv_data(data_file_path):
    print("Reading csv data...")

    X = pd.read_table(data_file_path, index_col=0, sep=',')

    X2 = X.T

    H = pd.read_table("Myeloid_NMF_Average_Gene_Spectra.txt", sep="\t", index_col=0)

    H2 = H.T

    X3 = X2.filter(items=H2.columns)

    H3 = H2.filter(items=X3.columns)

    H4 = H3.to_numpy()
    
    X4 = X3.values

    X5 = X4.astype(np.float64)

    test = non_negative_factorization(X5, W=None, H=H4, n_components=14, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=1000, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, random_state=None, verbose=0, shuffle=False)

    test2 = list(test)

    print("Calculation Completed...")

    processed = pd.DataFrame(test2[0], columns=H.columns, index=X3.index)

    row_sums = processed.sum(axis=1)

    processed_data = (processed.div(row_sums, axis=0) * 100)

    return processed_data

def process_myeloid_program_mtx_data(mtx_file_path, barcodes_file_path, features_file_path):
    print("MTX file path received:", mtx_file_path)
    print("Barcodes file path received:", barcodes_file_path)
    print("Features file path received:", features_file_path)

    with tempfile.TemporaryDirectory() as temp_dir:
        mtx_dest_path = os.path.join(temp_dir, 'matrix.mtx.gz')
        barcodes_dest_path = os.path.join(temp_dir, 'barcodes.tsv.gz')
        features_dest_path = os.path.join(temp_dir, 'features.tsv.gz')

        # Compress files if they are not already compressed
        if mtx_file_path.endswith('.gz'):
            shutil.copy(mtx_file_path, mtx_dest_path)
        else:
            compress_file(mtx_file_path, mtx_dest_path)

        if barcodes_file_path.endswith('.gz'):
            shutil.copy(barcodes_file_path, barcodes_dest_path)
        else:
            compress_file(barcodes_file_path, barcodes_dest_path)

        if features_file_path.endswith('.gz'):
            shutil.copy(features_file_path, features_dest_path)
        else:
            compress_file(features_file_path, features_dest_path)

        print("Reading features file to check format...")
        features = pd.read_csv(features_dest_path, sep='\t', header=None)
        print(features.head())

        if features.shape[1] < 2:
            print("Features file does not have enough columns. Expected at least 2 columns.")
            raise ValueError("Invalid features file format")

        print("Reading mtx data from temporary directory...")
        X = sc.read_10x_mtx(temp_dir, var_names='gene_symbols', cache=True)
      
    X.var_names_make_unique()
    
    print("MTX data loaded successfully")

    X2 = X.X.toarray()

    X3 = pd.DataFrame(data=X2, columns=X.var_names, index=X.obs.index)

    H = pd.read_table("Myeloid_NMF_Average_Gene_Spectra.txt", sep="\t", index_col=0)

    H2 = H.T

    H3 = H2.filter(items=X3.columns)

    X4 = X3.filter(items=H3.columns)
    
    X5 = X4.values
    
    H4 = H3.to_numpy()

    X6 = X5.astype(np.float64)

    test = non_negative_factorization(X6, W=None, H=H4, n_components=14, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=1000, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, random_state=None, verbose=0, shuffle=False)

    test2 = list(test)

    print("Calculation Completed...")

    processed = pd.DataFrame(test2[0], columns=H.columns, index=X.obs.index)

    row_sums = processed.sum(axis=1)

    processed_data = (processed.div(row_sums, axis=0) * 100)

    return processed_data
