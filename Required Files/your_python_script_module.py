def process_data(data_file_path):

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
	
	print("Reading data...")
	X = sc.read_h5ad(data_file_path)

	X2 = X.X.toarray()

	X3 = pd.DataFrame(data=X2, columns = X.var_names , index = X.obs.index)

	H = pd.read_table("Myeloid_NMF_Average_Gene_Spectra.txt", sep="\t", index_col=0)

	H2 = H.T

	H3 = H2.filter(items = X3.columns)

	X4 = X3.filter(items = H3.columns)

	H4 = H3.to_numpy()

	X5 = X4.astype(np.float64)

	test = sklearn.decomposition.non_negative_factorization(X5, W=None, H=H4, n_components= 14, init='random', update_H=False, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=1000, alpha=0.0, alpha_W=0.0, alpha_H='same', l1_ratio=0.0, regularization=None, random_state=None, verbose=0, shuffle=False)

	test2 = list(test)

	print("Calculation Completed...")

	processed = pd.DataFrame(test2[0], columns= H.columns, index=X.obs.index)

	row_sums = processed.sum(axis=1)

	processed_data = (processed.div(row_sums, axis=0) * 100)

	return processed_data
