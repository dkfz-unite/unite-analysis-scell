import anndata
from scipy.sparse import csc_matrix
from numpy import ndarray, unique
import json
import pandas
import scanpy
import celltypist
import os
import sys
import options

# Set global settings for scanpy
scanpy.settings.verbosity = 0
# scanpy.settings.autosave = True

if len(sys.argv) < 1:
    raise Exception("Provide data directory path as an argument")

if not os.path.exists(sys.argv[1]):
    raise Exception("Provided directory '" + sys.argv[1] + "' was not found")


root_path = sys.argv[1]
print("Performing analysis " + os.path.basename(root_path))


def list_dirs(dir_path: str):
    dirs = []

    for name in os.listdir(dir_path):
        if os.path.isdir(dir_path + "/" + name):
            dirs.append(name)

    return dirs


def read_options(dir_path: str):
    file_path = dir_path + "/options.json"

    if not os.path.exists(file_path):
        raise Exception("options.json not found")
    
    with open(file_path, "r") as file:
        dat = json.load(file)
        opt = options.Options(**dat)
        return opt
    

def read_annotations(dir_path: str):
    file_path = dir_path + "/annotations.tsv"

    if not os.path.exists(file_path):
        raise Exception("annotations.tsv not found")
    
    return pandas.read_csv(file_path, sep = "\t", index_col = 0)


def read_metadata(dir_path: str):
    file_path = dir_path + "/metadata.tsv"

    if not os.path.exists(file_path):
        raise Exception("metadata.tsv not found")
    
    return pandas.read_csv(file_path, sep = "\t")    


def read_data(dir_path: str):
    if not os.path.exists(dir_path):
        raise Exception("directory " + dir_path + " not found")
    
    return scanpy.read_10x_mtx(dir_path)


def enreach_data(data: anndata.AnnData, meta: pandas.DataFrame):
    for key in meta.keys():
        data.obs[key] = str(meta[key].values[0])
    return data


opt = read_options(root_path)
meta = read_metadata(root_path)
data = []
dirs = list_dirs(root_path)

for dir_name in dirs:
    dir_path = root_path + "/" + dir_name
    sample_meta = meta.loc[meta["sample_id"] == dir_name,:]
    sample_data = read_data(dir_path)
    enreach_data(sample_data, sample_meta)
    data.append(sample_data)

adata = anndata.concat(data)

has_cell_types = False

# Adding custom annotations
if opt.meta:
    print("Adding custom annotations")
    annotations = read_annotations(root_path)
    annotations = annotations[~annotations.index.duplicated(keep='first')]
    indices = annotations.index.intersection(adata.obs.index)
    if (len(indices) > 0):
        for column in annotations.columns:
            adata.obs[column] = annotations.loc[indices, column]
            if column == "cell_type":
                has_cell_types = True

adata.obs_names_make_unique()
adata.var_names_make_unique()

# Quality control (QC)
if (opt.qc):
    print("Calculating QC metrics")
    scanpy.pp.calculate_qc_metrics(adata, inplace = True)

# Sparse matrix
if (opt.sparse):
    print("Converting to sparse matrix")
    if (type(adata.X) is ndarray):
        adata.X = csc_matrix(adata.X)

# Preprocessing (PP)
if opt.pp == "seurat":
    # Performs total normalization and log transformation, but not scaling
    print("Running Seurat preprocessing")
    scanpy.pp.recipe_seurat(adata)
elif opt.pp == "zheng17":
    # Performs total normalization, log transformation and scaling
    print("Running Zheng17 preprocessing")
    scanpy.pp.recipe_zheng17(adata)
else:
    print("Running default preprocessing")
    scanpy.pp.filter_cells(adata, min_genes = opt.genes)
    scanpy.pp.filter_genes(adata, min_cells = opt.cells)
    scanpy.pp.normalize_total(adata, target_sum = 1e4)
    scanpy.pp.log1p(adata)

# Cell type annotation
if opt.model is not None:
    print("Predicting cell types")
    celltypist.models.download_models(model = opt.model)
    predictions = celltypist.annotate(adata, model = opt.model)
    if has_cell_types:
        adata.obs["cell_type_predicted"] = predictions.predicted_labels
    else:
        adata.obs["cell_type"] = predictions.predicted_labels

# Scaling the data
if opt.pp != "seurat":
    print("Scaling the data")
    scanpy.pp.scale(adata, zero_center = not opt.sparse)

# Principal Component Analysis (PCA)
if opt.pca:
    print("Running PCA")
    scanpy.pp.pca(adata, svd_solver="arpack", zero_center = not opt.sparse)

# Neighbors
if opt.neighbors:
    print("Calculating neighbors")
    scanpy.pp.neighbors(adata)

# Clustering
if opt.clustering:
    if opt.clustering == "louvain":
        print("Running Louvain clustering")
        scanpy.tl.louvain(adata, key_added = "cluster")
    elif opt.clustering == "leiden":
        print("Running Leiden clustering")
        scanpy.tl.leiden(adata, key_added = "cluster", flavor = "igraph", n_iterations = 2, directed = False)

# Embedding
if opt.embedding:
    if len(unique(adata.obs["cluster"].values)) < 10:
        palette = "tab10"
    else:
        palette = "tab20"

    if "umap" in opt.embedding:
        print("Creating UMAP embedding")
        scanpy.tl.umap(adata)
        # scanpy.pl.umap(adata, color="cluster", palette=palette, show=False, save="umap.svg")
    if "tsne" in opt.embedding:
        print("Creating t-SNE embedding")
        scanpy.tl.tsne(adata)
        # scanpy.pl.tsne(adata, color="cluster", palette=palette, show=False, save="tsne.svg")

# Save result data
print("Saving results")
adata.write(root_path + "/result.data.h5ad")

# Save result metadata
data = {
    "cells_number": adata.n_obs,
    "genes_number": adata.n_vars
}

with open(root_path + "/result.meta.json", "w") as file:
    json.dump(data, file)

# Cleaning up the data
for dir_name in dirs:
    dir_path = root_path + "/" + dir_name
    os.system("rm -r " + dir_path)

print("Done")
