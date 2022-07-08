# In this homework, we will implement PCA in three ways, and we will test the implementation on peripheral blood mononuclear cells(PBMCs) single-cell gene expression dataset 
# The dataset is stored as a matrix, where each row records the the gene expression feature of an individual cell. There a totally 3000 cells, with 3000 feature for each cell.

# TODO: Implement both function `pca', `svd' and `pca_sklearn', run the code and submit your result, including figures and matrix output. Note that we have loaded all necessary packages for the implementation, 
# the use additional packages is not allowed

# Note: 
# 1. To run the code, you need to install one additional package: [sklearn](https://scikit-learn.org/stable/)
# 2. We have implemented the data normalization in the functions using `StandardScaler()'.
# 4. plotting function is already provided as `plot_latent'. 

# Additional references:
# 1. sklearn Pipeline API: https://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html
# 2. sklearn PCA API: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
# 3. sklearn StandardScaler API: https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html#sklearn.preprocessing.StandardScaler

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


def pca(X, n_pcs = 2):
    """\
    Implementing `PCA' from scratch, using covariance matrix
    Parameters:
    -----------
        X: gene expression matrix, of the shape (n_cells, n_features)
        n_pcs: number of reduced dimensions
    Returns:
    -----------
        X_pca: pca representation of X, of the shape (n_cells, n_pcs).
    """
    # Data normalization  
    X = StandardScaler().fit_transform(X)
    
    # implement your code here
    pass



def svd(X, n_pcs = 2):
    """\
    Implementing `PCA' from scratch, using the data matrix directly with SVD
    Parameters:
    -----------
        X: gene expression matrix, of the shape (n_cells, n_features)
        n_pcs: number of reduced dimensions
    Returns:
    -----------
        X_pca: pca representation of X, of the shape (n_cells, n_pcs).
    """
    # Data normalization 
    X = StandardScaler().fit_transform(X)
    
    # implement your code here
    pass
    

def pca_sklearn(X, n_pcs = 2):
    """\
    Instead of implementing `PCA' step by step, we can also call `PCA' directly using python package `sklearn'.
    Implement `PCA' using sklearn, should include only 1-3 lines of code.
    Parameters:
    -----------
        X: gene expression matrix, of the shape (n_cells, n_features)
        n_pcs: number of reduced dimensions
    Returns:
    -----------
        X_pca: pca representation of X, of the shape (n_cells, n_pcs).
    """
    # Data normalization  
    X = StandardScaler().fit_transform(X)
    
    # implement your code here
    pass
    

# plotting function, no modification required
def plot_latent(z, anno, save = None, figsize = (10,10), axis_label = "Latent", **kwargs):
    _kwargs = {
        "s": 10,
        "alpha": 0.9,
    }
    _kwargs.update(kwargs)

    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot()
    cluster_types = set([x for x in np.unique(anno)])
    colormap = plt.cm.get_cmap("tab20", len(cluster_types))

    for i, cluster_type in enumerate(cluster_types):
        index = np.where(anno == cluster_type)[0]
        ax.scatter(z[index,0], z[index,1], color = colormap(i), label = cluster_type, **_kwargs)
    
    ax.legend(loc='upper left', prop={'size': 15}, frameon = False, ncol = 1, bbox_to_anchor=(1.04, 1))
    
    ax.tick_params(axis = "both", which = "major", labelsize = 15)

    ax.set_xlabel(axis_label + " 1", fontsize = 19)
    ax.set_ylabel(axis_label + " 2", fontsize = 19)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)  
    
    if save:
        fig.savefig(save, bbox_inches = "tight")
    
    print(save)


# read data matrix
expr_ctrl = pd.read_csv("counts_PBMC.csv", sep = ",", index_col = 0).values
anno_ctrl = pd.read_csv("celltypes_PBMC.txt", sep = "\t", header = None)

X_pca = pca(expr_ctrl)
plot_latent(X_pca, anno_ctrl, axis_label = "PCA", save = "PCA.pdf")
np.save(file = "PCA.npy", arr = X_pca)

X_svd = svd(expr_ctrl)
plot_latent(X_svd, anno_ctrl, axis_label = "PCA", save = "SVD.pdf")
np.save(file = "SVD.npy", arr = X_svd)

X_pca_skl = pca_sklearn(expr_ctrl)
plot_latent(X_pca_skl, anno_ctrl, axis_label = "PCA", save = "PCA_skl.pdf")
np.save(file = "PCA_skl.npy", arr = X_pca_skl)

