# -*- coding: utf-8 -*-
"""hw3_1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1m7ak7C0zccJbclealX6gHmkenGpcGAbX

# HW3-1: Comparison of different clustering algorithms
In our lecture, we talked about **K-means clustering**, **spectral clustering** and **Leiden clustering(a modified version of louvain clustering algorithm)**. 

In this homework, we will play around with those three algorithms on different datasets, and see how those algorithms perform in different circumstances.

#### Note

For all two notebooks, please insert your code under comments like this:
```
# =========================================
# instruction on what to implement
# REPLACE `pass` WITH YOUR CODE or INSERT YOUR CODE BELOW
# =========================================
```

**You can add additional lines of code when necessary.**
"""

# download louvain package
!pip install louvain python-igraph leidenalg graphtools

# Commented out IPython magic to ensure Python compatibility.
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline
plt.rc('font', size=14)

import sklearn
import sklearn.cluster
import sklearn.datasets

import louvain
import graphtools as gt

import warnings
warnings.filterwarnings("ignore")

"""## Dataset generation
First, we need to generate different datasets for clustering. We will use three datasets in $\mathbb{R}^2$: (**Note that the real data is usually high-dimensional, we use 2 dimensional dataset here just for better visualization**).

* **Double-circles**: Dataset of the shape of two circles with the same center but different radius.

* **Regular blobs**: several regular blobs, regular cluster shape.

* **Uniform distribution**: Data generated from uniform distribution.

You don't need to do anything at this step, since we have already implemented it for you.
"""

np.random.seed(0)
# ============
# Generate datasets. We choose the size big enough to see the scalability
# of the algorithms, but not too big to avoid too long running times
# ============
n_samples = 1500

# Circles
noisy_circles = sklearn.datasets.make_circles(
    n_samples=n_samples, 
    # Scale factor between inner and outer circle
    factor=.5,
    # Gaussian noise added to each point
    noise=.05)

# Uniform square
no_structure = (np.random.uniform(size=(n_samples, 2)), None)

# blobs with varied variances
varied = sklearn.datasets.make_blobs(n_samples=n_samples,
                                     cluster_std=[1.0, 2.5, 0.5],
                                     random_state=8)

# ============
# Associate each dataset with the correct # of clusters
# ============

default_base = {'n_clusters': 3}

generated_datasets = [
    (noisy_circles, {'n_clusters': 2}),
    (varied,      {}),
    (no_structure, {})]
  
fig, axes = plt.subplots(1,3,figsize=(12,4))

for i, (dataset, _) in enumerate(generated_datasets):
    ax = axes[i]
    X, y = dataset
    
    # normalize dataset for easier parameter selection
    X = sklearn.preprocessing.StandardScaler().fit_transform(X)
    ax.scatter(X[:,0], X[:,1], c=y, s = 0.5)
    
fig.tight_layout()

"""## K-means algorithm
You need to fill in the blank area of the k-means algorithm below. You will implement this part using `sklearn`.

**Useful tutorial**:
* https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html

"""

def kmeans_skl(X, n_clusts):
    """\
    Kmeans algorithm

    Parameters:
    -----------
    X:
      dataset, numpy array of the shape (n_samples, n_features)
    n_clusts:
      number of clusters
    
    Return:
    ----------
    group:
      group identity of the clustering result, numpy array of the shape (n_samples,). 
      e.g, if there are three data samples clustered into two groups, the return should be something like, numpy.array([0,1,1])
    """
    # =========================================
    # implement kmeans algorithm using sklearn
    # REPLACE `pass` WITH YOUR CODE
    # =========================================
    kmeans = sklearn.cluster.KMeans(n_clusters = n_clusts, random_state = 0).fit_predict(X)
    return kmeans

"""## Spectral clustering algorithm
You need to fill in the blank area of the spectral clustering algorithm below. You will implement this part using `sklearn`.

**Useful tutorial**:
* https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html

**Note**:
* Make sure to set the `affinity` in spectral clustering function to `"nearest_neighbors"`. The default setting in sklearn was `"rbf"`, which will not produce a very good result.
"""

def spectral_skl(X, n_clusts):
    """\
    spectral clustering algorithm

    Parameters:
    -----------
    X:
      dataset, numpy array of the shape (n_samples, n_features)
    n_clusts:
      number of clusters
    
    Return:
    ----------
    group:
      group identity of the clustering result, numpy array of the shape (n_samples,). 
      e.g, if there are three data samples clustered into two groups, the return should be something like, numpy.array([0,1,1])
    """
    # =========================================
    # implement spectral clustering algorithm using sklearn
    # REPLACE `pass` WITH YOUR CODE
    # =========================================
    spectral = sklearn.cluster.SpectralClustering(n_clusters=n_clusts, affinity="nearest_neighbors").fit_predict(X)
    
    return spectral

"""## Leiden algorithm
Leiden algorithm is a graph based algorithm, instead of specifying the number of clusters, leiden algorithm use resolution to control the size of clusters.

You don't need to do anything here, we have already implemented leiden algorithm for you, but you can play around with the resolution term to see how this will affect the number of clusters
"""

def leiden(X, resolution=0.01):
    """\
    leiden clustering algorithm.

    Parameters:
    -----------
    X:
      dataset, numpy array of the shape (n_samples, n_features)
    resolution:
      resolution term
    
    Return:
    ----------
    group:
      group identity of the clustering result, numpy array of the shape (n_samples,). 
      e.g, if there are three data samples clustered into two groups, the return should be something like, numpy.array([0,1,1])
    """ 
    G = gt.Graph(X)
    G_igraph = G.to_igraph()
    part = louvain.find_partition(G_igraph, louvain.RBConfigurationVertexPartition, 
                                  weights="weight", resolution_parameter=resolution)
    groups = np.array(part.membership)
    return groups

"""## Test the algorithms in different datasets
We will test the algorithms in different datasets. Run the code below and see how the algorithms perform under different circumstances. Submit the plots that you generate. Discuss with your collaborator and briefly discribe what you find in you submission.

"""

fig, axes = plt.subplots(3,3, figsize=(15, 10))
plot_title = True

titles = ['KMeans', 'Spectral', 'Leiden']

# loop through all three datasets
for i_dataset, (dataset, cluster_params) in enumerate(generated_datasets):
    
    params = default_base.copy()
    params.update(cluster_params)

    X, y = dataset
    
    # normalize dataset
    X = sklearn.preprocessing.StandardScaler().fit_transform(X)
    
    # ============
    # Run clustering algorithms
    # ============
    clusters = []
    
    # KMeans
    clusters.append(kmeans_skl(X, n_clusts = params['n_clusters']))
    
    # Spectral Clustering
    clusters.append(spectral_skl(X, n_clusts = params['n_clusters']))
    
    # Leiden
    clusters.append(leiden(X, resolution = 0.01))
 
    # ============
    # Plot clustering results
    # ============
    row_axes = axes[i_dataset]
    
    for i, ax in enumerate(row_axes.flatten()):
        curr_cluster = clusters[i]
        if plot_title:
            curr_title = '{}'.format(titles[i])
        else:
            curr_title = None
            
        ax.scatter(X[:,0], X[:,1], c=curr_cluster, s = 0.5)
        ax.set_title(curr_title)

    plot_title=False
fig.tight_layout()

"""## Running time
Next we will test the running time of the algorithms with datasets of different sizes. We will use the **Regular blobs** shaped dataset with different number of samples: 1000, 5000, 10000 and see how the running time changes. Report your running times. Discuss with your collaborator about what you see, briefly analyze the result in your submission.

"""

np.random.seed(0)

# generated datasets
generated_datasets = [(sklearn.datasets.make_blobs(n_samples=n_samples, random_state=8, cluster_std=1), n_samples) for n_samples in [1000, 5000, 10000]]
  
fig, axes = plt.subplots(1,3,figsize=(12,4))

for i, (dataset, _) in enumerate(generated_datasets):
    ax = axes[i]
    X, y = dataset
    
    # normalize dataset for easier parameter selection
    X = sklearn.preprocessing.StandardScaler().fit_transform(X)
    ax.scatter(X[:,0], X[:,1], c=y, s = 0.5)
    
fig.tight_layout()

fig, axes = plt.subplots(3,3, figsize=(15, 10))
plot_title = True


titles = ['KMeans', 'Spectral', 'Leiden']

# loop through all three datasets
for i_dataset, (dataset, n_samples) in enumerate(generated_datasets):
    print("dataset: "+str(i_dataset)+", number of samples: " + str(n_samples))
    
    params = default_base.copy()
    params.update(cluster_params)

    X, y = dataset
    
    # normalize dataset
    X = sklearn.preprocessing.StandardScaler().fit_transform(X)
    
    # ============
    # Run clustering algorithms
    # ============
    clusters = []
    
    
    # KMeans
    tic = time.time()
    clusters.append(kmeans_skl(X, n_clusts = 3))
    print("\t running time for k-means(sec):", time.time() - tic)
    
    # Spectral Clustering
    tic = time.time()
    clusters.append(spectral_skl(X, n_clusts = 3))
    print("\t running time for spectral clustering(sec):", time.time() - tic)
    
    # Leiden
    tic = time.time()
    clusters.append(leiden(X, resolution = 0.01))
    print("\t running time for leiden algorithm(sec):", time.time() - tic)

    # ============
    # Plot clustering results
    # ============
    row_axes = axes[i_dataset]
    
    for i, ax in enumerate(row_axes.flatten()):
        curr_cluster = clusters[i]
        if plot_title:
            curr_title = '{}'.format(titles[i])
        else:
            curr_title = None
            
        ax.scatter(X[:,0], X[:,1], c=curr_cluster, s = 0.5)
        ax.set_title(curr_title)
    


    plot_title=False
fig.tight_layout()

"""## Additional reading:
* **SpectralNet** is a fast implementation of spectral clustering algorithm with deep neural network. It has the advantage of spectral clustering and also scales well to large datasets: https://arxiv.org/abs/1801.01587  

* **Quora discussion** about the pros and cons of spectral clustering and k-means: https://www.quora.com/What-are-the-advantages-of-spectral-clustering-over-k-means-clustering

"""