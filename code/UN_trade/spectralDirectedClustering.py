#!/usr/bin/env python
#coding:utf-8
"""
  Purpose:  Spectral clustering algorithm
  Created:  12/02/2017
"""

import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

#----------------------------------------------------------------------
def DirectedSpectralClustering(n_clusters, matrix, matrix_name, n_eigenspace):
    """"""
    eigvals, eigvects = np.linalg.eig(matrix) # eigvects[:,i] is the eigenvector corresponding to the eigenvalue eigvals[i]

    if matrix_name in ["UnnormalizedLaplacianMatrix", "LaplacianMatrix"]:
#        print(eigvals)
        indices = np.abs(np.real(eigvals)).argsort()[-n_eigenspace:] # find the 'n_clusters' isolated eigenvalues
#        print(eigvals.argsort()[:n_eigenspace])
#        print(eigvals.argsort()[-n_eigenspace:])
    else:
        raise ValueError("Unknown matrix name")
    
    W = eigvects[:,indices]
    W = np.concatenate((W.real, W.imag), axis = 1)
    kmeans = KMeans(n_clusters=n_clusters).fit(W) # kmeans
    return kmeans.labels_, eigvals, eigvects, W
