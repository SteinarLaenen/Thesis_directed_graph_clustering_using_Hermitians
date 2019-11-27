from __future__ import division
import csv
import json
import math
import copy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import permutation_generator as pg
from spectralDirectedClustering import DirectedSpectralClustering
from spectralClustering import SpectralClustering
from matrices import *
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from tqdm import tqdm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


rcParams = {
    'axes.labelsize': 35,
#    'suptitle.fontsize': ,
    'font.size': 32,
    'lines.linewidth': 5,
    'lines.marker': '.',
    'lines.markeredgewidth': 20,
    'legend.fontsize': 50,
    'xtick.labelsize': 35,
    'ytick.labelsize': 35,
    'figure.figsize': [50, 23],
    'grid.linewidth': 0.0,
    'figure.subplot.hspace': 0.08,
    'figure.subplot.wspace': 0.08
}

plt.rcParams.update(rcParams)

def track_average(M, runs, top_pairs, n_clusters, spectral_matrix, matrixtype='norm', n_eigvects=2):
    vol_average = np.zeros((runs, top_pairs))
    size_average = np.zeros((runs, top_pairs))
    CI_average = np.zeros((runs, top_pairs))
    
    for run in tqdm(range(runs)):
        spectral_labels = spectral_clustering(spectral_matrix,
                                                   n_clusters, matrixtype, n_eigvects=n_eigvects)

        CI_average[run,:], size_average[run,:], vol_average[run,:] =  \
                        compute_top_pair(M, spectral_labels, top_pairs)

    
    return spectral_labels, CI_average, size_average, vol_average

def plot_topkth_pair(M, spectral_labels, top_pairs, spectral_matrix, to_plot, kth_pair, axes, axes_pos, n_clusters, matrix_name):
    top_CI, top_size, top_vol = compute_top_pair(M,
                                                 spectral_labels, top_pairs,
                                                 top_pair=True)


    print(axes_pos)
    if to_plot == "size":
        kth_best = top_size[kth_pair]
    elif to_plot == "vol":
        kth_best = top_vol[kth_pair]
    elif to_plot == "CI":
        kth_best = top_CI[kth_pair]
    else:
        raise ValueError("unkown metric to sort the top clusters by")

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    #world = gpd.read_file(shape_file)
    # exclude antartica
    world = world[(world.name!="Antarctica")]
#    world.plot(ax=axes[axes_pos[0],axes_pos[1]], color='w')
    world.plot(ax=axes, color='w')
    
    world = plot_clusters(spectral_labels, n_clusters,
                          cl_to_visualize=[kth_best[1], kth_best[2]],
                          sourcedest=kth_best[3])


    worldplot = world.plot(column='cluster',# color=('blue', 'yellow'),
#                           ax=axes[axes_pos[0],axes_pos[1]])#, color=('blue', 'yellow', 'grey'))
                           ax=axes)#, color=('blue', 'yellow', 'grey'))
    if axes_pos[0]==0 and axes_pos[1]==0:
        legend_elements = [Patch(facecolor='red', edgecolor='w',
                                 label='destination'),
                           Patch(facecolor='grey', edgecolor='w',
                                 label='source')]
        
#        axes[axes_pos[0],axes_pos[1]].legend(handles=legend_elements)
#        axes.legend(handles=legend_elements)




    # print(to_plot, kth_best[0][0])
    # print matrix_name + "Trade Flow (TF) = %.4E" % kth_best[0][0]
    # $\mathsf{TF}^{\mathsf{norm}}$ =
    axes.text(-190, -40, "%.2E" % kth_best[0][0]\
                    #   + "\nCI_vol = " + str(round(kth_best[0][2], 4))) \
                       + "\n%.2f" % kth_best[0][1], fontsize=160)
#    axes.set_title(matrix_name)
    axes.get_xaxis().set_visible(False)
    axes.get_yaxis().set_visible(False)

def sort_countries_net_flow(adjacency_matrix):
    
    net_flows = np.sum(adjacency_matrix - adjacency_matrix.T, axis=1)
    country_flow = [(flow, id_to_country_dic[matrixid_to_id[i]]) for i,flow in enumerate(net_flows)]
    country_flow = sorted(country_flow, key=lambda x: x[0])

    cluster_above_0 = []
    cluster_below_0 = []
    for country in country_flow:
        if country[0] > 0:
            cluster_above_0.append(country[1])
        else:
            cluster_below_0.append(country[1])

    spectral_labels = np.ones(len(net_flows))*(net_flows > 0)


    cluster1_idx = [j for j in range(len(spectral_labels)) if spectral_labels[j]==0]
    cluster2_idx = [j for j in range(len(spectral_labels)) if spectral_labels[j]==1]

    W_cl1_to_cl2 = 0
    W_cl2_to_cl1 = 0

    for i in cluster1_idx:
        for j in cluster2_idx:
            W_cl1_to_cl2 += adjacency_matrix[i,j]
            W_cl2_to_cl1 += adjacency_matrix[j,i]

    imbalance = W_cl1_to_cl2 - W_cl2_to_cl1

    if int(W_cl1_to_cl2) == 0 and int(W_cl2_to_cl1) == 0:
        cut_imbalance = 0
    else:
        cut_imbalance = 0.5*abs((imbalance/(W_cl1_to_cl2+W_cl2_to_cl1)))

    if imbalance >= 0:
        source_destination = ("source", "destination")
    elif imbalance < 0:
        source_destination = ("destination", "source")

    # compute volume 1 (sum of all edge values going in and out)
    vol_1 = 0
    for i in cluster1_idx:
        vol_1 += np.sum((adjacency_matrix + adjacency_matrix.T)[:, i])

    # compute volume 2
    vol_2 = 0
    for i in cluster2_idx:
        vol_2 += np.sum((adjacency_matrix + adjacency_matrix.T)[:, i])

    cut_imbalance_vol = cut_imbalance * min([vol_1, vol_2])
#            cut_imbalance_vol = cut_imbalance * min([vol_1/len(cluster1_idx), vol_2/len(cluster2_idx)])

    cut_imbalance_size = cut_imbalance * min([len(cluster1_idx), len(cluster2_idx)])

    cut_imbalance_CI = cut_imbalance / (len(cluster1_idx) + len(cluster2_idx))

    print 'CI:', cut_imbalance_CI
    print 'CI_vol:', cut_imbalance_vol
    print 'CI_size:', cut_imbalance_size


    return cluster_above_0, cluster_below_0, spectral_labels
    
    
def plot_top_clusters():
    """
    """
    for n_clus_vect in [(20,10)]:
        runs = 1
        n_clusters = n_clus_vect[0]
        n_eigvects = n_clus_vect[1]
        top_pairs = 1

        M = create_adjacency_matrix(import_export_dic, id_to_matrixid,
                                    matrixid_to_id, len(country_codes_withdata), matrixtype='assym')

        cluster1, cluster2, spectral_labels_sort = sort_countries_net_flow(M)
        spectral_labels_sort = np.array([int(i) for i in spectral_labels_sort])

        DISG_L = normalize(np.matmul(M.T, M))
        DISG_R = normalize(np.matmul(M, M.T))
        BI_SYM = np.matmul(M.T, M) + np.matmul(M, M.T)
        DD_SYM = normalize(BI_SYM)


        Herm = create_adjacency_matrix(import_export_dic, id_to_matrixid,
                                       matrixid_to_id, len(country_codes_withdata))
        Herm_RW = normalize(Herm)

        spectral_labels_Herm, CI_average_Herm, size_average_Herm, vol_average_Herm =  \
        track_average(M, runs, top_pairs, n_clusters, Herm, matrixtype='herm', n_eigvects=n_eigvects)

        # Herm_RW
        spectral_labels_Herm_RW, CI_average_Herm_RW, size_average_Herm_RW, vol_average_Herm_RW =  \
        track_average(M, runs, top_pairs, n_clusters, Herm_RW, matrixtype='herm', n_eigvects=n_eigvects)
        # DISG_R
        spectral_labels_DISG_R, CI_average_DISG_R, size_average_DISG_R, vol_average_DISG_R =  \
        track_average(M, runs, top_pairs, n_clusters, DISG_R, matrixtype='norm', n_eigvects=n_eigvects)
        # DISG_L
        spectral_labels_DISG_L, CI_average_DISG_L, size_average_DISG_L, vol_average_DISG_L =  \
        track_average(M, runs, top_pairs, n_clusters, DISG_L, matrixtype='norm', n_eigvects=n_eigvects)

        # BI_SYM
        spectral_labels_BI_SYM, CI_average_BI_SYM, size_average_BI_SYM, vol_average_BI_SYM =  \
        track_average(M, runs, top_pairs, n_clusters, BI_SYM, matrixtype='norm', n_eigvects=n_eigvects)

        # DD_SYM
        spectral_labels_DD_SYM, CI_average_DD_SYM, size_average_DD_SYM, vol_average_DD_SYM =  \
        track_average(M, runs, top_pairs, n_clusters, DD_SYM, matrixtype='norm', n_eigvects=n_eigvects)


        # plt.plot(np.mean(vol_average_Herm, axis=0), label="Herm")
        # plt.plot(np.mean(vol_average_Herm_RW, axis=0), label="Herm_RW")
        # plt.plot(np.mean(vol_average_DISG_R, axis=0), label="DISG_R")
        # plt.plot(np.mean(vol_average_DISG_L, axis=0), label="DISG_L")
        # plt.plot(np.mean(vol_average_BI_SYM, axis=0), label="BI_SYM")
        # plt.plot(np.mean(vol_average_DD_SYM, axis=0), label="DD_SYM")

        # plt.legend()
        # plt.show()

        # plt.plot(np.mean(size_average_Herm, axis=0), label="Herm")
        # plt.plot(np.mean(size_average_Herm_RW, axis=0), label="Herm_RW")
        # plt.plot(np.mean(size_average_DISG_R, axis=0), label="DISG_R")
        # plt.plot(np.mean(size_average_DISG_L, axis=0), label="DISG_L")
        # plt.plot(np.mean(size_average_BI_SYM, axis=0), label="BI_SYM")
        # plt.plot(np.mean(size_average_DD_SYM, axis=0), label="DD_SYM")

        # plt.legend()
        # plt.show()

        to_plot = "CI"

        for top_i in tqdm([-1, -2, -3]):
    #        plt.title("legend graph")
            fig, axes = plt.subplots(ncols=1, nrows=1, sharex=True, sharey=True)
            fig.subplots_adjust(0,0,1,1)
            #fig.suptitle('Trade Flow (TF) between the k=2 clusters after clustering on world wood trade network in year = ' + year, fontsize=30)
#            fig.suptitle('top pair #' + str(-1*top_i) +  ' sorted by TF_norm'  # + to_plot \
                         # + ', with k=' + str(n_clusters) + ' (mineral fuels HS 27 trade in year = ' + year + ')'
                         # , fontsize=30)
            
            plot_topkth_pair(M, spectral_labels_Herm_RW, top_pairs,
                             Herm_RW, to_plot, top_i, axes,
                             [0,0], n_clusters, "Herm_RW")

            # plot_topkth_pair(M, spectral_labels_Herm, top_pairs,
            #                  Herm, to_plot, top_i, axes,
            #                  [1,0], n_clusters, "Herm")

            # plot_topkth_pair(M, spectral_labels_DISG_R, top_pairs,
            #                  Herm_RW, to_plot, top_i, axes,
            #                  [0,1], n_clusters, "DISG_R")

            # plot_topkth_pair(M, spectral_labels_DISG_L, top_pairs,
            #                  Herm_RW, to_plot, top_i, axes,
            #                  [1,1], n_clusters, "DISG_L")

            # plot_topkth_pair(M, spectral_labels_BI_SYM, top_pairs,
            #                  Herm_RW, to_plot, top_i, axes,
            #                  [0,2], n_clusters, "BI_SYM")

            # plot_topkth_pair(M, spectral_labels_DD_SYM, top_pairs,
            #                  Herm_RW, to_plot, top_i, axes,
            #                  [0,1], n_clusters, "DD_SYM")
            
            # plot_topkth_pair(M, spectral_labels_sort, top_pairs,
            #                  Herm_RW, to_plot, top_i, axes,
            #                  [2,1], n_clusters, "Sort Algo")

#            axes[2,2].set_visible(False)
#            axes[2,0].set_visible(False)

#            plt.show()
#            exit()

#            plt.savefig('k_2_sorting_algos.pdf')

            plt.savefig('mineral27_net_trade_flow_per_symHerm_RW_k='+str(n_clusters)+'_' + year + '_top_pair=' + str(-1*top_i) + '.pdf')
            plt.close('all')

def track_average_klist(k_list, n_eigvects_list, M, runs, top_pairs, n_clusters, spectral_matrix, matrixtype='norm'):

    CI_k_average = np.zeros((len(k_list), top_pairs))
    vol_k_average = np.zeros((len(k_list), top_pairs))
    size_k_average = np.zeros((len(k_list), top_pairs))

    for i, k_num in enumerate(k_list):

        spectral_labels, CI_average, size_average, vol_average =  \
        track_average(M, runs, top_pairs, k_num, spectral_matrix, matrixtype=matrixtype, n_eigvects=n_eigvects_list[i])

        CI_k_average[i,:] = np.mean(CI_average, axis=0)
        vol_k_average[i,:] = np.mean(vol_average, axis=0)
        size_k_average[i,:] = np.mean(size_average, axis=0)


    return CI_k_average, vol_k_average, size_k_average
    
        
def plot_cut_imbalance_plots():
    """
    """
    runs = 100
    n_clusters = 20
    n_eigvects = 10
    top_pairs = 10

    M = create_adjacency_matrix(import_export_dic, id_to_matrixid,
                                matrixid_to_id, len(country_codes_withdata), matrixtype='assym')

    DISG_L = normalize(np.matmul(M.T, M))
    DISG_R = normalize(np.matmul(M, M.T))
    BI_SYM = np.matmul(M.T, M) + np.matmul(M, M.T)
    DD_SYM = normalize(BI_SYM)


    #### HERMITIAN MATRIX CONSTRUCTION
    Herm = create_adjacency_matrix(import_export_dic, id_to_matrixid,
                                   matrixid_to_id, len(country_codes_withdata))
    Herm_RW = normalize(Herm)

    k_list = [5, 10, 20]
    n_eigvects_list = [2, 5, 10]



    CI_k_average_Herm, vol_k_average_Herm, size_k_average_Herm = \
                                        track_average_klist(k_list, n_eigvects_list,
                                                            M, runs, top_pairs, n_clusters,
                                                            Herm, matrixtype='herm')

    CI_k_average_Herm_RW, vol_k_average_Herm_RW, size_k_average_Herm_RW = \
                                        track_average_klist(k_list, n_eigvects_list,
                                                            M, runs, top_pairs, n_clusters,
                                                            Herm_RW, matrixtype='herm')

    CI_k_average_DISG_R, vol_k_average_DISG_R, size_k_average_DISG_R = \
                                        track_average_klist(k_list, n_eigvects_list,
                                                            M, runs, top_pairs, n_clusters,
                                                            DISG_R, matrixtype='norm')

    CI_k_average_DISG_L, vol_k_average_DISG_L, size_k_average_DISG_L = \
                                        track_average_klist(k_list, n_eigvects_list,
                                                            M, runs, top_pairs, n_clusters,
                                                            DISG_L, matrixtype='norm')

    CI_k_average_BI_SYM, vol_k_average_BI_SYM, size_k_average_BI_SYM = \
                                        track_average_klist(k_list, n_eigvects_list,
                                                            M, runs, top_pairs, n_clusters,
                                                            BI_SYM, matrixtype='norm')

    CI_k_average_DD_SYM, vol_k_average_DD_SYM, size_k_average_DD_SYM = \
                                        track_average_klist(k_list, n_eigvects_list,
                                                            M, runs, top_pairs, n_clusters,
                                                            DD_SYM, matrixtype='norm')


    for j in range(3):
        if k_list[j] == 2:
            plt.subplot(2, 3, j+1)
            y_pos = np.arange(6)
            colors = ('b', 'g', 'r', 'c', 'm', 'y')
            performance = [CI_k_average_Herm[j,:][0],
                           CI_k_average_Herm_RW[j,:][0],
                           CI_k_average_DISG_R[j,:][0],
                           CI_k_average_DISG_L[j,:][0],
                           CI_k_average_BI_SYM[j,:][0],
                           CI_k_average_DD_SYM[j,:][0]]
            plt.bar(y_pos, performance, align='center', color = colors)
            plt.xticks(y_pos, ('Herm', 'Herm_RW', 'DISG_R', 'DISG_L', 'BI_SYM', 'DD_SYM'))
            plt.ylabel("TF_norm")
    #        plt.xlabel("top pairs")
            plt.title(str("K = " +  str(k_list[j])))
        else:
            plt.subplot(2, 3, j+1)
            plt.plot(CI_k_average_Herm[j,:], label="Herm")
            plt.plot(CI_k_average_Herm_RW[j,:], label="Herm_RW")
            plt.plot(CI_k_average_DISG_R[j,:], label="DISG_R")
            plt.plot(CI_k_average_DISG_L[j,:], label="DISG_L")
            plt.plot(CI_k_average_BI_SYM[j,:], label="BI_SYM")
            plt.plot(CI_k_average_DD_SYM[j,:], label="DD_SYM")
            plt.ylabel("TF_norm")
    #        plt.xlabel("top pairs")
            plt.title(str("K = " +  str(k_list[j])))
            
    plt.legend()


    # for j in range(4):
    #     if k_list[j] == 2:
    #         plt.subplot(3, 4, 4+j+1)
    #         y_pos = np.arange(6)
    #         colors = ('b', 'g', 'r', 'c', 'm', 'y')
    #         performance = [vol_k_average_Herm[j,:][0],
    #                        vol_k_average_Herm_RW[j,:][0],
    #                        vol_k_average_DISG_R[j,:][0],
    #                        vol_k_average_DISG_L[j,:][0],
    #                        vol_k_average_BI_SYM[j,:][0],
    #                        vol_k_average_DD_SYM[j,:][0]]
    #         plt.bar(y_pos, performance, align='center', color=colors)
    #         plt.xticks(y_pos, ('Herm', 'Herm_RW', 'DISG_R', 'DISG_L', 'BI_SYM', 'DD_SYM'))
    #         plt.ylabel("CI_vol")
    # #        plt.xlabel("top pairs")
    #         plt.title(str("K = " +  str(k_list[j])))
    #     else:
    #         plt.subplot(3, 4, 4+j+1)
    #         plt.plot(vol_k_average_Herm[j,:], label="Herm")
    #         plt.plot(vol_k_average_Herm_RW[j,:], label="Herm_RW")
    #         plt.plot(vol_k_average_DISG_R[j,:], label="DISG_R")
    #         plt.plot(vol_k_average_DISG_L[j,:], label="DISG_L")
    #         plt.plot(vol_k_average_BI_SYM[j,:], label="BI_SYM")
    #         plt.plot(vol_k_average_DD_SYM[j,:], label="DD_SYM")
    #         plt.ylabel("CI_vol")
    #        plt.xlabel("top pairs")
    #        plt.title(str("K = " +  str(k_list[j])))
#    plt.legend()

    for j in range(3):
        if k_list[j] == 2:
            plt.subplot(2, 3, 3+j+1)
            y_pos = np.arange(6)
            colors = ('b', 'g', 'r', 'c', 'm', 'y')
            performance = [size_k_average_Herm[j,:][0],
                           size_k_average_Herm_RW[j,:][0],
                           size_k_average_DISG_R[j,:][0],
                           size_k_average_DISG_L[j,:][0],
                           size_k_average_BI_SYM[j,:][0],
                           size_k_average_DD_SYM[j,:][0]]
            plt.bar(y_pos, performance, align='center', color=colors)
            plt.xticks(y_pos, ('Herm', 'Herm_RW', 'DISG_R', 'DISG_L', 'BI_SYM', 'DD_SYM'))
            plt.ylabel("TI_norm")
    #        plt.xlabel("top pairs")
            plt.title(str("K = " +  str(k_list[j])))
        else:
            plt.subplot(2, 3, 3+j+1)
            plt.plot(size_k_average_Herm[j,:], label="Herm")
            plt.plot(size_k_average_Herm_RW[j,:], label="Herm_RW")
            plt.plot(size_k_average_DISG_R[j,:], label="DISG_R")
            plt.plot(size_k_average_DISG_L[j,:], label="DISG_L")
            plt.plot(size_k_average_BI_SYM[j,:], label="BI_SYM")
            plt.plot(size_k_average_DD_SYM[j,:], label="DD_SYM")
            plt.ylabel("TI_norm")
            plt.xlabel("top pairs")
    #        plt.title(str("K = " +  str(k_list[j])))

#    plt.legend()
    plt.show()

            
def create_adjacency_matrix(impexp_dic, id_to_matrixid, matrixid_to_id, n, matrixtype='herm'):
    """creates a Hermitian adjacency matrix for the UN GlobalCom_data"""

    adj_matrix = np.zeros((n,n))
    
#     for element in impexp_dic:
#         country1_id = id_to_matrixid[element[0]]
#         country2_id = id_to_matrixid[element[1]]
#         # positive if export is higher than import
# #        difference = impexp_dic[element][1] - impexp_dic[element][0]
#         if max(impexp_dic[element]) == impexp_dic[element][1]:
#             difference = impexp_dic[element][1]
#         else:
#             difference = -1*impexp_dic[element][0]

        
#         adj_matrix[country1_id, country2_id] = difference

    # pick the maximum out of the two export percentages, or use the difference
    for element in impexp_dic:
        country1_id = id_to_matrixid[element[0]]
        country2_id = id_to_matrixid[element[1]]

        export_1_to_2 = impexp_dic[(element[0], element[1])][1]
        export_2_to_1 = impexp_dic[(element[1], element[0])][1]

        difference = export_1_to_2 - export_2_to_1
        
        # if max([export_1_to_2, export_2_to_1]) == export_1_to_2:
        #     difference = export_1_to_2
        # else:
        #     difference = -1*export_2_to_1
        
        adj_matrix[country1_id, country2_id] = difference

#    adj_matrix = np.exp(adj_matrix)
#    print(np.where(adj_matrix == 'NaN'))

    hermitian = np.ones((n, n), dtype=complex)*adj_matrix*1j
    
    if matrixtype == 'herm':
        return hermitian
    
    elif matrixtype == 'assym':
        return adj_matrix*(adj_matrix > 0)
    else:
        raise ValueError("unknown matrix type to create adjacency matrix")

def spectral_clustering(A, n_clusters, case, n_eigvects=2, lamda=0):
    """ performs spectral clustering and returns some values
    """
    # if A is hermitian, perform the complex valued case
    if case=="herm":
        spectral_labels, eigvals, \
        eigvects, W = DirectedSpectralClustering(n_clusters,
                                                 A,
                                                 "UnnormalizedLaplacianMatrix"
                                                 ,n_eigvects)

    else:
        spectral_labels, eigvals_normal, \
        eigvects_normal, W_normal = SpectralClustering(n_clusters,
                                                       A, "AdjacencyMatrix")

    return spectral_labels

def normalize(adjacency_matrix):
    """
    Return the symmetric normalized matrix
    :param matrix: np.array symetric square matrix
    :return: np.array square matrix
    """
    d = np.sum(np.absolute(adjacency_matrix), axis=0)
    D = np.identity(len(d))*d
    return np.matmul(linalg.fractional_matrix_power(D, -1), adjacency_matrix)

def normalize_and_regularize(adjacency_matrix):
    """
    Return the symmetric unnnormalized Laplacian matrix of a given graph
    :param matrix: np.array symetric square matrix
    :return: np.array square matrix
    """
    d = np.sum(np.absolute(adjacency_matrix), axis=0)
    D = np.identity(len(d))*d
    tau = sum(d)/len(d)
    regularize_matrix = np.identity(len(d))*tau
    adjacency_matrix = adjacency_matrix + regularize_matrix
    return np.matmul(linalg.fractional_matrix_power(D, -1), adjacency_matrix)


# make a list of lists for the clusters
def plot_clusters(spectral_labels, n_clusters, print_groups=False, \
                  cl_to_visualize=-1, sourcedest = (" ", " ")):
    

    # change the visualization if explicit labels are needed
    spectral_labels_cache = copy.deepcopy(spectral_labels)
    if cl_to_visualize == -1:
        # show all clusters
        spectral_labels_cache = spectral_labels
    else:
        # else we just move the other clusters
        for i, label in enumerate(spectral_labels):
            if label not in cl_to_visualize:
                spectral_labels_cache[i] = n_clusters+1


    country_clusters = [[] for i in range(n_clusters + 2)]

    # add countries to the list of lists based on what group they belong to
    for ID, group in enumerate(spectral_labels_cache):
        country_id = matrixid_to_id[ID]
        country_name = id_to_country_dic[country_id]

        country_clusters[group].append(country_to_iso[country_name])

    # rint the groups
    if print_groups==True:
        for group in country_clusters:
            print("NUMBER OF COUNTRIES IN CLUSTER:", len(group), "\n")
            print(group, '\n')

    # plot the world
    shape_file = '/home/steinar/Documents/edinburgh/data/naturalearth/ne_10m_admin_0_countries.shp'

    world = gpd.read_file(shape_file)
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    #world = gpd.read_file(shape_file)
    
    # exclude antartica
    world = world[(world.name!="Antarctica")]

#    clusterlist = [-1]*len(world['ADM0_A3'])
    clusterlist = [-1]*len(world['iso_a3'])
    
#    print(world['ADM0_A3'])

    # go through all countries in the world dataset from gdp

    for i, country in enumerate(zip(world['iso_a3'], world['name'])):
#    for i, country in enumerate(world["iso_a3"]):

        # get ISO and name of country
        ISO = country[0]
        name = country[1]

        #some ISO's don't seem to be correct for some reason (Norway and France)
        if ISO == '-99':
            try:
                ISO = country_to_iso[name]
            except:
                clusterlist[i] = "NO DATA"

        # get the cluster to which the country belongs to in the clustering
        cluster = [j for j in range(n_clusters) if ISO in country_clusters[j]]

        if len(cluster) == 1:
            value = cluster[0]
            clusterlist[i] = sourcedest[cl_to_visualize.index(value)]
#        elif country_to_id_dic[iso_to_country[ISO]] not in country_codes_withdata:
#            clusterlist[i] = "NO DATA"
        else:
            # try:
            #     if country_to_id_dic[iso_to_country[ISO]] not in country_codes_withdata:
            #         clusterlist[i] = "zno data available"
            #     else:
            #         clusterlist[i] = "zrest of clusters"
            # except:
            clusterlist[i] = "zrest of clusters"


    world['cluster'] = clusterlist
    world = world[(world.cluster!="zrest of clusters")]

    return world

def compute_top_pair(adj_matrix, spectral_labels, k, top_pair=False):
    """
    Function that computes the cut imbalance between all pairs of cluster. 
    return both CI_size and CI_vol of the top k pairs
    """

    top_k_vol = []
    top_k_size = []
    top_k_CI = []

    source_destination = ('source', 'destination')

    n_clusters = max(spectral_labels)+1
    
    for cluster1 in range(n_clusters):
        for cluster2 in range(cluster1, n_clusters):
            cluster1_idx = [j for j in range(len(spectral_labels)) if spectral_labels[j]==cluster1]
            cluster2_idx = [j for j in range(len(spectral_labels)) if spectral_labels[j]==cluster2]

            W_cl1_to_cl2 = 0
            W_cl2_to_cl1 = 0

            for i in cluster1_idx:
                for j in cluster2_idx:
                    W_cl1_to_cl2 += adj_matrix[i,j]
                    W_cl2_to_cl1 += adj_matrix[j,i]

            imbalance = W_cl1_to_cl2 - W_cl2_to_cl1
                        
            if int(W_cl1_to_cl2) == 0 and int(W_cl2_to_cl1) == 0:
                cut_imbalance = 0
            else:
                cut_imbalance = 0.5*abs((imbalance/(W_cl1_to_cl2+W_cl2_to_cl1)))
            
            if imbalance >= 0:
                source_destination = ("source", "destination")
            elif imbalance < 0:
                source_destination = ("destination", "source")
                
            # compute volume 1 (sum of all edge values going in and out)
            vol_1 = 0
            for i in cluster1_idx:
                vol_1 += np.sum((adj_matrix + adj_matrix.T)[:, i])

            # compute volume 2
            vol_2 = 0
            for i in cluster2_idx:
                vol_2 += np.sum((adj_matrix + adj_matrix.T)[:, i])

            cut_imbalance_vol = cut_imbalance * min([vol_1, vol_2])
#            cut_imbalance_vol = cut_imbalance * min([vol_1/len(cluster1_idx), vol_2/len(cluster2_idx)])

            cut_imbalance_size = cut_imbalance * min([len(cluster1_idx), len(cluster2_idx)])

#            cut_imbalance_CI = abs(imbalance)
            cut_imbalance_CI = abs(imbalance)/(len(cluster1_idx) + len(cluster2_idx))

            top_k_CI.append(((cut_imbalance_CI, cut_imbalance_size, cut_imbalance_vol),
                             cluster1, cluster2, source_destination))
            top_k_size.append(((cut_imbalance_CI, cut_imbalance_size, cut_imbalance_vol),
                               cluster1, cluster2, source_destination))
            top_k_vol.append(((cut_imbalance_CI, cut_imbalance_size, cut_imbalance_vol),
                              cluster1, cluster2, source_destination))

    top_k_CI_pair = sorted(top_k_CI, key=lambda x: x[0][0])
    top_k_vol_pair = sorted(top_k_vol, key=lambda x: x[0][2])
    top_k_size_pair = sorted(top_k_size, key=lambda x: x[0][1])

    if top_pair:
        return top_k_CI_pair, top_k_size_pair, top_k_vol_pair
    else:
    # add missing data for plotting
        to_add = 0
        if len(top_k_CI_pair) < k:
            to_add = k - len(top_k_CI)

            
        top_k_CI = np.array([0]*to_add + [pair[0][0] for pair in top_k_CI_pair])
        top_k_size = np.array([0]*to_add + [pair[0][1] for pair in top_k_size_pair])
        top_k_vol = np.array([0]*to_add + [pair[0][2] for pair in top_k_vol_pair])

            
        top_k_CI.sort()
        top_k_size.sort()
        top_k_vol.sort()

        return top_k_CI[-k:][::-1],top_k_size[-k:][::-1], top_k_vol[-k:][::-1]
        

if __name__ == "__main__":
    years = range(2001, 2004)# + range(2007, 2010)

    print(years)
    for year_int in tqdm(years):
        ###----------------------DATA PREPROCESSING -----------------------------------###
        shape_file = '/home/steinar/Documents/edinburgh/data/naturalearth/ne_10m_admin_0_countries.shp'
        #continent = 'europe'

        #### save bilateral trade distance in a dictionary ###
        bilateral_dist_dic = {}
        bil_data = pd.io.stata.read_stata('./process_comtrade/dist_cepii.dta')
        average_distance = 0
        n = 0
        max_distance = 0

        for origin, destination, dist in zip(bil_data['iso_o'], bil_data[ 'iso_d'], bil_data['dist']):
            n+=1
            bilateral_dist_dic[(origin, destination)] = int(dist)
            average_distance += dist

            if dist > max_distance:
                max_distance = dist


        average_distance = average_distance/n

        world = gpd.read_file(shape_file)
        world_info = [iso for iso in world["geometry"]]
        world_names = [iso for iso in world["NAME"]]

        #print([i for i in world['ADM0_A3']])
        world_isos = [iso for iso in world["ADM0_A3"]]
        #world_isos = []

        ISO_to_population = {}
        ISO_to_GDP = {}

        for ISO, POP_EST, GDP, POLYGON in zip(world["ISO_A3"], world["POP_EST"], world["GDP_MD_EST"], world["geometry"]):
            ISO_to_population[ISO] = int(POP_EST)
            ISO_to_GDP[ISO] = int(GDP) + 1


        #world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        #world = world[(world.continent == 'Europe')]
        #world_isos = [i for i in world]
        #print(world_isos)


        countries = []
        countries_idx = []
        id_to_country_dic = {}
        country_to_id_dic = {}
        country_codes_withdata = []

        # dictionary mapping country codes to indices in adjacency matrix
        id_to_matrixid = {}
        matrixid_to_id = {}

        # we need to create a list of all country codes that have data available to create a denser
        # adjacency matrix
        country_to_iso = {}
        iso_to_country = {}
        import_export_dic = {}
        country_world_imp_dic = {}
        country_world_exp_dic = {}
        reporters = []
        partners = []


        year = str(year_int)
        data_file = './process_comtrade/oil_27/processed_data_' + year +  '_27.csv'#

        with open(data_file, mode='rb') as clean_file:
            reader = csv.DictReader(clean_file)

            # we read some of the data from the csv file
            for row in reader:
                reporter = row['Reporter']        
                reporter_id = int(row['Reporter Code'])
                reporter_ISO = row['Reporter ISO']

                partner = row['Partner']
                partner_id = int(row['Partner Code'])
                partner_ISO = row['Partner ISO']

                trade_flow_code = int(row['Trade Flow Code'])
                money = int(row['Trade Value (US$)'])

                if partner_id == 0:
                    if trade_flow_code == 1:
                        country_world_imp_dic[reporter_id] = money
                    if trade_flow_code == 2:
                        country_world_exp_dic[reporter_id] = money

                if not(reporter_id == 0 or partner_id == 0):

                    if reporter_id not in country_codes_withdata and (reporter_ISO in world_isos):
                        country_codes_withdata.append(reporter_id)

                    if partner_id not in country_codes_withdata and (partner_ISO in world_isos):
                        country_codes_withdata.append(partner_id)


                    if reporter_id not in reporters:
                        reporters.append(reporter_id)
                    if partner_id not in partners:
                        partners.append(partner_id)

                    # we create our dictionaries 
                    country_to_iso[reporter] = reporter_ISO
                    country_to_iso[partner] = partner_ISO

                    iso_to_country[reporter_ISO] = reporter
                    iso_to_country[partner_ISO] = partner

                    country_to_id_dic[partner] = int(partner_id)
                    country_to_id_dic[reporter] = int(reporter_id)

                    id_to_country_dic[int(partner_id)] = partner
                    id_to_country_dic[int(reporter_id)] = reporter


#        print(len(country_codes_withdata))
#        print("number of partners:", len(partners))
#        print("number of reporters:", len(reporters))

        # we number the country codes from 1 to n, so we can use these new indices
        # for the adjacency matrix

        for i, ids in enumerate(country_codes_withdata):
            id_to_matrixid[ids] = i
            matrixid_to_id[i] = ids

        # dictionary that notes import and export value between each pair of countries,
        # for each reporting country. import is first value, export is second value
        for i in country_codes_withdata:
            for m in country_codes_withdata:
                import_export_dic[(i, m)] = [0,0]

        # plug in values for all the import and export countries
        with open(data_file, mode='rb') as clean_file:
            reader = csv.DictReader(clean_file)

            for row in reader:
                reporter_id = int(row['Reporter Code'])
                partner_id = int(row['Partner Code'])
                trade_flow_code = int(row['Trade Flow Code'])
                money = int(row['Trade Value (US$)'])

                reporter_ISO = row['Reporter ISO']
                partner_ISO = row['Partner ISO']


                if not (reporter_id==0 or partner_id==0) and ((reporter_ISO in world_isos) and (partner_ISO in world_isos)):
                    # import is the first index, export is second
                    if trade_flow_code == 1:
                        import_export_dic[(reporter_id, partner_id)][0] = money#/country_world_imp_dic[reporter_id]
                    if trade_flow_code == 2:
                        import_export_dic[(reporter_id, partner_id)][1] = money#/country_world_exp_dic[reporter_id]


        for pair in import_export_dic:
            original = import_export_dic[(pair[0], pair[1])]
            reverse = import_export_dic[(pair[1], pair[0])]

            new_original = [0,0]
            new_reverse = [0,0]

            # compare import and export of a pair
            if original[0] != 0 and reverse[1] != 0:
                average = (original[0] + reverse[1])/2
                new_original[0] = average
                new_reverse[1] = average
            elif original[0] != 0 and reverse[1] == 0:
                value = original[0]
                new_original[0] = value
                new_reverse[1] = value
            elif original[0] == 0 and reverse[1] != 0:
                value = reverse[1]
                new_original[0] = value
                new_reverse[1] = value
            else:
                if not pair[0] == pair[1]:
                    value = np.random.rand()
                    new_original[0] = value
                    new_reverse[1] = value

            # other way around
            if original[1] != 0 and reverse[0] != 0:
                average = (original[1] + reverse[0])/2
                new_original[1] = average
                new_reverse[0] = average
            elif original[1] != 0 and reverse[0] == 0:
                value = original[1]
                new_original[1] = value
                new_reverse[0] = value
            elif original[1] == 0 and reverse[0] != 0:
                value = reverse[0]
                new_original[1] = value
                new_reverse[0] = value
            else:
                if not pair[0] == pair[1]:
                    value = np.random.rand()*10000
                    new_original[1] = value
                    new_reverse[0] = value

            import_export_dic[(pair[0], pair[1])] = new_original
            import_export_dic[(pair[1], pair[0])] = new_reverse

        sum_total_export_dic = {}
        sum_total_import_dic = {}

        # get total trade after adding random noise to data
        for i in country_codes_withdata:
            sum_total_export_dic[i] = 0
            sum_total_import_dic[i] = 0

        for pair in import_export_dic:
            sum_total_export_dic[pair[0]] += import_export_dic[pair][1]
            sum_total_import_dic[pair[0]] += import_export_dic[pair][0]

        #normalize all the values
        # for pair in import_export_dic:
        #     try:
        #         new_export = (import_export_dic[pair][1] / ISO_to_GDP[country_to_iso[id_to_country_dic[pair[0]]]])
        #         new_import = (import_export_dic[pair][0] / ISO_to_GDP[country_to_iso[id_to_country_dic[pair[1]]]])
        #     except:
        #         new_export = import_export_dic[pair][1]
        #         new_import = import_export_dic[pair][0]
                
        #     import_export_dic[pair] = [new_import, new_export]

        ###----------------------DATA PREPROCESSING FINISHED-------------------###
    
        parser = argparse.ArgumentParser(description='choose what to plot')
        parser.add_argument('--to_plot', type=str, default='top_clusters', dest='toplot')
        args = parser.parse_args()
        if args.toplot == 'top_clusters':
            plot_top_clusters()
        elif args.toplot == 'CI_plots':
            plot_cut_imbalance_plots()
        else:
            raise ValueError("not a valid argument for --to_plot")
    
