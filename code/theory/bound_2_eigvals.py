from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg

def main():
    # define omega as third square root of unity
    omega = (-1 + np.sqrt(3)*1j)/2

    # number of vertices in the graph


    # probability of edge between two sets
    p = 1
    
    # direction
    q = 0.8

    # size of each cluster (all equal size)
    n_cluster = 100

    n_cluster_2 = 1000

    n = n_cluster*2 + n_cluster_2

    # cluster labels
    cluster_labels = [0]*n_cluster + [1]*n_cluster + [2]*n_cluster_2

    Herm_A = create_herm_adjacency_matrix(n, p, q, cluster_labels, 1j)

    D = np.identity(n)*np.sum(np.abs(Herm_A), axis=0)

    norm = linalg.fractional_matrix_power(D, -0.5)
    
    Laplacian = D - Herm_A

    norm_Laplacian = np.matmul(np.matmul(norm, Laplacian), norm)

    # eigvects[:,i] is the eigenvector corresponding to the eigenvalue eigvals[i]
    eigvals, eigvects = np.linalg.eig(norm_Laplacian)

#    plt.hist(eigvals, bins=np.arange(min(eigvals), max(eigvals), 5,dtype=int))
#    plt.show()

    indices = eigvals.argsort() # find the 'n_clusters' isolated eigenvalues
    print(np.round(np.real(np.sort(eigvals)[:5]),decimals=2))

    lamda1 = np.real(eigvals[indices[0]])
    lamda2 = np.real(eigvals[indices[1]])
#    print('second eigenvalue', lamda2)

#    print("lowerbound (lambda1 + lambda2) = ", lamda1 + lamda2)
#    print("lowerbound (lambda1) = ", lamda1)

    # compute the net flow
    net_flow_ST = compute_net_flow(Herm_A, 0, 1, n, cluster_labels)
    prop_vol = 2 - (2*net_flow_ST)/(sum(sum(D)))
    print("value proportional to net flow (2 - 2*(NF(S,T)/vol(V))):", prop_vol)#2 - 2*net_flow_ST/sum(sum(D)))
    print("1 over gamma", lamda2/prop_vol)#2 - 2*net_flow_ST/sum(sum(D)))
    print((lamda1+lamda2)/lamda2)


    # plotting eigenvectors
    d_norm = linalg.fractional_matrix_power(D, -0.5)
    d_norm = np.sum(d_norm, axis=0)
    
    top1_eigvec = d_norm*eigvects[:,indices[0]]

    x = np.real(top1_eigvec)
    y = np.imag(top1_eigvec)

    plt.scatter(x, y, s=150, color='r', label='1st eigvector')
#    plt.show()

    # create and plot the indicator vectors

    nf = 1/np.sqrt(2)
    # create first indicator vector
    indic1 = np.ones(n, dtype=complex)
    for i in range(n):
        if cluster_labels[i] == 0:
            indic1[i] = nf*(1 - 1j)
        elif cluster_labels[i] == 1:
            indic1[i] = nf*(-1 -1j)
        elif cluster_labels[i] == 2:
            indic1[i] = 0


    # create second indicator vector
    indic2 = np.ones(n, dtype=complex)
    for i in range(n):
        if cluster_labels[i] == 0:
            indic2[i] = nf*(-1 -1j)
        elif cluster_labels[i] == 1:
            indic2[i] = 0
        elif cluster_labels[i] == 2:
            indic2[i] = nf*(1 - 1j)

    indic3 = np.ones(n, dtype=complex)
    for i in range(n):
        if cluster_labels[i] == 0:
            indic2[i] = 0
        elif cluster_labels[i] == 1:
            indic2[i] = nf*(-1 -1j)
        elif cluster_labels[i] == 2:
            indic2[i] = nf*(-1 - 1j)

    # construct root of unity vector, whose quadratic sum should be zero
    quadform = np.ones(n, dtype=complex)
    for i in range(n):
        if cluster_labels[i] == 0:
            # third root of unity
            quadform[i] = omega
        elif cluster_labels[i] == 1:
            # just one
            quadform[i] = 1
        elif cluster_labels[i] == 2:
            # omega squared
            quadform[i] = omega**2

#    quadform = quadform/np.sqrt(np.dot(np.conj(quadform), quadform))
    

    summand_quad = 0
    summand_eigval = 0
    for i in range(n):
        for j in range(n):
#            if cluster_labels[i] != cluster_labels[j]:
#                value = quadform[i] - omega*quadform[j]
                value_quad = np.conj(quadform)[i]*Laplacian[i,j]*quadform[j]
                value_eigval = np.conj(top1_eigvec)[i]*Laplacian[i,j]*top1_eigvec[j]
                summand_quad += value_quad
                summand_eigval += value_eigval

#    print('manual quad', np.real(summand_quad))
#    print('manual eigval', np.real(summand_eigval))
            
#    print('quadratic sum:', np.real(quad_sum))
#    print('firsteigval sum:', np.real(firsteigval_sum))

    d_norm = linalg.fractional_matrix_power(D, 0.5)
    d_norm = np.sum(d_norm, axis=0)



    indic3 = (d_norm*indic3)/np.sum(d_norm)
    indic2 = (d_norm*indic2)/np.sum(d_norm)
    indic1 = (d_norm*indic1)/np.sum(d_norm)
    quadform = (d_norm*quadform)/np.sum(d_norm)

#    print(np.matmul(np.conj(indic1.T), indic1))

    quadform_x = np.real(quadform)
    quadform_y = np.imag(quadform)

    indic_1_x = np.real(indic1)
    indic_1_y = np.imag(indic1)
    
    indic_2_x = np.real(indic2)
    indic_2_y = np.imag(indic2)

    indic_3_x = np.real(indic3)
    indic_3_y = np.imag(indic3)

    # plt.scatter(indic_1_x, indic_1_y, s=130, color='b', label='indic1')
    
    # plt.scatter(indic_2_x, indic_2_y, s=150, color='g', label='indic2')
    # plt.scatter(indic_3_x, indic_3_y, s=170, color='y', label='indic3')
#    plt.scatter(omega, quadform_y, s=150, color='purple', label='quadform')
    plt.scatter(quadform_x, quadform_y, s=150, color='purple', label='quadform')
    plt.scatter(quadform_x, quadform_y, s=150, color='purple', label='quadform')
    plt.legend()
    plt.show()
            
def compute_net_flow(Herm, source, destination, n, cluster_labels):
    W_S_D = 0
    for i in np.arange(n):
        for j in np.arange(i,n):
            if cluster_labels[i] == source and cluster_labels[j] == destination:
                W_S_D += np.imag(Herm[i,j])

    W_D_S = 0
    for i in np.arange(n):
        for j in np.arange(i,n):
            if cluster_labels[i] == destination and cluster_labels[j] == source:
                W_D_S += np.imag(Herm[i,j])

    return W_S_D - W_D_S

def create_herm_adjacency_matrix(n, p, q, cluster_labels, omega):
    real_adjacency_matrix = np.zeros((n, n), dtype=complex)

    # construct 
    for i, clusnum_1 in zip(np.arange(n), cluster_labels):
        for j, clusnum_2 in zip(np.arange(i,n), cluster_labels[i:]):
              if np.random.rand() < p and not i == j:
                if clusnum_1 == clusnum_2 and np.random.rand() < 0.5:
#                   if edges are in the same cluster, with uniform random probability set edges
                    if np.random.rand() < 0.5:
                        real_adjacency_matrix[i, j] = 1
                    else:
                        real_adjacency_matrix[j, i] = 1

#               edges between clusters 
                if clusnum_1 == 0 and clusnum_2 == 1:
                    if np.random.rand() < q:
                        real_adjacency_matrix[i, j] = 1
                    else:
                        real_adjacency_matrix[j, i] = 1

#               edges between cluster 0 and 2
                if clusnum_1 == 0 and clusnum_2 == 2:
                    if np.random.rand() < 1-q:
                        real_adjacency_matrix[i, j] = 1
                    else:
                        real_adjacency_matrix[j, i] = 1

                # edges between cluster 1 and 2
                if clusnum_1 == 1 and clusnum_2 == 2:
                    if np.random.rand() < q:
                        real_adjacency_matrix[i, j] = 1
                    else:
                        real_adjacency_matrix[j, i] = 1

    imaginary_adj_matrix = real_adjacency_matrix*omega
    return imaginary_adj_matrix + np.conj(imaginary_adj_matrix.T) #
#    return real_adjacency_matrix*1j - real_adjacency_matrix.T*1j 
    




if __name__ == '__main__':
    main()
            

