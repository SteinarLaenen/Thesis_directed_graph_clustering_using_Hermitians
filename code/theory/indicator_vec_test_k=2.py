from __future__ import division
import numpy as np
from tqdm import tqdm

corresponding = ["S", "T", "O"]
combinations = [(i, j) for i in corresponding for j in corresponding]
# possible values for the indicator vector
nf = 1/np.sqrt(2)
options = np.array([1, -1, 1j, -1j ,
                    -1 + 1j, -1 - 1j, 1 - 1j, 1 + 1j],
#                    nf*(-1 + 1j), nf*(-1 - 1j), nf*(1 - 1j), nf*(1 + 1j)],
                   dtype=complex)


all_zeros_tuple = (0,0,0,0,0,0,0,0,0)


# every single possible combination of indicator vectors
possible_indicator_vecs = [(m, n, p) for m in options for n in options for p in options]

good_tuples = [(0.0, 4.0, 0.0, -4.0, 0.0, 4.0, 0.0, -4.0, 0.0),
                (0.0, 4.0, 0.0, -4.0, 0.0, -4.0, 0.0, 4.0, 0.0),
                (0.0, -4.0, 0.0, 4.0, 0.0, 4.0, 0.0, -4.0, 0.0),
                (0.0, -4.0, 0.0, 4.0, 0.0, -4.0, 0.0, 4.0, 0.0),
                (0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0),
                (0.0, 2.0, 0.0, -2.0, 0.0, -2.0, 0.0, 2.0, 0.0),
                (0.0, -2.0, 0.0, 2.0, 0.0, 2.0, 0.0, -2.0, 0.0),
                (0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0),
                (0.0, 4.0, 4.0, -4.0, 0.0, 0.0, -4.0, 0.0, 0.0), ]


good_vectors = []

for indicator in possible_indicator_vecs:
    # set the indicator vector
    indicator_vector = np.array([indicator[0], indicator[1], indicator[2]],
                                dtype=complex)

    products = []
    for cluster1, num1 in zip(corresponding, indicator_vector):
        for cluster2, num2 in zip(corresponding, indicator_vector):
            diff = (np.conj(num1)*num2 - np.conj(num2)*num1)*1j
            products.append((int(np.real(diff))))

            

    good_vectors.append((indicator, np.array(products)))

print(len(good_vectors))

q = 0

for indic1, product1 in tqdm(good_vectors):
    for indic2, product2 in tqdm(good_vectors):
        lin_product1 = indic1[0]*np.conj(indic2[0]) + indic1[1]*np.conj(indic2[1]) + indic1[2]*np.conj(indic2[2])
        addition = product1 + product2
        if lin_product1 == 0 and\
        addition[1] == -1*addition[3] and addition[2] == -1*addition[6] and addition[5] == -1* addition[7] \
           and tuple(addition) != all_zeros_tuple:
           # indic1[0]*np.conj(indic1[0]) == indic2[0]*np.conj(indic2[0]) and \
           # indic1[1]*np.conj(indic1[1]) == indic2[1]*np.conj(indic2[1]) and \
           # indic1[2]*np.conj(indic1[2]) == indic2[2]*np.conj(indic2[2]):
            # check if norm is equal for each entry such that we get equal volume

#            print(indic1[0]*np.conj(indic1[0]), indic2[0]*np.conj(indic2[0]))
#            print(indic1[1]*np.conj(indic1[1]), indic2[1]*np.conj(indic2[1]))

#            print("found")
#            print('add', addition)
            print(indic1)
            print(indic2)
            print(indic3)
            print(addition)
#            print(product1)
#            print(product2)
    
        
