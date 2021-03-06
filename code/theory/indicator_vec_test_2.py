from __future__ import division
import numpy as np
from tqdm import tqdm

omega = (-1 + np.sqrt(3)*1j)/2
omega5 = np.cos((2*np.pi)/5) + 1j*np.sin((2*np.pi)/5)
omega6 = np.cos((2*np.pi)/6) + 1j*np.sin((2*np.pi)/6)
omega7 = np.cos((2*np.pi)/7) + 1j*np.sin((2*np.pi)/7)
omega8 = np.cos((2*np.pi)/8) + 1j*np.sin((2*np.pi)/8)
omega9 = np.cos((2*np.pi)/9) + 1j*np.sin((2*np.pi)/9)
omega10 = np.cos((2*np.pi)/10) + 1j*np.sin((2*np.pi)/10)


# possible values for the indicator vector
nf = 1/np.sqrt(2)
options = np.array([1, #omega, omega**2,
#                    1j, -1j, -1,
                    nf*(1 + 1j), nf*(1 - 1j), nf*(-1 -1j), nf*(-1 +1j), 0],
#                    (1 + 1j), (1 - 1j), (-1 -1j), (-1 +1j)],
#                    1 + 1j, 1 - 1j, -1 -1j, -1 + 1j],
                    # omega5, omega5**2, omega5**3, omega5**4,
                    # omega6, omega6**2, omega6**3, omega6**4, omega6**5],
                    # omega7, omega7**2, omega7**3, omega7**4, omega7**5, omega7**6,
                    # omega8, omega8**2, omega8**3, omega8**4, omega8**5, omega8**6, omega8**7,
                    # omega9, omega9**2, omega9**3, omega9**4, omega9**5, omega9**6, omega9**7, omega9**8,
                    # omega10, omega10**2, omega10**3, omega10**4, omega10**5, omega10**6, omega10**7, omega10**8, omega10**9],
                   dtype=complex)




all_zeros_tuple = (0,0,0,0,0,0,0,0,0)
#                   0,0,0,0,0,0,0,0,0)


# every single possible combination of indicator vectors
possible_indicator_vecs = [(m, n, p) for m in options for n in options for p in options for q in options]

good_vectors = []

num1 = nf*(-1+1j)
num2 = omega

print((np.conj(num1)*num2 - num1*np.conj(num2))*1j)


for indicator in tqdm(possible_indicator_vecs):
    # set the indicator vector
    indicator_vector = np.array([indicator[0], indicator[1], indicator[2]],
                                dtype=complex)

    products = []
    for num1 in indicator_vector:
        for num2 in indicator_vector:
#            diff = (np.conj(num1)*num2*((-1 + np.sqrt(3)*1j)/2) + np.conj(num2)*num1*((-1 - np.sqrt(3)*1j)/2))
            diff = (np.conj(num1)*num2 - np.conj(num2)*num1)*1j
            products.append((np.real(diff)))

            

    good_vectors.append((indicator, np.array(products)))


#indicator_vector_lam0 = np.array([-1j, nf*(-1 + 1j), omega**2],
#                                dtype=complex)




indicator_vector_lam0 = np.array([omega6, omega6**2, omega6**3, omega6**4, omega6**5, omega6**6],
                                dtype=complex)
print(indicator_vector_lam0)
product_lam0 = []

for num1 in indicator_vector_lam0:
    for num2 in indicator_vector_lam0:
#        diff = (np.conj(num1)*num2*omega + np.conj(num2)*num1*omega)
        diff = 1j*(np.conj(num1)*num2 - np.conj(num2)*num1)

        print(np.sqrt(3))
        print(diff)
        product_lam0.append((np.real(diff)))

print(product_lam0)
q = 0

exit()

for indic1, product1 in good_vectors:#[(indicator_vector_lam0, product_lam0)]:
#    for indic2, product2 in tqdm(good_vectors):
#        for indic3, product3 in tqdm(good_vectors):
            # lin_product1 = indic1[0]*np.conj(indic2[0])
            # lin_product2 = indic1[1]*np.conj(indic2[1])
            # lin_product3 = indic1[2]*np.conj(indic2[2])
            
            # lin_product_sum1 = lin_product1 + lin_product2 + lin_product3

            addition = product1

#            print(indic1)
            print(addition)



            if addition[1] == -1*addition[3] and tuple(addition) != all_zeros_tuple and addition[1] == 4 and addition[2] == 0 and addition[5] == 0:

    #            print("found")
    #            print('add', addition)
                print(indic1)
#                print(indic2)

                print(addition)
#                print(product1)
#                print(product2)
