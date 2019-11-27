from __future__ import division
import numpy as np
import scipy.linalg as linalg

B = np.array([[ 1 +1j,  1 -1j, 0],
              [ 1 +1j, 0,  1 -1j],
              [0,  1 -1j,  1 +1j]], dtype=complex)


A = np.array([[ 0, 1, 0],
              [ 0, 0,  0],
              [1,  1,  0]], dtype=complex)

# total degree i.e. inner+outer
D = np.identity(4)*3#*np.sum(np.abs(A + A.T), axis=0)

Lap4 = np.array([[3,  2j, -1j, 1j],
                [-2j, 3,  -1j, 1j],
                [1j,  1j,  3, 1j],
                [-1j, -1j, -1j, 3]], dtype=complex)

Lap3 = np.array([[2,  1j, -1j],
                [-1j, 2,  -1j],
                [1j,  1j,  2]], dtype=complex)

D_inv_sqrt = linalg.fractional_matrix_power(D, -0.5)
D_sqrt = linalg.fractional_matrix_power(D, 0.5)
norm_Lap = np.matmul(np.matmul(D_inv_sqrt, Lap4), D_inv_sqrt)
#norm_Lap = Lap4

sqrt2 = np.sqrt(2)
x = (1/sqrt2)*np.array([[-1 + 1j, -1 + 1j, -1 + 1j, 1 - 1j]],
             dtype=complex)
x = np.matmul(x, D_sqrt)
print('x1:',x)
x_conj_T = np.conj(x).T
normalization_factor = np.sqrt(np.matmul(x, x_conj_T))
x = x / normalization_factor
print('x2',x)
x_conj_T = x_conj_T / normalization_factor
print("should be 12 bro:",np.matmul(x, x_conj_T))
Rayleigh = np.matmul(np.matmul(x, norm_Lap), x_conj_T)/np.matmul(x, x_conj_T)


print(x)
print(Rayleigh)
print("potential:", Rayleigh*12 - 12) #- np.matmul(x, x_conj_T))



