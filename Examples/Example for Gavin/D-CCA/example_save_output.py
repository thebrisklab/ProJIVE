# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 09:42:51 2017

@author: HShu
"""

import numpy as np
import dcca

Y_1 = np.loadtxt('Y_1.txt')
Y_2 = np.loadtxt('Y_2.txt')

#Y_k =  X_k + E_k = C_k + D_k + E_k for k=1,2
#where C_1 and C_2 share the same latent factors, but D_1 and D_2 have uncorrelated latent factors.

r_1, r_2, r_12 = 3, 5, 1 # r_k=rank{cov(x_k)}, r_12=rank{cov(x_1,x_2)}; Or see the definition of r_1, r_2 and r_12 in the JASA paper

#Use user-specified r_1, r_2 and r_12
X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor_hat, theta_hat = dcca.dCCA(Y_1, Y_2, r_1=r_1, r_2=r_2, r_12=r_12)

np.savetxt("X1_hat.csv", X_1_hat, delimiter=",")
np.savetxt("X2_hat.csv", X_2_hat, delimiter=",")

np.savetxt("C1_hat.csv", C_1_hat, delimiter=",")
np.savetxt("C2_hat.csv", C_2_hat, delimiter=",")

np.savetxt("D1_hat.csv", D_1_hat, delimiter=",")
np.savetxt("D2_hat.csv", D_2_hat, delimiter=",")

