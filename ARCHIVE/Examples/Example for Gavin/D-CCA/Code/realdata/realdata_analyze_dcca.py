# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 09:42:51 2017

@author: HShu
"""

import numpy as np
import dcca

#Methylation with SWISS <= 90
Y_1 = np.loadtxt('Expression.txt')
Y_2 = np.loadtxt('Methylation_below90.txt')

X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor, theta  = dcca.dCCA(Y_1, Y_2,method='ED')

ccor, theta,r_1_hat,r_2_hat

np.savetxt('X_1_hat_ours_below90.txt', X_1_hat, delimiter=' ')
np.savetxt('X_2_hat_ours_below90.txt', X_2_hat, delimiter=' ')
np.savetxt('C_1_hat_ours_below90.txt', C_1_hat, delimiter=' ')
np.savetxt('C_2_hat_ours_below90.txt', C_2_hat, delimiter=' ')
np.savetxt('D_1_hat_ours_below90.txt', D_1_hat, delimiter=' ')
np.savetxt('D_2_hat_ours_below90.txt', D_2_hat, delimiter=' ')

#Methylation with SWISS > 90
Y_2 = np.loadtxt('Methylation_above90.txt')

X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor, theta  = dcca.dCCA(Y_1, Y_2,method='ED')

ccor, theta,r_1_hat,r_2_hat

np.savetxt('X_1_hat_ours_above90.txt', X_1_hat, delimiter=' ')
np.savetxt('X_2_hat_ours_above90.txt', X_2_hat, delimiter=' ')
np.savetxt('C_1_hat_ours_above90.txt', C_1_hat, delimiter=' ')
np.savetxt('C_2_hat_ours_above90.txt', C_2_hat, delimiter=' ')
np.savetxt('D_1_hat_ours_above90.txt', D_1_hat, delimiter=' ')
np.savetxt('D_2_hat_ours_above90.txt', D_2_hat, delimiter=' ')

