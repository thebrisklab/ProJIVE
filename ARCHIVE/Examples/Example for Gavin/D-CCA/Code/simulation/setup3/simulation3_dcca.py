# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 09:42:51 2017

@author: HShu
"""

import numpy as np
import math
#import sys
from numpy.linalg import norm
import dcca
from rank_nullspace import nullspace


'''
Begin: Generate simulated data
'''

n =300

# noise standard deviation

p_1 = 900
p_2 = 300    
r_1, r_2, r_12 = 3, 5, 1
###################
# Define V_1, V_2. 
rng0 = np.random.RandomState(0)#(int(sys.argv[1]))

V_1 = np.ones(shape=[p_1,r_1])
V_1[:,0] = V_1[:,0]/norm(V_1[:,0])
for i in range(1,r_1):
    b_v = rng0.normal(0, 1, p_1-i) # coefficient for null space basis
    V_1_i_null = nullspace(V_1[:,:i].T)
    V_1[:,i] = V_1_i_null @ b_v
    V_1[:,i] = V_1[:,i]/norm(V_1[:,i])
    
V_2 = np.ones(shape=[p_2,r_2])
V_2[:,0] = V_2[:,0]/norm(V_2[:,0])
for i in range(1,r_2):
    b_v = rng0.normal(0, 1, p_2-i) # coefficient for null space basis
    V_2_i_null = nullspace(V_2[:,:i].T)
    V_2[:,i] = V_2_i_null @ b_v
    V_2[:,i] = V_2[:,i]/norm(V_2[:,i])
###################
######################################
# Define Lambda matrix
Lambda_1 = np.zeros(shape=(r_1,r_1))
Lambda_1[range(r_1),range(r_1)]=[500,300,100]
#for i in range(r_1):
#    Lambda_1[i,i] = 20*(r_1-i)
    
Lambda_2 = np.zeros(shape=(r_2,r_2))
Lambda_2[range(r_2),range(r_2)]=[500,400,300,200,100]
#for i in range(r_2):
#    Lambda_2[i,i] = 50*(r_2-i)
#####################################
std_noise = 1
###### noise matrices
E_1plus2 = rng0.normal(0, std_noise, (p_1+p_2)*n).reshape([p_1+p_2,n])            
E_1 = E_1plus2[:p_1,:]
E_2 = E_1plus2[p_1:,:]

angle_deg = 45      
 
# angles for each pair of canonical variables
angle =np.array([angle_deg])/180.0*math.pi 

D_theta = np.zeros(shape=(r_1,r_2)) # Let D_theta be a r_1*r_2 matrix
for i in range(r_12):
    D_theta[i,i] = math.cos(angle[i]) 
#####################################


#c_1: the common unit variable of the 1st pair of canonical variables
# discrete uniformly distributed variables with range {-2^0.5,-2^-0.5,0,2^-0.5,2^0.5}
#d_1: independent to c_1. Let it be standard Gaussian.

#rng1 = np.random.RandomState(1001)#(int(sys.argv[1]))
c_1 = rng0.randint(low=-2, high=3, size=n)
c_1_order=np.hstack([*np.where(c_1==-2),*np.where(c_1==-1),*np.where(c_1==0),*np.where(c_1==1),*np.where(c_1==2)])
c_1=c_1[c_1_order]
c_1=c_1/(2**0.5)

d_1 = rng0.normal(0, 1, n)

Z_1_star = rng0.normal(0, 1, r_1*n).reshape([r_1,n])
Z_2_star = rng0.normal(0, 1, r_2*n).reshape([r_2,n])

Z_1_star[0,:n] = (c_1 + d_1*math.tan(angle/2))/((1+math.tan(angle/2)**2)**0.5)
Z_2_star[0,:n] = (c_1 - d_1*math.tan(angle/2))/((1+math.tan(angle/2)**2)**0.5)


X_1 = V_1 @ Lambda_1**0.5 @ Z_1_star
X_2 = V_2 @ Lambda_2**0.5 @ Z_2_star


Y_1 = X_1 + E_1
Y_2 = X_2 + E_2

Sigma_1 = V_1 @ Lambda_1 @ V_1.T # Covariance matrix of X_1
Sigma_2 = V_2 @ Lambda_2 @ V_2.T # Covariance matrix of X_2


Lambda_1_inv_half = np.zeros(shape=(r_1,r_1))
Lambda_1_inv_half[range(r_1),range(r_1)] = Lambda_1.diagonal()**-0.5

Lambda_2_inv_half = np.zeros(shape=(r_2,r_2))
Lambda_2_inv_half[range(r_2),range(r_2)] = Lambda_2.diagonal()**-0.5

U_theta = np.identity(r_1)
V_theta = np.identity(r_2)


Gamma_1 = V_1 @ Lambda_1_inv_half @ U_theta[:,:r_12]
Gamma_2 = V_2 @ Lambda_2_inv_half @ V_theta[:,:r_12]

D_theta_diag = D_theta.diagonal()[:r_12]  #nonzero diagonal elements of D_theta

A_mat = np.diag( 0.5 * (1-  ((1-D_theta_diag)/(1+D_theta_diag))**0.5 ) )  

C_1 = Sigma_1 @ Gamma_1 @ A_mat @ Gamma_1.T @ X_1 + Sigma_1 @ Gamma_1 @ A_mat @ Gamma_2.T @ X_2    
C_2 = Sigma_2 @ Gamma_2 @ A_mat @ Gamma_1.T @ X_1 + Sigma_2 @ Gamma_2 @ A_mat @ Gamma_2.T @ X_2

D_1 = X_1 - C_1
D_2 = X_2 - C_2

'''
End: Generate simulated data
'''
np.savetxt('Y_1.txt', Y_1, delimiter=' ')
np.savetxt('Y_2.txt', Y_2, delimiter=' ')
np.savetxt('X_1.txt', X_1, delimiter=' ')
np.savetxt('X_2.txt', X_2, delimiter=' ')
np.savetxt('C_1.txt', C_1, delimiter=' ')
np.savetxt('C_2.txt', C_2, delimiter=' ')
np.savetxt('D_1.txt', D_1, delimiter=' ')
np.savetxt('D_2.txt', D_2, delimiter=' ')
np.savetxt('E_1.txt', E_1, delimiter=' ')
np.savetxt('E_2.txt', E_2, delimiter=' ')



X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor, theta= dcca.dCCA(Y_1, Y_2, r_1=r_1, r_2=r_2, r_12=r_12)
np.max(np.abs(D_1_hat@D_2_hat.T))

np.savetxt('X_1_hat_oracle_rank.txt', X_1_hat, delimiter=' ')
np.savetxt('X_2_hat_oracle_rank.txt', X_2_hat, delimiter=' ')
np.savetxt('C_1_hat_oracle_rank.txt', C_1_hat, delimiter=' ')
np.savetxt('C_2_hat_oracle_rank.txt', C_2_hat, delimiter=' ')
np.savetxt('D_1_hat_oracle_rank.txt', D_1_hat, delimiter=' ')
np.savetxt('D_2_hat_oracle_rank.txt', D_2_hat, delimiter=' ')


errF_X_1=norm(X_1-X_1_hat)/norm(X_1)
err2_X_1=norm(X_1-X_1_hat,2)/norm(X_1,2)
errF_X_2=norm(X_2-X_2_hat)/norm(X_2)
err2_X_2=norm(X_2-X_2_hat,2)/norm(X_2,2)  

 
errF_C_1=norm(C_1-C_1_hat)/norm(C_1)
err2_C_1=norm(C_1-C_1_hat,2)/norm(C_1,2)
errF_C_2=norm(C_2-C_2_hat)/norm(C_2)
err2_C_2=norm(C_2-C_2_hat,2)/norm(C_2,2)  

errF_D_1=norm(D_1-D_1_hat)/norm(D_1)
err2_D_1=norm(D_1-D_1_hat,2)/norm(D_1,2)
errF_D_2=norm(D_2-D_2_hat)/norm(D_2)
err2_D_2=norm(D_2-D_2_hat,2)/norm(D_2,2)  

out1=[errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2]


X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor,theta  = dcca.dCCA(Y_1, Y_2, method='ED')
np.max(np.abs(D_1_hat@D_2_hat.T))


np.savetxt('X_1_hat_est_rank.txt', X_1_hat, delimiter=' ')
np.savetxt('X_2_hat_est_rank.txt', X_2_hat, delimiter=' ')
np.savetxt('C_1_hat_est_rank.txt', C_1_hat, delimiter=' ')
np.savetxt('C_2_hat_est_rank.txt', C_2_hat, delimiter=' ')
np.savetxt('D_1_hat_est_rank.txt', D_1_hat, delimiter=' ')
np.savetxt('D_2_hat_est_rank.txt', D_2_hat, delimiter=' ')

errF_X_1=norm(X_1-X_1_hat)/norm(X_1)
err2_X_1=norm(X_1-X_1_hat,2)/norm(X_1,2)
errF_X_2=norm(X_2-X_2_hat)/norm(X_2)
err2_X_2=norm(X_2-X_2_hat,2)/norm(X_2,2)  


errF_C_1=norm(C_1-C_1_hat)/norm(C_1)
err2_C_1=norm(C_1-C_1_hat,2)/norm(C_1,2)
errF_C_2=norm(C_2-C_2_hat)/norm(C_2)
err2_C_2=norm(C_2-C_2_hat,2)/norm(C_2,2)  

errF_D_1=norm(D_1-D_1_hat)/norm(D_1)
err2_D_1=norm(D_1-D_1_hat,2)/norm(D_1,2)
errF_D_2=norm(D_2-D_2_hat)/norm(D_2)
err2_D_2=norm(D_2-D_2_hat,2)/norm(D_2,2)  

out2=[errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2]

propF_C_1=norm(C_1)/norm(X_1) # Frobenius norm
prop2_C_1=norm(C_1,2)/norm(X_1,2) # spectral norm
propF_C_2=norm(C_2)/norm(X_2)
prop2_C_2=norm(C_2,2)/norm(X_2,2)

propF_D_1=norm(D_1)/norm(X_1)
prop2_D_1=norm(D_1,2)/norm(X_1,2)
propF_D_2=norm(D_2)/norm(X_2)
prop2_D_2=norm(D_2,2)/norm(X_2,2)


out3=[propF_C_1, prop2_C_1, propF_C_2, prop2_C_2, propF_D_1, prop2_D_1, propF_D_2, prop2_D_2]

out4=[norm(X_1), norm(X_1,2), norm(X_2), norm(X_2,2), norm(X_1)/norm(X_2), norm(X_1,2)/norm(X_2,2), norm(X_1)/norm(X_2), norm(X_1,2)/norm(X_2,2)]

#np.savetxt('siml1_seed' + sys.argv[1] + '_p' + str(p_1) + '_angle' + str(angle_deg) + '_nstd' + str(std_noise)  +'.txt', np.array([out1,out2,out3,out4]), delimiter=' ')

                        
    
            
