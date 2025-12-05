# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 15:01:27 2017

"""

import numpy as np
import scipy as sp
import math

def sPOET(Y, K=None, K_max=None, K_min=None, method=None):
#assume Y is mean-zero p*n matrix, and K is the rank of its covariance matrix.
    p,n = Y.shape
        
    U, S, V_t = np.linalg.svd(Y, full_matrices=False)
    
    S = np.diag(S)
    
    Lambda = S**2/n # whose diagonal are eigenvalues of the sample covariance matrix
    
    #### select K
    if K is None:
        if method=='GR':
            K = select_K_by_GR(Lambda.diagonal(),K_max, K_min)
        else: 
            K = select_K_by_ED(Lambda.diagonal(),K_max, K_min)
    ####
                    
    c_hat = sum( Lambda.diagonal()[K:] ) / (p - K - p*K/n)
        
    Lambda_S = np.maximum(Lambda[:K,:K]-c_hat*p/n,0)
    
    X_hat = U[:,:K] @ np.sqrt(Lambda_S*n) @ V_t[:K,:]
    
    return X_hat, Lambda_S, U[:,:K], K, V_t # U is the Gamma matrix in Fan's Annals paper
 
        
def dCCA(Y_1,Y_2, r_1=None, r_2=None, r_12=None, method=None):

    p_1, n = Y_1.shape
    p_2, _ = Y_2.shape
    
   
    X_1_hat, Lambda_1, V_1, r_1, V_y1_t = sPOET(Y=Y_1, K=r_1, method=method)
    X_2_hat, Lambda_2, V_2, r_2, V_y2_t = sPOET(Y=Y_2, K=r_2, method=method)
    
        
    #Sigma_1 = 1.0/n * X_1_hat @ X_1_hat.T
    #Sigma_2 = 1.0/n * X_2_hat @ X_2_hat.T        
    #Sigma_12 = 1.0/n * X_1_hat @ X_2_hat.T


        
    Lambda_1_inv_half = np.zeros(Lambda_1.shape)
    index_nonzero_diag_Lambda_1 = np.nonzero(Lambda_1.diagonal()>0)
    Lambda_1_inv_half[index_nonzero_diag_Lambda_1,index_nonzero_diag_Lambda_1] = Lambda_1.diagonal()[index_nonzero_diag_Lambda_1]**-0.5
    
    Lambda_2_inv_half = np.zeros(Lambda_2.shape)
    index_nonzero_diag_Lambda_2 = np.nonzero(Lambda_2.diagonal()>0)
    Lambda_2_inv_half[index_nonzero_diag_Lambda_2,index_nonzero_diag_Lambda_2] = Lambda_2.diagonal()[index_nonzero_diag_Lambda_2]**-0.5
    
    
    Theta = (Lambda_1_inv_half @ V_1.T @ X_1_hat) @ (X_2_hat.T @ V_2 @ Lambda_2_inv_half)/n
            
    U_theta, D_theta, V_theta_t= np.linalg.svd(Theta, full_matrices=True) # D_theta is a vector
    
            
    Gamma_1 = V_1 @ Lambda_1_inv_half @ U_theta
    
    Gamma_2 = V_2 @ Lambda_2_inv_half @ V_theta_t.T
        
    D_theta= np.minimum(D_theta, 1) # modified to <=1
    
    r_theta = sum(D_theta>0) # the rank of Theta estimate
    
    
    if r_12 is None or r_12 < 1:        
        ccor_hat = np.linalg.svdvals(V_y1_t[:r_1,:] @ V_y2_t.T[:,:r_2])
        r_12 = select_r12_by_MDLIC(ccor_hat, r_1, r_2, n)
 

            
    if r_12 > r_theta:
        r_12 = r_theta
        
    
    A_mat_C = np.diag( 0.5 * (1-  ((1-D_theta[:r_theta])/(1+D_theta[:r_theta]))**0.5 ) )
        
    B_1 = X_1_hat @ (X_1_hat.T @ Gamma_1[:,:r_12])/n
    B_2 = X_2_hat @ (X_2_hat.T @ Gamma_2[:,:r_12])/n
    C_base = A_mat_C[:r_12,:r_12] @ ( Gamma_1[:,:r_12].T @ X_1_hat + Gamma_2[:,:r_12].T @ X_2_hat )
    
            
    C_1_hat = B_1 @ C_base
    C_2_hat = B_2 @ C_base
    
    B_1_rtheta = X_1_hat @ (X_1_hat.T @ Gamma_1[:,:r_theta])/n
    B_2_rtheta = X_2_hat @ (X_2_hat.T @ Gamma_2[:,:r_theta])/n
    C_base_rtheta = A_mat_C[:r_theta,:r_theta] @ ( Gamma_1[:,:r_theta].T @ X_1_hat + Gamma_2[:,:r_theta].T @ X_2_hat )
 
    
    D_1_hat = X_1_hat - B_1_rtheta @ C_base_rtheta        
    D_2_hat = X_2_hat - B_2_rtheta @ C_base_rtheta
     
    X_1_hat = C_1_hat + D_1_hat
    X_2_hat = C_2_hat + D_2_hat
                    
    return X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1, r_2, r_12, D_theta[:r_12], np.arccos(D_theta[:r_12])/math.pi*180
    
  
    
    
    
    
    



