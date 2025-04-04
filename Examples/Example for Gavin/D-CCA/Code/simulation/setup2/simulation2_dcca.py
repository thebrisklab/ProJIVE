# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 09:42:51 2017

"""


import numpy as np
import math
from scipy.stats import ortho_group
import sys
from numpy.linalg import norm
import dcca

######################This is the MKL parameter part #########################
#### MKL enalbes multithreads matrix computation, dramatically increase the run time efficiency
import ctypes
mkl_rt = ctypes.CDLL('libmkl_rt.so')
mkl_get_max_threads = mkl_rt.mkl_get_max_threads
mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(1)))

'''
Begin: Generate simulated data
'''

n =300

p_1_list = [[100,300,600,900,1200,1500],[900],[900]]
std_noise_list = [[1],[1],[0.1,0.5,1,2,3,4]]
angle_deg_list = [[45],[0,15,30,45,60,75],[45]]

for change in range(3):    
    for p_1 in p_1_list[change]:
        p_2 = 300    
        r_1, r_2, r_12 = 3, 5, 1
        ###################
        # Define V_1, V_2. 
        V_1 = ortho_group.rvs(dim=p_1,random_state=0)[:,:r_1]  
        V_2 = ortho_group.rvs(dim=p_2,random_state=0)[:,:r_2]  
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
        for std_noise in std_noise_list[change]:
            ###### noise matrices
            rng1 = np.random.RandomState(int(sys.argv[1]))
            E_1plus2 = rng1.normal(0, std_noise, (p_1+p_2)*n).reshape([p_1+p_2,n])            
            E_1 = E_1plus2[:p_1,:]
            E_2 = E_1plus2[p_1:,:]
            for angle_deg in angle_deg_list[change]:         
     
                # angles for each pair of canonical variables
                angle =np.array([angle_deg])/180.0*math.pi 
                #####################################
                
                D_theta = np.zeros(shape=(r_1,r_2)) # Let D_theta be a r_1*r_2 matrix
                
                for i in range(r_12):
                    D_theta[i,i] = math.cos(angle[i]) 
                
                
                Sigma_Z_1plus2 = np.bmat([ [np.identity(r_1), D_theta  ],  [D_theta.T, np.identity(r_2)] ]) #Covariance matrix of [Z_1,Z_2]
                
                rng2 = np.random.RandomState(int(sys.argv[1])+1001)  # random seed
                
                Z_1plus2 = rng2.multivariate_normal( np.zeros( shape=(r_1 + r_2,) ), Sigma_Z_1plus2, n ).T
                                                 
                Z_1 = Z_1plus2[:r_1,:]     
                Z_2 = Z_1plus2[r_1:,:] 
                
                #########################
                # Define  U_theta (size=r_1*r_1), V_theta (size=r_2*r_2). Let the full SVD of Theta matrix be: Theta = U_theta @ D_theta @ V_theta.T
                U_theta = np.identity(r_1)
                V_theta = np.identity(r_2)
                #####################
                
                
                Z_1_star = U_theta @ Z_1
                Z_2_star = V_theta @ Z_2
                
                ###################
                # Define V_1, V_2. 
               # V_1 = ortho_group.rvs(dim=p_1,random_state=0)[:,:r_1]  
               # V_2 = ortho_group.rvs(dim=p_2,random_state=0)[:,:r_2]  
                ###################
                
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
                
                ##select rank by Oracle
                X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor_Oracle, theta_Oracle  = dcca.dCCA(Y_1, Y_2, r_1=r_1, r_2=r_2, r_12=r_12)
                
                errF_X_1=norm(X_1-X_1_hat)/norm(X_1)
                err2_X_1=norm(X_1-X_1_hat,2)/norm(X_1,2)
                errF_X_2=norm(X_2-X_2_hat)/norm(X_2)
                err2_X_2=norm(X_2-X_2_hat,2)/norm(X_2,2)  
    
                
                errF_C_1 = norm(C_1-C_1_hat)/norm(C_1)
                err2_C_1 = norm(C_1-C_1_hat,2)/norm(C_1,2)
                errF_C_2 = norm(C_2-C_2_hat)/norm(C_2)
                err2_C_2 = norm(C_2-C_2_hat,2)/norm(C_2,2)  
                
                errF_D_1 = norm(D_1-D_1_hat)/norm(D_1)
                err2_D_1 = norm(D_1-D_1_hat,2)/norm(D_1,2)
                errF_D_2 = norm(D_2-D_2_hat)/norm(D_2)
                err2_D_2 = norm(D_2-D_2_hat,2)/norm(D_2,2)  
                
                out1 = [errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2, np.max(np.abs(D_1_hat@D_2_hat.T))]
                out1_rank = [r_1_hat, r_2_hat, r_12_hat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                
                
                ##select rank by ED
                X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor_ED, theta_ED  = dcca.dCCA(Y_1, Y_2, method='ED')
                
                errF_X_1=norm(X_1-X_1_hat)/norm(X_1)
                err2_X_1=norm(X_1-X_1_hat,2)/norm(X_1,2)
                errF_X_2=norm(X_2-X_2_hat)/norm(X_2)
                err2_X_2=norm(X_2-X_2_hat,2)/norm(X_2,2)  
    
                
                errF_C_1 = norm(C_1-C_1_hat)/norm(C_1)
                err2_C_1 = norm(C_1-C_1_hat,2)/norm(C_1,2)
                errF_C_2 = norm(C_2-C_2_hat)/norm(C_2)
                err2_C_2 = norm(C_2-C_2_hat,2)/norm(C_2,2)  
                
                errF_D_1 = norm(D_1-D_1_hat)/norm(D_1)
                err2_D_1 = norm(D_1-D_1_hat,2)/norm(D_1,2)
                errF_D_2 = norm(D_2-D_2_hat)/norm(D_2)
                err2_D_2 = norm(D_2-D_2_hat,2)/norm(D_2,2)  
                
                out2 = [errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2, np.max(np.abs(D_1_hat@D_2_hat.T))]
                out2_rank = [r_1_hat, r_2_hat, r_12_hat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                
                ##select rank by GR
                X_1_hat, X_2_hat, C_1_hat, C_2_hat, D_1_hat, D_2_hat, r_1_hat, r_2_hat, r_12_hat, ccor_GR, theta_GR  = dcca.dCCA(Y_1, Y_2, method='GR')
                
                errF_X_1=norm(X_1-X_1_hat)/norm(X_1)
                err2_X_1=norm(X_1-X_1_hat,2)/norm(X_1,2)
                errF_X_2=norm(X_2-X_2_hat)/norm(X_2)
                err2_X_2=norm(X_2-X_2_hat,2)/norm(X_2,2)  
    
                
                errF_C_1 = norm(C_1-C_1_hat)/norm(C_1)
                err2_C_1 = norm(C_1-C_1_hat,2)/norm(C_1,2)
                errF_C_2 = norm(C_2-C_2_hat)/norm(C_2)
                err2_C_2 = norm(C_2-C_2_hat,2)/norm(C_2,2)  
                
                errF_D_1 = norm(D_1-D_1_hat)/norm(D_1)
                err2_D_1 = norm(D_1-D_1_hat,2)/norm(D_1,2)
                errF_D_2 = norm(D_2-D_2_hat)/norm(D_2)
                err2_D_2 = norm(D_2-D_2_hat,2)/norm(D_2,2)  
                
                out3 = [errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2, np.max(np.abs(D_1_hat@D_2_hat.T))]
                out3_rank = [r_1_hat, r_2_hat, r_12_hat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                
                
                
                propF_C_1 = norm(C_1)/norm(X_1)
                prop2_C_1 = norm(C_1,2)/norm(X_1,2)
                propF_C_2 = norm(C_2)/norm(X_2)
                prop2_C_2 = norm(C_2,2)/norm(X_2,2)
                
                propF_D_1 = norm(D_1)/norm(X_1)
                prop2_D_1 = norm(D_1,2)/norm(X_1,2)
                propF_D_2 = norm(D_2)/norm(X_2)
                prop2_D_2 = norm(D_2,2)/norm(X_2,2)
                
                
                out4 = [propF_C_1, prop2_C_1, propF_C_2, prop2_C_2, propF_D_1, prop2_D_1, propF_D_2, prop2_D_2, 0, 0, 0, 0, 0]
                
                out5 = [norm(X_1), norm(X_1,2), norm(X_2), norm(X_2,2), norm(X_1)/norm(X_2), norm(X_1,2)/norm(X_2,2), 0, 0, 0, 0, 0, 0, 0]
                
                np.savetxt('siml2_seed' + sys.argv[1] + '_p' + str(p_1) + '_angle' + str(angle_deg) + '_nstd' + str(std_noise)  +'.txt', np.array([out1, out1_rank, out2, out2_rank, out3, out3_rank, out4, out5]), delimiter=' ')
                
                np.savetxt('siml2_theta_Oracle_seed' + sys.argv[1] + '_p' + str(p_1) + '_angle' + str(angle_deg) + '_nstd' + str(std_noise)  +'.txt', np.array([ccor_Oracle, theta_Oracle]), delimiter=' ')
                np.savetxt('siml2_theta_ED_seed' + sys.argv[1] + '_p' + str(p_1) + '_angle' + str(angle_deg) + '_nstd' + str(std_noise)  +'.txt', np.array([ccor_ED, theta_ED]), delimiter=' ')
                np.savetxt('siml2_theta_GR_seed' + sys.argv[1] + '_p' + str(p_1) + '_angle' + str(angle_deg) + '_nstd' + str(std_noise)  +'.txt', np.array([ccor_GR, theta_GR]), delimiter=' ')
                  
                                        
                    
                    
