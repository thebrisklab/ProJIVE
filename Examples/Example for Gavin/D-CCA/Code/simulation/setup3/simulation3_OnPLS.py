
import numpy as np
import sys
sys.path.append('/scratch/hshu/paper1_revision/OnPLS-master/')
import OnPLS
from test_zero_corr import test_zero_corr

######################This is the MKL parameter part #########################
#### MKL enalbes multithreads matrix computation, dramatically increase the run time efficiency
import ctypes
mkl_rt = ctypes.CDLL('libmkl_rt.so')
mkl_get_max_threads = mkl_rt.mkl_get_max_threads
mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(1)))


Y_1 = np.loadtxt('Y_1.txt')
Y_2 = np.loadtxt('Y_2.txt')

Y_1 = Y_1.T
Y_2 = Y_2.T

cv_scores = np.zeros([3])

np.random.seed(0)
for jointrank in range(1,4):
    # Create the estimator
    # Define the connections between blocks
    predComp = [[0, jointrank], [jointrank, 0]]
    # Define the numbers of non-global components
    orthComp = [3-jointrank, 5-jointrank]
    onpls = OnPLS.estimators.OnPLS(predComp, orthComp,numReps=20)# default numReps=1, to avoid local minima
    cv_scores_i = OnPLS.resampling.cross_validation(onpls, [Y_1, Y_2], cv_rounds=5, random_state=0)
    cv_scores[jointrank-1] = np.mean(cv_scores_i) 


opti_jointrank=cv_scores.argmax()+1

# Define the connections between blocks
predComp = [[0, opti_jointrank], [opti_jointrank, 0]]
# Define the numbers of non-global components
orthComp = [3-opti_jointrank, 5-opti_jointrank]
  
onpls = OnPLS.estimators.OnPLS(predComp, orthComp, numReps=20)# default numReps=1, to avoid local minima

# Fit a model
rlt = onpls.fit([Y_1, Y_2])

C_1_OnPLS=np.transpose( np.dot(rlt.T[0],np.transpose(rlt.P[0])) )
C_2_OnPLS=np.transpose( np.dot(rlt.T[1],np.transpose(rlt.P[1])) )

D_1_OnPLS=np.transpose( np.dot(rlt.To[0],np.transpose(rlt.Po[0])) )
D_2_OnPLS=np.transpose( np.dot(rlt.To[1],np.transpose(rlt.Po[1])) )

X_1_OnPLS=C_1_OnPLS+D_1_OnPLS
X_2_OnPLS=C_2_OnPLS+D_2_OnPLS

np.savetxt('C_1_OnPLS.txt', C_1_OnPLS, delimiter=' ')
np.savetxt('C_2_OnPLS.txt', C_2_OnPLS, delimiter=' ')

np.savetxt('D_1_OnPLS.txt', D_1_OnPLS, delimiter=' ')
np.savetxt('D_2_OnPLS.txt', D_2_OnPLS, delimiter=' ')

np.savetxt('X_1_OnPLS.txt', X_1_OnPLS, delimiter=' ')
np.savetxt('X_2_OnPLS.txt', X_2_OnPLS, delimiter=' ')


#test correlations between D_1_onpls and D_2_onpls  
nlogP = [-np.log(test_zero_corr(D_1_OnPLS[i,:], D_2_OnPLS[j,:]))  for i in range(D_1_OnPLS.shape[0]) for j in range(D_2_OnPLS.shape[0])]

np.savetxt('siml3_OnPLS_nlogP.txt', nlogP, delimiter=' ')















