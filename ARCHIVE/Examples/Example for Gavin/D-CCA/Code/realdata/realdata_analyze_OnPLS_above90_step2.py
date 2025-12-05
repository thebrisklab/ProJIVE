
import numpy as np
import sys
sys.path.append('/scratch/hshu/paper1_revision/OnPLS-master/')
import OnPLS


######################This is the MKL parameter part #########################
#### MKL enalbes multithreads matrix computation, dramatically increase the run time efficiency
#import ctypes
#mkl_rt = ctypes.CDLL('libmkl_rt.so')
#mkl_get_max_threads = mkl_rt.mkl_get_max_threads
#mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(1)))


Y_1 = np.loadtxt('Expression.txt')
Y_2 = np.loadtxt('Methylation_above90.txt')

Y_1 = Y_1.T
Y_2 = Y_2.T

cv_scores=np.zeros(shape=[10,11,11])

for i in np.arange(1,11):
    for j in range(11):
        for k in range(11):
            cv_scores[i-1,j,k] = np.loadtxt('OnPLS_above90_CVscore_common' + str(i) + '_distinct' + str(j) + 'and' + str(k) + '.txt')

            
loc_max=np.where(cv_scores==np.amax(cv_scores))

np.random.seed(0)
# Define the connections between blocks
predComp = [[0, loc_max[0][0]], [loc_max[0][0], 0]]
# Define the numbers of non-global components
orthComp = [loc_max[1][0], loc_max[2][0]]
  
onpls = OnPLS.estimators.OnPLS(predComp, orthComp, numReps=20)# default numReps=1, to avoid local minima

# Fit a model
rlt = onpls.fit([Y_1, Y_2])

C_1_onpls=np.transpose( np.dot(rlt.T[0],np.transpose(rlt.P[0])) )
C_2_onpls=np.transpose( np.dot(rlt.T[1],np.transpose(rlt.P[1])) )

D_1_onpls=np.transpose( np.dot(rlt.To[0],np.transpose(rlt.Po[0])) )
D_2_onpls=np.transpose( np.dot(rlt.To[1],np.transpose(rlt.Po[1])) )

X_1_onpls=C_1_onpls+D_1_onpls
X_2_onpls=C_2_onpls+D_2_onpls

np.savetxt('C_1_hat_OnPLS_above90.txt', C_1_onpls, delimiter=' ')
np.savetxt('C_2_hat_OnPLS_above90.txt', C_2_onpls, delimiter=' ')

np.savetxt('D_1_hat_OnPLS_above90.txt', D_1_onpls, delimiter=' ')
np.savetxt('D_2_hat_OnPLS_above90.txt', D_2_onpls, delimiter=' ')

np.savetxt('X_1_hat_OnPLS_above90.txt', X_1_onpls, delimiter=' ')
np.savetxt('X_2_hat_OnPLS_above90.txt', X_2_onpls, delimiter=' ')