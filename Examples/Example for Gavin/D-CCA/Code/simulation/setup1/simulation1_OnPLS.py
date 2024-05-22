
import numpy as np
import sys
sys.path.append('/scratch/hshu/paper1_revision/OnPLS-master/')
import OnPLS
from numpy.linalg import norm
from test_zero_corr import test_zero_corr

######################This is the MKL parameter part #########################
#### MKL enalbes multithreads matrix computation, dramatically increase the run time efficiency
import ctypes
mkl_rt = ctypes.CDLL('libmkl_rt.so')
mkl_get_max_threads = mkl_rt.mkl_get_max_threads
mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(1)))

seed = sys.argv[1]

Y_1 = np.loadtxt('/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed'+seed+'_p900_angle45_nstd1_Y_1.txt')
Y_2 = np.loadtxt('/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed'+seed+'_p900_angle45_nstd1_Y_2.txt')

Y_1 = Y_1.T
Y_2 = Y_2.T

cv_scores = np.zeros([3])

np.random.seed(0)
for jointrank in range(1,4):
    # Create the estimator
    # Define the connections between blocks
    predComp = [[0, jointrank], [jointrank, 0]]
    # Define the numbers of non-global components
    orthComp = [3-jointrank, 3-jointrank]
    onpls = OnPLS.estimators.OnPLS(predComp, orthComp,numReps=20)# default numReps=1, to avoid local minima
    cv_scores_i = OnPLS.resampling.cross_validation(onpls, [Y_1, Y_2], cv_rounds=5, random_state=0)
    cv_scores[jointrank-1] = np.mean(cv_scores_i) 


opti_jointrank=cv_scores.argmax()+1

# Define the connections between blocks
predComp = [[0, opti_jointrank], [opti_jointrank, 0]]
# Define the numbers of non-global components
orthComp = [3-opti_jointrank, 3-opti_jointrank]
  
onpls = OnPLS.estimators.OnPLS(predComp, orthComp, numReps=20)# default numReps=1, to avoid local minima

# Fit a model
rlt = onpls.fit([Y_1, Y_2])

C_1_onpls=np.transpose( np.dot(rlt.T[0],np.transpose(rlt.P[0])) )
C_2_onpls=np.transpose( np.dot(rlt.T[1],np.transpose(rlt.P[1])) )

D_1_onpls=np.transpose( np.dot(rlt.To[0],np.transpose(rlt.Po[0])) )
D_2_onpls=np.transpose( np.dot(rlt.To[1],np.transpose(rlt.Po[1])) )

X_1_onpls=C_1_onpls+D_1_onpls
X_2_onpls=C_2_onpls+D_2_onpls


X_1 = np.loadtxt('/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed'+seed+'_p900_angle45_nstd1_X_1.txt')
X_2 = np.loadtxt('/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed'+seed+'_p900_angle45_nstd1_X_2.txt')

C_1 = np.loadtxt('/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed'+seed+'_p900_angle45_nstd1_C_1.txt')
C_2 = np.loadtxt('/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed'+seed+'_p900_angle45_nstd1_C_2.txt')

D_1 = X_1 - C_1
D_2 = X_2 - C_2

errF_X_1=norm(X_1-X_1_onpls)/norm(X_1)
err2_X_1=norm(X_1-X_1_onpls,2)/norm(X_1,2)
errF_X_2=norm(X_2-X_2_onpls)/norm(X_2)
err2_X_2=norm(X_2-X_2_onpls,2)/norm(X_2,2)  


errF_C_1 = norm(C_1-C_1_onpls)/norm(C_1)
err2_C_1 = norm(C_1-C_1_onpls,2)/norm(C_1,2)
errF_C_2 = norm(C_2-C_2_onpls)/norm(C_2)
err2_C_2 = norm(C_2-C_2_onpls,2)/norm(C_2,2)  

errF_D_1 = norm(D_1-D_1_onpls)/norm(D_1)
err2_D_1 = norm(D_1-D_1_onpls,2)/norm(D_1,2)
errF_D_2 = norm(D_2-D_2_onpls)/norm(D_2)
err2_D_2 = norm(D_2-D_2_onpls,2)/norm(D_2,2)  

output = [errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2]
   
np.savetxt('siml1_seed'+seed+'_p900_angle45_nstd1_OnPLS.txt', output, delimiter=' ')
  

#test correlations between D_1_onpls and D_2_onpls
nlogP = [-np.log(test_zero_corr(D_1_onpls[i,:], D_2_onpls[j,:]))  for i in range(D_1_onpls.shape[0]) for j in range(D_2_onpls.shape[0])]

np.savetxt('siml1_seed'+seed+'_p900_angle45_nstd1_OnPLS_nlogP.txt', nlogP, delimiter=' ')











