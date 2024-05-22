
import numpy as np
import sys
sys.path.append('/scratch/hshu/paper1_revision/OnPLS-master/')
import OnPLS


######################This is the MKL parameter part #########################
#### MKL enalbes multithreads matrix computation, dramatically increase the run time efficiency
import ctypes
mkl_rt = ctypes.CDLL('libmkl_rt.so')
mkl_get_max_threads = mkl_rt.mkl_get_max_threads
mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(1)))


Y_1 = np.loadtxt('Expression.txt')
Y_2 = np.loadtxt('Methylation_below90.txt')

Y_1 = Y_1.T
Y_2 = Y_2.T


# Create the estimator
# Define the connections between blocks
predComp = [[0, int(sys.argv[1])], [int(sys.argv[1]), 0]]
# Define the numbers of non-global components
orthComp = [int(sys.argv[2]), int(sys.argv[3])]
onpls = OnPLS.estimators.OnPLS(predComp, orthComp,numReps=20)# default numReps=1, to avoid local minima
cv_scores = OnPLS.resampling.cross_validation(onpls, [Y_1, Y_2], cv_rounds=5, random_state=0)
  


np.savetxt('OnPLS_below90_CVscore_common' + sys.argv[1] + '_distinct' + sys.argv[2] + 'and' + sys.argv[3] + '.txt', [np.mean(cv_scores)])


