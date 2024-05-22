#Python 2.7.11
import numpy as np
import sys
sys.path.append('G:/Hai Windows/work/softwares/python_package/OnPLS-master/')
import OnPLS

X1=np.loadtxt('G:\Hai Windows\work\PVBD_3sets_app\Expression_hard_centered.txt')
X2=np.loadtxt('G:\Hai Windows\work\PVBD_3sets_app\Methylation_hard_centered.txt')
X3=np.loadtxt('G:\Hai Windows\work\PVBD_3sets_app\miRNA_hard_centered.txt')

X1=np.transpose(X1)
X2=np.transpose(X2)
X3=np.transpose(X3)

##warning if using model setting from the result of JIVE1
# Define the connections between blocks
#predComp = [[0, 2, 2], [2, 0, 2], [2, 2, 0]]
# Define the numbers of non-global components
#orthComp = [25, 14, 20]


## use model setting from the result of JIVE2
# Define the connections between blocks
predComp = [[0, 2, 2], [2, 0, 2], [2, 2, 0]]
# Define the numbers of non-global components
orthComp = [25, 14, 20]

# Create the estimator
np.random.seed(0)
onpls = OnPLS.estimators.OnPLS(predComp, orthComp,numReps=10)# default numReps=1


# Fit a model
rlt=onpls.fit([X1, X2, X3])


C1_onpls=np.transpose( np.dot(rlt.T[0],np.transpose(rlt.P[0])) )
C2_onpls=np.transpose( np.dot(rlt.T[1],np.transpose(rlt.P[1])) )
C3_onpls=np.transpose( np.dot(rlt.T[2],np.transpose(rlt.P[2])) )


D1_onpls=np.transpose( np.dot(rlt.To[0],np.transpose(rlt.Po[0])) )
D2_onpls=np.transpose( np.dot(rlt.To[1],np.transpose(rlt.Po[1])) )
D3_onpls=np.transpose( np.dot(rlt.To[2],np.transpose(rlt.Po[2])) )

X1_onpls=C1_onpls+D1_onpls
X2_onpls=C2_onpls+D2_onpls
X3_onpls=C3_onpls+D3_onpls

np.savetxt('Expression_C1_OnPLS.txt', C1_onpls, delimiter=' ')
np.savetxt('Methylation_C2_OnPLS.txt', C2_onpls, delimiter=' ')
np.savetxt('miRNA_C3_OnPLS.txt', C3_onpls, delimiter=' ')

np.savetxt('Expression_D1_OnPLS.txt', D1_onpls, delimiter=' ')
np.savetxt('Methylation_D2_OnPLS.txt', D2_onpls, delimiter=' ')
np.savetxt('miRNA_D3_OnPLS.txt', D3_onpls, delimiter=' ')

np.savetxt('Expression_X1_OnPLS.txt', X1_onpls, delimiter=' ')
np.savetxt('Methylation_X2_OnPLS.txt', X2_onpls, delimiter=' ')
np.savetxt('miRNA_X3_OnPLS.txt', X3_onpls, delimiter=' ')

