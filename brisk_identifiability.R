######################
#######################
########################
#########################

# Example that shows
#    i) the joint versus individual variation is not identifiable
#    ii) the number of individual components is not identifiable


# for dataset 1, these can be anything
# if we make T_1 have singular values<=1
# a sufficient condition is that 
# (I - T_1 T_1^\top) is psd
#WJ1 = matrix(rnorm(10),nrow=5)
#WI1 = matrix(rnorm(10),nrow=5)
WJ1 = diag(5)[,c(1,2)]
WI1 = diag(5)[,c(2,3)]


# for dataset 2, we will 
# create a non-zero overlap in the
# col(WJ2) and col(WI2)
WJ2 = diag(5)[,c(1,2)]
WI2 = diag(5)[,c(1,3)]
# then (I-T_2 T_2^\top) can be negative
# in the overlapping subspace. 
# i.e., we can construct a valid
# individual components model 
# in which the class of T_1 is larger than 
# the class of orthogonal matrices

T1 = diag(c(sqrt(0.5+0.1),1)) # eigenvalue <1 for T1, so this will result in 
# an alternative decomposition of the individual components

T2 = t(solve(T1)) # the key here is that 
# the negative eigenvalue lie in the subspace shared by W_{J2} and W_{I2}
# and that the eigenvalue is less than sqrt(2)
# check that T2 has elements less than sqrt(2):
T2


WJ1star = WJ1%*%T1
WJ2star = WJ2%*%T2

# individual components covariance matrix:
covWI1star = WI1%*%t(WI1) + WJ1%*%(diag(ncol(WJ1))-T1%*%t(T1))%*%t(WJ1)
eig.WI1star = eigen(covWI1star)
# this example shows the number of individual components is not identifiable:
WI1star = eig.WI1star$vectors[,1:3]%*%diag(sqrt(eig.WI1star$values)[1:3])


covWI2star = WI2%*%t(WI2) + WJ2%*%(diag(ncol(WJ2))-T2%*%t(T2))%*%t(WJ2)
eig.WI2star = eigen(covWI2star)
eig.WI2star$values
# this is a WI2star matrix:
WI2star = eig.WI2star$vectors[,1:2]%*%diag(sqrt(eig.WI2star$values)[1:2])



# check covariance matrices are equivalent:
cov1 = WJ1%*%t(WJ1)+WI1%*%t(WI1)
cov2 = WJ2%*%t(WJ2)+WI2%*%t(WI2)

cov1star = WJ1star%*%t(WJ1star)+WI1star%*%t(WI1star)
cov2star = WJ2star%*%t(WJ2star)+WI2star%*%t(WI2star)

norm(cov1-cov1star,type = 'F')
norm(cov2-cov2star,type = 'F')


#########################
###########################
##############################
## E.g. in which W_{J2} and W_{I2} are linearly independent:
library(steadyICA)

WJ1 = diag(5)[,c(1,2)]
WI1 = diag(5)[,c(3,4)]

WJ2 = matrix(rnorm(10),5)
#WI2 = WJ2%*%theta2W(pi/4)
WI2 = WJ2%*%matrix(4:1,2)
#WI2 = 11*WJ2
U.WJ2 = svd(WJ2)$u
U.WI2 = svd(WI2)$u
svd(t(U.WI2)%*%U.WJ2)$d


T1 = diag(c(sqrt(0.5+0.4),sqrt(0.5+0.4))) # eigenvalue <1 for T1, so this will result in 
# an alternative decomposition of the individual components

T2 = t(solve(T1)) # the key here is that 
# the negative eigenvalue lie in the subspace shared by W_{J2} and W_{I2}
# and that the eigenvalue is less than sqrt(2)
# check that T2 has elements less than sqrt(2):
T2


WJ1star = WJ1%*%T1
WJ2star = WJ2%*%T2

# individual components covariance matrix:
covWI1star = WI1%*%t(WI1) + WJ1%*%(diag(ncol(WJ1))-T1%*%t(T1))%*%t(WJ1)
eigen(covWI1star)$values

covWI2star = WI2%*%t(WI2) + WJ2%*%(diag(ncol(WJ2))-T2%*%t(T2))%*%t(WJ2)
eig.WI2star = eigen(covWI2star)

# if these are non-negative, then we can construct a WI2star:
eig.WI2star$values
# this is a WI2star matrix:
WI2star = eig.WI2star$vectors[,1:2]%*%diag(sqrt(eig.WI2star$values[1:2]))

# check covariance matrices are equivalent:
cov1 = WJ1%*%t(WJ1)+WI1%*%t(WI1)
cov2 = WJ2%*%t(WJ2)+WI2%*%t(WI2)

cov1star = WJ1star%*%t(WJ1star)+covWI1star
cov2star = WJ2star%*%t(WJ2star)+WI2star%*%t(WI2star)

norm(cov1-cov1star,type = 'F')
norm(cov2-cov2star,type = 'F')


