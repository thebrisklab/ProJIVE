#############################################################################################################################
### simulate a simple model:                                                                #################################
### Author:  Benjamin Risk                                                             #################################
#############################################################################################################################
# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_to_SimulateData.R'
prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
library(singR); library(CJIVE); library(reticulate)

W = matrix(c(1,0,0,0,2,0,0,0,0,3,0,0,0,0,4,0,0,0),6,3,byrow = TRUE)
sigmasq1 = 1
sigmasq2 = 0.5
nobs=10000

theta = matrix(rnorm(n=nobs*3),nrow=nobs)
errors = cbind(matrix(sqrt(sigmasq1)*rnorm(n=nobs*3),nrow=nobs),matrix(sqrt(sigmasq2)*rnorm(n=nobs*3),nrow=nobs))
simdata = theta%*%t(W)+errors

################
PJIVE.res = ProJIVE_EM(Y=simdata, P=c(3,3), Q=c(1,1,1), Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE", init.loads = "AJIVE")
round(PJIVE.res$LoadingMatrix,2)
# loading matrix off by /sqrt(nobs)

PJIVE.res$ErrorVariances
# Errors off by /nobs
sqrt(nobs)


# initialize from the truth:
PJIVE.res = ProJIVE_EM(Y=simdata, P=c(3,3), Q=c(1,1,1), Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE", init.loads = list(list(matrix(W[1:3,1],3),matrix(W[4:6,1],3)),list(matrix(W[1:3,2],3),matrix(W[4:6,3],3))))

round(PJIVE.res$LoadingMatrix,2)
# weird results. 
PJIVE.res$ErrorVariances
# error variances are not correct
PJIVE.res$ErrorVariances/100














############################
##########################3
#####################3
CJIVE.res = cc.jive(list(simdata[,1:3],simdata[,4:6]), signal.ranks=c(2,2),joint.rank=1,perm.test = FALSE)
names(CJIVE.res)
names(CJIVE.res$CanCorRes)
CJIVE.res$CanCorRes$Canonical_Correlations
# I found CJIVE a bit confusing. sJIVE is where the results are. 
# why are there two canonical correlations if joint.rank=1? maybe i need 
# to remind myself what's being done here
# this correlation should be close to one in this example I think?
J.hat = CJIVE.res$sJIVE$joint_matrices
I.hat = CJIVE.res$sJIVE$indiv_matrices

# CJIVE loading estimates
(WJ=lapply(J.hat, function(x) x[['v']]))
t(WJ[[1]])%*%WJ[[1]]
# the columns have been scaled to have norm equal to one
# is there a way to recover the correct values?
chord.norm.diff(WJ[[1]],W[1:3,1])

lapply(I.hat, function(x) x[['v']])
