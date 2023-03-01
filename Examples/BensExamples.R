#############################################################################################################################
### simulate a simple model:                                                                #################################
### Author:  Benjamin Risk                                                             #################################
#############################################################################################################################
# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_for_PJIVE.R'
source("Functions_for_PJIVE.R")
library(singR); library(CJIVE); library(reticulate)

W = matrix(c(1,0,0,0,2,0,0,0,0,3,0,0,0,0,4,0,0,0),6,3,byrow = TRUE)
sigmasq1 = 1
sigmasq2 = 0.5
nobs=10000

theta = matrix(rnorm(n=nobs*3),nrow=nobs)
errors = cbind(matrix(sqrt(sigmasq1)*rnorm(n=nobs*3),nrow=nobs),matrix(sqrt(sigmasq2)*rnorm(n=nobs*3),nrow=nobs))
simdata = theta%*%t(10*W)+errors

################
PJIVE.res = ProJIVE_EM(Y=simdata, P=c(3,3), Q=c(1,1,1), Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE", init.loads = "AJIVE")
round(PJIVE.res$LoadingMatrix,2)
# loading matrix off by /sqrt(nobs)

PJIVE.res$ErrorVariances
# Errors off by /nobs
sqrt(nobs)

# Check Asymptotic Varinace calc
PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:3,1:2]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:3),c(1,3)]
PJIVE.err.var = PJIVE.res$ErrorVariances
ProJIVE.info = ProJIVE_AsymVar(list(PJIVE.loads.X, PJIVE.loads.Y), PJIVE.err.var, PJIVE.scores, 1, simdata)
show.image.2(ProJIVE.info$ObservedEmpericalInformationMatrix)
show.image.2(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)
Asym.Vars = diag(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)

plot(Asym.Vars)
show.image.2(matrix(Asym.Vars[1:6], ncol = 2))
show.image.2(matrix(Asym.Vars[7:12], ncol = 2))

load.dat = data.frame(Name = c(paste0("X.Joint.Ld1.P", 1:3), paste0("X.Indiv.Ld1.P", 1:3)),
                      Loading = c(PJIVE.loads.X), Upper = c(PJIVE.loads.X)+sqrt(Asym.Vars[1:6])*1.96,
                      Lower = c(PJIVE.loads.X)-sqrt(Asym.Vars[1:6])*1.96)
  
ggplot(load.dat) +
  geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper), width=0.4, colour="orange") 

# initialize from the truth:
PJIVE.res = ProJIVE_EM(Y=simdata, P=c(3,3), Q=c(1,1,1), Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE",
                       init.loads = list(list(matrix(W[1:3,1],3),matrix(W[4:6,1],3)),list(matrix(W[1:3,2],3),matrix(W[4:6,3],3))))

round(PJIVE.res$LoadingMatrix,2)
# weird results. 
PJIVE.res$ErrorVariances

# Check Asymptotic Varinace calc
PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:3,1:2]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:3),c(1,3)]
PJIVE.err.var = PJIVE.res$ErrorVariances
ProJIVE.info = ProJIVE_AsymVar(list(PJIVE.loads.X, PJIVE.loads.Y), PJIVE.err.var, PJIVE.scores, 1, simdata)
show.image.2(ProJIVE.info$ObservedEmpericalInformationMatrix)
show.image.2(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)
Asym.Vars = diag(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)

plot(Asym.Vars)
show.image.2(matrix(Asym.Vars[1:6], ncol = 2))
show.image.2(matrix(Asym.Vars[7:12], ncol = 2))

load.dat = data.frame(Name = c(paste0("X.Joint.Ld1.P", 1:3), paste0("X.Indiv.Ld1.P", 1:3)),
                      Loading = c(PJIVE.loads.X), Upper = c(PJIVE.loads.X)+sqrt(Asym.Vars[1:6])*1.96,
                      Lower = c(PJIVE.loads.X)-sqrt(Asym.Vars[1:6])*1.96)

ggplot(load.dat) +
  geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper), width=0.4, colour="orange") 














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
