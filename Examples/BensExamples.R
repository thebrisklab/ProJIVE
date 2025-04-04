#############################################################################################################################
### simulate a simple model:                                                                #################################
### Author:  Benjamin Risk                                                             #################################
#############################################################################################################################
# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_for_PJIVE.R'
source("Functions_for_PJIVE.R")
library(singR); library(CJIVE); library(reticulate)

W = matrix(c(1,0,0,0,2,0,0,0,0,
             3,0,0,0,0,0,0,0,4),6,3,byrow = TRUE)
# W = matrix(c(2,4,0,1,0,0,0,2,0,
#              3,0,2,0,0,1,4,0,0),6,3,byrow = TRUE)
sigmasq1 = 1
sigmasq2 = 2
nobs=100

theta = matrix(rnorm(n=nobs*3),nrow=nobs)
errors = cbind(matrix(sqrt(sigmasq1)*rnorm(n=nobs*3),nrow=nobs),matrix(sqrt(sigmasq2)*rnorm(n=nobs*3),nrow=nobs))
simdata = theta%*%t(W)+errors

################
PJIVE.res = ProJIVE_EM(Y=simdata, P=c(3,3), Q=c(1,1,1), Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE", init.loads = "AJIVE")
round(PJIVE.res$LoadingMatrix,2)
# loading matrix off by /sqrt(nobs)

PJIVE.res$ErrorVariances
sqrt(nobs)

# Check Asymptotic Varinace calc
PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:3,1:2]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:3),c(1,3)]
# PJIVE.loads.X.scaled = PJIVE.res$ErrorVariances[1]^(-.5)*PJIVE.loads.X
# PJIVE.loads.Y.scaled = PJIVE.res$ErrorVariances[2]^(-.5)*PJIVE.loads.Y
PJIVE.err.var = PJIVE.res$ErrorVariances
ProJIVE.info = ProJIVE_AsymVar(list(PJIVE.loads.X, PJIVE.loads.Y), PJIVE.err.var, PJIVE.scores, 1, simdata)
# show.image.2(ProJIVE.info$ObservedEmpericalInformationMatrix)
# show.image.2(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)
Asym.Vars = diag(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)/sqrt(nobs)
plot(ProJIVE.info$MeanScoreVector)

plot(Asym.Vars)
# show.image.2(matrix(Asym.Vars[1:6], ncol = 2))
# show.image.2(matrix(Asym.Vars[7:12], ncol = 2))

load.dat.X = data.frame(Name = c(paste0("X.Joint.Ld1.P", 1:3), paste0("X.Indiv.Ld1.P", 1:3)),
                        Loading = c(PJIVE.loads.X), Upper = c(PJIVE.loads.X)+sqrt(Asym.Vars[1:6])*1.96,
                        Lower = c(PJIVE.loads.X)-sqrt(Asym.Vars[1:6])*1.96)

ggplot(load.dat.X) +
  geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper), width=0.4, colour="orange") 

load.dat.Y = data.frame(Name = c(paste0("Y.Joint.Ld1.P", 1:3), paste0("Y.Indiv.Ld1.P", 1:3)),
                        Loading = c(PJIVE.loads.Y), Upper = c(PJIVE.loads.Y)+sqrt(Asym.Vars[7:12])*1.96,
                        Lower = c(PJIVE.loads.Y)-sqrt(Asym.Vars[7:12])*1.96)

ggplot(load.dat.Y) +
  geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper), width=0.4, colour="orange") 

#############################      Bootstrap?      #############################
B = 25
bstraps = list()
bstrap.mean = 0
theta.hat = PJIVE.res$SubjectScoreMatrix
W.hat = PJIVE.res$LoadingMatrix

for(b in 1:B){
  errors = cbind(matrix(sqrt(PJIVE.res$ErrorVariances[1])*rnorm(n=nobs*3),nrow=nobs),
                 matrix(sqrt(PJIVE.res$ErrorVariances[2])*rnorm(n=nobs*3),nrow=nobs))
  tempdata = theta.hat%*%t(W.hat)+errors
  
  
  temp.res = ProJIVE_EM(Y=tempdata, P=c(3,3), Q=c(1,1,1), Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE",
                        init.loads = "AJIVE", plots = FALSE, verbose = FALSE)
  
  bstrap.mean = bstrap.mean + temp.res$LoadingMatrix/B
  bstraps[[b]] = temp.res$LoadingMatrix
}

bstrap.var = 0
for(b in 1:B){
  bstrap.var = bstrap.var + matrix(bstraps[[b]] - bstrap.mean, ncol = 1)%*%matrix(bstraps[[b]] - bstrap.mean, nrow = 1)/(B-1)
}
round(bstrap.var)

bstrap.se = sqrt(diag(bstrap.var))

load.dat.X = data.frame(Name = c(paste0("X.Joint.Ld1.P", 1:3), paste0("X.Indiv.Ld1.P", 1:3)),
                        Loading = c(PJIVE.loads.X), Upper = c(PJIVE.loads.X)+bstrap.se[1:6]*1.96,
                        Lower = c(PJIVE.loads.X)-bstrap.se[1:6]*1.96)

load.dat.Y = data.frame(Name = factor(c(paste0("Y.Joint.Ld1.P", 1:3), paste0("Y.Indiv.Ld1.P", 1:3))),
                        Loading = c(PJIVE.loads.Y), Upper = c(PJIVE.loads.Y)+bstrap.se[7:12]*1.96,
                        Lower = c(PJIVE.loads.Y)-bstrap.se[7:12]*1.96)

p1 = ggplot(load.dat.X) +
      geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
      geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper),
                    width=0.4, colour="orange") + ylim(c(-3,3)) 

p2 = ggplot(load.dat.Y) +
      geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
      geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper),
                    width=0.4, colour="orange") + ylim(c(-3,3))

cowplot::plot_grid(p1 + ylim(c(-3,3)), p2 + ylim(c(-3,3)))


#######################    initialize from the truth:    #######################
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
Asym.Vars = diag(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)/sqrt(nobs)

plot(Asym.Vars)
show.image.2(matrix(Asym.Vars[1:6], ncol = 2))
show.image.2(matrix(Asym.Vars[7:12], ncol = 2))

load.dat.X = data.frame(Name = c(paste0("X.Joint.Ld1.P", 1:3), paste0("X.Indiv.Ld1.P", 1:3)),
                      Loading = c(PJIVE.loads.X), Upper = c(PJIVE.loads.X)+sqrt(Asym.Vars[1:6])*1.96,
                      Lower = c(PJIVE.loads.X)-sqrt(Asym.Vars[1:6])*1.96)

ggplot(load.dat.X) +
  geom_bar(aes(x=Name, y=Loading), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(x=Name, ymin=Lower, ymax=Upper), width=0.4, colour="orange") 

load.dat.Y = data.frame(Name = c(paste0("Y.Joint.Ld1.P", 1:3), paste0("Y.Indiv.Ld1.P", 1:3)),
                        Loading = c(PJIVE.loads.Y), Upper = c(PJIVE.loads.Y)+sqrt(Asym.Vars[7:12])*1.96,
                        Lower = c(PJIVE.loads.Y)-sqrt(Asym.Vars[7:12])*1.96)

ggplot(load.dat.Y) +
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
