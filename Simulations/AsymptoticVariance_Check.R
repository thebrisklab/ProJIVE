#############################################################################################################################
### Generates 100 replicate pairs of data sets as outlined in the ProJIVE manuscript        #################################
###         replicates 2-100 use the same loading matrix as the first replicate             #################################
### Author: Raphiel J. Murden                                                               #################################
### Supervised by Benjamin Risk                                                             #################################
#############################################################################################################################
rep_number = 1
r.J = 3
r.I1 = 1
r.I2 = 1
#outdir = args[2]
n = 1000
p1 = 20
p2 = 20 ####Note that p1 and p2 differ when compared to values used in simulations
JntVarEx1 = 0.5
JntVarEx2 = 0.5
#files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25

# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_to_SimulateData.R'
prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
prog.gipca.dir = "H:/My Documents/P-JIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
# source(file.path(prog.dir, "Functions_for_CJIVE.R"))
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
library(CJIVE); library(singR)
# source(file.path(prog.dir, "Functions_for_PJIVE_DEBUG.R"))
gipca.files = list.files(prog.gipca.dir,full.names = TRUE)
for(file.nm in gipca.files){source(file.nm)}
doc.dir = "H:/My Documents/P-JIVE/Programs/Examples/Output"

#######################################################################
set.seed(rep_number) ##To ensure that any randomized processes give the same results each time

###Construct Datasets
true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true ranks of overall signals
ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), equal.eig = TRUE, 
                         IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J, ind_ranks = c(r.I1, r.I2),
                         JntVarAdj = T, SVD.plots = F, Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
## Proportions of groups for mixture
mix.probs = c(0.2, 0.5, 0.3)
diagnoses = factor(c(rep(1, each = n*mix.probs[1]),rep(2, each = n*mix.probs[2]),rep(3, each = n*mix.probs[3])))
blocks <- ToyDat[["Data Blocks"]]
##Setup input parameters to use Gavin's EM code
Y = do.call(cbind, blocks); P = c(p1, p2); Q = c(r.J,r.I1,r.I2)

# Theta
JntScores = ToyDat[['Scores']][['Joint']]
IndivScore.X = ToyDat[['Scores']][["Indiv"]][,1:Q[2]]
IndivScore.Y = ToyDat[['Scores']][["Indiv"]][,Q[2]+(1:Q[2])]

# WJ1 and WJ2
JntLd.X = t(ToyDat[['Loadings']][["Joint"]][[1]])
JntLd.Y = t(ToyDat[['Loadings']][["Joint"]][[2]])

# WI1 and WI2
IndivLd.X = t(ToyDat[['Loadings']][["Indiv"]][[1]])
IndivLd.Y = t(ToyDat[['Loadings']][["Indiv"]][[2]])

#############################################
############### ProJIVE  ####################
WJ.init = list(JntLd.X, JntLd.Y)
WI.init = list(IndivLd.X, IndivLd.Y)
init.loads = list(WJ.init, WI.init)
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, init.loads = init.loads, sig_hat = "MLE")

sig.hat.X = sig.hat.Y = NULL
sig.hat.X = c(PJIVE.res$ErrorVariances[1], sig.hat.X)
sig.hat.Y = c(PJIVE.res$ErrorVariances[2], sig.hat.Y)

PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+(1:r.I1))]

sim.loads.X = sim.loads.Y = list()
sim.loads.X[[1]] = PJIVE.loads.X
sim.loads.Y[[1]] = PJIVE.loads.Y

mse.sim.loads.X = mse.sim.loads.Y = list()
mse.sim.loads.X[[1]] = (cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)^2
mse.sim.loads.Y[[1]] = (cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)^2

pmse.sim.loads.X = pmse.sim.loads.Y = NULL
pmse.sim.loads.X = c(pmse.sim.loads.X, pmse(S1 = cbind(JntLd.X, IndivLd.X), S2 = PJIVE.loads.X))
pmse.sim.loads.Y = c(pmse.sim.loads.Y, pmse(S1 = cbind(JntLd.Y, IndivLd.Y), S2 = PJIVE.loads.Y))

abs.err.sim.loads.X = abs.err.sim.loads.Y = list()
abs.err.sim.loads.X[[1]] = abs(cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)
abs.err.sim.loads.X[[1]][abs(cbind(JntLd.X, IndivLd.X)) != 0] = (abs(cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)/abs(cbind(JntLd.X, IndivLd.X)))[abs(cbind(JntLd.X, IndivLd.X)) != 0]

abs.err.sim.loads.Y[[1]] = abs(cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)
abs.err.sim.loads.Y[[1]][abs(cbind(JntLd.Y, IndivLd.Y)) != 0] = (abs(cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)/abs(cbind(JntLd.Y, IndivLd.Y)))[abs(cbind(JntLd.Y, IndivLd.Y)) != 0]

num.sims = 100
for(l in 1:(num.sims-1)){
  set.seed(l)
  ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), 
                           IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J, ind_ranks = c(r.I1, r.I2),
                           JntVarAdj = T, SVD.plots = F, Error = T, print.cor = F, Loads = init.loads, Scores = "Gaussian")
  
  blocks <- ToyDat[["Data Blocks"]]
  Y = do.call(cbind, blocks)
  PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, sig_hat = "MLE", plots = FALSE)
  
  sig.hat.X = c(PJIVE.res$ErrorVariances[1], sig.hat.X)
  sig.hat.Y = c(PJIVE.res$ErrorVariances[2], sig.hat.Y)
  
  PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
  PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]
  
  Subj.Scores.X = PJIVE.res$SubjectScoreMatrix[,1:(r.J+r.I1)]
  Subj.Scores.Y = PJIVE.res$LoadingMatrix[,-r.J+(1:r.I1)]
  
  sim.loads.X[[l+1]] = PJIVE.loads.X
  sim.loads.Y[[l+1]] = PJIVE.loads.Y
  
  mse.sim.loads.X[[l+1]] = (cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)^2
  mse.sim.loads.Y[[l+1]] = (cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)^2

  pmse.sim.loads.X = c(pmse.sim.loads.X, pmse(S1 = cbind(JntLd.X, IndivLd.X), S2 = PJIVE.loads.X))
  pmse.sim.loads.Y = c(pmse.sim.loads.Y, pmse(S1 = cbind(JntLd.Y, IndivLd.Y), S2 = PJIVE.loads.Y))
  
  abs.err.sim.loads.X[[l+1]] = abs(cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)
  abs.err.sim.loads.X[[l+1]][abs(cbind(JntLd.X, IndivLd.X)) != 0] = (abs(cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)/abs(cbind(JntLd.X, IndivLd.X)))[abs(cbind(JntLd.X, IndivLd.X)) != 0]

  abs.err.sim.loads.Y[[l+1]] = abs(cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)
  abs.err.sim.loads.Y[[l+1]][abs(cbind(JntLd.Y, IndivLd.Y)) != 0] = (abs(cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)/abs(cbind(JntLd.Y, IndivLd.Y)))[abs(cbind(JntLd.Y, IndivLd.Y)) != 0]
}

#######################################################################################################################################
#######################################      Images of Standard Error Matrices     ####################################################
mean.loads.X = abs.loads.X = 0
for(l in 1:num.sims){
  mean.loads.X = mean.loads.X + (mse.sim.loads.X[[l]])/num.sims
  abs.loads.X = abs.loads.X + (abs.err.sim.loads.X[[l]])/num.sims
}

mean.loads.Y = abs.loads.Y = 0
for(l in 1:num.sims){
  mean.loads.Y = mean.loads.Y + (mse.sim.loads.Y[[l]])/num.sims
  abs.loads.Y = abs.loads.Y + (abs.err.sim.loads.Y[[l]])/num.sims
}

pdf(file = file.path(doc.dir, paste0("CheckAsympVariance_GaussianLoadings_", lubridate::today(), ".pdf")))
layout(matrix(1:8, nrow = 4))
show.image.2(cbind(JntLd.X, IndivLd.X), main = "True X1 Loading matrix")
show.image.2(cbind(JntLd.Y, IndivLd.Y), main = "True X2 Loading matrix")

show.image.2(mean.loads.X, main = "X1 Loadings MSE")
show.image.2(mean.loads.Y, main = "X2 Loadings MSE")

show.image.2(abs.loads.X, main = "X1 Loadings Abs. Rel. Error")
show.image.2(abs.loads.Y, main = "X2 Loadings Abs. Rel. Error")

hist(pmse.sim.loads.X, main = "X1 Loadings PMSE")
hist(pmse.sim.loads.Y, main = "X2 Loadings PMSE")
layout(1)
dev.off()

hist(sig.hat.X, main = "Hist. of Error Variances for X1", breaks = 20)
hist(sig.hat.Y, main = "Hist. of Error Variances for X2", breaks = 20)

###############################################################################################################################
#########################   Examine information based variance estimator from Meilijson    ####################################
W.mats = w_to_w_k(PJIVE.res$LoadingMatrix, P, Q)
err.sigs = PJIVE.res$ErrorVariances

W.mat1.inv = solve(W.mats[[1]]%*%t(W.mats[[1]])+err.sigs[1]*diag(P[1]))
W.mat2.inv = solve(W.mats[[2]]%*%t(W.mats[[2]])+err.sigs[2]*diag(P[2]))

sig.1.star = diag(rep(1,(r.J+r.I1))) - t(W.mats[[1]])%&%W.mat1.inv%*%W.mats[[1]]
sig.2.star = diag(rep(1,(r.J+r.I2))) - t(W.mats[[2]])%&%W.mat2.inv%*%W.mats[[2]]

Y1.star = t(W.mats[[1]])%*%W.mat1.inv%*%t(blocks[[1]])
Y2.star = t(W.mats[[2]])%*%W.mat1.inv%*%t(blocks[[2]])

I.W1.complete = err.sigs[1]*diag(P[1]) 
I.W2.complete = err.sigs[2]*diag(P[2])

Var.t.t.1 = 0
Var.t.t.2 = 0
for(i in 1:n){
  Var.t.t.1 = Var.t.t.1 + 2*(sig.1.star + Y1.star[,i]%*%t(Y1.star[,i]))%*%(sig.1.star + Y1.star[,i]%*%t(Y1.star[,i]))
                        + sum(Y1.star[,i]^2)*(sig.1.star - Y1.star[,i]%*%t(Y1.star[,i]))
                        + sum(diag(sig.1.star))*(sig.1.star+ Y1.star[,i]%*%t(Y1.star[,i])) - Y1.star[,i]%*%t(Y1.star[,i])
  Var.t.t.2 = Var.t.t.2 + 2*(sig.2.star + Y2.star[,i]%*%t(Y2.star[,i]))%*%(sig.2.star + Y2.star[,i]%*%t(Y2.star[,i]))
  + sum(Y2.star[,i]^2)*(sig.2.star - Y2.star[,i]%*%t(Y2.star[,i]))
  + sum(diag(sig.2.star))*(sig.2.star+ Y2.star[,i]%*%t(Y2.star[,i])) - Y2.star[,i]%*%t(Y2.star[,i])
}

Var.S.W1.cond = err.sigs[1]^(-2)*W.mats[[1]]%*%(Var.t.t.1)%*%t(W.mats[[1]])
Var.S.W2.cond = err.sigs[2]^(-2)*W.mats[[2]]%*%(Var.t.t.2)%*%t(W.mats[[2]])

I.y.W1 = I.W1.complete + Var.S.W1.cond
Var.W1 = solve(I.y.W1)
show.image.2(as.matrix(Var.W1))

I.y.W2 = I.W2.complete + Var.S.W2.cond
Var.W2 = solve(I.y.W2)
show.image.2(as.matrix(Var.W2))

###############################################################################################################################
#######################################      Histograms of Norms for X1    ####################################################
################## X1 Loading Matrix
pdf(file = file.path(doc.dir, paste0("CheckAsympVariance_HistogramsForX1_", lubridate::today(), ".pdf")))
layout(matrix(1:6, nrow =  2, byrow = TRUE))
chord.norms.X = sapply(sim.loads.X, function(X) chord.norm.diff(X, PJIVE.loads.X))
hist(chord.norms.X, main = "Histogram of X1 Loading Chordal Norms", breaks = 20)
F.norms.X = sapply(sim.loads.X, function(X) norm(X-PJIVE.loads.X,"F")/norm(PJIVE.loads.X,"F"))
hist(F.norms.X, main = "Histogram of X1 Loading Rel. F Norms", breaks = 20)

################## X1 Joint Loading Matrix
chord.norms.JX = sapply(sim.loads.X, function(X) chord.norm.diff(X[,1:r.J], JntLd.X))
hist(chord.norms.JX, main = "Histogram of X1 Joint Loading Chordal Norms", breaks = 20)
F.norms.JX = sapply(sim.loads.X, function(X) norm(X[,1:r.J]-JntLd.X,"F")/norm(JntLd.X,"F"))
hist(F.norms.JX, main = "Histogram of X1 Joint Loading Rel. F Norms", breaks = 20)

################## X1 Indvidual Loading Matrix
chord.norms.IX = sapply(sim.loads.X, function(X) chord.norm.diff(X[,-(1:r.J)], IndivLd.X))
hist(chord.norms.IX, main = "Histogram of X1 Indep. Loading Chordal Norms", breaks = 20)
F.norms.IX = sapply(sim.loads.X, function(X) norm(X[,-(1:r.J)]-IndivLd.X,"F")/norm(IndivLd.X,"F"))
hist(F.norms.IX, main = "Histogram of X1 Indep Loading Rel. F Norms", breaks = 20)
dev.off()
layout(1)

###############################################################################################################################
#######################################      Histograms of Norms for X2    ####################################################
################## X2 Loading MatriY
pdf(file = file.path(doc.dir, paste0("CheckAsympVariance_HistogramsForX2_", lubridate::today(), ".pdf")))
layout(matrix(1:6, nrow =  2, byrow = TRUE))
chord.norms.Y = sapply(sim.loads.Y, function(Y) chord.norm.diff(Y, PJIVE.loads.Y))
hist(chord.norms.Y, main = "Histogram of X2 Loading Chordal Norms", breaks = 20)
F.norms.Y = sapply(sim.loads.Y, function(Y) norm(Y-PJIVE.loads.Y,"F")/norm(PJIVE.loads.Y,"F"))
hist(F.norms.Y, main = "Histogram of X2 Loading Rel. F Norms", breaks = 20)

################## X2 Joint Loading MatriY
chord.norms.JY = sapply(sim.loads.Y, function(Y) chord.norm.diff(Y[,1:r.J], JntLd.Y))
hist(chord.norms.JY, main = "Histogram of X2 Joint Loading Chordal Norms", breaks = 20)
F.norms.JY = sapply(sim.loads.Y, function(Y) norm(Y[,1:r.J]-JntLd.Y,"F")/norm(JntLd.Y,"F"))
hist(F.norms.JY, main = "Histogram of X2 Joint Loading Rel. F Norms", breaks = 20)

################## X2 Indvidual Loading MatriY
chord.norms.IY = sapply(sim.loads.Y, function(Y) chord.norm.diff(Y[,-(1:r.J)], IndivLd.Y))
hist(chord.norms.IY, main = "Histogram of X2 Indep. Loading Chordal Norms", breaks = 20)
F.norms.IY = sapply(sim.loads.Y, function(Y) norm(Y[,-(1:r.J)]-IndivLd.Y,"F")/norm(IndivLd.Y,"F"))
hist(F.norms.IY, main = "Histogram of X2 Indep Loading Rel. F Norms", breaks = 20)
dev.off()
layout(1)
