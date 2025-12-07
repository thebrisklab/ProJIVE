####################################################################################################
### Generates a replicate pair of synthetic data blocks from the simulations study described in ####
###         the ProJIVE Manuscript. After data have been generated, apply  ProJIVE          ########
### Author: Raphiel J. Murden                                                               ########
### Supervised by Benjamin Risk                                                             ########
####################################################################################################
# devtools::install.packages("ProJIVE_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(singR); library(CJIVE); library(reticulate); library(cowplot);
library(r.jive); library(ProJIVE)

#Simulation parameters
rep_number = 1
r.J = 1
r.I1 = 2
r.I2 = 2
#outdir = args[2]
n = 200
p1 = 20
p2 = 200 ####Note that p1 and p2 differ when compared to values used in simulations
JntVarEx1 = 0.1
JntVarEx2 = 0.1
#files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25
nparams = p1*(r.J+r.I1)+p2*(r.J+r.I2)+2
prop = n/nparams

#######################   Generate Toy  Data   #################################
set.seed(rep_number) ##To ensure that any randomized processes give the same results each time

###Construct Datasets
true_signal_ranks = r.J + c(r.I1,r.I2)    ##true ranks of overall signals
ToyDat = ProJIVE::GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), #mix.probs = c(.5,.5),
                         equal.eig = F, IndVarEx = c(IndVarEx1, IndVarEx2),
                         jnt_rank = r.J, ind_ranks = c(r.I1, r.I2), JntVarAdj = F,
                         SVD.plots = T, Error = T, print.cor = TRUE, Loads = "Gaussian",
                         Scores = "Gaussian",
                         # error.variances = c(1,1),
                         error.variances = list(rep(1:2, p1/2), rep(1:2, p2/2))
)

blocks <- ToyDat[["Data Blocks"]]
## Setup data to input into ProJIVE
Y = do.call(cbind, blocks)
P = c(p1, p2)
Q = c(r.J,r.I1,r.I2)

# Theta
JntScores = ToyDat[['Scores']][['Joint']]
IndivScore.X = ToyDat[['Scores']][["Indiv"]][,1:Q[2], drop = FALSE]
IndivScore.Y = ToyDat[['Scores']][["Indiv"]][,Q[2]+(1:Q[3]), drop = FALSE]

# WJ1 and WJ2
JntLd.X = t(ToyDat[['Loadings']][["Joint"]][[1]])
JntLd.Y = t(ToyDat[['Loadings']][["Joint"]][[2]])

# WI1 and WI2
IndivLd.X = t(ToyDat[['Loadings']][["Indiv"]][[1]])
IndivLd.Y = t(ToyDat[['Loadings']][["Indiv"]][[2]])

# Signal matrices
JX = ToyDat$`Data Components`[["JointSignalMatrices"]][[1]]
JY = ToyDat$`Data Components`[["JointSignalMatrices"]][[2]]
IX = ToyDat$`Data Components`[["IndivSignalMatrices"]][[1]]
IY = ToyDat$`Data Components`[["IndivSignalMatrices"]][[2]]

# Noise Matrices
EX = ToyDat$`Data Components`[["NoiseMatrices"]][[1]]
# plot(apply(EX, 2, var))
EY = ToyDat$`Data Components`[["NoiseMatrices"]][[2]]
# plot(apply(EY, 2, var))

# Total Signal matrices
AX = JX + IX
AY = JY + IY

## Proportions of variation explained
JVE.X = MatVar(JX)/MatVar(blocks[[1]])
JVE.Y = MatVar(JY)/MatVar(blocks[[2]])

IVE.X = MatVar(IX)/MatVar(blocks[[1]])
IVE.Y = MatVar(IY)/MatVar(blocks[[2]])

TotVE.X = MatVar((JX + IX))/MatVar(blocks[[1]])
TotVE.Y = MatVar((JY + IY))/MatVar(blocks[[2]])

Jnt.all = cbind(JX, JY)
Indiv.all = cbind(IX, IY)
s2n = var(c(Jnt.all+Indiv.all))/var(c(cbind(EX,EY)))


############       ProJIVE with loadings and error variance initialized from the truth        ###########
WJ.init = list(JntLd.X, JntLd.Y)
WI.init = list(IndivLd.X, IndivLd.Y)
init.loads = list(WJ.init, WI.init)
PJIVE.res.RJM.all = ProJIVE(Y=Y, P=P, Q=Q, Max.iter = 5000, init.loads = init.loads, num.starts = 1,
                           sig_hat = "MLE", diff.tol = 1e-14, return.all.starts = TRUE,
                           plots = TRUE, isotropic.error = TRUE)
PJIVE.res.RJM = PJIVE.res.RJM.all$ProJIVE_Results[[1]]

PJIVE.Jnt.all = PJIVE.res.RJM$SubjectScoreMatrix[,1:r.J]%*%t(PJIVE.res.RJM$LoadingMatrix[,1:r.J])
PJIVE.Indiv.all = PJIVE.res.RJM$SubjectScoreMatrix[,-(1:r.J)]%*%t(PJIVE.res.RJM$LoadingMatrix[,-(1:r.J)])

RSE.Jnt = norm(Jnt.all - PJIVE.Jnt.all, type = "F")^2/norm(Jnt.all, type = "F")^2
RSE.Jnt - MatVar(Jnt.all - PJIVE.Jnt.all)/MatVar(Jnt.all)

RSE.Indiv = norm(Indiv.all - PJIVE.Indiv.all, type = "F")^2/norm(Indiv.all, type = "F")^2
RSE.Indiv - MatVar(Indiv.all - PJIVE.Indiv.all)/MatVar(Indiv.all)

PJIVE.scores = PJIVE.res.RJM$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res.RJM$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res.RJM$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

PJIVE.err.var = PJIVE.res.RJM$ErrorVariances
summary(PJIVE.err.var[[2]])

## Since the model calls for the latent scores to have variance = 1, we should
##    scale the estimated scores by the inverse of their standard deviation
##    and scale the loadings by the same scalar

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores",
     col = c(rep("orange",n), rep("green",n), rep("purple",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J), drop = FALSE], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1), drop = FALSE], xlab = "True Indiv X1 Scores", ylab = "ProJIVE Indiv X1 Scores",
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1), drop = FALSE], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv X2 Scores", ylab = "ProJIVE Indiv X2 Scores",
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings and error variance \ninitialized  from True Values)",
                  sub = bquote(sigma["1"]*"="*.(round(mean(PJIVE.err.var[[1]]),2))*"; "*sigma["2"]*"="*.(round(mean(PJIVE.err.var[[2]]),2))))
legend("left", paste("Comp.", 1:r.J), pch = 1, col  = c("orange", "green", "purple")[1:r.J], bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X1", ylab = "ProJIVE Joint Loadings X1",
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X2", ylab = "ProJIVE Joint Loadings X2",
     col = c(rep("orange",p2), rep("green",p2), rep("purple",p2)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X",
     col = c(rep("orange",p1), rep("green",p1)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings X2", ylab = "ProJIVE Individual Loadings X2",
     col = c(rep("orange",p2), rep("green",p2)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

########################      Feng Diagnostics      ############################
layout(matrix(1:4, nrow = 2))
##########       True Ranks       ##########
X1.svd = svd(blocks[[1]], nu=r.J+r.I1, nv=r.J+r.I1)
X2.svd = svd(blocks[[2]], nu=r.J+r.I2, nv=r.J+r.I2)

M = cbind(X1.svd$u, X2.svd$u)
M.svd = svd(M)
wedin.bounds = get_wedin_bound_samples(M, M.svd, signal_rank = 1)
angs = sort(wedin.bounds*180/pi)
ang.probs = ecdf(angs)(angs)

random.bounds = get_random_direction_bound(n_obs = n,
                                           dims = sapply(blocks, ncol),
                                           ranks = r.J + c(r.I1,r.I2))
rand.angs = sort(acos(random.bounds-1)*180/pi)
rand.angs.probs = ecdf(rand.angs)(rand.angs)

M.angs = t(X1.svd$u)%*%X2.svd$u
M.angs.svd = svd(M.angs)
real.angs = acos(M.angs.svd$d)*180/pi

plot(angs, 1-ang.probs, xlim = c(0,90), col = 'blue', pch = '+',
     xlab = "Principal Angle", ylab = "Probability",
     main = paste0("X1:", r.J+r.I1, " & X2:", r.J+r.I2),
     sub = "True ranks")
lines(angs, 1-ang.probs, xlim = c(0,90), col = 'blue')
abline(v=quantile(angs, 0.95),col='blue',lty = 3)
points(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
lines(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
abline(v=quantile(rand.angs, 0.05),col='red',lty = 4)
abline(h=0.05,col='darkgray',lty = 5)
segments(x0 = real.angs, y0 = 0.2, y1 = 0.8)

##########       One Under       ##########
X1.svd = svd(blocks[[1]], nu=r.J+r.I1, nv=r.J+r.I1)
X2.svd = svd(blocks[[2]], nu=r.J+r.I2-1, nv=r.J+r.I2-1)

M = cbind(X1.svd$u, X2.svd$u)
M.svd = svd(M)
wedin.bounds = get_wedin_bound_samples(M, M.svd, signal_rank = 1)
angs = sort(wedin.bounds*180/pi)
ang.probs = ecdf(angs)(angs)

random.bounds = get_random_direction_bound(n_obs = n,
                                           dims = sapply(blocks, ncol),
                                           ranks = r.J + c(r.I1,(r.I2-1)))
rand.angs = sort(acos(random.bounds-1)*180/pi)
rand.angs.probs = ecdf(rand.angs)(rand.angs)

M.angs = t(X1.svd$u)%*%X2.svd$u
M.angs.svd = svd(M.angs)
real.angs = acos(M.angs.svd$d)*180/pi

plot(angs, 1-ang.probs, xlim = c(0,90), col = 'blue', pch = '+',
     xlab = "Principal Angle", ylab = "Probability",
     main = paste0("X1:", r.J+r.I1, " & X2:", r.J+r.I2-1),
     sub = "rank(X2) under")
lines(angs, 1-ang.probs, xlim = c(0,90), col = 'blue')
abline(v=quantile(angs, 0.95),col='blue',lty = 3)
points(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
lines(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
abline(v=quantile(rand.angs, 0.05),col='red',lty = 4)
abline(h=0.05,col='darkgray',lty = 5)
segments(x0 = real.angs, y0 = 0.2, y1 = 0.8)

##########       One over       ##########
X1.svd = svd(blocks[[1]], nu=r.J+r.I1+1, nv=r.J+r.I1+1)
X2.svd = svd(blocks[[2]], nu=r.J+r.I2, nv=r.J+r.I2)

M = cbind(X1.svd$u, X2.svd$u)
M.svd = svd(M)
wedin.bounds = get_wedin_bound_samples(M, M.svd, signal_rank = 1)
angs = sort(wedin.bounds*180/pi)
ang.probs = ecdf(angs)(angs)

random.bounds = get_random_direction_bound(n_obs = n,
                                           dims = sapply(blocks, ncol),
                                           ranks = r.J + c(r.I1+1,(r.I2)))
rand.angs = sort(acos(random.bounds-1)*180/pi)
rand.angs.probs = ecdf(rand.angs)(rand.angs)

M.angs = t(X1.svd$u)%*%X2.svd$u
M.angs.svd = svd(M.angs)
real.angs = acos(M.angs.svd$d)*180/pi

plot(angs, 1-ang.probs, xlim = c(0,90), col = 'blue', pch = '+',
     xlab = "Principal Angle", ylab = "Probability",
     main = paste0("X1:", r.J+r.I1+1, " & X2:", r.J+r.I2),
     sub = "rank(X1) over")
lines(angs, 1-ang.probs, xlim = c(0,90), col = 'blue')
abline(v=quantile(angs, 0.95),col='blue',lty = 3)
points(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
lines(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
abline(v=quantile(rand.angs, 0.05),col='red',lty = 4)
abline(h=0.05,col='darkgray',lty = 5)
segments(x0 = real.angs, y0 = 0.2, y1 = 0.8)

##########       Both over       ##########
X1.svd = svd(blocks[[1]], nu=r.J+r.I1+1, nv=r.J+r.I1+1)
X2.svd = svd(blocks[[2]], nu=r.J+r.I2+1, nv=r.J+r.I2+1)

M = cbind(X1.svd$u, X2.svd$u)
M.svd = svd(M)
wedin.bounds = get_wedin_bound_samples(M, M.svd, signal_rank = 1)
angs = sort(wedin.bounds*180/pi)
ang.probs = ecdf(angs)(angs)

random.bounds = get_random_direction_bound(n_obs = n,
                                           dims = sapply(blocks, ncol),
                                           ranks = r.J + c(r.I1+1,(r.I2+1)))
rand.angs = sort(acos(random.bounds-1)*180/pi)
rand.angs.probs = ecdf(rand.angs)(rand.angs)

M.angs = t(X1.svd$u)%*%X2.svd$u
M.angs.svd = svd(M.angs)
real.angs = acos(M.angs.svd$d)*180/pi

plot(angs, 1-ang.probs, xlim = c(0,90), col = 'blue', pch = '+',
     xlab = "Principal Angle", ylab = "Probability",
     main = paste0("X1:", r.J+r.I1+1, " & X2:", r.J+r.I2+1),
     sub = "Both over")
lines(angs, 1-ang.probs, xlim = c(0,90), col = 'blue')
abline(v=quantile(angs, 0.95),col='blue',lty = 3)
points(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
lines(rand.angs, rand.angs.probs, xlim = c(0,90), col = 'red')
abline(v=quantile(rand.angs, 0.05),col='red',lty = 4)
abline(h=0.05,col='darkgray',lty = 5)
segments(x0 = real.angs, y0 = 0.2, y1 = 0.8)
layout(1)
