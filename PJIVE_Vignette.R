#############################################################################################################################
### Generates a pair of toy data blocks in the manner outlined within the CJIVE manuscript  #################################
###         After data have been generated                                                  #################################
### Author: Raphiel J. Murden                                                               #################################
### Supervised by Benjamin Risk                                                             #################################
#############################################################################################################################
rep_number = 0
r.J = 2
r.I1 = 2
r.I2 = 2
#outdir = args[2]
n = 200
p1 = 30
p2 = 50 ####Note that we have p2 = 1000 here as opposed to p2 = 10,000 in simulations
JntVarEx1 = 0.05
JntVarEx2 = 0.05
#files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25

# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_to_SimulateData.R'
prog.dir = "H:/My Documents/P-JIVE/Programs"
source(file.path(prog.dir, "Functions_for_CJIVE.R"))
source(file.path(prog.dir, "Functions_for_PJIVE.R"))

#######################################################################
set.seed(rep_number) ##To ensure that any randomized processes give the same results each time

###Construct Datasets
true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true ranks of overall signals
ToyDat = GenerateToyData(n = n, p1 = p1, p2 = p2, JntVarEx1 = JntVarEx1, JntVarEx2 = JntVarEx2, 
                         IndVarEx1 = IndVarEx1, IndVarEx2 =  IndVarEx2, jnt_rank = r.J,
                         equal.eig = F,ind_rank1 = r.I1, ind_rank2 = r.I2, JntVarAdj = T, SVD.plots = T,
                         Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
## Proportions of groups for mixture
mix.probs = c(0.2, 0.5, 0.3)
diagnoses = factor(rep(1:3, each = n*mix.probs))
blocks <- ToyDat[[2]]
##Setup input parameters to use Gavin's EM code
Y = do.call(cbind, blocks); P = c(p1, p2); Q = c(r.J,r.I1,r.I2)

# Theta
JntScores = ToyDat[['Scores']][['Joint']]
IndivScore.X = ToyDat[['Scores']][["Indiv_1"]]
IndivScore.Y = ToyDat[['Scores']][["Indiv_2"]]

# WJ1 and WJ2
JntLd.X = t(ToyDat[['Loadings']][["Joint_1"]])
JntLd.Y = t(ToyDat[['Loadings']][["Joint_2"]])

# WI1 and WI2
IndivLd.X =t(ToyDat[['Loadings']][["Indiv_1"]])
IndivLd.Y = t(ToyDat[['Loadings']][["Indiv_2"]])

# Signal matrices
JX = ToyDat[[1]]$J1
JY = ToyDat[[1]]$J2
IX = ToyDat[[1]]$I1
IY = ToyDat[[1]]$I2

# Noise Matrices
EX = ToyDat[[1]]$E1
EY = ToyDat[[1]]$E2

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

#### CJIVE
CJIVE.res = cc.jive(blocks, true_signal_ranks, r.J, perm.test = FALSE)

# CJIVE signal matrix estimates
J.hat = CJIVE.res$sJIVE$joint_matrices
I.hat = CJIVE.res$sJIVE$indiv_matrices

# CJIVE loading estimates
WJ = lapply(J.hat, function(x) x[['v']])
WI = lapply(I.hat, function(x) x[['v']])

# Plots of CJIVE estimates against true counterparts
layout(matrix(1:6,2, byrow = TRUE))
plot(JntScores, CJIVE.res$CanCorRes$Jnt_Scores, xlab = "True Joint Scores",ylab = "CJIVE Joint Scores", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, CJIVE.res$CanCorRes$Jnt_Scores), 3)))
plot(JntLd.X, WJ[[1]], xlab = "True Joint Loadings X", ylab = "CJIVE Joint Loadings X", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, WJ[[1]]), 3)))
plot(JntLd.Y, WJ[[2]], xlab = "True Joint Loadings Y", ylab = "CJIVE Joint Loadings Y", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, WJ[[2]]), 3)))
# plot(1,lwd=0,axes=F,xlab="",ylab="", "n")
plot.new(); legend("left", paste("Comp.", 1:2), pch = 1, col  = c("orange", "green"),bty = "n" )
plot(IndivLd.X, WI[[1]], xlab = "True Individual Loadings X", ylab = "CJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, WI[[1]]), 3)))
plot(IndivLd.Y, WI[[2]], xlab = "True Individual Loadings Y", ylab = "CJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, WI[[2]]), 3)))
layout(1)

##### ProJIVE
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

layout(matrix(1:6,2, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores), 3)))
plot(JntLd.X, PJIVE.loads.X[,1:r.J], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J]), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J]), 3)))

# plot(1,lwd=0,axes=F,xlab="",ylab="", "n")
plot.new(); legend("left", paste("Comp.", 1:2), pch = 1, col  = c("orange", "green"),bty = "n" )
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3)))
layout(1)

##### ProJIVE with sigma_hat[k] from pPCA MLE
sigma.hat = c(mean(svd(scale(blocks[[1]]))$d[-(1:(r.J+r.I1))]), 
              mean(svd(scale(blocks[[2]]))$d[-(1:(r.J+r.I2))]))
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5, sig_hat = sigma.hat)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

layout(matrix(1:6,2, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores), 3)))
plot(JntLd.X, PJIVE.loads.X[,1:r.J], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J]), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J]), 3)))

# plot(1,lwd=0,axes=F,xlab="",ylab="", "n")
plot.new(); legend("left", paste("Comp.", 1:2), pch = 1, col  = c("orange", "green"),bty = "n" )
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3)))
layout(1)


