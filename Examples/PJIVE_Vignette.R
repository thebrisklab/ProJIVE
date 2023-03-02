#############################################################################################################################
### Generates a pair of toy data blocks in the manner outlined within the CJIVE manuscript  #################################
###         After data have been generated, apply CJIVE, RJIVE, GIPCA, and ProJIVE          #################################
### Author: Raphiel J. Murden                                                               #################################
### Supervised by Benjamin Risk                                                             #################################
#############################################################################################################################
# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_to_SimulateData.R'
prog.dir = "C:/Users/rmurden/OneDrive - Emory University/Documents/GitHub/ProJIVE"
prog.gipca.dir = "H:/My Documents/P-JIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
library(singR); library(CJIVE); library(reticulate)
prog.dcca.dir = "H:/My Documents/P-JIVE/Programs/D-CCA"
source_python(file.path(prog.dcca.dir, "dcca_given.py"))
gipca.files = list.files(prog.gipca.dir,full.names = TRUE)
for(file.nm in gipca.files){source(file.nm)}
doc.dir = "H:/My Documents/P-JIVE/Programs/Examples/Output"

##Simulation parameters
rep_number = 1
r.J = 3
r.I1 = 2
r.I2 = 2
#outdir = args[2]
n = 1000
p1 = 20
p2 = 200 ####Note that p1 and p2 differ when compared to values used in simulations
JntVarEx1 = 0.5
JntVarEx2 = 0.5
#files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25

#######################################################################
set.seed(rep_number) ##To ensure that any randomized processes give the same results each time

###Construct Datasets
true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true ranks of overall signals
ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), equal.eig = FALSE,
                         IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J, ind_ranks = c(r.I1, r.I2),
                         JntVarAdj = T, SVD.plots = F, Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
## Proportions of groups for mixture
mix.probs = c(0.2, 0.5, 0.3)
diagnoses = factor(c(rep(1, each = n*mix.probs[1]),rep(2, each = n*mix.probs[2]),rep(3, each = n*mix.probs[3])))
blocks <- ToyDat[["Data Blocks"]]
##Setup input parameters to use Gavin's EM code
Y = do.call(cbind, blocks)
# Y = do.call(cbind, lapply(blocks, function(x) scale(x, scale = FALSE)))
P = c(p1, p2); Q = c(r.J,r.I1,r.I2)

# Theta
JntScores = ToyDat[['Scores']][['Joint']]
IndivScore.X = ToyDat[['Scores']][["Indiv"]][,1:Q[2], drop = FALSE]
IndivScore.Y = ToyDat[['Scores']][["Indiv"]][,Q[2]+(1:Q[2]), drop = FALSE]

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

##Looks like the variances are about equal to me ( +/- 10%)?
# plot(apply((blocks[[1]] - AX), 2, sd))
# plot(apply((blocks[[2]] - AY), 2, sd))

#########################################################################################################
############       ProJIVE with loadings and error variance initialized from the truth        ###########
WJ.init = list(JntLd.X, JntLd.Y)
WI.init = list(IndivLd.X, IndivLd.Y)
init.loads = list(WJ.init, WI.init)
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, init.loads = init.loads, sig_hat = c(1,1), plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]
PJIVE.err.var = PJIVE.res$ErrorVariances

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J), drop = FALSE], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1), drop = FALSE], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1), drop = FALSE], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings and error variance \ninitialized  from True Values)", 
                  sub = bquote(sigma["1"]*"="*.(round(PJIVE.err.var[1],2))*"; "*sigma["2"]*"="*.(round(PJIVE.err.var[2],2))))
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

ProJIVE.info = ProJIVE_AsymVar(list(PJIVE.loads.X, PJIVE.loads.Y), PJIVE.err.var, PJIVE.scores, r.J, Y)
show.image.2(ProJIVE.info$ObservedEmpericalInformationMatrix)
show.image.2(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)

show.image.2(pnorm(abs(PJIVE.loads.X/PJIVE.err.var[1]*n), lower.tail = FALSE))
show.image.2(pnorm(abs(PJIVE.loads.Y/PJIVE.err.var[2]*n), lower.tail = FALSE))

#########################################################################################################################
############       ProJIVE with loadings and error variance initialized from the truth and centering Y        ###########
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, init.loads = init.loads, sig_hat = c(1,1), plots = TRUE, center = TRUE)
PJIVE.res = ProJIVE(Y=Y, P=P, Q=Q, init.loads = init.loads, sig_hat = c(1,1), plots = TRUE, num.starts = 2,
                    center = TRUE, return.all.starts = TRUE)

PJIVE.scores = PJIVE.res$ProJIVE_Results[[1]]$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$ProJIVE_Results[[1]]$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$ProJIVE_Results[[1]]$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]
PJIVE.err.var = PJIVE.res$ProJIVE_Results[[1]]$ErrorVariances

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J), drop = FALSE], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1), drop = FALSE], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1), drop = FALSE], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings and error variance \ninitialized  from True Values)", 
                  sub = bquote(sigma["1"]*"="*.(round(PJIVE.err.var[1],2))*"; "*sigma["2"]*"="*.(round(PJIVE.err.var[2],2))))
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

ProJIVE.info = ProJIVE_AsymVar(list(PJIVE.loads.X, PJIVE.loads.Y), PJIVE.err.var, PJIVE.scores, r.J, Y)
show.image.2(ProJIVE.info$ObservedEmpericalInformationMatrix)
show.image.2(ProJIVE.info$Inverse_ObservedEmpericalInformationMatrix)

show.image.2(pnorm(abs(PJIVE.loads.X/PJIVE.err.var[1]*n), lower.tail = FALSE))
show.image.2(pnorm(abs(PJIVE.loads.Y/PJIVE.err.var[2]*n), lower.tail = FALSE))

#########################################################################################################
###########   ProJIVE with loadings initialized from the truth and randomly for error var   #############
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, init.loads = init.loads, sig_hat = "MLE")

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings initialized  \nfrom True Values)")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

#############################################
############### ProJIVE  ####################
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results with \nRandom initial values")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)


#############################################################
#####    ProJIVE with sigma_hat[k] from pPCA MLE   ##########
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5, sig_hat = "MLE", plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]


layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n(noise variance \ninitalized at pPCA MLE)")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)


#################################################################################################
##########   ProJIVE with loadings initialized from Right Singular Vectors solutions   ##########
svd.X1 = svd(blocks[[1]])
svd.X2 = svd(blocks[[2]])

WJ.init = list(svd.X1$v[,1:r.J, drop = FALSE], svd.X2$v[,1:r.J, drop = FALSE])
WI.init = list(svd.X1$v[,1:r.I1], svd.X2$v[,1:r.I2])
init.loads = list(WJ.init, WI.init)
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5, init.loads = init.loads, plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]


layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings initialized  \nat Right Singular Vectors)")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

#################################################################################################
##########   ProJIVE with loadings initialized from Right Singular Vectors solutions   ##########
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5, init.loads = init.loads, sig_hat = "MLE", plots = FALSE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]


layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings initialized  \nat Right Singular Vectors\n and noise variance\n at pPCA MLE)")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

#########################################################################
##########   ProJIVE with loadings initialized at CJIVE soln   ##########
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5, init.loads = "CJIVE", plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings initialized  \nfrom CJIVE solution")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)


##################################################################################################
##########   ProJIVE with loadings initialized at CJIVE soln and sig_hat from pPCA MLE  ##########
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-5, init.loads = "CJIVE", sig_hat = "MLE", plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]


layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("ProJIVE Results \n (loadings initialized  \nat CJIVE solution and \nnoise variance at pPCA MLE")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

####################################################################################################################
##########   ProJIVE with loadings initialized at CJIVE soln sig_hat from pPCA MLE and misspecified ranks ##########
PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q+1, diff.tol=1e-5, init.loads = "CJIVE", sig_hat = "MLE", plots = TRUE)

PJIVE.scores = PJIVE.res$SubjectScoreMatrix
PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+1:r.I1)]

r.J.wrong = r.J+1
r.I1.wrong = r.I1+1
r.I2.wrong = r.I2+1

layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, PJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, PJIVE.scores[,1:(r.J), drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = PJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
points(JntScores, PJIVE.scores[,2])
plot(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "ProJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, PJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = PJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
points(IndivScore.X, PJIVE.scores[,3])
plot(IndivScore.Y, PJIVE.scores[,5], xlab = "True Indiv 2 Scores", ylab = "ProJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, PJIVE.scores[,5]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = PJIVE.scores[,6], standardize = TRUE), 3)))
points(IndivScore.Y, PJIVE.scores[,6])
plot.new(); title("ProJIVE Results \n (loadings initialized  \nat CJIVE solution and \nnoise variance at pPCA MLE")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "ProJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, PJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = PJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "ProJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, PJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = PJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, PJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "ProJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, PJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = PJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "ProJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, PJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = PJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

r.jive::show.image.2(pnorm(abs(PJIVE.loads.X/PJIVE.err.var[1]*n), lower.tail = FALSE)<0.05)
show.image.2(pnorm(abs(PJIVE.loads.Y/PJIVE.err.var[1]*n), lower.tail = FALSE)<0.05)


#######################################
##############  CJIVE   ###############
CJIVE.res = cc.jive(blocks, true_signal_ranks, r.J, perm.test = FALSE)

# CJIVE loading estimates
CJIVE.scores = cbind(CJIVE.res$CanCorRes$Jnt_Scores, CJIVE.res$sJIVE$indiv_matrices[[1]]$u, CJIVE.res$sJIVE$indiv_matrices[[2]]$u)
CJIVE.loads.X = cbind(CJIVE.res$sJIVE$joint_matrices[[1]]$v, CJIVE.res$sJIVE$indiv_matrices[[1]]$v)
CJIVE.loads.Y = cbind(CJIVE.res$sJIVE$joint_matrices[[2]]$v, CJIVE.res$sJIVE$indiv_matrices[[2]]$v)

# Plots of CJIVE estimates against true counterparts
layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, CJIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "CJIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, CJIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = CJIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, CJIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "CJIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, CJIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = CJIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, CJIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "CJIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, CJIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = CJIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("CJIVE Results")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, CJIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "CJIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, CJIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = CJIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, CJIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "CJIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, CJIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = CJIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, CJIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "CJIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, CJIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = CJIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, CJIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "CJIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, CJIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = CJIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

########################################
##############  R.JIVE   ###############
blocks.t = lapply(blocks, t)
R.JIVE.res = jive(blocks.t, rankJ = r.J, rankA = true_signal_ranks - r.J, method = "given")

# R.JIVE signal matrix estimates
J.hat = lapply(R.JIVE.res$joint, svd, nu = r.J, nv = r.J)
I.hat = lapply(R.JIVE.res$individual, svd, nu = r.I1, r.I2)

R.JIVE.jnt.scores = svd(do.call(rbind, R.JIVE.res$joint), nu = r.J, nv = r.J)$v

# R.JIVE loading estimates
WJ = lapply(J.hat, function(x) x[['u']])
WI = lapply(I.hat, function(x) x[['u']])

R.JIVE.scores = cbind(R.JIVE.jnt.scores, I.hat[[1]]$v, I.hat[[2]]$v)
R.JIVE.loads.X = cbind(WJ[[1]],WI[[1]])
R.JIVE.loads.Y = cbind(WJ[[2]],WI[[2]])

# Plots of R.JIVE estimates against true counterparts
layout(matrix(c(1:6,4,7,8),3, byrow = TRUE))
plot(JntScores, R.JIVE.scores[,1:r.J, drop = FALSE], xlab = "True Joint Scores", ylab = "R.JIVE Joint Scores", 
     col = c(rep("orange",n), rep("green",n), rep("purple",n)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, R.JIVE.scores[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntScores, S2 = R.JIVE.scores[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivScore.X, R.JIVE.scores[,r.J+(1:r.I1)], xlab = "True Indiv 1 Scores", ylab = "R.JIVE Indiv 1 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.X, R.JIVE.scores[,r.J+(1:r.I1)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.X, S2 = R.JIVE.scores[,r.J+(1:r.I1)], standardize = TRUE), 3)))
plot(IndivScore.Y, R.JIVE.scores[,(r.J+r.I1)+(1:r.I2)], xlab = "True Indiv 2 Scores", ylab = "R.JIVE Indiv 2 Scores", 
     col = c(rep("orange",n), rep("green",n)),
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivScore.Y, R.JIVE.scores[,(r.J+r.I1)+(1:r.I2)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivScore.Y, S2 = R.JIVE.scores[,(r.J+r.I1)+(1:r.I2)], standardize = TRUE), 3)))
plot.new(); title("R.JIVE Results")
legend("left", paste("Comp.", 1:3), pch = 1, col  = c("orange", "green", "purple"),bty = "n" )
plot(JntLd.X, R.JIVE.loads.X[,1:r.J, drop = FALSE], xlab = "True Joint Loadings X", ylab = "R.JIVE Joint Loadings X", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, R.JIVE.loads.X[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.X, S2 = R.JIVE.loads.X[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(JntLd.Y, R.JIVE.loads.Y[,1:r.J, drop = FALSE], xlab = "True Joint Loadings Y", ylab = "R.JIVE Joint Loadings Y", 
     col = c(rep("orange",p1), rep("green",p1), rep("purple",p1)),  
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, R.JIVE.loads.Y[,1:r.J, drop = FALSE]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = JntLd.Y, S2 = R.JIVE.loads.Y[,1:r.J, drop = FALSE], standardize = TRUE), 3)))
plot(IndivLd.X, R.JIVE.loads.X[,-(1:r.J)], xlab = "True Individual Loadings X", ylab = "R.JIVE Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, R.JIVE.loads.X[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.X, S2 = R.JIVE.loads.X[,-(1:r.J)], standardize = TRUE), 3)))
plot(IndivLd.Y, R.JIVE.loads.Y[,-(1:r.J)], xlab = "True Individual Loadings Y", ylab = "R.JIVE Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, R.JIVE.loads.Y[,-(1:r.J)]), 3), ":",
                  " PMSE = ", round(pmse.2(S1 = IndivLd.Y, S2 = R.JIVE.loads.Y[,-(1:r.J)], standardize = TRUE), 3)))
layout(1)

################################################################
######################     GIPCA      ##########################
GIPCA.res = EPCAJIVEMissbio(blocks, r.J, c(r.I1,r.I2), D = P, family = rep("gaussian",2), tol = 1E-5)
# GIPCABIC(blocks, ranka = c(r.I1,r.I2), rankj = r.J)
GIPCA.joint.scores = GIPCA.res$U.joint[,-1, drop = FALSE]
GIPCA.joint.loads.X = t(GIPCA.res$V.joint[-1,1:P[1], drop = FALSE])
GIPCA.joint.loads.Y = t(GIPCA.res$V.joint[-1,-(1:P[1]), drop = FALSE])

GIPCA.indiv.loads.X = t(GIPCA.res$V.ind[[1]])
GIPCA.indiv.loads.Y = t(GIPCA.res$V.ind[[2]])

layout(matrix(1:6,2, byrow = TRUE))
plot(JntScores, GIPCA.joint.scores, xlab = "True Joint Scores", ylab = "GIPCA Joint Scores", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, GIPCA.joint.scores), 3)))
plot(JntLd.X, GIPCA.joint.loads.X, xlab = "True Joint Loadings X", ylab = "GIPCA Joint Loadings X", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, GIPCA.joint.loads.X), 3)))
plot(JntLd.Y, GIPCA.joint.loads.Y, xlab = "True Joint Loadings Y", ylab = "GIPCA Joint Loadings Y", 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, GIPCA.joint.loads.Y), 3)))

# plot(1,lwd=0,axes=F,xlab="",ylab="", "n")
plot.new(); title("GIPCA Results") 
legend("left", paste("Comp.", 1:2), pch = 1, col  = c("orange", "green"),bty = "n" )
plot(IndivLd.X, GIPCA.indiv.loads.X, xlab = "True Individual Loadings X", ylab = "GIPCA Individual Loadings X", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, GIPCA.indiv.loads.X), 3)))
plot(IndivLd.Y, GIPCA.indiv.loads.Y, xlab = "True Individual Loadings Y", ylab = "GIPCA Individual Loadings Y", 
     col = c(rep("orange",p1), rep("green",p2)), 
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, GIPCA.indiv.loads.Y), 3)))
layout(1)


##################################################################
######################     ARCHIVE      ##########################
relnormF.X1 = norm(AX - PJIVE.res$SubjectScoreMatrix[,1:(r.J+r.I1)]%*%t(PJIVE.res$LoadingMatrix[1:p1,1:(r.J+r.I1)]), tpye = "F")/norm(AX, tpye = "F")
relnormF.X2 = norm(AY - PJIVE.res$SubjectScoreMatrix[,-c(r.J+(1:r.I1))]%*%t(PJIVE.res$LoadingMatrix[-(1:p1),-c(r.J+(1:r.I1))]), tpye = "F")/norm(AY, tpye = "F")

relnormF.J.X1 = norm(JX - PJIVE.res$SubjectScoreMatrix[,1:(r.J)]%*%t(PJIVE.res$LoadingMatrix[1:p1,1:r.J, drop = FALSE]), tpye = "F")/norm(JX, tpye = "F")
relnormF.J.X2 = norm(JY - PJIVE.res$SubjectScoreMatrix[,1:(r.J)]%*%t(PJIVE.res$LoadingMatrix[-(1:p1),1:r.J, drop = FALSE]), tpye = "F")/norm(JY, tpye = "F")

relnormF.I.X1 = norm(IX - PJIVE.res$SubjectScoreMatrix[,r.J+(1:r.I1)]%*%t(PJIVE.res$LoadingMatrix[1:p1,r.J+(1:r.I1)]), tpye = "F")/norm(IX, tpye = "F")
relnormF.I.X2 = norm(IY - PJIVE.res$SubjectScoreMatrix[,(r.J+r.I1)+(1:r.I2)]%*%t(PJIVE.res$LoadingMatrix[-(1:p1),(r.J+r.I1)+(1:r.I2)]), tpye = "F")/norm(IY, tpye = "F")

relnormS.X1 = norm(AX - PJIVE.res$SubjectScoreMatrix[,1:(r.J+r.I1)]%*%t(PJIVE.res$LoadingMatrix[1:p1,1:(r.J+r.I1)]), tpye = "2")/norm(AX, tpye = "2")
relnormS.X2 = norm(AY - PJIVE.res$SubjectScoreMatrix[,-c(r.J+(1:r.I1))]%*%t(PJIVE.res$LoadingMatrix[-(1:p1),-c(r.J+(1:r.I1))]), tpye = "2")/norm(AY, tpye = "2")

relnormS.J.X1 = norm(JX - PJIVE.res$SubjectScoreMatrix[,1:(r.J)]%*%t(PJIVE.res$LoadingMatrix[1:p1,1:r.J, drop = FALSE]), tpye = "2")/norm(JX, tpye = "2")
relnormS.J.X2 = norm(JY - PJIVE.res$SubjectScoreMatrix[,1:(r.J)]%*%t(PJIVE.res$LoadingMatrix[-(1:p1),1:r.J, drop = FALSE]), tpye = "2")/norm(JY, tpye = "2")

relnormS.I.X1 = norm(IX - PJIVE.res$SubjectScoreMatrix[,r.J+(1:r.I1)]%*%t(PJIVE.res$LoadingMatrix[1:p1,r.J+(1:r.I1)]), tpye = "2")/norm(IX, tpye = "2")
relnormS.I.X2 = norm(IY - PJIVE.res$SubjectScoreMatrix[,(r.J+r.I1)+(1:r.I2)]%*%t(PJIVE.res$LoadingMatrix[-(1:p1),(r.J+r.I1)+(1:r.I2)]), tpye = "2")/norm(IY, tpye = "2")

