#############################################################################################################################
### Generates a pair of toy data blocks in the manner outlined within the CJIVE manuscript  #################################
###         After data have been generated, apply CJIVE, RJIVE, GIPCA, and ProJIVE          #################################
### Author: Raphiel J. Murden                                                               #################################
### Supervised by Benjamin Risk                                                             #################################
#############################################################################################################################
# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_to_SimulateData.R'
prog.dir = "H:/My Documents/ProJIVE/Programs/Functions"
prog.gipca.dir = "H:/My Documents/ProJIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
# source(file.path(prog.dir, "Functions_for_CJIVE.R"))
# source(file.path(prog.dir, "Functions_for_PJIVE.R"))
source(file.path(prog.dir, "ProJIVE_Kge2.R"))
gipca.files = list.files(prog.gipca.dir,full.names = TRUE)
lapply(gipca.files, source)
require(stringr); require(CJIVE)
ajive.dir = "H:/My Documents/Applications2/r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

w_to_w_k=function(w, P, Q){
  
  P=c(0,P)
  K=length(P[-1])
  w_k=list()
  for(k in 1:K){
    w_k[[k]] = w[sum(P[1:k])+(1:P[k+1]),c(1:Q[1],sum(Q[1:k])+1:Q[k+1])]
  }
  
  
  return(w_k)
}

rep_number = 0
r.J = 1
r.I = c(2, 2, 2, 3)
n = 150
p = c(20, 50, 60, 70) 
K = length(p)
JntVarEx = c(0.5, 0.15, 0.15, 0.05)
IndVarEx = rep(0.25, 4)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb.cols = cbPalette[c(8,2,3,5)]
#######################################################################
set.seed(rep_number) ##To ensure that any randomized processes give the same results each time

###Construct Datasets
true_signal_ranks = r.J + r.I                          ##true ranks of overall signals
ToyDat = GenerateToyData_Kge2(n = n, p = p, JntVarEx = JntVarEx, IndVarEx = IndVarEx,
                         jnt_rank = r.J, equal.eig = F, ind_ranks = r.I,
                         JntVarAdj = T, SVD.plots = T, 
                         Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
blocks <- ToyDat[[2]]
##Setup input parameters to use Gavin's EM code

# Theta
JntScores = ToyDat[['Scores']][['Joint']]
IndScores = ToyDat[['Scores']][['Indiv']]

X1.loads = t(rbind(ToyDat$Loadings$Joint[[1]], ToyDat$Loadings$Indiv[[1]]))
X2.loads = t(rbind(ToyDat$Loadings$Joint[[2]], ToyDat$Loadings$Indiv[[2]]))
X3.loads = t(rbind(ToyDat$Loadings$Joint[[3]], ToyDat$Loadings$Indiv[[3]]))

Y = do.call(cbind, blocks);
P = p; Q = c(r.J,r.I)
#############################################
############### ProJIVE  ####################
PJIVE.res = ProJIVE_EM_Kge2(Y=Y, P=P, Q=Q, diff.tol=1e-7, sig_hat = "MLE", init.loads = "AJIVE")
# PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, diff.tol=1e-7, sig_hat = "MLE", init.loads = "AJIVE")
JntScores = ToyDat$Scores$Joint
ProJIVE.JntScores = PJIVE.res$SubjectScoreMatrix[,1:r.J]
PJIVE.IndScores = PJIVE.res$SubjectScoreMatrix[,-(1:r.J)]
ProJIVE.loads = w_to_w_k(PJIVE.res$LoadingMatrix, P, Q)
ProJIVE.X1.loads = ProJIVE.loads[[1]]
ProJIVE.X2.loads = ProJIVE.loads[[2]]
ProJIVE.X3.loads = ProJIVE.loads[[3]]

plot(JntScores,ProJIVE.JntScores, main = "Joint Scores", xlab = "True Joint Scores", ylab = "ProJIVE Joint Scores",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores,ProJIVE.JntScores),2)))

chord.val = chord.norm.diff(ToyDat$Scores$Indiv[,str_which(colnames(ToyDat$Scores$Indiv), "Indiv X1")], 
                            PJIVE.IndScores[,str_which(colnames(PJIVE.IndScores), "Data1")])
plot(ToyDat$Scores$Indiv[,str_which(colnames(ToyDat$Scores$Indiv), "Indiv X1")], 
     PJIVE.IndScores[,str_which(colnames(PJIVE.IndScores), "Data1")],
     sub = paste0("Chordal Norm = ", round(chord.val,2)),
     xlab = "True X1 Indiv Scores", ylab = "ProJIVE X1 Indiv Scores", 
     col = rep(cb.cols[1:Q[3]],rep(P[2],Q[3])))

chord.val = chord.norm.diff(ToyDat$Scores$Indiv[,str_which(colnames(ToyDat$Scores$Indiv), "Indiv X2")], 
                            PJIVE.IndScores[,str_which(colnames(PJIVE.IndScores), "Data2")])
plot(ToyDat$Scores$Indiv[,str_which(colnames(ToyDat$Scores$Indiv), "Indiv X2")], 
     PJIVE.IndScores[,str_which(colnames(PJIVE.IndScores), "Data2")],
     sub = paste0("Chordal Norm = ", round(chord.val,2)),
     xlab = "True X2 Indiv Scores", ylab = "ProJIVE X2 Indiv Scores", 
     col = rep(cb.cols[1:Q[3]],rep(P[2],Q[3])))
chord.val = chord.norm.diff(ToyDat$Scores$Indiv[,str_which(colnames(ToyDat$Scores$Indiv), "Indiv X3")], 
                            PJIVE.IndScores[,str_which(colnames(PJIVE.IndScores), "Data3")])
plot(ToyDat$Scores$Indiv[,str_which(colnames(ToyDat$Scores$Indiv), "Indiv X3")], 
     PJIVE.IndScores[,str_which(colnames(PJIVE.IndScores), "Data3")],
     sub = paste0("Chordal Norm = ", round(chord.val,2)),
     xlab = "True X3 Indiv Scores", ylab = "ProJIVE X3 Indiv Scores", 
     col = rep(cb.cols[1:Q[3]],rep(P[2],Q[3])))


plot(X1.loads[,1],ProJIVE.X1.loads[,1], main = "X1 Joint Loadings", xlab = "True Joint Loadings", ylab = "ProJIVE Joint Loadings",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(X1.loads[,1],ProJIVE.X1.loads[,1]),2)))
plot(X2.loads[,1],ProJIVE.X2.loads[,1], main = "X2 Joint Loadings", xlab = "True Joint Loadings", ylab = "ProJIVE Joint Loadings",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(X2.loads[,1],ProJIVE.X2.loads[,1]),2)))
plot(X3.loads[,1],ProJIVE.X3.loads[,1], main = "X3 Joint Loadings", xlab = "True Joint Loadings", ylab = "ProJIVE Joint Loadings",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(X3.loads[,1],ProJIVE.X3.loads[,1]),2)))

plot(X1.loads[,-1],ProJIVE.X1.loads[,-1], main = "X1 Indiv Loadings", xlab = "True Indiv Loadings", ylab = "ProJIVE Indiv Loadings",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(X1.loads[,-1],ProJIVE.X1.loads[,-1]),2)),
     col = rep(cb.cols[2:3],c(P[1],P[1])))
plot(X2.loads[,-1],ProJIVE.X2.loads[,-1], main = "X2 Indiv Loadings", xlab = "True Indiv Loadings", ylab = "ProJIVE Indiv Loadings",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(X2.loads[,-1],ProJIVE.X2.loads[,-1]),2)),
     col = rep(cb.cols[2:3],rep(P[2],2)))
plot(X3.loads[,-1],ProJIVE.X3.loads[,-1], main = "X3 Indiv Loadings", xlab = "True Indiv Loadings", ylab = "ProJIVE Indiv Loadings",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(X3.loads[,-1],ProJIVE.X3.loads[,-1]),2)),
     col = rep(cb.cols[2:3],c(P[3],P[3])))

# # layout(matrix(1:6, nrow = 2))
# for(k in 1:ncol(IndScores)){
#   plot(IndScores[,k], PJIVE.IndScores[,k], main = paste0("Indiv Scores: Data set #", round(k/2), " Component ", ifelse(k%%2 ==0, 2,1)),
#        sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndScores[,k], PJIVE.IndScores[,k]),2)))
# }
# for(k in 1:3){
#   print(paste0("Chordal Norm = ", round(chord.norm.diff(IndScores[,(2*k-1):(2*k)], PJIVE.IndScores[,(2*k-1):(2*k)]),2)))
# }
# layout(1)


#######################################################################################
###############################     R.JIVE      #######################################
blocks.t = lapply(blocks, t)
R.JIVE.res = jive(blocks.t, rankJ = r.J, rankA = Q[-1], method = "given")

R.JntScores =  svd(do.call(rbind, R.JIVE.res$joint), nu = r.J, nv = r.J)$v
R.IndScores = do.call(cbind, lapply(R.JIVE.res$individual, function(x) svd(x, nu = 2, nv = 2)$v))
R.J.svd = R.I.svd = list()
for(k in 1:length(P)){
  R.J.svd[[k]] = svd(R.JIVE.res$joint[[k]], nu = Q[1], nv = Q[1])
  R.I.svd[[k]] = svd(R.JIVE.res$individual[[k]], nu = Q[k+1], nv = Q[k+1])
}

plot(JntScores, R.JntScores, main = "Joint Scores",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, R.JntScores),2)))
# layout(matrix(1:6, nrow = 2))

chord.val = chord.norm.diff(IndScores[,str_which(colnames(IndScores), "Indiv X1")],
                            R.I.svd[[1]]$v)
plot(IndScores[,str_which(colnames(IndScores), "Indiv X1")],
     R.I.svd[[1]]$v,
     main = "Indiv X1 Scores",
     sub = paste0("Chordal Norm = ", round(chord.val,2)),
     xlab = "True X1 Indiv Scores", ylab = "R.JIVE X1 Indiv Scores", 
     col = rep(cb.cols[1:Q[2]],rep(P[1],Q[2])))

chord.val = chord.norm.diff(IndScores[,str_which(colnames(IndScores), "Indiv X2")],
                            R.I.svd[[2]]$v)
plot(IndScores[,str_which(colnames(IndScores), "Indiv X2")],
     R.I.svd[[2]]$v,
     main = "Indiv X2 Scores",
     sub = paste0("Chordal Norm = ", round(chord.val,2)),
     xlab = "True X2 Indiv Scores", ylab = "R.JIVE X2 Indiv Scores", 
     col = rep(cb.cols[1:Q[3]],rep(P[2],Q[3])))

chord.val = chord.norm.diff(IndScores[,str_which(colnames(IndScores), "Indiv X3")],
                            R.I.svd[[3]]$v)
plot(IndScores[,str_which(colnames(IndScores), "Indiv X3")],
     R.I.svd[[3]]$v,
     main = "Indiv X3 Scores",
     sub = paste0("Chordal Norm = ", round(chord.val,2)),
     xlab = "True X3 Indiv Scores", ylab = "R.JIVE X3 Indiv Scores", 
     col = rep(cb.cols[1:Q[3]],rep(P[2],Q[3])))

#######################################################################################
###############################     AJIVE      #######################################
AJIVE.res<-ajive(blocks,initial_signal_ranks = Q[1]+Q[-1], joint_rank = r.J)

A.JntScores = AJIVE.res$joint_scores
A.IndScores = do.call(cbind, lapply(AJIVE.res$block_decomps, function(x) x$individual$u))

plot(JntScores, A.JntScores, main = "Joint Scores",
     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntScores, A.JntScores),2)))
layout(matrix(1:6, nrow = 2))
for(k in 1:ncol(IndScores)){
  plot(IndScores[,k], A.IndScores[,k], main = paste0("Indiv Scores: Data set #", round(k/2), " Component ", ifelse(k%%2 ==0, 2,1)),
       sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndScores[,k], A.IndScores[,k]),2)))
}
layout(1)

for(k in 1:3){
  print(paste0("Chordal Norm = ", round(chord.norm.diff(IndScores[,(2*k-1):(2*k)], A.IndScores[,(2*k-1):(2*k)]),2)))
}

IndScores = cbind(ToyDat$Scores$Indiv_1, ToyDat$Scores$Indiv_1, ToyDat$Scores$Indiv_2)
PJIVE.IndScores = PJIVE.res$SubjectScoreMatrix[,-(1:r.J)]

layout(matrix(1:6, nrow = 2))
for(k in 1:ncol(IndScores)){
  plot(IndScores[,k], PJIVE.IndScores[,k], main = paste0("Indiv Scores: Data set #", round(k/2), " Component ", ifelse(k%%2 ==0, 2,1)),
       sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndScores[,k], PJIVE.IndScores[,k]),2)))
}




