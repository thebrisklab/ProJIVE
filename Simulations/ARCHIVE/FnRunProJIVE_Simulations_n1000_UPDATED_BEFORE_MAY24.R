#############################################################################################################################
### Compare methods of implementing JIVE analysis on two different pairs of toy datasets ####################################
### Author: Raphiel J. Murden                                                            ####################################
### Supervised by Benjamin Risk                                                          ####################################
#############################################################################################################################

#####################################################      NOTES      #######################################################
### 11MAR2022: While looking into the issue of how different results were obtained in the mixed variability simulations, 
###                I discovered that the solutions are remarkably consistent and symmetric if CJIVE (i.e. AJIVE) loadings are
###                 used to initiate the EM algorithm. I will run simulations again with, but this time initiating loadings at  
###                 at the CJIVE/AJIVE solution and error variances at their pPCA MLEs
###
### 29NOV2022: Updating simulation study to include dCCA as a competing method. Note: this version of dCCA has been modified 
###                to NOT include all of the rank choice methods discussed in their paper. Instead, we force the joint rank to 
###                it's true value and the total signal rank to it's true value
### 06DEC2022: R version of dCCA (programed by RJM) did not seem to work properly. Updating this script to use the Python version
###                programmed by the authors article and available on their GitHub page (https://github.com/shu-hai/D-CCA/blob/master/dcca.py)
###                via the 'reticulate' R package. 
### 25JAN2023: Including PMSE (Permutation-invariant mean squared error) as a measure of estimation accuracy.See Risk, Matterson, and Rupert 2019
### 31JAN2023: Added "drop = FALSE" for subsetting columns of score and loading matrices with ranks is 1. 
###             This helps to ensure that PMSE is calculated and doesn't cause the entire simulation to fail
#############################################################################################################################
args = commandArgs(trailingOnly=TRUE)
# args = args[-1]
args
rep_number = as.numeric(args[1])
outdir = args[2]
n = as.numeric(args[3])
p1 = as.numeric(args[4])
p2 = as.numeric(args[5])
JntVarEx1 = as.numeric(args[6])
JntVarEx2 = as.numeric(args[7])
files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25

r.J = as.numeric(args[8])
r.I1 = 2
r.I2 = 2

nm = paste("sim_", JntVarEx1, "_", JntVarEx2, rep_number, ".csv", sep="")
if (!(nm %in% files)){
  
  library(r.jive); library(reticulate); library(CJIVE); library(singR)
  prog.dir = "/home/rmurden/PJIVE/Programs"
  prog.gipca.dir = "/home/rmurden/PJIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
  prog.dcca.dir = "/projects/guo_lab/cbis/users/raphiel/PJIVE/Programs/D-CCA"
  files= list.files(file.path(prog.dir, "R"))
  for (i in files) source(file.path(prog.dir, 'R', i))
  source(file.path(prog.dir, "Functions_for_PJIVE.R"))
  source_python(file.path(prog.dcca.dir, "dcca_given.py"))
  
  gipca.files = list.files(prog.gipca.dir,full.names = TRUE)
  for(i in gipca.files) source(i)
  
  TIME = Sys.time()
  
  #######################################################################
  print(paste0("Is rep_number=",rep_number," a double? ",is.double(rep_number)))
  set.seed(round(rep_number)) ##To ensure that any randomized processes give the same results each time
  ###Construct Datasets
  #####Add noise to already generated Joint and individual signal matrices for X, Y data blocks:
  
  #######JIVE Implementations for Toy Data No 1################################################################################
  ###Setup input parameters for JIVE implementations
  true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true jranks of overall signals
  ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), 
                           IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J,
                           equal.eig = F, ind_ranks = c(r.I1, r.I2), JntVarAdj = T, SVD.plots = F,
                           Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
  
  blocks <- ToyDat[[2]]
  rnd.smp = sample(n, n/2)
  blocks.sub1 = lapply(blocks, function(x){x[rnd.smp,]})
  blocks.sub2 = lapply(blocks, function(x){x[-rnd.smp,]})
  blocks.t <- lapply(blocks,t)   ##iJIVE and dCCA requires datasets organized as p-by-n matrices
  JntScores = as.matrix(ToyDat[['Scores']][['Joint']], nrow = n)
  IndivScore.X = ToyDat[['Scores']]$Indiv[,1:r.I1, drop = FALSE]
  IndivScore.Y = ToyDat[['Scores']]$Indiv[,r.I1+(1:r.I2), drop = FALSE]
  
  JntLd.X = t(ToyDat[['Loadings']][["Joint"]][[1]])
  JntLd.Y = t(ToyDat[['Loadings']][["Joint"]][[2]])
  IndivLd.X =t(ToyDat[['Loadings']][["Indiv"]][[1]])
  IndivLd.Y = t(ToyDat[['Loadings']][["Indiv"]][[2]])
  
  JX = ToyDat$`Data Components`$JointSignalMatrices[[1]]
  JY = ToyDat$`Data Components`$JointSignalMatrices[[2]]
  IX = ToyDat$`Data Components`$IndivSignalMatrices[[1]]
  IY = ToyDat$`Data Components`$IndivSignalMatrices[[2]]
  EX = ToyDat$`Data Components`$NoiseMatrices[[1]]
  EY = ToyDat$`Data Components`$NoiseMatrices[[2]]
  
  AX = JX + IX
  AY = JY + IY
  JntBlock = cbind(JX,JY)
  
  TotVar.X = MatVar(ToyDat[[2]][[1]])
  
  TotVar.Y = MatVar(ToyDat[[2]][[2]])
  
  JVE.X = MatVar(JX)/TotVar.X
  JVE.Y = MatVar(JY)/TotVar.Y
  
  IVE.X = MatVar(IX)/TotVar.X
  IVE.Y = MatVar(IY)/TotVar.Y
  
  TotVE.X = MatVar((JX + IX))/TotVar.X
  TotVE.Y = MatVar((JY + IY))/TotVar.Y
  
  #### Apply R.JIVE methodology
  i.time = system.time({jive.unknown.perm<-jive(blocks.t,rankJ = r.J, rankA = c(true_signal_ranks-r.J), method="given")})
  print(paste("R.JIVE done in", round(i.time['elapsed'], 3), "sec."))
  
  #####Retrieve restuls from iJIVE package with ranks not specified: use permutation method to find them
  i.rJ = jive.unknown.perm[["rankJ"]]
  i.rI = jive.unknown.perm[["rankA"]]
  
  i.JX.hat = jive.unknown.perm$joint[[1]]
  i.JY.hat = jive.unknown.perm$joint[[2]]
  i.JntBlock = rbind(i.JX.hat, i.JY.hat)
  i.JntScores.hat= svd(i.JntBlock, nu = i.rJ, nv = i.rJ)[['v']]
  
  if(i.rJ == 0){
    i.JntScores.hat = rep(0,n)
    i.JX.hat = matrix(0, nrow = n, ncol = p1)
    i.JY.hat = matrix(0, nrow = n, ncol = p2)
    
    i.jnt.subnorm = chord.norm.diff(JntScores, i.JntScores.hat)
    i.jnt.loadX.norm = chord.norm.diff(JntLd.X, t(i.JX.hat))
    i.jnt.loadY.norm = chord.norm.diff(JntLd.Y, t(i.JY.hat))
    
    i.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = i.JntScores.hat)
    i.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = t(i.JX.hat))
    i.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = t(i.JY.hat))
  } else {
    i.JX.hat = jive.unknown.perm$joint[[1]]
    i.JY.hat = jive.unknown.perm$joint[[2]]
    i.jnt.loadX = svd(i.JX.hat, nu = i.rJ, nv = i.rJ)[['u']]
    i.jnt.loadY = svd(i.JY.hat, nu = i.rJ, nv = i.rJ)[['u']]
    
    i.jnt.subnorm = chord.norm.diff(JntScores, i.JntScores.hat)
    i.jnt.loadX.norm = chord.norm.diff(JntLd.X, i.jnt.loadX)
    i.jnt.loadY.norm = chord.norm.diff(JntLd.Y, i.jnt.loadY)
    
    i.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = i.JntScores.hat)
    i.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = i.jnt.loadX)
    i.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = i.jnt.loadY)
  } 
  
  i.IX.hat = jive.unknown.perm$individual[[1]]
  i.IY.hat = jive.unknown.perm$individual[[2]]
  
  i.IX.svd = svd(i.IX.hat, nu = i.rI[1], nv = i.rI[1])
  i.IY.svd = svd(i.IY.hat, nu = i.rI[2], nv = i.rI[2])
  i.bX.hat = i.IX.svd[['v']]; i.bY.hat = i.IY.svd[['v']]
  i.W.IX.hat = i.IX.svd[['u']]; i.W.IY.hat = i.IY.svd[['u']]
  
  i.indiv.X.subnorm = chord.norm.diff(IndivScore.X, i.bX.hat)
  i.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, i.bY.hat)
  i.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, i.W.IX.hat)
  i.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, i.W.IY.hat)

  i.indiv.X.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.X, S2 = i.bX.hat)
  i.indiv.Y.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.Y, S2 = i.bY.hat)
  i.indiv.X.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.X, S2 = i.W.IX.hat)
  i.indiv.Y.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.Y, S2 = i.W.IY.hat)
  
  i.TotVar.X = sum(svd(jive.unknown.perm$data[[1]])[['d']]^2)
  i.TotVar.Y = sum(svd(jive.unknown.perm$data[[2]])[['d']]^2)
  
  i.jnt.VarEx.X = sum(svd(i.JX.hat)[['d']]^2)/i.TotVar.X
  i.jnt.VarEx.Y = sum(svd(i.JY.hat)[['d']]^2)/i.TotVar.Y
  i.ind.VarEx.X = sum(svd(i.IX.hat)[['d']]^2)/i.TotVar.X
  i.ind.VarEx.Y = sum(svd(i.IY.hat)[['d']]^2)/i.TotVar.Y
  
  #### Apply AJIVE-Oracle
  a.c.time = system.time({ajive.oracle<-cc.jive(dat.blocks = blocks, signal.ranks = true_signal_ranks, joint.rank = r.J)})
  print(paste("AJIVE done in", round(a.c.time['elapsed'], 3), "sec."))
  
  #####Retrieve results from AJIVE_Oracle with ranks not specified: use permutation method to find them
  a.c.JntScores.hat = ajive.oracle$CanCorRes$Jnt_Scores
  a.c.jnt.loadX = ajive.oracle$sJIVE$joint_matrices[[1]]$v
  a.c.jnt.loadY = ajive.oracle$sJIVE$joint_matrices[[2]]$v
  
  a.c.jnt.subnorm = chord.norm.diff(JntScores, a.c.JntScores.hat)
  a.c.jnt.loadX.norm = chord.norm.diff(JntLd.X, a.c.jnt.loadX)
  a.c.jnt.loadY.norm = chord.norm.diff(JntLd.Y, a.c.jnt.loadY)
  
  a.c.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = a.c.JntScores.hat)
  a.c.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = a.c.jnt.loadX)
  a.c.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = a.c.jnt.loadY)
  
  a.c.W.IX.hat = ajive.oracle$sJIVE$indiv_matrices[[1]]$v
  a.c.W.IY.hat = ajive.oracle$sJIVE$indiv_matrices[[2]]$v
  a.c.bX.hat = ajive.oracle$sJIVE$indiv_matrices[[1]]$u
  a.c.bY.hat = ajive.oracle$sJIVE$indiv_matrices[[2]]$u
  
  a.c.indiv.X.subnorm = chord.norm.diff(IndivScore.X, a.c.bX.hat)
  a.c.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, a.c.bY.hat)
  a.c.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, a.c.W.IX.hat)
  a.c.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, a.c.W.IY.hat)
  
  a.c.indiv.X.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.X, S2 = a.c.bX.hat)
  a.c.indiv.Y.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.Y, S2 = a.c.bY.hat)
  a.c.indiv.X.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.X, S2 = a.c.W.IX.hat)
  a.c.indiv.Y.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.Y, S2 = a.c.W.IY.hat)
  
  a.c.jnt.VarEx.X = sum(ajive.oracle$sJIVE$joint_matrices[[1]]$d^2)/TotVar.X
  a.c.jnt.VarEx.Y = sum(ajive.oracle$sJIVE$joint_matrices[[2]]$d^2)/TotVar.Y
  a.c.ind.VarEx.X = sum(ajive.oracle$sJIVE$indiv_matrices[[1]]$d^2)/TotVar.X
  a.c.ind.VarEx.Y = sum(ajive.oracle$sJIVE$indiv_matrices[[2]]$d^2)/TotVar.Y
  
  #### Apply ProJIVE-Oracle
  Y=do.call(cbind, blocks); P=c(p1, p2); Q=c(r.J,(true_signal_ranks-r.J))
  WJ.init = list(JntLd.X, JntLd.Y)
  WI.init = list(IndivLd.X, IndivLd.Y)
  init.loads = list(WJ.init, WI.init)
  pro.oracle.time = system.time({pro.oracle.jive.res.all = ProJIVE(Y=Y, P=P, Q=Q, plots = TRUE, sig_hat = "MLE", init.loads = init.loads, num.starts = 10,
                                                               center = TRUE, return.all.starts = FALSE)})
  print(paste("ProJIVE with loadings initiated from the truth done in", round(pro.oracle.time['elapsed'], 3), "sec."))
  
  pro.oracle.jive.res = pro.oracle.jive.res.all[[1]]
  
  #####Retrieve results from pro_Oracle with ranks pre-specified, and inital loadings from CJIVE
  pro.oracle.rJ = pro.oracle.jive.res$Ranks["Joint"]
  pro.oracle.r1 = pro.oracle.jive.res$Ranks["Indiv_1"]
  pro.oracle.r2 = pro.oracle.jive.res$Ranks["Indiv_2"]
  pro.oracle.JntScores.hat = pro.oracle.jive.res$SubjectScoreMatrix[,1:pro.oracle.rJ, drop = FALSE]
  pro.oracle.jnt.subnorm = chord.norm.diff(JntScores, pro.oracle.JntScores.hat)
  pro.oracle.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = pro.oracle.JntScores.hat)
  
  pro.oracle.load.mats = w_to_w_k(pro.oracle.jive.res$LoadingMatrix, P, Q)
  pro.oracle.jnt.loadX = pro.oracle.load.mats[[1]][,1:pro.oracle.rJ, drop = FALSE]
  pro.oracle.jnt.loadY = pro.oracle.load.mats[[2]][,1:pro.oracle.rJ, drop = FALSE]
  pro.oracle.jnt.loadX.norm = chord.norm.diff(JntLd.X, pro.oracle.jnt.loadX)
  pro.oracle.jnt.loadY.norm = chord.norm.diff(JntLd.Y, pro.oracle.jnt.loadY)
  pro.oracle.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = pro.oracle.jnt.loadX)
  pro.oracle.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = pro.oracle.jnt.loadY)
  
  pro.oracle.IX.load = pro.oracle.load.mats[[1]][,-(1:pro.oracle.rJ), drop = FALSE]
  pro.oracle.IY.load = pro.oracle.load.mats[[2]][,-(1:pro.oracle.rJ), drop = FALSE]
  pro.oracle.bX.hat = pro.oracle.jive.res$SubjectScoreMatrix[,pro.oracle.rJ+(1:pro.oracle.r1), drop = FALSE]
  pro.oracle.bY.hat = pro.oracle.jive.res$SubjectScoreMatrix[,(pro.oracle.rJ+pro.oracle.r1+(1:pro.oracle.r2)), drop = FALSE]
  
  pro.oracle.indiv.X.subnorm = chord.norm.diff(IndivScore.X, pro.oracle.bX.hat)
  pro.oracle.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, pro.oracle.bY.hat)
  pro.oracle.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, pro.oracle.IX.load)
  pro.oracle.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, pro.oracle.IY.load)
  
  pro.oracle.indiv.X.sub.pmse = pmse.2(standardize = TRUE, S1 =IndivScore.X, S2 = pro.oracle.bX.hat)
  pro.oracle.indiv.Y.sub.pmse = pmse.2(standardize = TRUE, S1 =IndivScore.Y, S2 = pro.oracle.bY.hat)
  pro.oracle.indiv.X.load.pmse = pmse.2(standardize = TRUE, S1 =IndivLd.X, S2 = pro.oracle.IX.load)
  pro.oracle.indiv.Y.load.pmse = pmse.2(standardize = TRUE, S1 =IndivLd.Y, S2 = pro.oracle.IY.load)
  
  pro.oracle.jnt.VarEx.X = pro.oracle.jive.res$VarianceExplained[[1]][1]
  pro.oracle.jnt.VarEx.Y = pro.oracle.jive.res$VarianceExplained[[2]][1]
  pro.oracle.ind.VarEx.X = pro.oracle.jive.res$VarianceExplained[[1]][2]
  pro.oracle.ind.VarEx.Y = pro.oracle.jive.res$VarianceExplained[[2]][2]
  
  #### Apply ProJIVE-Init Oracle
  pro.time = system.time({pro.jive.res = ProJIVE_EM(Y=Y, P=P, Q=Q, plots = TRUE, sig_hat = "MLE", init.loads = "CJIVE")})
  print(paste("ProJIVE (with AJIVE initiated Loadings) done in", round(pro.time['elapsed'], 3), "sec."))
  
  #####Retrieve results from pro_Oracle with ranks not specified: use permutation method to find them
  pro.rJ = pro.jive.res$Ranks["Joint"]
  pro.r1 = pro.jive.res$Ranks["Indiv_1"]
  pro.r2 = pro.jive.res$Ranks["Indiv_2"]
  pro.JntScores.hat = pro.jive.res$SubjectScoreMatrix[,1:pro.rJ, drop = FALSE]
  pro.jnt.subnorm = chord.norm.diff(JntScores, pro.JntScores.hat)
  pro.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = pro.JntScores.hat)
  
  pro.load.mats = w_to_w_k(pro.jive.res$LoadingMatrix, P, Q)
  pro.jnt.loadX = pro.load.mats[[1]][,1:pro.rJ, drop = FALSE]
  pro.jnt.loadY = pro.load.mats[[2]][,1:pro.rJ, drop = FALSE]
  pro.jnt.loadX.norm = chord.norm.diff(JntLd.X, pro.jnt.loadX)
  pro.jnt.loadY.norm = chord.norm.diff(JntLd.Y, pro.jnt.loadY)
  pro.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = pro.jnt.loadX)
  pro.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = pro.jnt.loadY)
  
  pro.IX.load = pro.load.mats[[1]][,-(1:pro.rJ), drop = FALSE]
  pro.IY.load = pro.load.mats[[2]][,-(1:pro.rJ), drop = FALSE]
  pro.bX.hat = pro.jive.res$SubjectScoreMatrix[,pro.rJ+(1:pro.r1), drop = FALSE]
  pro.bY.hat = pro.jive.res$SubjectScoreMatrix[,(pro.rJ+pro.r1+(1:pro.r2)), drop = FALSE]
  
  pro.indiv.X.subnorm = chord.norm.diff(IndivScore.X, pro.bX.hat)
  pro.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, pro.bY.hat)
  pro.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, pro.IX.load)
  pro.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, pro.IY.load)
  
  pro.indiv.X.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.X, S2 = pro.bX.hat)
  pro.indiv.Y.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.Y, S2 = pro.bY.hat)
  pro.indiv.X.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.X, S2 = pro.IX.load)
  pro.indiv.Y.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.Y, S2 = pro.IY.load)
  
  pro.jnt.VarEx.X = pro.jive.res$VarianceExplained[[1]][1]
  pro.jnt.VarEx.Y = pro.jive.res$VarianceExplained[[2]][1]
  pro.ind.VarEx.X = pro.jive.res$VarianceExplained[[1]][2]
  pro.ind.VarEx.Y = pro.jive.res$VarianceExplained[[2]][2]
  
  #################    Apply GIPCA      
  GIPCA.time = system.time({GIPCA.res = EPCAJIVEMissbio(blocks, r.J, c(r.I1,r.I2), D = P, family = rep("gaussian",2), tol = 15)})
  print(paste("GIPCA done in", round(GIPCA.time['elapsed'], 3), "sec."))
  
  GIPCA.joint.scores = GIPCA.res$U.joint[,-1, drop = FALSE]
  GI.jnt.subnorm = chord.norm.diff(JntScores, GIPCA.joint.scores)
  GI.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = GIPCA.joint.scores)
  
  GIPCA.joint.loads.X = t(GIPCA.res$V.joint[-1,1:P[1], drop = FALSE])
  GIPCA.joint.loads.Y = t(GIPCA.res$V.joint[-1,-(1:P[1]), drop = FALSE])
  GI.jnt.loadX.norm = chord.norm.diff(JntLd.X, GIPCA.joint.loads.X)
  GI.jnt.loadY.norm = chord.norm.diff(JntLd.Y, GIPCA.joint.loads.Y)
  GI.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = GIPCA.joint.loads.X)
  GI.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = GIPCA.joint.loads.Y)
  
  GIPCA.bX.hat = GIPCA.res$U.ind[[1]]
  GIPCA.bY.hat = GIPCA.res$U.ind[[2]]
  GIPCA.indiv.loads.X = t(GIPCA.res$V.ind[[1]])
  GIPCA.indiv.loads.Y = t(GIPCA.res$V.ind[[2]])

  GI.indiv.X.subnorm = chord.norm.diff(IndivScore.X, GIPCA.bX.hat)
  GI.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, GIPCA.bY.hat)
  GI.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, GIPCA.indiv.loads.X)
  GI.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, GIPCA.indiv.loads.Y)
  
  GI.indiv.X.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.X, S2 = GIPCA.bX.hat)
  GI.indiv.Y.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.Y, S2 = GIPCA.bY.hat)
  GI.indiv.X.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.X, S2 = GIPCA.indiv.loads.X)
  GI.indiv.Y.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.Y, S2 = GIPCA.indiv.loads.Y)

  # GIPCA.jnt.VarEx.X = protrue.jive.res$VarianceExplained$Joint[[1]][1]
  # GIPCA.jnt.VarEx.Y = protrue.jive.res$VarianceExplained$Joint[[2]][1]
  # GIPCA.ind.VarEx.X = protrue.jive.res$VarianceExplained$Indiv[[1]][2]
  # GIPCA.ind.VarEx.Y = protrue.jive.res$VarianceExplained$Indiv[[2]][2]
  
  #################    Apply dCCA      
  dCCA.time = system.time({dCCA.res = dCCA(blocks.t[[1]], blocks.t[[2]], r_1 = as.integer(true_signal_ranks[1]),  r_2 = as.integer(true_signal_ranks[2]), r_12 = as.integer(r.J))})
  print(paste("dCCA (with reticulate) done in", round(dCCA.time['elapsed'], 3), "sec."))
  
  dCCA.CX.hat = dCCA.res[[3]]; svd.CX.dcca = svd(dCCA.CX.hat, nu = r.J, nv = r.J)
  dCCA.CY.hat = dCCA.res[[4]]; svd.CY.dcca = svd(dCCA.CY.hat, nu = r.J, nv = r.J)
  dCCA.DX.hat = dCCA.res[[5]]; svd.DX.dcca = svd(dCCA.DX.hat, nu = r.I1, nv = r.I1)
  dCCA.DY.hat = dCCA.res[[6]]; svd.DY.dcca = svd(dCCA.DY.hat, nu = r.I2, nv = r.I2)
  
  # dCCA.joint.scores = svd(cbind(svd.CX.dcca$v, svd.CY.dcca$v), nv = r.J, nu = r.J)$u
  temp.diag = diag((1-sqrt(1 - dCCA.res[[10]])/sqrt(1 + dCCA.res[[10]])), nrow = r.J, ncol = r.J)
  dCCA.joint.scores = 0.5*(svd.CX.dcca$v+svd.CY.dcca$v)%*%temp.diag; rm(temp.diag)
  
  dCCA.jnt.subnorm = chord.norm.diff(JntScores, dCCA.joint.scores)
  dCCA.jnt.sub.pmse = pmse.2(standardize = TRUE, S1 = JntScores, S2 = dCCA.joint.scores)
  
  dCCA.jnt.loadX.norm = chord.norm.diff(JntLd.X, svd.CX.dcca[['u']])
  dCCA.jnt.loadY.norm = chord.norm.diff(JntLd.Y, svd.CY.dcca[['u']])
  dCCA.jnt.loadX.pmse = pmse.2(standardize = TRUE, S1 = JntLd.X, S2 = svd.CX.dcca[['u']])
  dCCA.jnt.loadY.pmse = pmse.2(standardize = TRUE, S1 = JntLd.Y, S2 = svd.CY.dcca[['u']])
  
  dCCA.indiv.X.subnorm = chord.norm.diff(IndivScore.X, svd.DX.dcca[['v']])
  dCCA.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, svd.DY.dcca[['v']])
  dCCA.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, svd.DX.dcca[['u']])
  dCCA.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, svd.DY.dcca[['u']])
  
  dCCA.indiv.X.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.X, S2 = svd.DX.dcca[['v']])
  dCCA.indiv.Y.sub.pmse = pmse.2(standardize = TRUE, S1 = IndivScore.Y, S2 = svd.DY.dcca[['v']])
  dCCA.indiv.X.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.X, S2 = svd.DX.dcca[['u']])
  dCCA.indiv.Y.load.pmse = pmse.2(standardize = TRUE, S1 = IndivLd.Y, S2 = svd.DY.dcca[['u']])

  dCCA.TotVE.X = MatVar(dCCA.res[[1]])/TotVar.X
  dCCA.TotVE.Y = MatVar(dCCA.res[[2]])/TotVar.Y
  dCCA.jnt.VarEx.X = MatVar(dCCA.CX.hat)/TotVar.X
  dCCA.jnt.VarEx.Y = MatVar(dCCA.CY.hat)/TotVar.Y
  dCCA.ind.VarEx.X = MatVar(dCCA.DX.hat)/TotVar.X
  dCCA.ind.VarEx.Y = MatVar(dCCA.DY.hat)/TotVar.Y
  
  
  CompEvals = c(pro.oracle.jnt.subnorm, pro.oracle.jnt.loadX.norm, pro.oracle.jnt.loadY.norm, pro.oracle.indiv.X.subnorm, pro.oracle.indiv.Y.subnorm, pro.oracle.indiv.X.loadnorm, 
                pro.oracle.indiv.Y.loadnorm,
                
                pro.jnt.subnorm, pro.jnt.loadX.norm, pro.jnt.loadY.norm, pro.indiv.X.subnorm, pro.indiv.Y.subnorm, pro.indiv.X.loadnorm, 
                pro.indiv.Y.loadnorm,
                
                a.c.jnt.subnorm, a.c.jnt.loadX.norm, a.c.jnt.loadY.norm, a.c.indiv.X.subnorm, a.c.indiv.Y.subnorm, a.c.indiv.X.loadnorm, 
                a.c.indiv.Y.loadnorm,
                
                i.jnt.subnorm, i.jnt.loadX.norm, i.jnt.loadY.norm, i.indiv.X.subnorm, i.indiv.Y.subnorm, i.indiv.X.loadnorm, 
                i.indiv.Y.loadnorm,
                
                GI.jnt.subnorm, GI.jnt.loadX.norm, GI.jnt.loadY.norm, GI.indiv.X.subnorm, GI.indiv.Y.subnorm, GI.indiv.X.loadnorm, 
                GI.indiv.Y.loadnorm,
                
                dCCA.jnt.subnorm, dCCA.jnt.loadX.norm, dCCA.jnt.loadY.norm, dCCA.indiv.X.subnorm, dCCA.indiv.Y.subnorm, dCCA.indiv.X.loadnorm, 
                dCCA.indiv.Y.loadnorm, 
                
                pro.oracle.jnt.sub.pmse, pro.oracle.jnt.loadX.pmse, pro.oracle.jnt.loadY.pmse, pro.oracle.indiv.X.sub.pmse, pro.oracle.indiv.Y.sub.pmse, pro.oracle.indiv.X.load.pmse, 
                pro.oracle.indiv.Y.load.pmse,
                 
                pro.jnt.sub.pmse, pro.jnt.loadX.pmse, pro.jnt.loadY.pmse, pro.indiv.X.sub.pmse, pro.indiv.Y.sub.pmse, pro.indiv.X.load.pmse, 
                pro.indiv.Y.load.pmse,
                 
                a.c.jnt.sub.pmse, a.c.jnt.loadX.pmse, a.c.jnt.loadY.pmse, a.c.indiv.X.sub.pmse, a.c.indiv.Y.sub.pmse, a.c.indiv.X.load.pmse, 
                a.c.indiv.Y.load.pmse,
                 
                i.jnt.sub.pmse, i.jnt.loadX.pmse, i.jnt.loadY.pmse, i.indiv.X.sub.pmse, i.indiv.Y.sub.pmse, i.indiv.X.load.pmse, 
                i.indiv.Y.load.pmse,
                 
                GI.jnt.sub.pmse, GI.jnt.loadX.pmse, GI.jnt.loadY.pmse, GI.indiv.X.sub.pmse, GI.indiv.Y.sub.pmse, GI.indiv.X.load.pmse, 
                GI.indiv.Y.load.pmse,
                 
                dCCA.jnt.sub.pmse, dCCA.jnt.loadX.pmse, dCCA.jnt.loadY.pmse, dCCA.indiv.X.sub.pmse, dCCA.indiv.Y.sub.pmse, dCCA.indiv.X.load.pmse, 
                dCCA.indiv.Y.load.pmse)

  names(CompEvals) = c("ProJIVE - Oracle Joint Subj Scores- Norm", "ProJIVE - Oracle Joint Loads X- Norm", "ProJIVE - Oracle Joint Loads Y- Norm", "ProJIVE - Oracle Indiv Subj Scores X- Norm",
                       "ProJIVE - Oracle Indiv Subj Scores Y- Norm", "ProJIVE - Oracle Indiv Loads X- Norm", "ProJIVE - Oracle Indiv Loads Y- Norm",
                       
                       "ProJIVE Joint Subj Scores- Norm", "ProJIVE Joint Loads X- Norm", "ProJIVE Joint Loads Y- Norm", "ProJIVE Indiv Subj Scores X- Norm",
                       "ProJIVE Indiv Subj Scores Y- Norm", "ProJIVE Indiv Loads X- Norm", "ProJIVE Indiv Loads Y- Norm",
                       
                       "AJIVE Joint Subj Scores- Norm", "AJIVE Joint Loads X- Norm", "AJIVE Joint Loads Y- Norm", "AJIVE Indiv Subj Scores X- Norm",
                       "AJIVE Indiv Subj Scores Y- Norm", "AJIVE Indiv Loads X- Norm", "AJIVE Indiv Loads Y- Norm",
                       
                       "R.JIVE Joint Subj Scores- Norm", "R.JIVE Joint Loads X- Norm", "R.JIVE Joint Loads Y- Norm", "R.JIVE Indiv Subj Scores X- Norm",
                       "R.JIVE Indiv Subj Scores Y- Norm", "R.JIVE Indiv Loads X- Norm", "R.JIVE Indiv Loads Y- Norm",
                       
                       "GIPCA Joint Subj Scores- Norm", "GIPCA Joint Loads X- Norm", "GIPCA Joint Loads Y- Norm", "GIPCA Indiv Subj Scores X- Norm",
                       "GIPCA Indiv Subj Scores Y- Norm", "GIPCA Indiv Loads X- Norm", "GIPCA Indiv Loads Y- Norm",
                       
                       "dCCA Joint Subj Scores- Norm", "dCCA Joint Loads X- Norm", "dCCA Joint Loads Y- Norm", "dCCA Indiv Subj Scores X- Norm",
                       "dCCA Indiv Subj Scores Y- Norm", "dCCA Indiv Loads X- Norm", "dCCA Indiv Loads Y- Norm", 
                       
                       "ProJIVE - Oracle Joint Subj Scores- PMSE", "ProJIVE - Oracle Joint Loads X- PMSE", "ProJIVE - Oracle Joint Loads Y- PMSE", "ProJIVE - Oracle Indiv Subj Scores X- PMSE",
                       "ProJIVE - Oracle Indiv Subj Scores Y- PMSE", "ProJIVE - Oracle Indiv Loads X- PMSE", "ProJIVE - Oracle Indiv Loads Y- PMSE",
                       
                       "ProJIVE Joint Subj Scores- PMSE", "ProJIVE Joint Loads X- PMSE", "ProJIVE Joint Loads Y- PMSE", "ProJIVE Indiv Subj Scores X- PMSE",
                       "ProJIVE Indiv Subj Scores Y- PMSE", "ProJIVE Indiv Loads X- PMSE", "ProJIVE Indiv Loads Y- PMSE",
                       
                       "AJIVE Joint Subj Scores- PMSE", "AJIVE Joint Loads X- PMSE", "AJIVE Joint Loads Y- PMSE", "AJIVE Indiv Subj Scores X- PMSE",
                       "AJIVE Indiv Subj Scores Y- PMSE", "AJIVE Indiv Loads X- PMSE", "AJIVE Indiv Loads Y- PMSE",
                       
                       "R.JIVE Joint Subj Scores- PMSE", "R.JIVE Joint Loads X- PMSE", "R.JIVE Joint Loads Y- PMSE", "R.JIVE Indiv Subj Scores X- PMSE",
                       "R.JIVE Indiv Subj Scores Y- PMSE", "R.JIVE Indiv Loads X- PMSE", "R.JIVE Indiv Loads Y- PMSE",
                       
                       "GIPCA Joint Subj Scores- PMSE", "GIPCA Joint Loads X- PMSE", "GIPCA Joint Loads Y- PMSE", "GIPCA Indiv Subj Scores X- PMSE",
                       "GIPCA Indiv Subj Scores Y- PMSE", "GIPCA Indiv Loads X- PMSE", "GIPCA Indiv Loads Y- PMSE",
                       
                       "dCCA Joint Subj Scores- PMSE", "dCCA Joint Loads X- PMSE", "dCCA Joint Loads Y- PMSE", "dCCA Indiv Subj Scores X- PMSE",
                       "dCCA Indiv Subj Scores Y- PMSE", "dCCA Indiv Loads X- PMSE", "dCCA Indiv Loads Y- PMSE")
  
  VarEx.true = c(JVE.X, JVE.Y, IVE.X, IVE.Y, TotVE.X, TotVE.Y)
  VarEx.ajive = c(a.c.jnt.VarEx.X, a.c.jnt.VarEx.Y, a.c.ind.VarEx.X, a.c.ind.VarEx.Y)
  VarEx.rjive = c(i.jnt.VarEx.X, i.jnt.VarEx.Y, i.ind.VarEx.X, i.ind.VarEx.Y)
  VarEx.projive = c(pro.jnt.VarEx.X, pro.jnt.VarEx.Y, pro.ind.VarEx.X, pro.ind.VarEx.Y)
  VarEx.pro.oracle = c(pro.oracle.jnt.VarEx.X, pro.oracle.jnt.VarEx.Y, pro.oracle.ind.VarEx.X, pro.oracle.ind.VarEx.Y)
  VarEx.dCCA = c(dCCA.jnt.VarEx.X, dCCA.jnt.VarEx.Y, dCCA.ind.VarEx.X, dCCA.ind.VarEx.Y)
  
  names(VarEx.true) = c("JointVarExp_X", "JointVarExp_Y", "IndivVarExp_X", "IndivVarExp_Y", 
                        "TotalVarExp_X",  "TotalVarExp_Y")
  
  VarEx = c(VarEx.true, VarEx.ajive, VarEx.rjive, VarEx.projive, VarEx.pro.oracle, VarEx.dCCA)
  names(VarEx) = c(paste0("True_", names(VarEx.true)), 
                   paste0("AJIVE_", names(VarEx.true)[1:4]),
                   paste0("RJIVE_", names(VarEx.true)[1:4]),
                   paste0("ProJIVE_", names(VarEx.true)[1:4]),
                   paste0("ProJIVE_Oracle_", names(VarEx.true)[1:4]),
                   paste0("dCCA", names(VarEx.true)[1:4]))
  
  total_time = Sys.time() - TIME
  
  out = c(VarEx, CompEvals, total_time['elapsed'], i.time['elapsed'], a.c.time['elapsed'], pro.time['elapsed'], pro.oracle.time['elapsed'], 
          GIPCA.time['elapsed'], dCCA.time['elapsed'])
  
  names(out) = c(names(VarEx), names(CompEvals),"Total_Time", "R.JIVE_Time", 
                 "AJIVE_Time", "ProJIVE_Time", "OracleProJIVE_Time", "GIPCA_Time", "dCCA_Time")
  
  fname = file.path(outdir, paste("sim_", JntVarEx1, "_", JntVarEx2, rep_number, ".csv", sep=""))

  write.csv(file=fname, out, row.names = T)
  print(paste("Output saved. Total time:", round(total_time, 3), "minutes"))
}

