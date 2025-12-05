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
  
  library(r.jive) 
  prog.dir = "/home/rmurden/PJIVE/Programs"
  prog.gipca.dir = "/home/rmurden/PJIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
  files= list.files(file.path(prog.dir, "R"))
  for (i in files) source(file.path(prog.dir, 'R', i))
  source(file.path(prog.dir, "Functions_for_PJIVE.R"))
  source(file.path(prog.dir, "Functions_for_CJIVE.R"))
  gipca.files = list.files(prog.gipca.dir,full.names = TRUE)
  for(i in gipca.files) source(i)
  
  TIME = Sys.time()
  
  #######################################################################
  print(paste0("Is rep_number=",rep_number," an double? ",is.double(rep_number)))
  set.seed(round(rep_number)) ##To ensure that any randomized processes give the same results each time
  ###Construct Datasets
  #####Add noise to already generated Joint and individual signal matrices for X, Y data blocks:
  
  #######JIVE Implementations for Toy Data No 1################################################################################
  ###Setup input parameters for JIVE implementations
    true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true jranks of overall signals
    ToyDat = GenerateToyData(n = n, p1 = p1, p2 = p2, JntVarEx1 = JntVarEx1, JntVarEx2 = JntVarEx2, 
                             IndVarEx1 = IndVarEx1, IndVarEx2 =  IndVarEx2, jnt_rank = r.J,
                             equal.eig = F,ind_rank1 = r.I1, ind_rank2 = r.I2, JntVarAdj = T, SVD.plots = F,
                             Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
    
    blocks <- lapply(ToyDat[[2]], scale)
    rnd.smp = sample(n, n/2)
    blocks.sub1 = lapply(blocks, function(x){x[rnd.smp,]})
    blocks.sub2 = lapply(blocks, function(x){x[-rnd.smp,]})
    blocks.t <- lapply(blocks,t)   ##iJIVE requires datasets organized as p-by-n matrices
    JntScores = ToyDat[['Scores']][['Joint']]
    IndivScore.X = ToyDat[['Scores']][["Indiv_1"]]
    IndivScore.Y = ToyDat[['Scores']][["Indiv_2"]]
    
    JntLd.X = t(ToyDat[['Loadings']][["Joint_1"]])
    JntLd.Y = t(ToyDat[['Loadings']][["Joint_2"]])
    IndivLd.X =t(ToyDat[['Loadings']][["Indiv_1"]])
    IndivLd.Y = t(ToyDat[['Loadings']][["Indiv_2"]])
    
    JX = ToyDat[[1]]$J1
    JY = ToyDat[[1]]$J2
    IX = ToyDat[[1]]$I1
    IY = ToyDat[[1]]$I2
    EX = ToyDat[[1]]$E1
    EY = ToyDat[[1]]$E2
    
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
  } else {
    i.JX.hat = jive.unknown.perm$joint[[1]]
    i.JY.hat = jive.unknown.perm$joint[[2]]
    i.jnt.loadX = svd(i.JX.hat, nu = i.rJ, nv = i.rJ)[['u']]
    i.jnt.loadY = svd(i.JY.hat, nu = i.rJ, nv = i.rJ)[['u']]
    
    i.jnt.subnorm = chord.norm.diff(JntScores, i.JntScores.hat)
    i.jnt.loadX.norm = chord.norm.diff(JntLd.X, i.jnt.loadX)
    i.jnt.loadY.norm = chord.norm.diff(JntLd.Y, i.jnt.loadY)
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
  
  i.TotVar.X = sum(svd(jive.unknown.perm$data[[1]])[['d']]^2)
  i.TotVar.Y = sum(svd(jive.unknown.perm$data[[2]])[['d']]^2)
  
  i.jnt.VarEx.X = sum(svd(i.JX.hat)[['d']]^2)/i.TotVar.X
  i.jnt.VarEx.Y = sum(svd(i.JY.hat)[['d']]^2)/i.TotVar.Y
  i.ind.VarEx.X = sum(svd(i.IX.hat)[['d']]^2)/i.TotVar.X
  i.ind.VarEx.Y = sum(svd(i.IY.hat)[['d']]^2)/i.TotVar.Y
  
  #### Apply AJIVE-Oracle
  a.c.time = system.time({ajive.oracle<-ajive(blocks,initial_signal_ranks = true_signal_ranks, joint_rank = r.J)})
  print(paste("AJIVE done in", round(a.c.time['elapsed'], 3), "sec."))
  
  #####Retrieve results from AJIVE_Oracle with ranks not specified: use permutation method to find them
  a.c.JntScores.hat = ajive.oracle$joint_scores
  a.c.jnt.loadX = ajive.oracle$block_decomps[[1]]$joint$v
  a.c.jnt.loadY = ajive.oracle$block_decomps[[2]]$joint$v
  
  a.c.jnt.subnorm = chord.norm.diff(JntScores, a.c.JntScores.hat)
  a.c.jnt.loadX.norm = chord.norm.diff(JntLd.X, a.c.jnt.loadX)
  a.c.jnt.loadY.norm = chord.norm.diff(JntLd.Y, a.c.jnt.loadY)
  
  a.c.W.IY.hat = ajive.oracle$block_decomps[[2]][["individual"]][["v"]]
  a.c.W.IX.hat = ajive.oracle$block_decomps[[1]][["individual"]][["v"]]
  a.c.bX.hat = ajive.oracle$block_decomps[[1]][["individual"]][["u"]]
  a.c.bY.hat = ajive.oracle$block_decomps[[2]][["individual"]][["u"]]
  
  a.c.indiv.X.subnorm = chord.norm.diff(IndivScore.X, a.c.bX.hat)
  a.c.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, a.c.bY.hat)
  a.c.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, a.c.W.IX.hat)
  a.c.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, a.c.W.IY.hat)
  
  a.c.jnt.VarEx.X = sum(svd(ajive.oracle[["block_decomps"]][[1]]$joint[["full"]])[['d']]^2)/TotVar.X
  a.c.jnt.VarEx.Y = sum(svd(ajive.oracle[["block_decomps"]][[2]]$joint[["full"]])[['d']]^2)/TotVar.Y
  a.c.ind.VarEx.X = sum(svd(ajive.oracle[["block_decomps"]][[1]]$individual[["full"]])[['d']]^2)/TotVar.X
  a.c.ind.VarEx.Y = sum(svd(ajive.oracle[["block_decomps"]][[2]]$individual[["full"]])[['d']]^2)/TotVar.Y
  
  #### Apply ProJIVE-Oracle
  Y=do.call(cbind, blocks); P=c(p1, p2); Q=c(r.J,(true_signal_ranks-r.J))
  pro.true.time = system.time({protrue.jive.res = ProJIVE_EM(Y=Y, P=P, Q=Q, plots = TRUE, sig_hat = "MLE", init.loads = "CJIVE")})
  print(paste("ProJIVE done in", round(pro.true.time['elapsed'], 3), "sec."))
  
  #####Retrieve results from pro_Oracle with ranks not specified: use permutation method to find them
  pro.rJ = protrue.jive.res$Ranks["Joint"]
  pro.r1 = protrue.jive.res$Ranks["Indiv_1"]
  pro.r2 = protrue.jive.res$Ranks["Indiv_2"]
  pro.JntScores.hat = protrue.jive.res$SubjectScoreMatrix[,1:pro.rJ]
  pro.jnt.subnorm = chord.norm.diff(JntScores, pro.JntScores.hat)
  
  pro.load.mats = w_to_w_k(protrue.jive.res$LoadingMatrix, P, Q)
  pro.jnt.loadX = pro.load.mats[[1]][,1:pro.rJ]
  pro.jnt.loadY = pro.load.mats[[2]][,1:pro.rJ]
  pro.jnt.loadX.norm = chord.norm.diff(JntLd.X, pro.jnt.loadX)
  pro.jnt.loadY.norm = chord.norm.diff(JntLd.Y, pro.jnt.loadY)
  
  pro.IX.load = pro.load.mats[[1]][,-(1:pro.rJ)]
  pro.IY.load = pro.load.mats[[2]][,-(1:pro.rJ)]
  pro.bX.hat = protrue.jive.res$SubjectScoreMatrix[,pro.rJ+(1:pro.r1)]
  pro.bY.hat = protrue.jive.res$SubjectScoreMatrix[,(pro.rJ+pro.r1+(1:pro.r2))]
  
  pro.indiv.X.subnorm = chord.norm.diff(IndivScore.X, pro.bX.hat)
  pro.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, pro.bY.hat)
  pro.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, pro.IX.load)
  pro.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, pro.IY.load)
  
  pro.jnt.VarEx.X = protrue.jive.res$VarianceExplained$Joint[1]
  pro.jnt.VarEx.Y = protrue.jive.res$VarianceExplained$Joint[2]
  pro.ind.VarEx.X = protrue.jive.res$VarianceExplained$Indiv[1]
  pro.ind.VarEx.Y = protrue.jive.res$VarianceExplained$Indiv[2]
  
  #################    Apply GIPCA      
  GIPCA.time = system.time({GIPCA.res = EPCAJIVEMissbio(blocks, r.J, c(r.I1,r.I2), D = P, family = rep("gaussian",2), tol = 1E-5)})
  print(paste("GIPCA done in", round(GIPCA.time['elapsed'], 3), "sec."))
  
  GIPCA.joint.scores = GIPCA.res$U.joint[,-1, drop = FALSE]
  GI.jnt.subnorm = chord.norm.diff(JntScores, GIPCA.joint.scores)
  
  GIPCA.joint.loads.X = t(GIPCA.res$V.joint[-1,1:P[1], drop = FALSE])
  GIPCA.joint.loads.Y = t(GIPCA.res$V.joint[-1,-(1:P[1]), drop = FALSE])
  GI.jnt.loadX.norm = chord.norm.diff(JntLd.X, GIPCA.joint.loads.X)
  GI.jnt.loadY.norm = chord.norm.diff(JntLd.Y, GIPCA.joint.loads.Y)
  
  GIPCA.bX.hat = GIPCA.res$U.ind[[1]]
  GIPCA.bY.hat = GIPCA.res$U.ind[[2]]
  GIPCA.indiv.loads.X = t(GIPCA.res$V.ind[[1]])
  GIPCA.indiv.loads.Y = t(GIPCA.res$V.ind[[2]])

  GI.indiv.X.subnorm = chord.norm.diff(IndivScore.X, GIPCA.bX.hat)
  GI.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, GIPCA.bY.hat)
  GI.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, GIPCA.indiv.loads.X)
  GI.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, GIPCA.indiv.loads.Y)
  
  # GI.jnt.VarEx.X = sum(svd(ajive.oracle[["block_decomps"]][[1]]$joint[["full"]])[['d']]^2)/TotVar.X
  # GI.jnt.VarEx.Y = sum(svd(ajive.oracle[["block_decomps"]][[2]]$joint[["full"]])[['d']]^2)/TotVar.Y
  # GI.ind.VarEx.X = sum(svd(ajive.oracle[["block_decomps"]][[1]]$individual[["full"]])[['d']]^2)/TotVar.X
  # GI.ind.VarEx.Y = sum(svd(ajive.oracle[["block_decomps"]][[2]]$individual[["full"]])[['d']]^2)/TotVar.Y
  
  # #### Apply CC.Oracle JIVE
# cc.true.time = system.time({cctrue.jive.res = cc.jive(blocks, signal.ranks = true_signal_ranks, center = F)})
# print(paste("CJIVE done in", round(cc.true.time['elapsed'], 3), "sec."))
# 
# #####Retrieve results from CC_Oracle with ranks not specified: use permutation method to find them
# cc.JntScores.hat = cctrue.jive.res$CanCorRes$Jnt_Scores
# cc.jnt.subnorm = chord.norm.diff(JntScores, cc.JntScores.hat)
# 
# cc.W.JX = cctrue.jive.res$sJIVE$joint_matrices[[1]]$v
# cc.W.JY = cctrue.jive.res$sJIVE$joint_matrices[[2]]$v
# cc.jnt.loadX.norm = chord.norm.diff(JntLd.X, cc.W.JX)
# cc.jnt.loadY.norm = chord.norm.diff(JntLd.Y, cc.W.JY)
# 
# cc.bX.hat = cctrue.jive.res$sJIVE$indiv_matrices[[1]]$u
# cc.bY.hat = cctrue.jive.res$sJIVE$indiv_matrices[[2]]$u
# cc.W.IX.hat = cctrue.jive.res$sJIVE$indiv_matrices[[1]]$v
# cc.W.IY.hat = cctrue.jive.res$sJIVE$indiv_matrices[[2]]$v
# 
# cc.indiv.X.subnorm = chord.norm.diff(IndivScore.X, cc.bX.hat)
# cc.indiv.Y.subnorm = chord.norm.diff(IndivScore.Y, cc.bY.hat)
# cc.indiv.X.loadnorm = chord.norm.diff(IndivLd.X, cc.W.IX.hat)
# cc.indiv.Y.loadnorm = chord.norm.diff(IndivLd.Y, cc.W.IY.hat)
# 
# ######    Predictions     #####
# #AJIVE 
# AJIVE.res.sub1 = ajive(blocks.sub1, true_signal_ranks, joint_rank = r.J)
# a.jnt_score.sub1 = AJIVE.res.sub1$joint_scores
# JL.1_sub1 = t(AJIVE.res.sub1$block_decomps[[1]]$joint$full)%*%a.jnt_score.sub1
# JL.2_sub1 = t(AJIVE.res.sub1$block_decomps[[2]]$joint$full)%*%a.jnt_score.sub1
# JL.1_sub1.ginv = ginv(JL.1_sub1)
# JL.2_sub1.ginv = ginv(JL.2_sub1)
# a.pred.aJIVE.sub2 = (blocks.sub2[[1]]%*%t(JL.1_sub1.ginv) + blocks.sub2[[2]]%*%t(JL.2_sub1.ginv))/2
# 
# a.c.pred.cor = diag(cor(a.pred.aJIVE.sub2, JntScores[-rnd.smp,]))
# print(paste("AJIVE Prediction done"))
# 
# #ProJIVE
# ProJIVE.res.sub1 = ProJIVE_EM(do.call(cbind, blocks.sub1), P, Q, plots = FALSE)
# Load.Mats = w_to_w_k(ProJIVE.res.sub1$LoadingMatrix, P, Q)
# LoadMat.1_sub1 = Load.Mats[[1]]
# LoadMat.2_sub1 = Load.Mats[[2]]
# JL.1_sub1.ginv = ginv(matrix(LoadMat.1_sub1[,1:r.J], ncol=r.J))
# JL.2_sub1.ginv = ginv(matrix(LoadMat.2_sub1[,1:r.J], ncol=r.J))
# naive.pred.ProJIVE.sub2 = (blocks.sub2[[1]]%*%t(JL.1_sub1.ginv) + blocks.sub2[[2]]%*%t(JL.2_sub1.ginv))/2
# 
# pro.pred.cor = diag(cor(naive.pred.ProJIVE.sub2, JntScores[-rnd.smp,]))
# print(paste("ProJIVE Prediction done"))

#  PredEvals = c(a.c.pred.cor, pro.pred.cor)
#  names(PredEvals) = c(paste("AJIVE_Pred_Corr", 1:r.J, sep = ""), paste("ProJIVE_Pred_Corr", 1:r.J, sep = ""))
  
#  print(paste("Corr. of ProJIVE v AJIVE Scores:" ,
#              round(diag(cor(ProJIVE.res.sub1$SubjectScoreMatrix[,1:r.J], AJIVE.res.sub1$joint_scores)),3)))
  
  CompEvals = c(pro.jnt.subnorm, pro.jnt.loadX.norm, pro.jnt.loadY.norm, pro.indiv.X.subnorm, pro.indiv.Y.subnorm, pro.indiv.X.loadnorm, 
                pro.indiv.Y.loadnorm,
                
                a.c.jnt.subnorm, a.c.jnt.loadX.norm, a.c.jnt.loadY.norm, a.c.indiv.X.subnorm, a.c.indiv.Y.subnorm, a.c.indiv.X.loadnorm, 
                a.c.indiv.Y.loadnorm,
                
                i.jnt.subnorm, i.jnt.loadX.norm, i.jnt.loadY.norm, i.indiv.X.subnorm, i.indiv.Y.subnorm, i.indiv.X.loadnorm, 
                i.indiv.Y.loadnorm,
                
                GI.jnt.subnorm, GI.jnt.loadX.norm, GI.jnt.loadY.norm, GI.indiv.X.subnorm, GI.indiv.Y.subnorm, GI.indiv.X.loadnorm, 
                GI.indiv.Y.loadnorm)
  
  names(CompEvals) = c("ProJIVE Joint Subj Scores", "ProJIVE Joint Loads X", "ProJIVE Joint Loads Y", "ProJIVE Indiv Subj Scores X",
                       "ProJIVE Indiv Subj Scores Y", "ProJIVE Indiv Loads X", "ProJIVE Indiv Loads Y",
                       
                       "AJIVE-Oracle Joint Subj Scores", "AJIVE-Oracle Joint Loads X", "AJIVE-Oracle Joint Loads Y", "AJIVE-Oracle Indiv Subj Scores X",
                       "AJIVE-Oracle Indiv Subj Scores Y", "AJIVE-Oracle Indiv Loads X", "AJIVE-Oracle Indiv Loads Y",
                       
                       "R.JIVE Joint Subj Scores", "R.JIVE Joint Loads X", "R.JIVE Joint Loads Y", "R.JIVE Indiv Subj Scores X",
                       "R.JIVE Indiv Subj Scores Y", "R.JIVE Indiv Loads X", "R.JIVE Indiv Loads Y",
                       
                       "GIPCA Joint Subj Scores", "GIPCA Joint Loads X", "GIPCA Joint Loads Y", "GIPCA Indiv Subj Scores X",
                       "GIPCA Indiv Subj Scores Y", "GIPCA Indiv Loads X", "GIPCA Indiv Loads Y")
  
  VarEx.true = c(JVE.X, JVE.Y, IVE.X, IVE.Y, TotVE.X, TotVE.Y)
  VarEx.ajive = c(a.c.jnt.VarEx.X, a.c.jnt.VarEx.Y, a.c.ind.VarEx.X, a.c.ind.VarEx.Y)
  VarEx.rjive = c(i.jnt.VarEx.X, i.jnt.VarEx.Y, i.ind.VarEx.X, i.ind.VarEx.Y)
  VarEx.projive = c(pro.jnt.VarEx.X, pro.jnt.VarEx.Y, pro.ind.VarEx.X, pro.ind.VarEx.Y)
  
  names(VarEx.true) = c("JointVarExp_X", "JointVarExp_Y", "IndivVarExp_X", "IndivVarExp_Y", 
                        "TotalVarExp_X",  "TotalVarExp_Y")
  
  VarEx = c(VarEx.true, VarEx.ajive, VarEx.rjive, VarEx.projive)
  names(VarEx) = c(paste0("True_", names(VarEx.true)), 
                   paste0("AJIVE_", names(VarEx.true)[1:4]),
                   paste0("RJIVE_", names(VarEx.true)[1:4]),
                   paste0("ProJIVE_", names(VarEx.true)[1:4]))
  total_time = Sys.time() - TIME
  
  # out = c(VarEx, CompEvals, PredEvals, total_time['elapsed'],i.time['elapsed'],a.c.time['elapsed'], pro.true.time['elapsed'])
  out = c(VarEx, CompEvals, total_time['elapsed'],i.time['elapsed'],a.c.time['elapsed'], pro.true.time['elapsed'], GIPCA.time['elapsed'])
  # names(out) = c(names(VarEx), names(CompEvals), names(PredEvals), "Total_Time", "R.JIVE_Time", 
  #                "CC_Oracle_Time", "ProJIVE_Time", "GIPCA_Time")
  names(out) = c(names(VarEx), names(CompEvals),"Total_Time", "R.JIVE_Time", 
                 "AJIVE_Time", "ProJIVE_Time", "GIPCA_Time")
  
  fname = file.path(outdir, paste("sim_", JntVarEx1, "_", JntVarEx2, rep_number, ".csv", sep=""))
  
  write.csv(file=fname, out, row.names = T)
  print(paste("Output saved. Total time:", total_time['elapsed']))
}

