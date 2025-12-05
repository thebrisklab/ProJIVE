#############################################################################################################################
### Generates a pair of toy data blocks in the manner outlined within the CJIVE manuscript  #################################
###         After data have been generated, apply CJIVE, RJIVE, GIPCA, and ProJIVE          #################################
### Author: Raphiel J. Murden                                                               #################################
### Supervised by Benjamin Risk                                                             #################################
#############################################################################################################################
args = commandArgs(trailingOnly=TRUE)
rep_number = as.numeric(args[1])
outdir = args[2]
n = as.numeric(args[3])
p1 = as.numeric(args[4])
p2 = as.numeric(args[5])
JntVarEx1 = as.numeric(args[6])
JntVarEx2 = as.numeric(args[7])
files = list.files(outdir, pattern = ".csv")
files.loads = list.files(outdir, pattern = ".RData")
IndVarEx1 = 0.25
IndVarEx2 = 0.25
r.J = as.numeric(args[8])
r.I1 = as.numeric(args[9])
r.I2 = as.numeric(args[9])

nm = paste("ConsistencyCheck_n", n, "_", rep_number, ".csv", sep="")
nm.loads = file.path(outdir, paste0("SimLoadings_rJ",r.J,"_rI",r.I1, ".RData"))
if (!(nm %in% files)){
  prog.dir = "/home/rmurden/PJIVE/Programs"
  prog.gipca.dir = "/home/rmurden/PJIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
  prog.dcca.dir = "/projects/guo_lab/cbis/users/raphiel/PJIVE/Programs/D-CCA"
  files= list.files(file.path(prog.dir, "R"))
  for (i in files) source(file.path(prog.dir, 'R', i))
  source(file.path(prog.dir, "Functions_for_PJIVE.R"))
  library(CJIVE); library(singR); library(lubridate)
  source(file.path(prog.dir, "Functions_for_CJIVE.R"))
  
  #######################################################################
  set.seed(rep_number) ##To ensure that any randomized processes give the same results each time
  
  ###Construct Datasets
  true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true ranks of overall signals
  if (!(nm.loads %in% files.loads)){
    
    ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), equal.eig = TRUE,
                             IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J, ind_ranks = c(r.I1, r.I2),
                             JntVarAdj = T, SVD.plots = F, Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
    # WJ1 and WJ2
    JntLd.X = t(ToyDat[['Loadings']][["Joint"]][[1]])
    JntLd.Y = t(ToyDat[['Loadings']][["Joint"]][[2]])
    WJ.init = list(JntLd.X, JntLd.Y)
    
    # WI1 and WI2
    IndivLd.X = t(ToyDat[['Loadings']][["Indiv"]][[1]])
    IndivLd.Y = t(ToyDat[['Loadings']][["Indiv"]][[2]])
    WI.init = list(IndivLd.X, IndivLd.Y)
    
    init.loads = list(WJ.init, WI.init)
    
    save(init.loads, file = nm.loads)
    
  } else {
    load(nm.loads)
    ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), 
                             IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J, ind_ranks = c(r.I1, r.I2),
                             JntVarAdj = T, SVD.plots = F, Error = T, print.cor = F, Loads = init.loads, Scores = "Gaussian")
    # WJ1 and WJ2
    JntLd.X = t(ToyDat[['Loadings']][["Joint"]][[1]])
    JntLd.Y = t(ToyDat[['Loadings']][["Joint"]][[2]])
    WJ.init = list(JntLd.X, JntLd.Y)
    
    # WI1 and WI2
    IndivLd.X = t(ToyDat[['Loadings']][["Indiv"]][[1]])
    IndivLd.Y = t(ToyDat[['Loadings']][["Indiv"]][[2]])
    WI.init = list(IndivLd.X, IndivLd.Y)
  }
  
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
  
  #############################################
  ############### ProJIVE  ####################
  PJIVE.res = ProJIVE_EM(Y=Y, P=P, Q=Q, init.loads = init.loads, sig_hat = "MLE",plots = FALSE)
  
  sig.hat.X = PJIVE.res$ErrorVariances[1]
  sig.hat.Y = PJIVE.res$ErrorVariances[2]
  
  PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
  PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(r.J+(1:r.I1))]
  
  sim.loads.X = PJIVE.loads.X
  sim.loads.Y = PJIVE.loads.Y
  
  mse.sim.loads.X = (cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)^2
  mse.sim.loads.Y = (cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)^2
  
  pmse.sim.loads.X = pmse(S1 = cbind(JntLd.X, IndivLd.X), S2 = PJIVE.loads.X, standardize = TRUE)
  pmse.sim.loads.Y = pmse(S1 = cbind(JntLd.Y, IndivLd.Y), S2 = PJIVE.loads.Y, standardize = TRUE)
  
  norm.sim.loads.X = chord.norm.diff(cbind(JntLd.X, IndivLd.X), PJIVE.loads.X)
  norm.sim.loads.Y = chord.norm.diff(cbind(JntLd.Y, IndivLd.Y), PJIVE.loads.Y)
  
  abs.err.sim.loads.X = abs(cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)
  abs.err.sim.loads.X[abs(cbind(JntLd.X, IndivLd.X)) != 0] = (abs(cbind(JntLd.X, IndivLd.X) - PJIVE.loads.X)/abs(cbind(JntLd.X, IndivLd.X)))[abs(cbind(JntLd.X, IndivLd.X)) != 0]
  
  abs.err.sim.loads.Y = abs(cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)
  abs.err.sim.loads.Y[abs(cbind(JntLd.Y, IndivLd.Y)) != 0] = (abs(cbind(JntLd.Y, IndivLd.Y) - PJIVE.loads.Y)/abs(cbind(JntLd.Y, IndivLd.Y)))[abs(cbind(JntLd.Y, IndivLd.Y)) != 0]
  
  out = list("PJIVE.Loads.X1" = sim.loads.X, "PJIVE.Loads.X2" = sim.loads.Y,
             "MSE.PJIVE.Loads.X1" = mse.sim.loads.X, "MSE.PJIVE.Loads.X2" = mse.sim.loads.Y,
             "Abs.Err.PJIVE.Loads.X1" = abs.err.sim.loads.X, "Abs.Err.PJIVE.Loads.X2" = abs.err.sim.loads.Y,
             "Norm.PJIVE.Loads.X1" = norm.sim.loads.X, "Norm.PJIVE.Loads.X2" = norm.sim.loads.Y,
             "PMSE.PJIVE.Loads.X1" = pmse.sim.loads.X, "PMSE.PJIVE.Loads.X2" = pmse.sim.loads.Y)
  
  save(out, file = file.path(outdir, nm))
  if(file.exists(file.path(outdir, nm))) print("Results saved.")
}
