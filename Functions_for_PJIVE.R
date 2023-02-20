###########################################################################################################################
#######                 Functions for JIVE-PNC Project(s)                              ####################################
#######                 Author: Raphiel J. Murden                                      ####################################
#######                 Supervised by Benjamin Risk                                    ####################################
###########################################################################################################################
# require(rootSolve); require(optimx); require(gplots); require(r.jive); 
ajive.dir = "H:/My Documents/Applications2/r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

require(Matrix); require(ggplot2); require(reshape2); require(fields); require(mvtnorm)
require(dplyr); require(xtable);  require(MASS); require(extraDistr); require(stringr)

######################################################################################################################
###########   Generates K Simulated Datasets that follow JIVE Model using binary subject scores   ####################
######################################################################################################################
GenerateToyData <- function(n, p, JntVarEx, IndVarEx, jnt_rank = 1, equal.eig = FALSE, ind_ranks, JntVarAdj = TRUE, mix.probs = NULL,
                                 SVD.plots = TRUE, Error = TRUE, print.cor = TRUE, Loads = "Rademacher", Scores = "Gaussian_Mixture", 
                            error.variances = NULL){
  
  # Write out both joint and indiv subject scores for both data sets first
  r.J = jnt_rank
  r.I = ind_ranks
  K = length(p)
  Scores.text = ifelse(is.numeric(Scores), "as defined by user", paste("from", Scores, "distributions"))
  Loads.text = ifelse(is.list(Loads), "as defined by user.", paste("from", Loads, "distributions."))
  cat(paste0("Generating Scores ", Scores.text, " and Loadings ", Loads.text, ". \n"))
  
  # Generate scores
  if(Scores=="Binomial"){
      JntScores = matrix(rbinom(n*r.J, size=1, prob=0.2), nrow = n, ncol = r.J)
      
      b = rbinom(n*sum(r.I), size=1, prob=0.4)
      b = 1 - 2*b
      IndivScores = matrix(b, nrow = n, ncol = sum(r.I))
  } else if (Scores=="Gaussian_Mixture"){
      if(is.null(mix.probs)){mix.probs = c(0.2, 0.5, 0.3)}
        n.groups = length(mix.probs)
        n.vals = (n.groups-1)/2 
      if(is.integer(length(mix.probs)/2)){
        group.means = c(-4*(n.vals:1),4*(1:n.vals))
      } else {
        group.means = c(-4*(n.vals:1),0,4*(1:n.vals))
      }
      JointScores.groups = list()
      for(l in 1:length(mix.probs)){JointScores.groups[[l]] = matrix(rnorm(n*mix.probs[l]*r.J, mean = group.means[l]), ncol = r.J)}
      JntScores = do.call(rbind, JointScores.groups)
      IndivScores = matrix(rnorm(n*sum(r.I)), nrow = n, ncol = sum(r.I))
  } else if(Scores=="Gaussian"){
      JntScores = matrix(rnorm(n*r.J), ncol = r.J)
      IndivScores = matrix(rnorm(n*sum(r.I)), nrow = n, ncol = sum(r.I))
  } else if (is.numeric(Scores)){
      if(nrow(Scores == n) & ncol(Scores = r.J+sum(r.I))){
        JntScores = Scores[,1:r.J]
        IndivScores = Scores[,-(1:r.J)]
      } else{
        cat("Matrix of subject scores must have n rows and jnt_rank+sum(ind_ranks) columns.")
      }
  } else {
    message("Please use one of the following three options to generate data with corresponding subject scores: 1) 'Gaussian' 2) 'Gaussian_Mixture' 3) Binomial.\n")
    message("Alternatively, the user can enter a matrix of subject scores with n rows and jnt_rank+sum(ind_ranks) columns.\n")
  }
  
  colnames(JntScores) = paste("Jnt Score", 1:r.J)
  IndivScores.names = NULL
  
  for(k in 1:K){
    IndivScores.names = c(IndivScores.names, paste0("Indiv X", k, " Score ", 1:r.I[k]))
  }
  colnames(IndivScores) = IndivScores.names
  
  if(print.cor){cat("The correlation between subject scores is given by"); print(round(cor(cbind(JntScores, IndivScores)),4))}
  
  # Generate Loadings
  Jnt.Loads.All = list()
  Indiv.Loads.All = list()
  D.J = D.I = list()
  Noise = Joint.Sigs = Indiv.Sigs = list()
  Sig.Mats =  Data.Mats = list()
  temp.fcn = function(x){x[sample(round(length(x)/2))] = 1; x}
  
  if(is.null(error.variances)){
    error.variances = rep(1,K)
  } else if(is.numeric(error.variances) & length(error.variances) != K){
    cat("Number of entries in 'error.variances' must equal number of data sets and thus have the same length as 'p'.\n")
  }
  
  for(k in 1:K){
    if(length(Loads)==1){
      if(Loads == "Gaussian"){
        Jnt.Loads.All[[k]] = matrix(rnorm(r.J*p[k]), nrow = r.J, ncol = p[k])
      } else if (Loads == "Fixed"){
        Jnt.Loads.All[[k]] = matrix(apply(matrix(0, nrow = r.J, ncol = p[k]), 1, temp.fcn), nrow = r.J)
      } else if (Loads == "Double_Exp"){
        Jnt.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.J)), nrow = r.J)
      } else if (Loads == "Rademacher"){
        Jnt.Loads.All[[k]] = matrix(rsign(p[k]*(r.J)), nrow = r.J)
      }
    } else if (is.list(Loads)){
      Jnt.Loads.All[[k]] = t(Loads[[1]][[k]])
    } 
    D.J[[k]] = (1 - equal.eig)*diag(r.J:1) + equal.eig*diag(rep(1,r.J))
    
    Joint.Sigs[[k]] = JntScores%*%sqrt(D.J[[k]])%*%Jnt.Loads.All[[k]]
    
    if(SVD.plots){
      plot(svd(Joint.Sigs[[k]])$d, ylab = "Singular Values")
      title(paste0("SVD of Joint Signal from X", k))
    }
    temp.IndScores = IndivScores[,(k>1)*sum(r.I[1:(k-1)])+(1:r.I[k])]
    
    if(length(Loads)==1){
      if(Loads == "Gaussian"){
        Indiv.Loads.All[[k]] = matrix(rnorm(n = p[k]*r.I[k]), nrow = r.I[k], ncol = p[k])
      } else if (Loads == "Fixed"){
        temp.fcn = function(x){x[sample(round(length(x)/4))] = 1; x}
        Indiv.Loads.All[[k]] = matrix(apply(matrix(-1, nrow = r.I[k], ncol = p[k]), 1, temp.fcn), nrow = r.I[k])
      } else if (Loads == "Double_Exp"){
        Indiv.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.I[k])), nrow = r.I[k])
      } else if (Loads == "Rademacher"){
        Indiv.Loads.All[[k]] = matrix(rsign(p[k]*(r.I[k])), nrow = r.I[k])
      }
    } else if (is.list(Loads)){
      Indiv.Loads.All[[k]] = t(Loads[[2]][[k]])
    } 
    
    D.I[[k]] = (1 - equal.eig)*diag(r.I[k]:1) + equal.eig*diag(rep(1,r.I[k]))
    
    Indiv.Sigs[[k]] = temp.IndScores%*%sqrt(D.I[[k]])%*%Indiv.Loads.All[[k]]
    
    if(SVD.plots){
      plot(svd(Indiv.Sigs[[k]])$d, ylab = "Singular Values")
      title(paste0("SVD of Individual Signal from X", k))
    }
    
    Sig.Mats[[k]] = Joint.Sigs[[k]] + Indiv.Sigs[[k]]
    
    Noise[[k]] = matrix(rnorm(n*p[k], sd = sqrt(error.variances[k])), nrow  = n)
    
    Data.Mats[[k]] = AdjSigVarExp(Joint.Sigs[[k]], Indiv.Sigs[[k]], Noise[[k]],
                                  JntVarEx[k], IndVarEx[k])
    
    Joint.Sigs[[k]] = Data.Mats[[k]]$J
    Indiv.Sigs[[k]] = Data.Mats[[k]]$I
  }
  
  ##############################Define Datasets##############################
  Blocks = lapply(Data.Mats, function(x) x[["Data"]])
  
  Dat.Comps = list(Joint.Sigs, Indiv.Sigs, Noise)
  names(Dat.Comps) = c("JointSignalMatrices", "IndivSignalMatrices", "NoiseMatrices")
  
  Scores = list(JntScores, IndivScores)
  names(Scores) = c("Joint", "Indiv")
  
  Loadings = list(list(),list())
  names(Loadings) = c("Joint", "Indiv")
  for(k in 1:K){
    Loadings[["Joint"]][[k]] = D.J[[k]]%*%Jnt.Loads.All[[k]]
    Loadings[["Indiv"]][[k]] = D.I[[k]]%*%Indiv.Loads.All[[k]]
  }
  
  out = list(Dat.Comps, Blocks, Scores, Loadings)
  names(out) = c("Data Components", "Data Blocks", "Scores", "Loadings")
  
  out
}


#######################################################################################################################
###########   Generates 2 Simulated Datasets that follow JIVE Model using binary subject scores    ####################
###########   NOTE: This version sohuld be archived once confirmed that it coresponds with newer   ####################
###########           version above built for K>2                                                  ####################
#######################################################################################################################
GenerateToyData_OLD <- function(n, p1, p2, JntVarEx1, JntVarEx2, IndVarEx1, IndVarEx2, jnt_rank = 1, equal.eig = FALSE,
                             ind_rank1 = 2, ind_rank2 = 2, JntVarAdj = TRUE, SVD.plots = TRUE, Error = TRUE, print.cor = TRUE,
                             Loads = "Rademacher", Scores = "Gaussian_Mixture"){
  
  # Write out both joint and indiv subject scores for both data sets first
  r.J = jnt_rank
  r.I1 = ind_rank1
  r.I2 = ind_rank2
  
  print(paste("Generating Scores from", Scores, "distributions and Loadings from", Loads, "distributions."))
  
  # Generate Joint scores
  if(Scores=="Binomial"){
    JntScores = matrix(rbinom(n*r.J, size=1, prob=0.2), nrow = n, ncol = r.J)
    
    b = rbinom(n*(r.I1 + r.I2), size=1, prob=0.4)
    b = 1 - 2*b
    IndivScores = matrix(b, nrow = n, ncol = (r.I1 + r.I2))
  } else if (Scores=="Gaussian_Mixture"){
    mix.probs = c(0.2, 0.5, 0.3)
    JointScores.g1 = matrix(rnorm(n*mix.probs[1]*r.J), ncol = r.J)
    JointScores.g2 = matrix(rnorm(n*mix.probs[2]*r.J, -4, 1), ncol = r.J)
    JointScores.g3 = matrix(rnorm(n*mix.probs[3]*r.J, 4, 1), ncol = r.J)
    JntScores = rbind(JointScores.g1, JointScores.g2, JointScores.g3)
    
    IndivScores = matrix(rnorm(n*(r.I1+r.I2)), nrow = n, ncol = (r.I1 + r.I2))
  } else if(Scores=="Gaussian"){
    JntScores = matrix(rnorm(n*r.J), ncol = r.J)
    IndivScores = matrix(rnorm(n*(r.I1+r.I2)), nrow = n, ncol = (r.I1 + r.I2))
  }
    
  colnames(JntScores) = paste("Jnt Score", 1:r.J)
  colnames(IndivScores) = c(paste("Ind X Score", 1:r.I1), paste("Ind Y Score", 1:r.I2))
  
  if(print.cor){print("The correlation between subject scores is given by"); print(round(cor(cbind(JntScores, IndivScores)),4))}
  
  ##############################Define X Dataset##############################
  if(Loads == "Gaussian"){
    AdjJntLoad.X = matrix(rnorm(r.J*p1), nrow = r.J, ncol = p1)
  } else if (Loads == "Fixed"){
    temp.fcn = function(x){x[sample(round(length(x)/2))] = 1; x}
    AdjJntLoad.X = matrix(apply(matrix(0, nrow = r.J, ncol = p1), 1, temp.fcn), nrow = r.J)
  } else if (Loads == "Double_Exp"){
    AdjLoad.X = matrix(extraDistr::rlaplace(p1*(r.J+r.I1)), nrow = r.J+r.I1)
    AdjJntLoad.X = AdjLoad.X[1:r.J, , drop = FALSEALSE]
  } else if (Loads == "Rademacher"){
    AdjLoad.X = matrix(extraDistr::rsign(p1*(r.J+r.I1)), nrow = r.J+r.I1)
    AdjJntLoad.X = AdjLoad.X[1:r.J, , drop = FALSEALSE]
  } 
  
  #change relavent scaling of joint components
  D.J = (equal.eig + 1)*diag(r.J:1) + equal.eig*diag(rep(1,r.J))
  JX = JntScores%*%sqrt(D.J)%*%AdjJntLoad.X
  
  if(SVD.plots){
    plot(svd(JX)$d, ylab = "Singular Values")
    title("SVD of Joint Signal from X")
  }
  
  ##Note that the individual signal has rank = 2 as well
  IndScores.X = IndivScores[,1:r.I1]
  IndLoad.X = matrix(rnorm(n = p1*r.I1), nrow = r.I1, ncol = p1)
  if(Loads == "Gaussian"){
    IndLoad.X = matrix(rnorm(n = p1*r.I1), nrow = r.I1, ncol = p1)
  } else if (Loads == "Fixed"){
    temp.fcn = function(x){x[sample(round(length(x)/4))] = 1; x}
    IndLoad.X = matrix(apply(matrix(-1, nrow = r.I1, ncol = p1), 1, temp.fcn), nrow = r.I1)
  } else if (Loads == "Double_Exp"){
    IndLoad.X = AdjLoad.X[-(1:r.J), , drop = FALSEALSE]
  } else if (Loads == "Rademacher"){
    IndLoad.X = AdjLoad.X[-(1:r.J), ,drop = FALSEALSE]
  } 
  D.IX = (equal.eig + 1)*diag((r.I1:1)) + equal.eig*diag(rep(1,r.I1))
  
  IX = IndScores.X%*%sqrt(D.IX)%*%IndLoad.X
  
  if(SVD.plots){
    plot(svd(IX)$d, ylab = "Singular Values")
    title("SVD of Individual Signal from X")
  }
  
  AX = JX + IX
  
  ##############################Define Y Dataset##############################
  if(Loads == "Gaussian"){
    AdjJntLoad.Y = matrix(rnorm(r.J*p2), nrow = r.J, ncol = p2)
  } else if (Loads == "Fixed"){
    temp.fcn = function(x){x[sample(round(length(x)/2))] = 2; x}
    AdjJntLoad.Y = matrix(apply(matrix(1, nrow = r.J, ncol = p2), 1, temp.fcn), nrow = r.J)
  } else if (Loads == "Double_Exp"){
    AdjLoad.Y = matrix(extraDistr::rlaplace(p2*(r.J+r.I2)), nrow = r.J+r.I2)
    AdjJntLoad.Y = AdjLoad.Y[1:r.J, , drop = FALSEALSE]
  } else if (Loads == "Rademacher"){
    AdjLoad.Y = matrix(extraDistr::rsign(p2*(r.J+r.I2)), nrow = r.J+r.I2)
    AdjJntLoad.Y = AdjLoad.Y[1:r.J, , drop = FALSEALSE]
  } 
  
  ##Note that the joint signal has rank = 3
  JY = JntScores%*%sqrt(D.J)%*%AdjJntLoad.Y
  
  if(SVD.plots){
    plot(svd(JY)$d, ylab = "Singular Values")
    title("SVD of Joint Signal from Y")
  }
  
  IndScores.Y = IndivScores[,(r.I1 + 1:r.I2)]
  IndLoad.Y = matrix(rnorm(n = p2*r.I1), nrow = r.I1, ncol = p2)
  if(Loads == "Gaussian"){
    IndLoad.Y = matrix(rnorm(n = p2*r.I1), nrow = r.I1, ncol = p2)
  } else if (Loads == "Fixed"){
    temp.fcn = function(x){x[sample(round(length(x)/4))] = 1; x}
    IndLoad.Y = matrix(apply(matrix(-1, nrow = r.I2, ncol = p2), 1, temp.fcn), nrow = r.I2)
  } else if (Loads == "Double_Exp"){
    IndLoad.Y = AdjLoad.Y[-(1:r.J), , drop = FALSEALSE]
  } else if (Loads == "Rademacher"){
    IndLoad.Y = AdjLoad.Y[-(1:r.J), ,drop = FALSEALSE]
  } 
  D.IY = (equal.eig + 1)*diag((r.I2:1)) + equal.eig*diag(rep(1,r.I2)) 
  
  ##Note that the individual signal has rank=2
  IY = IndScores.Y%*%sqrt(D.IY)%*%IndLoad.Y
  
  if(SVD.plots){
    plot(svd(IY)$d, ylab = "Singular Values")
    title("SVD of Individual Signal from Y")
  }
  
  ##Error matrix
  EX = matrix(rnorm(n*p1), nrow=n, ncol=p1)*Error
  
  ##Error matrix
  EY = matrix(rnorm(n*p2), nrow=n, ncol=p2)*Error
  
  Dat.X = AdjSigVarExp(JX, IX, EX, JntVarEx1, IndVarEx1)
  JX = Dat.X$J
  IX = Dat.X$I
  
  Dat.Y = AdjSigVarExp(JY, IY, EY, JntVarEx2, IndVarEx2)
  JY = Dat.Y$J
  IY = Dat.Y$I
  
  Blocks = list(Dat.X[["Data"]], Dat.Y[["Data"]])
  
  Dat.Comps = list(JX, JY, IX, IY, EX, EY)
  names(Dat.Comps) = c("J1", "J2", "I1", "I2", "E1", "E2" )
  
  Scores = list(JntScores, IndScores.X, IndScores.Y)
  names(Scores) = c("Joint", "Indiv_1", "Indiv_2")
  
  Loadings = list(sqrt(D.J)%*%AdjJntLoad.X, sqrt(D.IX)%*%IndLoad.X, 
                  sqrt(D.J)%*%AdjJntLoad.Y, sqrt(D.IY)%*%IndLoad.Y)
  names(Loadings) = c("Joint_1", "Indiv_1", "Joint_2", "Indiv_2")
  
  out = list(Dat.Comps, Blocks, Scores, Loadings)
  names(out) = c("Data Components", "Data Blocks", "Scores", "Loadings")
  
  out
}

####################################################################################
###########   Adjust Dataset Components to get Desired R^2 Values  #################
####################################################################################
AdjSigVarExp <-function(J, I, E, JntVarEx, IndVarEx){
  simul.quads = function(x, parms){
    JJ = parms[1]
    II = parms[2]
    EE = parms[3]
    JE = parms[4]
    IE = parms[5]
    R_J = parms[6]
    R_I = parms[7]
    
    y1 = x[1]^2*II*(1 - R_I) - 2*x[1]*IE*R_I - R_I*(x[2]^2*JJ + 2*x[2]*JE + EE)
    y2 = x[2]^2*JJ*(1 - R_J) - 2*x[2]*JE*R_J - R_J*(x[1]^2*II + 2*x[1]*IE + EE)
    
    y = c(y1,y2)
    return(y)
  }
  
  JJ = MatVar(J)
  II = MatVar(I)
  EE = MatVar(E)
  JE = sum(diag(J%*%t(E)))
  IE = sum(diag(I%*%t(E)))
  R_J = JntVarEx
  R_I = IndVarEx
  
  parms = c(JJ, II, EE, JE, IE, R_J, R_I)
  
  A = J + I
  AA = MatVar(A)
  EE = MatVar(E)
  AE = sum(diag(A%*%t(E)))
  
  ##Desired Total Variance Explained
  d0 = IndVarEx + JntVarEx
  a = AA*(1 - d0)
  b = -2*AE
  c = -d0*EE
  
  d.A = (-b+sqrt(b^2 - 4*a*c))/(2*a)
  
  start = c(0.5, 0.5)*d.A
  
  roots = rootSolve::multiroot(simul.quads, start, parms = parms)
  c = roots$root[1]
  d = roots$root[2]
  
  J.d = d*J; I.c = c*I;
  Dat = J.d + I.c + E
  
  res = list(J.d, I.c, Dat)
  names(res) = c("J", "I", "Data")
  return(res)
}


#######################################################################################################################
#####################   Convert ProJIVE  Simulation Results to a form that will allow ggplot       ####################
#######################################################################################################################
ConvSims_gg_ProJIVE<-function(AllSims){
  sim.names = colnames(AllSims)
  #CompEvals.names = c(sim.names[grep("Scores", sim.names)], sim.names[grep("Loads", sim.names)])
  num.sims = dim(AllSims)[1]
  
  p1 = unique(AllSims$p1)
  p2 = unique(AllSims$p2)
  
  JVE_1 =  c(as.matrix(AllSims[,grep("JntVarEx1",sim.names)]))
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))
  JVE_1.raw = factor(JVE_1, labels = JVE1.labs, levels = c(0.05, 0.5))
  
  JVE_2 = c(as.matrix(AllSims[,grep("JntVarEx2",sim.names)]))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)), bquote("R"[J2]^2*"=0.5, p"[2]*"="*.(p2)) )
  JVE_2.raw = factor(JVE_2, labels = JVE2.labs, levels = c(0.05, 0.5))
  
  IVE_1 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.X", sim.names)]))
  IVE_1 = round(as.numeric(IVE_1), 2)
  
  IVE_2 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.Y", sim.names)]))
  IVE_2 = round(as.numeric(IVE_2), 2)
  
  Norm = c(as.numeric(as.matrix(AllSims[,grep("Scores", sim.names)])))
  
  Method1 = sim.names[grep("Scores",sim.names)]
  Method2 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Method3 = str_replace_all(Method2, "\\.", " ")
  Method4 = str_replace_all(Method3, " [:digit:]", "")
  Method5 = str_replace_all(Method4, " Joint.+", "")
  Method6 = str_replace_all(Method5, " Indiv.+", "")
  Method7 = str_replace_all(Method6, "R JIVE", "R.JIVE")
  Method = rep(factor(Method7, levels = unique(Method7)), each = num.sims)
  
  Type1 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Type2 = str_replace_all(Type1, "\\.", " ")
  Type3 = str_replace_all(Type2, " [:digit:]", "")
  Type4 = str_replace_all(Type3, "X", "X[1]")
  Type5 = str_replace_all(Type4, "Y", "X[2]")
  Type6 = str_replace_all(Type5, "^.+Joint", "Joint")
  Type7 = str_replace_all(Type6, "^.+Indiv", "Indiv")
  Type = rep(as.factor(Type7), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  sim.score.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1.raw, each = n.levs), 
                                  JVE_2 = rep(JVE_2.raw, each = n.levs),
                                  IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  Norm = c(as.numeric(as.matrix(AllSims[,grep("Loads", sim.names)])))
  
  Method1 = sim.names[grep("Loads",sim.names)]
  Method2 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Method3 = str_replace_all(Method2, "\\.", " ")
  Method4 = str_replace_all(Method3, " [:digit:]", "")
  Method5 = str_replace_all(Method4, " Joint.+", "")
  Method6 = str_replace_all(Method5, " Indiv.+", "")
  Method7 = str_replace_all(Method6, "R JIVE", "R.JIVE")
  Method = rep(factor(Method7, levels = unique(Method7)), each = num.sims)
  
  Type1 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Type2 = str_replace_all(Type1, "\\.", " ")
  Type3 = str_replace_all(Type2, " [:digit:]", "")
  Type4 = str_replace_all(Type3, "X", "X[1]")
  Type5 = str_replace_all(Type4, "Y", "X[2]")
  Type6 = str_replace_all(Type5, "^.+Joint", "Joint")
  Type7 = str_replace_all(Type6, "^.+Indiv", "Indiv")
  Type = rep(as.factor(Type7), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  JVE_1 = rep(JVE_1.raw, each = n.levs)
  JVE_2 = rep(JVE_2.raw, each = n.levs)
  
  sim.load.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = JVE_1, JVE_2 = JVE_2, 
                                 IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))

  sim.all.norms.gg = rbind(sim.score.norms.gg, sim.load.norms.gg)
  out = sim.all.norms.gg 
  out
}

#######################################################################################################################
#####################   Convert ProJIVE  Simulation Results to a form that will allow ggplot       ####################
#######################################################################################################################
ConvSims_gg_ProJIVE2<-function(AllSims, n){
  sim.names = colnames(AllSims)
  #CompEvals.names = c(sim.names[grep("Scores", sim.names)], sim.names[grep("Loads", sim.names)])
  num.sims = dim(AllSims)[1]
  
  p1 = unique(AllSims$p1)
  p2 = unique(AllSims$p2)
  
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)*", n="*.(n)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))# *", n="*.(n)))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)*", n="*.(n)), bquote("R"[J2]^2*"=0.5, p"[2]*"="*.(p2)))# *", n="*.(n)))
  
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))# *", n="*.(n)))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)), bquote("R"[J2]^2*"=0.5, p"[2]*"="*.(p2)))# *", n="*.(n)))

  JVE_1 =  c(as.matrix(AllSims[,grep("JntVarEx1",sim.names)]))
  JVE_1.raw = factor(JVE_1, labels = JVE1.labs, levels = c(0.05, 0.5))
  
  JVE_2 = c(as.matrix(AllSims[,grep("JntVarEx2",sim.names)]))
  JVE_2.raw = factor(JVE_2, labels = JVE2.labs, levels = c(0.05, 0.5))
  
  IVE_1 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.X", sim.names)]))
  IVE_1 = round(as.numeric(IVE_1), 2)
  
  IVE_2 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.Y", sim.names)]))
  IVE_2 = round(as.numeric(IVE_2), 2)
  
  Norm = c(as.numeric(
    as.matrix(
      AllSims[,sort(c(grep("Scores..Norm",sim.names), grep("Scores.X..Norm",sim.names), grep("Scores.Y..Norm",sim.names)))])))
  
  Method1 = sim.names[sort(c(grep("Scores..Norm",sim.names), grep("Scores.X..Norm",sim.names), grep("Scores.Y..Norm",sim.names)))]
  Method2 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Method3 = str_replace_all(Method2, "\\.", " ")
  Method4 = str_replace_all(Method3, " [:digit:]", "")
  Method5 = str_replace_all(Method4, " Joint.+", "")
  Method6 = str_replace_all(Method5, " Indiv.+", "")
  Method7 = str_replace_all(Method6, "R JIVE", "R.JIVE")
  Method = rep(factor(Method7, levels = unique(Method7)), each = num.sims)
  
  Type1 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Type2 = str_replace_all(Type1, "\\.", " ")
  Type3 = str_replace_all(Type2, " [:digit:]", "")
  Type4 = str_replace_all(Type3, "X", "X[1]")
  Type5 = str_replace_all(Type4, "Y", "X[2]")
  Type6 = str_replace_all(Type5, "^.+Joint", "Joint")
  Type7 = str_replace_all(Type6, "^.+Indiv", "Indiv")
  Type = rep(as.factor(Type7), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  sim.score.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1.raw, each = n.levs), 
                                  JVE_2 = rep(JVE_2.raw, each = n.levs),
                                  IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  Norm = c(as.numeric(as.matrix(AllSims[,sort(c(grep("Loads.X..Norm", sim.names), grep("Loads.Y..Norm", sim.names)))])))
  
  Method1 = sim.names[sort(c(grep("Loads.X..Norm", sim.names), grep("Loads.Y..Norm", sim.names)))]
  Method2 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Method3 = str_replace_all(Method2, "\\.", " ")
  Method4 = str_replace_all(Method3, " [:digit:]", "")
  Method5 = str_replace_all(Method4, " Joint.+", "")
  Method6 = str_replace_all(Method5, " Indiv.+", "")
  Method7 = str_replace_all(Method6, "R JIVE", "R.JIVE")
  Method = rep(factor(Method7, levels = unique(Method7)), each = num.sims)
  
  Type1 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Type2 = str_replace_all(Type1, "\\.", " ")
  Type3 = str_replace_all(Type2, " [:digit:]", "")
  Type4 = str_replace_all(Type3, "X", "X[1]")
  Type5 = str_replace_all(Type4, "Y", "X[2]")
  Type6 = str_replace_all(Type5, "^.+Joint", "Joint")
  Type7 = str_replace_all(Type6, "^.+Indiv", "Indiv")
  Type8 = str_replace_all(Type7, "Loads", "Loadings")
  Type = rep(as.factor(Type8), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  JVE_1 = rep(JVE_1.raw, each = n.levs)
  JVE_2 = rep(JVE_2.raw, each = n.levs)
  
  sim.load.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = JVE_1, JVE_2 = JVE_2, 
                                 IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  sim.all.norms.gg = rbind(sim.score.norms.gg, sim.load.norms.gg)
  
  PMSE = c(as.numeric(
    as.matrix(
      AllSims[,sort(c(grep("Scores..PMSE",sim.names), grep("Scores.X..PMSE",sim.names), grep("Scores.Y..PMSE",sim.names)))])))
  
  Method1 = sim.names[sort(c(grep("Scores..PMSE",sim.names), grep("Scores.X..PMSE",sim.names), grep("Scores.Y..PMSE",sim.names)))]
  Method2 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Method3 = str_replace_all(Method2, "\\.", " ")
  Method4 = str_replace_all(Method3, " [:digit:]", "")
  Method5 = str_replace_all(Method4, " Joint.+", "")
  Method6 = str_replace_all(Method5, " Indiv.+", "")
  Method7 = str_replace_all(Method6, "R JIVE", "R.JIVE")
  Method = rep(factor(Method7, levels = unique(Method7)), each = num.sims)
  
  Type1 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Type2 = str_replace_all(Type1, "\\.", " ")
  Type3 = str_replace_all(Type2, " [:digit:]", "")
  Type4 = str_replace_all(Type3, "X", "X[1]")
  Type5 = str_replace_all(Type4, "Y", "X[2]")
  Type6 = str_replace_all(Type5, "^.+Joint", "Joint")
  Type7 = str_replace_all(Type6, "^.+Indiv", "Indiv")
  Type = rep(as.factor(Type7), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  sim.score.pmse.gg = data.frame(PMSE = PMSE, Method = Method, Type = Type, JVE_1 = rep(JVE_1.raw, each = n.levs), 
                                  JVE_2 = rep(JVE_2.raw, each = n.levs),
                                  IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  PMSE = c(as.numeric(as.matrix(AllSims[,sort(c(grep("Loads.X..PMSE", sim.names), grep("Loads.Y..PMSE", sim.names)))])))
  
  Method1 = sim.names[sort(c(grep("Loads.X..PMSE", sim.names), grep("Loads.Y..PMSE", sim.names)))]
  Method2 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Method3 = str_replace_all(Method2, "\\.", " ")
  Method4 = str_replace_all(Method3, " [:digit:]", "")
  Method5 = str_replace_all(Method4, " Joint.+", "")
  Method6 = str_replace_all(Method5, " Indiv.+", "")
  Method7 = str_replace_all(Method6, "R JIVE", "R.JIVE")
  Method = rep(factor(Method7, levels = unique(Method7)), each = num.sims)
  
  Type1 = str_replace_all(Method1, "\\.\\.\\.", " ") 
  Type2 = str_replace_all(Type1, "\\.", " ")
  Type3 = str_replace_all(Type2, " [:digit:]", "")
  Type4 = str_replace_all(Type3, "X", "X[1]")
  Type5 = str_replace_all(Type4, "Y", "X[2]")
  Type6 = str_replace_all(Type5, "^.+Joint", "Joint")
  Type7 = str_replace_all(Type6, "^.+Indiv", "Indiv")
  Type8 = str_replace_all(Type7, "Loads", "Loadings")
  Type = rep(as.factor(Type8), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  JVE_1 = rep(JVE_1.raw, each = n.levs)
  JVE_2 = rep(JVE_2.raw, each = n.levs)
  
  sim.load.pmse.gg = data.frame(PMSE = PMSE, Method = Method, Type = Type, JVE_1 = JVE_1, JVE_2 = JVE_2, 
                                 IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  sim.all.pmse.gg = rbind(sim.score.pmse.gg, sim.load.pmse.gg)
  
  temp = as.matrix(AllSims[,grep("VarExp_", sim.names)])
  temp = temp[,str_which(colnames(temp), "True", negate = T)]
  temp.names = colnames(temp)
  VarEx = c(as.numeric(temp))
  
  Method1 = rep(temp.names[grep("VarExp_",temp.names)], each = num.sims)
  Method2 = str_replace_all(Method1, "JIVE_Oracle", "JIVE-Oracle")
  Method3 = str_replace(Method2, "\\.[:digit:]", "") #
  Method4 = str_replace_all(Method3, "_.+", "") 
  Method5 = str_replace_all(Method4, "Joint.+", "")
  Method6 = str_replace_all(Method5, "Indiv.+", "")
  Method = factor(Method6, levels = unique(Method6)[c(4,3,1,2,5)])#, levels = unique(Method6)[c(4,3,1,2,5)])
  #table(Method, Method6)
  
  Type1 = rep(temp.names[grep("VarExp_", temp.names)], num.sims)
  Type2 = str_replace(Type1, "\\.[:digit:]", "")
  Type3 = str_replace(Type2, "_X", "_X\\[1\\]")
  Type4 = str_replace(Type3, "_Y", "_X\\[2\\]")
  Type5 = str_replace(Type4, ".+Joint", "Joint")
  Type6 = str_replace(Type5, ".+Indiv", "Indiv")
  Type = as.factor(Type6)
  
  n.levs = nlevels(Type)*nlevels(Method)
  JVE_1 = rep(JVE_1.raw, each = n.levs)
  JVE_2 = rep(JVE_2.raw, each = n.levs)
  
  VarEx.out = data.frame(VarEx = VarEx, Type = Type, Method = Method, JVE_1 = JVE_1.raw, JVE_2 = JVE_2.raw, 
                         IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  out = list(sim.all.norms.gg, sim.all.pmse.gg, VarEx.out)
  names(out) = c("Norms", "PMSEs", "VarEx")
  out
}

#################################################################################################
#####################         Plot PJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.score.pmse.plot<-function(pmse.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(pmse.dat$Type)[grep("Score",levels(pmse.dat$Type))][c(3,1,2)]
  labs.ex = c("Joint Subj Scores", expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]))
  ggplot(data = pmse.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs, labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}

#################################################################################################
#####################         Plot PJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.load.pmse.plot<-function(pmse.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(pmse.dat$Type)[grep("Load",levels(pmse.dat$Type))][c(3,4,1,2)]
  labs.ex = c(expression("Joint Variable Loadings"*"X"[1]),expression("Joint Variable Loadings"*"X"[2]),
              expression("Indiv Variable Loadings"*"X"[1]), expression("Indiv Variable Loadings"*"X"[2]))
  ggplot(data = pmse.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs, labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}

########################################################################
####################START OF PRE-DEFINED FUNCS##########################
########################################################################
## FUNC to Generate the corresponding Ak and Bk:
generate_ab=function(m,n){
  return(matrix(rep(0,n*m),nrow=n,ncol=m))
}


## FUNC to generate list(G) of W matrix with right dimension
generate_w=function(k, wk_hat, P, Q){
  r=lapply(Q,generate_ab,P[k])
  r[[1]]=as.matrix(wk_hat[[k]][1:P[k],1:Q[1]], nrow=P[k], ncol=Q[1])
  r[[k+1]]=as.matrix(wk_hat[[k]][1:P[k],Q[1]+1:Q[k+1]], nrow=P[k], ncol=Q[k+1])
  return(do.call('cbind',r))
}

## FUNC to generate list(G) of W matrix from W_k matrices
wk_to_w=function(wk, P, Q){
  
  w=do.call('rbind',lapply(1:length(wk),function(k)generate_w(k, wk, P, Q)))
  
  return(w)
}

## FUNC to generate list of W_k matrices from W
w_to_w_k=function(w, P, Q){
  
  P=c(0,P)
  K=length(P[-1])
  w_k=list()
  for(k in 1:K){
    w_k[[k]] = w[sum(P[1:k])+(1:P[k+1]),c(1:Q[1],sum(Q[1:k])+1:Q[k+1])]
  }
  
  
  return(w_k)
}

## FUNC to generate list(G) of D matrix from sigma vectors
generate_d=function(sig_lst, p_vec){
  
  d=NULL
  for(k in 1:length(p_vec)){
    d = c(d, rep(sig_lst[k],p_vec[k]))
  }
  D=diag(d)
  
  return(D)
}

obs_LogLik<-function(Y, mu, w, d){
  ##############################################
  # input:    -mu    :a G list of d dimension vectors
  #           -w     :a G list of dxp matrices
  #           -d     :a G list of d length vector indicating the noise
  #           -Y     :a nxd data frame as the observations
  #           
  # output:   a real value of the log likelihood
  ##############################################
  
  N=dim(Y)[1]
  
  
  lik<-rep(0, N)
  
  
  c=w%*%t(w)+d
  s=t(Y)%*%Y/N
  
  LogLik=-N/2*(ncol(Y)*log(2*pi) + log(det(c)) + sum(diag(solve(c)%*%s)))
  
  # lik=dmvnorm(Y,mu,c,log = TRUE)
  # LogLik=sum(lik)
  
  return(LogLik)
}

complete_LogLik<-function(Y, theta, mu, w, d){
  ##############################################
  # input:    -mu    :a G list of d dimension vectors
  #           -w     :a G list of dxp matrices
  #           -d     :a G list of d length vector indicating the noise
  #           -Y     :a nxd data frame as the observations
  #           
  # output:   a real value of the log likelihood
  ##############################################
  
  N=dim(Y)[1]
  
  constant.term = -(N/2)*(sum(c(ncol(Y),ncol(w)))*log(2*pi) + log(det(d)))
  
  kernel.term = 0
  for(i in 1:N){
    Yi = matrix(Y[i,], ncol = 1)
    theta_i = matrix(theta[i,], ncol = 1)
    kernel.term = kernel.term + t(Yi - w%*%theta_i)%*%diag(diag(d)^-1)%*%(Yi - w%*%theta_i) + sum(theta_i^2)
  }
  
  LogLik=constant.term-0.5*kernel.term 
  
  return(LogLik)
}

## FUNC to evaluate convergence in the loop
eval_converge=function(vals_vec,all_obs.LogLik, diff.tol){
  len = length(all_obs.LogLik)
  if(length(vals_vec)==1) {
    return(TRUE)
  }else{
    if(vals_vec[length(vals_vec)]==-Inf){
      return(TRUE)
    }else{
      diff.ll = all_obs.LogLik[len]-all_obs.LogLik[len-1]
      return((diff.ll>=(diff.tol) & diff.ll>0))
    }
  }
}
########################################################################
####################END OF PRE-DEFINED FUNCS############################
########################################################################

###############################################################################################################
###########          pJIVE ML estimation of JIVE model that uses EM algorithm                ##################
###########       This version was originally built for K>2                                  ##################
###############################################################################################################
ProJIVE_EM=function(Y,P,Q,Max.iter=10000,diff.tol=1e-5,plots=TRUE,chord.tol=-1,sig_hat=NULL, init.loads = NULL){
  
  ## init.loads must be a list of two lists - first item contains a list of joint loading matrices
  ##                                          second item is a list of indiv loading matrices
  ##                                          each matrix must have dimension p_k-by-r_
  
  # Total sample size
  N=dim(Y)[1]
  
  # Total number feature blocks
  if(length(P)==(length(Q)-1)){
    K=length(P)
  } else{
    stop("Error: The number of feature blocks does not match the number of invididual signals")
  }
  Q.tot = Q[1] + Q[-1]
  
  # Selection matrices A_knd B_k
  A=list()
  B=list()
  for(k in 1:K){
    up=lapply(Q,generate_ab,Q[1])
    up[[1]]=diag(Q[1])
    down=lapply(Q,generate_ab,Q[(k+1)])
    down[[(k+1)]]=diag(Q[(k+1)])
    B[[k]]=rbind(do.call("cbind",up),
                 do.call("cbind",down))
    
    a=lapply(P,generate_ab,P[k])
    a[[k]]=diag(P[k])
    A[[k]]=do.call("cbind",a)
  }
  
  # get initial estimates of loadings matrices via cc.jive
  if (is.null(init.loads)){
    WJ = WI = list()
    for(k in 1:K){
      WJ[[k]] = matrix(rnorm(Q[1]*P[k]), nrow = P[k])
      WI[[k]] = matrix(rnorm(Q[k+1]*P[k]), nrow = P[k])
    }
  } else if (is.list(init.loads)){
    WJ = init.loads[[1]]
    WI = init.loads[[2]]
  } else if(init.loads == "AJIVE" | init.loads == "CJIVE"){
    dat.blocks = list()
    dat.blocks[[1]] = Y[,1:P[1]]
    for(k in 2:K){
      dat.blocks[[k]] = Y[,cumsum(P[k-1])+(1:P[k])]
    }
    ajive.solution = ajive(dat.blocks, initial_signal_ranks = Q[1]+Q[-1], joint_rank = Q[1])
    
    WJ = lapply(ajive.solution$block_decomps, function(x) x$joint$v)
    WI = lapply(ajive.solution$block_decomps, function(x) x$individual$v)
  }
  
  # Block specific loading matrices W_k
  wk_hat=wji_hat=list()
  for(k in 1:K){ 
    wji_hat[[k]]=cbind(WJ[[k]], WI[[k]])
  }
  wk_hat=wji_hat
  
  # Total loading matrices W
  w_hat=wk_to_w(wk_hat,P,Q)
  
  # Initializationf of sig_hat
  # sig_hat=list(rep(1,K))                       
  if(is.null(sig_hat)){
    sig_hat = rnorm(length(P))
  } else if (is.character(sig_hat) & sig_hat[1] == "MLE"){
    for(k in 1:K){
      temp.Y = Y[,(k>1)*sum(P[1:(k-1)])+(1:P[k])]
      sig_hat[k] = mean(svd(temp.Y)$d[-(1:Q.tot[k])]^2)
    }
    rm(temp.Y)
  }
  
  sig_hat = as.numeric(sig_hat)
  # Total noise matrices D
  d_hat=generate_d(sig_hat,P)
  
  
  # Initiate Chordal Distance
  chord.dist = NULL
  mu_hat=apply(as.matrix(Y),2,sum) / N   
  
  ## Store some values to save computation time
  Yc=sweep(Y, 2, mu_hat)
  S=t(Yc)%*%Yc
  
  Iq=diag(sum(Q))
  Ip=diag(sum(P))
  
  # Initiate LogLik
  # all_obs.LogLik=c(-Inf)
  # all_complete.LogLik = c(-Inf)
  c_solv=solve(Iq+t(w_hat)%*%solve(d_hat)%*%w_hat)
  exp.theta = U = Y%*%solve(d_hat)%*%w_hat%*%c_solv
  
  all_obs.LogLik=obs_LogLik(Y, mu_hat, w_hat, d_hat)
  all_complete.LogLik = complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat)
  
  # Set initial iteration number:
  iter=1
  
  #create flag to stop EM alrogithm if/when necessary
  flag=FALSE
  while (Max.iter>=iter 
         & eval_converge(all_obs.LogLik,all_obs.LogLik,N*diff.tol)
         & (ifelse(is.null(chord.dist), 1, chord.dist[iter-1]) > chord.tol)
         & eval_converge(all_complete.LogLik,all_complete.LogLik,N*diff.tol)
         & flag==FALSE) 
  {
    ################## START OF EM-ALGORITHM ######################
    w=w_hat
    d=d_hat
    
    c_solv=solve(Iq+t(w)%*%solve(d_hat)%*%w)
    c_inv = solve(d_hat)%*%(diag(sum(P))-w%*%c_solv%*%t(w)%*%solve(d_hat))
    c_inv.w = c_inv%*%w
    w.c_in.w = t(w)%*%c_inv.w
    U = N*diag(sum(Q)) - N*w.c_in.w + t(c_inv.w)%*%S%*%c_inv.w
    
    for(k in 1:K){
      temp1 = A[[k]]%*%S%*%c_inv.w
      temp2 = B[[k]]%*%U%*%t(B[[k]])
      temp3 = A[[k]]%*%w
      ## Update wk_hat 
      wk_hat[[k]] = temp1%*%t(B[[k]])%*%solve(temp2)
      
      ## Update sigma_hat
      sig_hat[k]=mean(diag(A[[k]]%*%S%*%t(A[[k]]) + temp3%*%U%*%t(temp3) - 2*temp3%*%t(temp1)))/N
      
      rm(temp1, temp2, temp3)
    }
    
    w_hat=wk_to_w(wk_hat, P, Q)
    chord.dist = c(chord.dist, chord.norm.diff(w, w_hat))
    
    ## Update d_hat 
    d_hat=generate_d(sig_hat,P)
    
    ################## End of EM-ALGORITHM ######################
    
    # Compute subject scores
    exp.theta = U = Y%*%solve(d_hat)%*%w%*%c_solv
    
    all_obs.LogLik=append(all_obs.LogLik, obs_LogLik(Y, mu_hat, w_hat, d_hat))  
    all_complete.LogLik=append(all_complete.LogLik, complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat))  
    
    len = length(all_obs.LogLik)
    if(all_obs.LogLik[len]<=all_obs.LogLik[len-1] | all_complete.LogLik[len]<=all_complete.LogLik[len-1]){
      d_hat = d; w_hat = w; c_solv=solve(Iq+t(w)%*%solve(d_hat)%*%w)
      exp.theta = U = Y%*%solve(d_hat)%*%w%*%c_solv
      all_obs.LogLik = all_obs.LogLik[-len]
      all_complete.LogLik = all_complete.LogLik[-len]
      flag = TRUE
    }
    ## Update iter
    iter=iter + 1
  }
  
  obs_BIC = (sum(P*(Q[1]+Q[-1]))+2)*log(N)-2*all_obs.LogLik[length(all_obs.LogLik)]
  obs_AIC = (sum(P*(Q[1]+Q[-1]))+2)*2-2*all_obs.LogLik[length(all_obs.LogLik)]
  
  cat(paste0("Total Iterations = ",toString(iter), " \n",
             "Observed Data Likelihood = ",toString(round(all_obs.LogLik[length(all_obs.LogLik)],4)), " \n", 
             "Complete Data Likelihood = ",toString(round(all_complete.LogLik[length(all_complete.LogLik)],4)), " \n",
             "BIC = ",toString(round(obs_BIC,4)), " \n",
             "AIC = ",toString(round(obs_AIC,4)), " \n",
             "Chordal Norm = ",toString(chord.dist[length(chord.dist)]), "\n"))
  
  if(plots){
    if(is.finite(min(all_obs.LogLik))) {
      layout(matrix(1:3, nrow = 1))
      plot(all_obs.LogLik, ylab = "Log-Likelihood", main = "Observed Data Log-Likelihood")
    } else {
      layout(matrix(1:2, nrow = 1))
    }
    if(is.finite(min(all_complete.LogLik))) {plot(all_complete.LogLik, ylab = "Log-Likelihood", 
                                                      main = "Complete Data Log-Likelihood")}
    plot(chord.dist, ylab = "Chordal Norm", main = "Distance between consecutive \n estimates of 'W'")
    layout(1)
  }
  theta.names = paste0("Joint_Score_", 1:Q[1])
  for(k in 1+1:K){theta.names = c(theta.names, paste0("Individual_Data", k-1, "_Score_", 1:Q[k]))}
  colnames(exp.theta) = theta.names
  names(Q) = c("Joint", paste0("Indiv_", 1:K))
  
  J.hat = I.hat = list()
  dat.svd = list()
  
  for(k in 1:K){
    J.hat[[k]] = exp.theta[,1:Q[1]]%*%t(wk_hat[[k]][,1:Q[1]])
    I.hat[[k]] = exp.theta%*%t(B[[k]])[,-(1:Q[1])]%*%t(wk_hat[[k]][,-(1:Q[1])])
    dat.svd[[k]] = svd(Y%*%t(A[[k]]))
  }
  
  signal_matrices = list(J.hat, I.hat)
  names(signal_matrices) = c("Joint", "Individual")
  tot.var = sapply(dat.svd, function(x) sum(x$d^2))
  Joint.signal.var = sapply(signal_matrices[["Joint"]], function(x) sum(svd(x)$d^2))
  Individual.signal.var = sapply(signal_matrices[["Individual"]], function(x) sum(svd(x)$d^2))
  
  VarEx = list()
  for(k in 1:K){temp = c(Joint.signal.var[k]/tot.var[k], Individual.signal.var[k]/tot.var[k])
                names(temp) = c("Joint", "Individual")
                VarEx[[k]] = temp}
  names(VarEx) = paste0("Data_Block", 1:K)
  
  out = list(exp.theta, w_hat, Q, VarEx, sig_hat, chord.dist, all_complete.LogLik, all_obs.LogLik,obs_BIC,obs_AIC)
  names(out) = c("SubjectScoreMatrix", "LoadingMatrix", "Ranks", "VarianceExplained", "ErrorVariances",
                 "ChordalDistances","Complete-Data-Log-Likelihood", "Observed-Data-Log-Likelihood", "BIC", "AIC")
  
  return(out)
}

#############################################################################
##########      Constructs descriptive statistics table    ##################
#############################################################################
brain.descrip = function(dat,x,out.row.names = NULL){
  dat1 = dat[,-which(colnames(dat)==x)]
  dat.mean1 = apply(dat1, 2, function(y) aggregate(y~Covars_AllCog[,x],
                                                   FUN = function(x) round(mean(x),1)))
  dat.mean = t(sapply(dat.mean1, function(mean.list) mean.list$y))
  dat.mean = cbind(dat.mean, round(colMeans(dat1),1))
  rownames(dat.mean) = names(dat.mean1)
  
  dat.sd1 = apply(dat1, 2, function(y) aggregate(y~Covars_AllCog[,x],
                                                 FUN = function(x) round(sd(x), 2)))
  dat.sd = t(sapply(dat.sd1, function(mean.list) mean.list$y))
  rownames(dat.sd) = names(dat.sd1)
  dat.sd = cbind(dat.sd, apply(dat1, 2, function(y) round(sd(y),2)))
  
  dat.f.stats = apply(dat1, 2, function(y) summary(lm(y~Covars_AllCog[,x]))[["fstatistic"]])
  dat.pvals = round(apply(dat.f.stats, 2, function(f) 1-pf(f[1],f[2],f[3])),3)
  
  out.tab = matrix(paste0(dat.mean, " (", dat.sd, ")"), ncol = 4)
  out.tab = cbind(out.tab, dat.pvals)
  
  rownames(out.tab) = rownames(dat.mean)
  colnames(out.tab) = c(dat.mean1[[1]][,1], "Total", "p value")
  return(out.tab)
}


#############################################################################
###########   modification of rain cloud plots function   ##################
#############################################################################
raincloundplots_rjm <- function(data, x, y, cols.in = NULL, lab.inputs){
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  x = as.factor(data[,x]); 
  y = data[,y]
  levs = levels(x); 
  num.levs = length(levels(x))
  cols.in = ifelse(is.null(cols.in), okabe[1:num.levs])
  out.plot = ggplot(data = data) 
  for(i in 1:num.levs){
    dat.temp = data[x == levs[i],]
    out.plot = out.plot +
      geom_jitter(data = dat.temp, aes(x = x, y = y), color = cols.in[i], alpha = 0.6, 
                  show.legend = FALSEALSE, position = position_nudge(x = 0.2), size = 1) + 
      raincloudplots:::geom_flat_violin(data = dat.temp, aes(x = x, y = y), fill = cols.in[i],
                                        color = cols.in[i], position = position_nudge(x = 0.35), alpha = 0.6, show.legend = FALSEALSE) + 
      gghalves::geom_half_boxplot(data = dat.temp, aes(x = x, y = y), fill = cols.in[i], color = cols.in[i], 
                                  position = position_nudge(x = 0.2), side = "r", outlier.shape = NA, center = TRUE, errorbar.draw = FALSEALSE, 
                                  width = 0.2, alpha = 0.6, show.legend = FALSEALSE) + 
      labs(x = lab.inputs[1], y = lab.inputs[2], subtitle = lab.inputs[3])
  }
  out.plot
}

#############################################################################
###########   construct dataset to make the var ex ggplots   ###############
#############################################################################
MakeVarEx.data.gg <- function(AllSims.rows, p1, p2, n){
  
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)*", n="*.(n)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)*", n="*.(n)))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)*", n="*.(n)), bquote("R"[J2]^2*"=0.5, p"[2]*"="*.(p2)*", n="*.(n)))
  
  AllSims.rows$JntVarEx1.fac =  factor(AllSims.rows$JntVarEx1, labels = JVE1.labs, levels = c(0.05, 0.5))
  AllSims.rows$JntVarEx2.fac =  factor(AllSims.rows$JntVarEx2, labels = JVE2.labs, levels = c(0.05, 0.5))
  
  JntVar.Table_p220_rJ1.SD = aggregate(True_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(JntVar.Table_p220_rJ1.SD)[3:4] = c("Mean_Jnt_Var_X_1","SD_Jnt_Var_X_1")
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(JntVar.Table_p220_rJ1.SD)[5:6] = c("Mean_Jnt_Var_X_2","SD_Jnt_Var_X_2")
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(JntVar.Table_p220_rJ1.SD)[7:8] = c("Mean_Ind_Var_X_1","SD_Ind_Var_X_1")
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  JntVar.Table_p220_rJ1.SD = cbind(JntVar.Table_p220_rJ1.SD, aggregate(True_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(JntVar.Table_p220_rJ1.SD)[9:10] = c("Mean_Ind_Var_X_2","SD_Ind_Var_X_2")
  JntVar.Table_p220_rJ1.SD[,-(1:2)] = round(JntVar.Table_p220_rJ1.SD[,-(1:2)],3)
  
  ## AJIVE Variance Explained
  AJIVE_JntVar.Table_p220_rJ1.SD = aggregate(AJIVE_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(AJIVE_JntVar.Table_p220_rJ1.SD)[3:4] = c("Mean_Jnt_Var_X_1","SD_Jnt_Var_X_1")
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(AJIVE_JntVar.Table_p220_rJ1.SD)[5:6] = c("Mean_Jnt_Var_X_2","SD_Jnt_Var_X_2")
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(AJIVE_JntVar.Table_p220_rJ1.SD)[7:8] = c("Mean_Ind_Var_X_1","SD_Ind_Var_X_1")
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  AJIVE_JntVar.Table_p220_rJ1.SD = cbind(AJIVE_JntVar.Table_p220_rJ1.SD, aggregate(AJIVE_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(AJIVE_JntVar.Table_p220_rJ1.SD)[9:10] = c("Mean_Ind_Var_X_2","SD_Ind_Var_X_2")
  
  AJIVE_JntVar.Table_p220_rJ1.SD[,-(1:2)]  = round(AJIVE_JntVar.Table_p220_rJ1.SD[,-(1:2)],3)
  
  ## ProJIVE-Oracle Variance Explained
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = aggregate(ProJIVE_Oracle_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD)[3:4] = c("Mean_Jnt_Var_X_1","SD_Jnt_Var_X_1")
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD)[5:6] = c("Mean_Jnt_Var_X_2","SD_Jnt_Var_X_2")
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD)[7:8] = c("Mean_Ind_Var_X_1","SD_Ind_Var_X_1")
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_Oracle_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD)[9:10] = c("Mean_Ind_Var_X_2","SD_Ind_Var_X_2")
  
  ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,-(1:2)] = round(ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,-(1:2)],3)
  
  ## ProJIVE Variance Explained
  ProJIVE_JntVar.Table_p220_rJ1.SD = aggregate(ProJIVE_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_JntVar.Table_p220_rJ1.SD)[3:4] = c("Mean_Jnt_Var_X_1","SD_Jnt_Var_X_1")
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_JntVar.Table_p220_rJ1.SD)[5:6] = c("Mean_Jnt_Var_X_2","SD_Jnt_Var_X_2")
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_JntVar.Table_p220_rJ1.SD)[7:8] = c("Mean_Ind_Var_X_1","SD_Ind_Var_X_1")
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  ProJIVE_JntVar.Table_p220_rJ1.SD = cbind(ProJIVE_JntVar.Table_p220_rJ1.SD, aggregate(ProJIVE_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(ProJIVE_JntVar.Table_p220_rJ1.SD)[9:10] = c("Mean_Ind_Var_X_2","SD_Ind_Var_X_2")
  
  ProJIVE_JntVar.Table_p220_rJ1.SD[,-(1:2)] = round(ProJIVE_JntVar.Table_p220_rJ1.SD[,-(1:2)],3)
  
  ## RJIVE Variance Explained
  RJIVE_JntVar.Table_p220_rJ1.SD = aggregate(RJIVE_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_JointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(RJIVE_JntVar.Table_p220_rJ1.SD)[3:4] = c("Mean_Jnt_Var_X_1","SD_Jnt_Var_X_1")
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_JointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(RJIVE_JntVar.Table_p220_rJ1.SD)[5:6] = c("Mean_Jnt_Var_X_2","SD_Jnt_Var_X_2")
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_IndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(RJIVE_JntVar.Table_p220_rJ1.SD)[7:8] = c("Mean_Ind_Var_X_1","SD_Ind_Var_X_1")
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  RJIVE_JntVar.Table_p220_rJ1.SD = cbind(RJIVE_JntVar.Table_p220_rJ1.SD, aggregate(RJIVE_IndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(RJIVE_JntVar.Table_p220_rJ1.SD)[9:10] = c("Mean_Ind_Var_X_2","SD_Ind_Var_X_2")
  
  RJIVE_JntVar.Table_p220_rJ1.SD[,-(1:2)] = round(RJIVE_JntVar.Table_p220_rJ1.SD[,-(1:2)],3)
  
  ## dCCA Variance Explained
  dCCA_JntVar.Table_p220_rJ1.SD = aggregate(dCCAJointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAJointVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(dCCA_JntVar.Table_p220_rJ1.SD)[3:4] = c("Mean_Jnt_Var_X_1","SD_Jnt_Var_X_1")
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAJointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAJointVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(dCCA_JntVar.Table_p220_rJ1.SD)[5:6] = c("Mean_Jnt_Var_X_2","SD_Jnt_Var_X_2")
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAIndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAIndivVarExp_X ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(dCCA_JntVar.Table_p220_rJ1.SD)[7:8] = c("Mean_Ind_Var_X_1","SD_Ind_Var_X_1")
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAIndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = mean)[3])
  dCCA_JntVar.Table_p220_rJ1.SD = cbind(dCCA_JntVar.Table_p220_rJ1.SD, aggregate(dCCAIndivVarExp_Y ~ JntVarEx1.fac + JntVarEx2.fac, data = AllSims.rows, FUN = sd)[3])
  colnames(dCCA_JntVar.Table_p220_rJ1.SD)[9:10] = c("Mean_Ind_Var_X_2","SD_Ind_Var_X_2")
  
  dCCA_JntVar.Table_p220_rJ1.SD[,-(1:2)] = round(dCCA_JntVar.Table_p220_rJ1.SD[,-(1:2)],3)
  
  VarEx.dat.gg = data.frame(JntVarEx1.fac = JntVar.Table_p220_rJ1.SD[,"JntVarEx1.fac"], JntVarEx2.fac = JntVar.Table_p220_rJ1.SD[,"JntVarEx2.fac"],
                            Mean_EmpJntVarEx = c(JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_1"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_1"], 
                                                 ProJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_1"], AJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_1"],
                                                 RJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_1"], dCCA_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_1"], 
                                                 
                                                 JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_2"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_2"], 
                                                 ProJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_2"], AJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_2"],
                                                 RJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_2"], dCCA_JntVar.Table_p220_rJ1.SD[,"Mean_Jnt_Var_X_2"], 
                                                 
                                                 JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_1"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_1"], 
                                                 ProJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_1"], AJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_1"], 
                                                 RJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_1"], dCCA_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_1"], 
                                                 
                                                 JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_2"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_2"], 
                                                 ProJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_2"], AJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_2"],
                                                 RJIVE_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_2"], dCCA_JntVar.Table_p220_rJ1.SD[,"Mean_Ind_Var_X_2"]),
                            
                            SD_EmpJntVarEx = c(JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_1"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_1"], 
                                               ProJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_1"], AJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_1"],
                                               RJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_1"], dCCA_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_1"], 
                                               
                                               JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_2"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_2"], 
                                               ProJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_2"], AJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_2"],
                                               RJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_2"], dCCA_JntVar.Table_p220_rJ1.SD[,"SD_Jnt_Var_X_2"], 
                                               
                                               JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_1"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_1"], 
                                               ProJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_1"], AJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_1"],
                                               RJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_1"], dCCA_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_1"], 
                                               
                                               JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_2"], ProJIVE_Oracle_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_2"], 
                                               ProJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_2"], AJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_2"],
                                               RJIVE_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_2"], dCCA_JntVar.Table_p220_rJ1.SD[,"SD_Ind_Var_X_2"]),
                            Type = factor(rep(1:4, each = 24),labels = c("Joint_X_1", "Joint_X_2","Indiv_X_1", "Indiv_X_2")),
                            Method = factor(rep(c(1:6), each = 4), labels = c("True","ProJIVE-Oracle", "ProJIVE","AJIVE","R.JIVE", "dCCA")))

  VarEx.dat.gg
}


###############################################################################################################################
#####################         Plot ProJIVE Proportions of Variance Explained from Simulation Study      #######################
###############################################################################################################################
gg.varex.plot<-function(varex.dat, cols, show.legend = FALSE, text.size = 11, lty = 1, y.max = 1){
  varex.dat = varex.dat[str_which(as.character(varex.dat$Type), "TotalVarExp", T),]
  labs = levels(varex.dat$Type)
  x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
             expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))
  
  ggplot(data = varex.dat, aes(x = Type, y = VarEx)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = FALSE) +
    labs(y = "Variance Explained", x = "Source") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}


#################################################################################################
#####################         Plot PJIVE PMSEs from Simulation Study      #######################
#################################################################################################
gg.pmse.plot<-function(pmse.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(pmse.dat$Type)[c(3,6,7,1,2,4,5)]
  labs.ex = c("Joint Subj Scores", expression("Joint Loadings"*"X"[1]), expression("Joint Loadings"*"X"[2]), 
              expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]),
              expression("Indiv Loadings"*"X"[1]), expression("Indiv Loadings"*"X"[2]))
  ggplot(data = pmse.dat, aes(x = Type, y = PMSE)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.25, lwd = 0.25) +
    labs(y = "PMSE", x = "Type") +
    facet_grid(JVE_1 ~ JVE_2, labeller = label_parsed) +
    scale_x_discrete(limits = levels(pmse.dat$Type)[c(3,6,7,1,2,4,5)], labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 35, size = text.size-2),
          text = element_text(size = text.size), legend.position = "bottom", legend.direction = "horizontal")
}

###################################################################################################
#####################         Plot Chordal Norms from Simulation Study      #######################
###################################################################################################
gg.norm.plot.2<-function(norm.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)]
  labs.ex = c("Joint Subj Scores", expression("Joint Loadings"*"X"[1]), expression("Joint Loadings"*"X"[2]), 
              expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]),
              expression("Indiv Loadings"*"X"[1]), expression("Indiv Loadings"*"X"[2]))
  ggplot(data = norm.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.25, lwd = 0.25) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_1 ~ JVE_2, labeller = label_parsed) +
    scale_x_discrete(limits = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)], labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 35, size = text.size-2),
          text = element_text(size = text.size) ,legend.position = "bottom", legend.direction = "horizontal")
}

#################################################################################################
#####################         Plot CJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.score.norm.plot<-function(norm.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[grep("Score",levels(norm.dat$Type))][c(3,1,2)]
  labs.ex = c("Joint Subj Scores", expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]))
  ggplot(data = norm.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs, labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}

#################################################################################################
#####################         Plot CJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.load.norm.plot<-function(norm.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[grep("Load",levels(norm.dat$Type))][c(3,4,1,2)]
  labs.ex = c(expression("Joint Variable Loadings"*"X"[1]),expression("Joint Variable Loadings"*"X"[2]),
              expression("Indiv Variable Loadings"*"X"[1]), expression("Indiv Variable Loadings"*"X"[2]))
  ggplot(data = norm.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs, labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}


######################################################################################################
#####################       Adaptation of PMSE function from singR Package     #######################
######################################################################################################
pmse.2<-function (M1 = NULL, M2 = NULL, S1 = NULL, S2 = NULL, standardize = FALSE){
  tfun = function(x) all(x == 0)
  if (is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2)) 
    stop("need to supply either M1 and M2 or S1 and S2")
  if (!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if (!is.null(M1) && nrow(M1) > ncol(M1)){ 
    stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")
  }
  if (is.null(M1)) {
    if(is.vector(S1)) {S1 = matrix(S1, ncol = 1)}
    if(is.vector(S2)) {S2 = matrix(S2, ncol = 1)}
    
    nS = nrow(S1)
    if (nS != nrow(S2)) 
      stop("S1 and S2 must have the same number of rows")
    if (sum(apply(S1, 2, tfun)) + sum(apply(S2, 2, tfun))) 
      stop("pmse not defined when S1 or S2 has a column of all zeros")
    if (standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if (p < q) {
      S1 = cbind(S1, matrix(0, nS, (q - p)))
    }
    if (q < p) {
      S2 = cbind(S2, matrix(0, nS, (p - q)))
    }
    Stemp = matchICA.2(S = S1, template = S2)
    n.comp = max(q, p)
    indices = c(1:n.comp)[!(apply(Stemp, 2, tfun) | apply(S2, 
                                                          2, tfun))]
    return(sqrt(sum((Stemp[, indices] - S2[, indices])^2))/sqrt(nS * 
                                                                  min(p, q)))
  }
  else {
    if (sum(apply(M1, 1, tfun)) + sum(apply(M2, 1, tfun))) 
      stop("pmse not defined when M1 or M2 has a row of all zeros")
    if (standardize) {
      temp = diag((diag(M1 %*% t(M1)))^(-1/2))
      M1 = temp %*% M1
      temp = diag((diag(M2 %*% t(M2)))^(-1/2))
      M2 = temp %*% M2
    }
    p = ncol(M1)
    if (p != ncol(M2)) 
      stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp = max(d, q)
    if (n.comp > p) 
      warning("M should be d x p")
    if (d < q) {
      M1 = rbind(M1, matrix(0, (q - d), p))
    }
    if (q < d) {
      M2 = rbind(M2, matrix(0, (d - q), p))
    }
    l2.mat1 = l2.mat2 = matrix(NA, nrow = n.comp, ncol = n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        l2.mat1[i, j] = sum((M2[i, ] - M1[j, ])^2)
        l2.mat2[i, j] = sum((M2[i, ] + M1[j, ])^2)
      }
    }
    l2.mat1 = sqrt(l2.mat1)
    l2.mat2 = sqrt(l2.mat2)
    l2.mat = l2.mat1 * (l2.mat1 <= l2.mat2) + l2.mat2 * (l2.mat2 < 
                                                           l2.mat1)
    map = as.vector(solve_LSAP(l2.mat))
    l2.1 = ifelse(nrow(l2.mat1)==1 && ncol(l2.mat1)==1 && l2.mat1[1,1] == 0, 0, diag(l2.mat1[, map]))
    l2.2 = ifelse(nrow(l2.mat2)==1 && ncol(l2.mat2)==1 && l2.mat2[1,1] == 0, 0, diag(l2.mat2[, map]))
    sign.change = -1 * (l2.2 < l2.1) + 1 * (l2.1 <= l2.2)
    perm = diag(n.comp)[, map] %*% diag(sign.change)
    M.perm = t(perm) %*% M1
    indices = c(1:n.comp)[!(apply(M.perm, 1, tfun) | apply(M2, 
                                                           1, tfun))]
    return(sqrt(sum((M.perm[indices, ] - M2[indices, ])^2))/sqrt(p * 
                                                                   min(d, q)))
  }
}

matchICA.2<-function (S, template, M = NULL) {
  n.comp = ncol(S)
  n.obs = nrow(S)
  if (n.comp > n.obs) 
    warning("Input should be n x d")
  if (!all(dim(template) == dim(S))) 
    warning("Template should be n x d")
  S = t(S)
  template = t(template)
  l2.mat1 = matrix(NA, nrow = n.comp, ncol = n.comp)
  l2.mat2 = l2.mat1
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      l2.mat1[i, j] = sum((template[i, ] - S[j, ])^2)/n.obs
      l2.mat2[i, j] = sum((template[i, ] + S[j, ])^2)/n.obs
    }
  }
  l2.mat1 = sqrt(l2.mat1)
  l2.mat2 = sqrt(l2.mat2)
  l2.mat = l2.mat1 * (l2.mat1 <= l2.mat2) + l2.mat2 * (l2.mat2 < 
                                                         l2.mat1)
  map = as.vector(clue::solve_LSAP(l2.mat))
  if(nrow(l2.mat1)==1 && ncol(l2.mat1)==1 && l2.mat1[1,1] == 0){
    l2.1 = 0
  } else if (nrow(l2.mat1)==1 && ncol(l2.mat1)==1){
    l2.1 = as.matrix(l2.mat1[1,1])
  } else {
  l2.1 = diag(l2.mat1[, map])
  }
  if(nrow(l2.mat2)==1 && ncol(l2.mat2)==1 && l2.mat2[1,1] == 0){
    l2.2 = 0
  } else if (nrow(l2.mat2)==1 && ncol(l2.mat2)==1){
    l2.2 = as.matrix(l2.mat2[1,1])
  } else {
    l2.2 = diag(l2.mat2[, map])
  }
  sign.change = -1 * (l2.2 < l2.1) + 1 * (l2.1 <= l2.2)
  perm = diag(n.comp)[, map] %*% diag(sign.change)
  s.perm = t(perm) %*% S
  if (!is.null(M)) {
    M.perm = t(M) %*% perm
    return(list(S = t(s.perm), M = t(M.perm)))
  }
  else {
    t(s.perm)
  }
}
# ######################################################################################################################
# ###########   S-POET: Spiked covariance estimation (Wang and Fan 2017)                            ####################
# ######################################################################################################################
# sPOET <- function(Y, K = NULL, K_max = NULL, K_min = NULL, method = NULL){
#   n = ncol(Y); p = nrow(Y)
#   Y.svd = svd(Y); U = Y.svd[['u']]; S = diag(Y.svd[['d']]); V_t = t(Y.svd[['v']])
# 
#   Lambda = S^2 #eigenvalues of the sample covariance matrix
# 
#   if(is.null(K) & !is.null(method)){
#     K = ifelse(method=="GR", select_K_by_GR(diag(Lambda), K_max, K_min),
#                select_K_by_ED(diag(Lambda), K_max, K_min))
#   }
# 
#   c_hat = sum(Lambda[-(1:K)])/(p - K - p*K/n)
# 
#   Lambda_S = Lambda[1:K,1:K]
#   for(k in 1:K){Lambda_S[k,k] = max(Lambda[k,k] - c_hat*p/n,0)}
# 
#   X_hat = U[,1:K]%*%sqrt(Lambda_S)%*%V_t[1:K,]
# 
#   return(list(X_hat, Lambda_S, U[,1:K], K, V_t[1:K]))
# }
# 
# 
# ######################################################################################################################
# ###########                                         d-CCA                                         ####################
# ######################################################################################################################
# dCCA <- function(Y_1, Y_2, r_1=NULL, r_2 = NULL, r_12 = NULL, method = NULL){
#   n = ncol(Y_1); p_1 = nrow(Y_1); p_2 = nrow(Y_2)
# 
#   s.out1 = sPOET(Y = Y_1, K = r_1, method = method)
#   s.out2 = sPOET(Y = Y_2, K = r_2, method = method)
# 
#   X1_hat = s.out1[[1]]; Lambda_1 = s.out1[[2]]; V_1 = s.out1[[3]]; r_1 = s.out1[[4]]; V_y1_t = s.out1[[5]]
#   X2_hat = s.out2[[1]]; Lambda_2 = s.out2[[2]]; V_2 = s.out2[[3]]; r_2 = s.out2[[4]]; V_y2_t = s.out2[[5]]
# 
#   Lambda1_inv_half = Lambda_1*0
#   index_nonzero_diag_lambda1 = which(diag(Lambda_1)>0)
#   Lambda1_inv_half[index_nonzero_diag_lambda1, index_nonzero_diag_lambda1] = diag(Lambda_1)[index_nonzero_diag_lambda1]^(-0.5)
# 
#   Lambda2_inv_half = Lambda_2*0
#   index_nonzero_diag_lambda2 = which(diag(Lambda_2)>0)
#   Lambda2_inv_half[index_nonzero_diag_lambda2, index_nonzero_diag_lambda2] = diag(Lambda_2)[index_nonzero_diag_lambda2]^(-0.5)
# 
#   Theta = Lambda1_inv_half%*%t(V_1)%*%X1_hat%*%t(X2_hat)%*%V_2%*%Lambda2_inv_half/n
# 
#   svd.theta = svd(Theta)
#   U_theta = svd.theta[['u']]; D_theta = svd.theta[['d']]; V_theta_t = svd.theta[['v']]
# 
#   Gamma1 = V_1%*%Lambda1_inv_half%*%U_theta
#   Gamma2 = V_2%*%Lambda2_inv_half%*%V_theta_t
# 
#   D_theta_star = D_theta
#   D_theta_star[which(D_theta>1)] = 1
# 
#   r_theta = sum(D_theta_star>0)
# 
#   if(is.null(r_12) | r_12<1){
#     ccor_hat = svd(V_y1_t[1:r_1,]%*%t(V_y1_t[1:r_2,]))[['d']]
#     r_12 = select_r12_by_MDLIC(ccor_hat, r_1, r_2, n)
#   }
# 
#   # if(r_12>r_theta)
#     r_theta = r_12
# 
#   A_mat_C = diag(0.5*(1 - ((1 - D_theta_star[1:r_theta])/(1+D_theta_star[1:r_theta]))^0.5), nrow = r_theta, ncol = r_theta)
# 
#   B_1 = X1_hat%*%t(X1_hat)%*%Gamma1[,1:r_12]/n
#   B_2 = X2_hat%*%t(X2_hat)%*%Gamma2[,1:r_12]/n
#   C_base = A_mat_C[1:r_12, 1:r_12]%*%(t(Gamma1[,1:r_12])%*%X1_hat + t(Gamma2[,1:r_12])%*%X2_hat)
# 
#   C1_hat = B_1%*%C_base
#   C2_hat = B_2%*%C_base
# 
#   B_1_rtheta = X1_hat%*%t(X1_hat)%*%Gamma1[,1:r_theta]/n
#   B_2_rtheta = X2_hat%*%t(X2_hat)%*%Gamma2[,1:r_theta]/n
#   C_base_rtheta = A_mat_C[1:r_theta, 1:r_theta]%*%(t(Gamma1[,1:r_theta])%*%X1_hat + t(Gamma2[,1:r_theta])%*%X2_hat)
# 
#   D1_hat = X1_hat - B_1%*%C_base_rtheta
#   D2_hat = X2_hat - B_2%*%C_base_rtheta
# 
#   X1_hat = D1_hat + C1_hat
#   X2_hat = D2_hat + C2_hat
# 
#   out = list(X1_hat, X2_hat, C1_hat, C2_hat, D1_hat, D2_hat, r_1, r_2, r_12, D_theta_star[1:r_12], acos(D_theta_star[1:r_12])/(pi*180))
#   names(out) = c("SignalMatrix_X1","SignalMatrix_X2","Common_SignalMatrix_X1","Common_SignalMatrix_X2","Distinct_SignalMatrix_X1","Distinct_SignalMatrix_X2",
#                  "SignalRank_X1", "SignalRank_X2", "CommonSignalRank", "CanonicalCorrelations", "PrincipalAngles")
#   return(out)
# }
# 
# 
# 
# ######################################################################################################
# ###########          pJIVE ML estimation of JIVE model that uses EM algorithm       ##################
# ###########   This version should be archived after confirmation that it matches    ##################
# ###########     with the newer version above, which was built for K>2               ##################
# ######################################################################################################
# ProJIVE_EM_Kequal2=function(Y,P,Q,Max.iter=10000,diff.tol=1e-5,plots=TRUE,chord.tol=-1,sig_hat=NULL, init.loads = NULL){
#   
#   ## init.loads must be a list of two lists - first item contains a list of joint loading matrices
#   ##                                          second item is a list of indiv loading matrices
#   ##                                          each matrix must have dimension p_k-by-r_
#   
#   # Total sample size
#   N=dim(Y)[1]
#   
#   # Total number feature blocks
#   if(length(P)==(length(Q)-1)){
#     K=length(P)
#   } else{
#     stop("Error: The number of feature blocks do not match the number of invididual scores")
#   }
#   
#   # Selection matrices A_knd B_k
#   A=list()
#   B=list()
#   for(k in 1:K){
#     up=lapply(Q,generate_ab,Q[1])
#     up[[1]]=diag(Q[1])
#     down=lapply(Q,generate_ab,Q[(k+1)])
#     down[[(k+1)]]=diag(Q[(k+1)])
#     B[[k]]=rbind(do.call("cbind",up),
#                  do.call("cbind",down))
#     
#     a=lapply(P,generate_ab,P[k])
#     a[[k]]=diag(P[k])
#     A[[k]]=do.call("cbind",a)
#   }
#   
#   # get initial estimates of loadings matrices via cc.jive
#   if (is.null(init.loads)){
#     WJ = list(matrix(rnorm(Q[1]*P[1]), nrow = P[1]),matrix(rnorm(Q[1]*P[2]), nrow = P[2]))
#     WI = list(matrix(rnorm(Q[2]*P[1]), nrow = P[1]),matrix(rnorm(Q[3]*P[2]), nrow = P[2]))
#   } else if (is.list(init.loads)){
#     WJ = init.loads[[1]]
#     WI = init.loads[[2]]
#   } else if(init.loads == "AJIVE"){
#     dat.blocks = list()
#     dat.blocks[[1]] = Y[,1:P[1]]
#     for(k in 2:K){
#       dat.blocks[[k]] = Y[,cumsum(P[k-1])+(1:P[k])]
#     }
#     ajive.solution = ajive(dat.blocks, initial_signal_ranks = Q[1]+Q[-1], joint_rank = Q[1])
#     
#     WJ = lapply(ajive.solution$block_decomps, function(x) x$joint$v)
#     WI = lapply(ajive.solution$block_decomps, function(x) x$individual$v)
#   }
#   
#   
#   # Block specific loading matrices W_k
#   wk_hat=wji_hat=list()
#   for(k in 1:K){ 
#     wji_hat[[k]]=cbind(WJ[[k]], WI[[k]])
#   }
#   wk_hat=wji_hat
#   
#   # Total loading matrices W
#   w_hat=wk_to_w(wk_hat,P,Q)
#   
#   # Initialization of sig_hat
#   # sig_hat=list(rep(1,K))                       
#   if(is.null(sig_hat)){
#     sig_hat = rnorm(length(P))
#   } else if (is.character(sig_hat) & sig_hat[1] == "MLE"){
#     sig_hat=mean(svd(Y[,1:P[1]])$d[-(1:sum(Q[1:2]))]^2)
#     for(k in 2:K){sig_hat = c(sig_hat, mean(svd(Y[,sum(P[1:(k-1)])+(1:P[k])])$d[-(1:sum(Q[c(1,k+1)]))]^2))}
#   }
#   
#   # Total noise matrices D
#   d_hat=generate_d(sig_hat,P)
#   
#   # Initiate LogLik
#   all_obs.LogLik=c(-Inf)
#   all_complete.LogLik = c(-Inf)
#   
#   # Initiate Chordal Distance
#   chord.dist = 1
#   mu_hat=apply(as.matrix(Y),2,sum) / N   
#   
#   ## Store some values to save computation time
#   Yc=sweep(Y, 2, mu_hat)
#   S=t(Yc)%*%Yc
#   
#   Iq=diag(sum(Q))
#   Ip=diag(sum(P))
#   
#   # Set initial iteration number:
#   iter=0
#   
#   while (Max.iter>=iter & eval_converge(all_obs.LogLik,all_obs.LogLik,N*diff.tol) & (chord.dist[iter+1] >= chord.tol)) 
#   {
#     ################## START OF EM-ALGORITHM ######################
#     w=w_hat
#     
#     c_solv=solve(Iq+t(w)%*%solve(d_hat)%*%w)
#     c_inv = solve(d_hat)%*%(diag(sum(P))-w%*%c_solv%*%t(w)%*%solve(d_hat))
#     c_inv.w = c_inv%*%w
#     w.c_in.w = t(w)%*%c_inv.w
#     U = N*diag(sum(Q)) - N*w.c_in.w + t(c_inv.w)%*%S%*%c_inv.w
#     
#     for(k in 1:K){
#       temp1 = A[[k]]%*%S%*%c_inv.w
#       temp2 = B[[k]]%*%U%*%t(B[[k]])
#       temp3 = A[[k]]%*%w
#       ## Update wk_hat 
#       wk_hat[[k]] = temp1%*%t(B[[k]])%*%solve(temp2)
#       
#       ## Update sigma_hat
#       sig_hat[k]=mean(diag(A[[k]]%*%S%*%t(A[[k]]) + temp3%*%U%*%t(temp3) - 2*temp3%*%t(temp1)))/N
#       
#       rm(temp1, temp2, temp3)
#     }
#     
#     w_hat=wk_to_w(wk_hat, P, Q)
#     chord.dist = c(chord.dist, chord.norm.diff(w, w_hat))
#     
#     ## Update d_hat 
#     d_hat=generate_d(sig_hat,P)
#     
#     
#     
#     ################## End of EM-ALGORITHM ######################
#     
#     ## Update iter
#     iter=iter + 1
#     
#     # Compute subject scores
#     exp.theta = U = Y%*%solve(d_hat)%*%w%*%c_solv
#     
#     all_obs.LogLik=append(all_obs.LogLik,obs_LogLik(Y, mu_hat, w_hat, d_hat))  
#     all_complete.LogLik=append(all_complete.LogLik, complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat))  
#   }
#   
#   obs_BIC = (sum(P*(Q[1]+Q[-1]))+2)*log(N)-2*all_obs.LogLik[length(all_obs.LogLik)]
#   obs_AIC = (sum(P*(Q[1]+Q[-1]))+2)*2-2*all_obs.LogLik[length(all_obs.LogLik)]
#   
#   cat("Total Iterations = ",toString(iter), "\n",
#       "Observed Data Likelihood = ",toString(round(all_obs.LogLik[length(all_obs.LogLik)],4)), "\n",
#       "Complete Data Likelihood = ",toString(round(all_complete.LogLik[length(all_complete.LogLik)],4)), "\n",
#       "BIC = ",toString(round(obs_BIC,4)), "\n",
#       "AIC = ",toString(round(obs_AIC,4)), "\n",
#       "Chordal Norm = ",toString(chord.dist[length(chord.dist)]))
#   
#   if(plots){
#     if(is.finite(min(all_obs.LogLik[-1]))) {
#       layout(matrix(1:3, nrow = 1))
#       plot(all_obs.LogLik[-1], ylab = "Log-Likelihood", main = "Observed Data Log-Likelihood")
#     } else {
#       layout(matrix(1:2, nrow = 1))
#     }
#     if(is.finite(min(all_complete.LogLik[-1]))) {plot(all_complete.LogLik[-1], ylab = "Log-Likelihood", 
#                                                       main = "Complete Data Log-Likelihood")}
#     plot(chord.dist, ylab = "Chordal Norm", main = "Distance between consecutive \n estimates of 'W'")
#     layout(1)
#   }
#   theta.names = paste0("Joint_Score_", 1:Q[1])
#   for(k in 1+1:K){theta.names = c(theta.names, paste0("Individual_Data", k-1, "_Score_", 1:Q[k]))}
#   colnames(exp.theta) = theta.names
#   names(Q) = c("Joint", "Indiv_1", "Indiv_2")
#   
#   J1.hat = exp.theta[,1:Q[1]]%*%t(w_hat[1:P[1],1:Q[1]])
#   J2.hat = exp.theta[,1:Q[1]]%*%t(w_hat[-(1:P[1]),1:Q[1]])
#   
#   I1.hat = exp.theta[,(Q[1]+1):(Q[1]+Q[2])]%*%t(w_hat[1:P[1],(Q[1]+1):(Q[1]+Q[2])])
#   I2.hat = exp.theta[,(Q[1]+Q[2]):(Q[1]+Q[2]+Q[3])]%*%t(w_hat[-(1:P[1]),(Q[1]+Q[2]):(Q[1]+Q[2]+Q[3])])
#   
#   signal_matrices = list(J1.hat, J2.hat, I1.hat, I2.hat)
#   dat.svd = list(svd(Y[,1:P[1]]),svd(Y[,-(1:P[1])]))
#   tot.var = sapply(dat.svd, function(x) sum(x$d^2))
#   signal.var = sapply(signal_matrices, function(x) sum(svd(x)$d^2))
#   
#   VarEx = list(c(signal.var[1]/tot.var[1],signal.var[2]/tot.var[2]),
#                c(signal.var[3]/tot.var[1],signal.var[4]/tot.var[2]))
#   names(VarEx) = c("Joint", "Indiv")
#   
#   out = list(exp.theta, w_hat, Q, VarEx, sig_hat, chord.dist, all_complete.LogLik[-1], all_obs.LogLik[-1],obs_BIC,obs_AIC)
#   names(out) = c("SubjectScoreMatrix", "LoadingMatrix", "Ranks", "VarianceExplained", "ErrorVariances",
#                  "ChordalDistances","Complete-Data-Log-Likelihood", "Observed-Data-Log-Likelihood", "BIC", "AIC")
#   
#   return(out)
# }
# 
# #############################################################################
# ###########          pJIVE ML estimation of JIVE model     ##################
# ###########          that uses EM algorithm                ##################
# #############################################################################
# ProJIVE_EM_OLD=function(Y,P,Q,Max.iter=10000,diff.tol=1e-5,plots=TRUE,chord.tol=-1,sig_hat=NULL, init.loads = NULL){
#   
#   ## init.loads must be a list of two lists - first item contains a list of joint loading matrices
#   ##                                          second item is a list of indiv loading matrices
#   ##                                          each matrix must have dimension p_k-by-r_
#   
#   # Total sample size
#   N=dim(Y)[1]
#   
#   # Total number feature blocks
#   if(length(P)==(length(Q)-1)){
#     K=length(P)
#   } else{
#     stop("Error: The number of feature blocks do not match the number of invididual scores")
#   }
#   
#   # Selection matrices A_knd B_k
#   A=list()
#   B=list()
#   for(k in 1:K){
#     up=lapply(Q,generate_ab,Q[1])
#     up[[1]]=diag(Q[1])
#     down=lapply(Q,generate_ab,Q[(k+1)])
#     down[[(k+1)]]=diag(Q[(k+1)])
#     B[[k]]=rbind(do.call("cbind",up),
#                  do.call("cbind",down))
#     
#     a=lapply(P,generate_ab,P[k])
#     a[[k]]=diag(P[k])
#     A[[k]]=do.call("cbind",a)
#   }
#   
#   # get initial estimates of loadings matrices via cc.jive
#   if (is.null(init.loads)){
#     WJ = list(matrix(rnorm(Q[1]*P[1]), nrow = P[1]),matrix(rnorm(Q[1]*P[2]), nrow = P[2]))
#     WI = list(matrix(rnorm(Q[2]*P[1]), nrow = P[1]),matrix(rnorm(Q[3]*P[2]), nrow = P[2]))
#   } else if (is.list(init.loads)){
#     WJ = init.loads[[1]]
#     WI = init.loads[[2]]
#   } else if(init.loads == "CJIVE"){
#     dat.blocks=list(Y[,1:P[1]],Y[,P[1]+(1:P[2])])
#     cjive.solution = cc.jive(dat.blocks = dat.blocks, signal.ranks = Q[1]+Q[-1], 
#                              joint.rank = Q[1], perm.test = FALSEALSE)
#     
#     WJ = lapply(cjive.solution$sJIVE$joint_matrices, function(x) x$v)
#     WI =lapply(cjive.solution$sJIVE$indiv_matrices, function(x) x$v)
#   }
#   
#   # Block specific loading matrices W_k
#   wk_hat=wji_hat=list()
#   for(k in 1:K){ # Selecting a sub matrx of the Cholesky Decomp solution L
#     wji_hat[[k]]=cbind(WJ[[k]], WI[[k]])
#   }
#   wk_hat=wji_hat
#   
#   # Total loading matrices W
#   w_hat=wk_to_w(wk_hat,P,Q)
#   
#   # Initializationf of sig_hat
#   # sig_hat=list(rep(1,K))                       
#   if(is.null(sig_hat)){
#     sig_hat = rnorm(length(P))
#   } else if (is.character(sig_hat) & sig_hat[1] == "MLE"){
#     sig_hat=mean(svd(Y[,1:P[1]])$d[-(1:sum(Q[1:2]))]^2)
#     for(k in 2:K){sig_hat = c(sig_hat, mean(svd(Y[,sum(P[1:(k-1)])+(1:P[k])])$d[-(1:sum(Q[c(1,k+1)]))]^2))}
#   }
#   
#   # Total noise matrices D
#   d_hat=generate_d(sig_hat,P)
#   
#   # Initiate LogLik
#   all_obs.LogLik=c(-Inf)
#   all_complete.LogLik = c(-Inf)
#   
#   # Initiate Chordal Distance
#   chord.dist = 1
#   mu_hat=apply(as.matrix(Y),2,sum) / N   
#   
#   ## Store some values to save computation time
#   Yc=sweep(Y, 2, mu_hat)
#   S=t(Yc)%*%Yc
#   
#   Iq=diag(sum(Q))
#   Ip=diag(sum(P))
#   
#   # Set initial iteration number:
#   iter=0
#   
#   while (Max.iter>=iter & eval_converge(all_obs.LogLik,all_obs.LogLik,N*diff.tol) & (chord.dist[iter+1] >= chord.tol)) 
#   {
#     ################## START OF EM-ALGORITHM ######################
#     w=w_hat
#     
#     c_solv=solve(Iq+t(w)%*%solve(d_hat)%*%w)
#     
#     U=S%*%solve(d_hat)%*%w%*%c_solv
#     V=c_solv+c_solv%*%t(w)%*%solve(d_hat)%*%S%*%solve(d_hat)%*%w%*%c_solv
#     
#     ## Update d_tild
#     
#     d_tild=S-2*w%*%t(U)+w%*%V%*%t(w)   
#     for(k in 1:K){
#       
#       ## Update wk_hat 
#       wk_hat[[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))
#       
#       ## Update sigma_hat
#       sig_hat[k]=mean(diag(A[[k]]%*%diag(diag(d_tild))%*%t(A[[k]])))/N
#     }
#     
#     # w=w_hat
#     # 
#     # Cinv = solve(d_hat)%*%(Ip-w%*%solve(Iq+t(w)%*%solve(d_hat)%*%w)%*%t(w)%*%solve(d_hat))
#     # 
#     # # Compute subject scores
#     # exp.theta = U = Y%*%Cinv%*%w
#     # # Compute expected outer product
#     # cov.theta = Iq - t(w)%*%Cinv%*%w
#     # exp.theta2 = V = cov.theta + t(exp.theta)%*%exp.theta
#     # 
#     # ## Update d_tild
#     # 
#     # # d_tild=(S-2*w%*%t(U)%*%Y+w%*%V%*%t(w))/N   
#     # d_tild=(S + 2*t(Y)%*%U%*%t(w) - w%*%V%*%t(w))/N   
#     # for(k in 1:K){
#     #   U = t(Y)%*%exp.theta
#     #   
#     #   ## Update wk_hat 
#     #   wk_hat[[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))
#     #   
#     #   ## Update sigma_hat
#     #   sig_hat[k]=ifelse(is.null(sig_hat),mean(diag(A[[k]]%*%diag(diag(d_tild))%*%t(A[[k]]))),sig_hat[k])
#     # }
#     # 
#     ## Update w_hat 
#     w_hat=wk_to_w(wk_hat, P, Q)
#     chord.dist = c(chord.dist, chord.norm.diff(w, w_hat))
#     
#     ## Update d_hat 
#     d_hat=generate_d(sig_hat,P)
#     
#     
#     ################## End of EM-ALGORITHM ######################
#     
#     ## Update iter
#     iter=iter + 1
#     
#     # Compute subject scores
#     exp.theta = U = Y%*%solve(d_hat)%*%w%*%c_solv
#     
#     all_obs.LogLik=append(all_obs.LogLik,obs_LogLik(Y, mu_hat, w_hat, d_hat))  
#     all_complete.LogLik=append(all_complete.LogLik, complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat))  
#   }
#   
#   obs_BIC = sum(P*(Q[1]+Q[-1]))*log(N)-2*all_obs.LogLik[length(all_obs.LogLik)]
#   obs_AIC = sum(P*(Q[1]+Q[-1]))*2-2*all_obs.LogLik[length(all_obs.LogLik)]
#   
#   cat("Total Iterations = ",toString(iter), "\n",
#       "Observed Data Likelihood = ",toString(round(all_obs.LogLik[length(all_obs.LogLik)],4)), "\n",
#       "Complete Data Likelihood = ",toString(round(all_complete.LogLik[length(all_complete.LogLik)],4)), "\n",
#       "BIC = ",toString(round(obs_BIC,4)), "\n",
#       "AIC = ",toString(round(obs_AIC,4)), "\n",
#       "Chordal Norm = ",toString(chord.dist[length(chord.dist)]))
#   
#   if(plots){
#     if(is.finite(min(all_obs.LogLik[-1]))) {
#       layout(matrix(1:3, nrow = 1))
#       plot(all_obs.LogLik[-1], ylab = "Log-Likelihood", main = "Observed Data Log-Likelihood")
#     } else {
#       layout(matrix(1:2, nrow = 1))
#     }
#     if(is.finite(min(all_complete.LogLik[-1]))) {plot(all_complete.LogLik[-1], ylab = "Log-Likelihood", 
#                                                       main = "Complete Data Log-Likelihood")}
#     plot(chord.dist, ylab = "Chordal Norm", main = "Distance between consecutive \n estimates of 'W'")
#     layout(1)
#   }
#   theta.names = paste0("Joint_Score_", 1:Q[1])
#   for(k in 1+1:K){theta.names = c(theta.names, paste0("Individual_Data", k-1, "_Score_", 1:Q[k]))}
#   colnames(exp.theta) = theta.names
#   names(Q) = c("Joint", "Indiv_1", "Indiv_2")
#   
#   J1.hat = exp.theta[,1:Q[1]]%*%t(w_hat[1:P[1],1:Q[1]])
#   J2.hat = exp.theta[,1:Q[1]]%*%t(w_hat[-(1:P[1]),1:Q[1]])
#   
#   I1.hat = exp.theta[,(Q[1]+1):(Q[1]+Q[2])]%*%t(w_hat[1:P[1],(Q[1]+1):(Q[1]+Q[2])])
#   I2.hat = exp.theta[,(Q[1]+Q[2]):(Q[1]+Q[2]+Q[3])]%*%t(w_hat[-(1:P[1]),(Q[1]+Q[2]):(Q[1]+Q[2]+Q[3])])
#   
#   signal_matrices = list(J1.hat, J2.hat, I1.hat, I2.hat)
#   dat.svd = list(svd(Y[,1:P[1]]),svd(Y[,-(1:P[1])]))
#   tot.var = sapply(dat.svd, function(x) sum(x$d^2))
#   signal.var = sapply(signal_matrices, function(x) sum(svd(x)$d^2))
#   
#   VarEx = list(c(signal.var[1]/tot.var[1],signal.var[2]/tot.var[2]),
#                c(signal.var[3]/tot.var[1],signal.var[4]/tot.var[2]))
#   names(VarEx) = c("Joint", "Indiv")
#   
#   out = list(exp.theta, w_hat, Q, VarEx, sig_hat, chord.dist, all_complete.LogLik[-1], all_obs.LogLik[-1],obs_BIC,obs_AIC)
#   names(out) = c("SubjectScoreMatrix", "LoadingMatrix", "Ranks", "VarianceExplained", "ErrorVariances",
#                  "ChordalDistances","Complete-Data-Log-Likelihood", "Observed-Data-Log-Likelihood", "BIC", "AIC")
#   
#   # out = list(exp.theta, w_hat, Q, sig_hat, chord.dist, all_complete.LogLik[-1], all_obs.LogLik[-1])
#   # names(out) = c("SubjectScoreMatrix", "LoadingMatrix", "Ranks", "ErrorVariances",
#   #                "ChordalDistances","Complete-Data-Log-Likelihood", "Observed-Data-Log-Likelihood")
#   return(out)
# }
# 
