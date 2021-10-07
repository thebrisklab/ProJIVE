source("H:/My Documents/P-JIVE/Programs/Functions/Functions_for_CJIVE.R")
source("H:/My Documents/P-JIVE/Programs/Functions/Functions_for_PJIVE.R")

######################################################################################################################
###########   Generates K Simulated Datasets that follow JIVE Model using binary subject scores   ####################
######################################################################################################################
GenerateToyData_Kge2 <- function(n, p, JntVarEx, IndVarEx, jnt_rank = 1, equal.eig = F, ind_ranks, JntVarAdj = T,
                             SVD.plots = T, Error = T, print.cor = T, Loads = "Rademacher", Scores = "Gaussian_Mixture"){
  
  # Write out both joint and indiv subject scores for both data sets first
  r.J = jnt_rank
  r.I = ind_ranks
  K = length(p)

  print(paste("Generating Scores from", Scores, "distributions and Loadings from", Loads, "distributions."))
  
  # Generate Joint scores
  if(Scores=="Binomial"){
    JntScores = matrix(rbinom(n*r.J, size=1, prob=0.2), nrow = n, ncol = r.J)
    
    b = rbinom(n*sum(r.I), size=1, prob=0.4)
    b = 1 - 2*b
    IndivScores = matrix(b, nrow = n, ncol = sum(r.I))
  } else if (Scores=="Gaussian_Mixture"){
    mix.probs = c(0.2, 0.5, 0.3)
    JointScores.g1 = matrix(rnorm(n*mix.probs[1]*r.J), ncol = r.J)
    JointScores.g2 = matrix(rnorm(n*mix.probs[2]*r.J, -4, 1), ncol = r.J)
    JointScores.g3 = matrix(rnorm(n*mix.probs[3]*r.J, 4, 1), ncol = r.J)
    JntScores = rbind(JointScores.g1, JointScores.g2, JointScores.g3)
    
    IndivScores = matrix(rnorm(n*sum(r.I)), nrow = n, ncol = sum(r.I))
  } else if(Scores=="Gaussian"){
    JntScores = matrix(rnorm(n*r.J), ncol = r.J)
    IndivScores = matrix(rnorm(n*sum(r.I)), nrow = n, ncol = sum(r.I))
  }

    colnames(JntScores) = paste("Jnt Score", 1:r.J)
    IndivScores.names = NULL
      
  for(k in 1:K){
    IndivScores.names = c(IndivScores.names, paste0("Indiv X", k, " Score ", 1:r.I[k]))
  }
  colnames(IndivScores) = IndivScores.names
    
  if(print.cor){print("The correlation between subject scores is given by"); print(round(cor(cbind(JntScores, IndivScores)),4))}
  
  ##############################Define Datasets##############################
  Jnt.Loads.All = list()
  Indiv.Loads.All = list()
  D.J = D.I = list()
  Noise = Joint.Sigs = Indiv.Sigs = list()
  Sig.Mats =  Data.Mats = list()
  temp.fcn = function(x){x[sample(round(length(x)/2))] = 1; x}

  for(k in 1:K){
    if(Loads == "Gaussian"){
      Jnt.Loads.All[[k]] = matrix(rnorm(r.J*p[k]), nrow = r.J, ncol = p[k])
    } else if (Loads == "Fixed"){
      Jnt.Loads.All[[k]] = matrix(apply(matrix(0, nrow = r.J, ncol = p[k]), 1, temp.fcn), nrow = r.J)
    } else if (Loads == "Double_Exp"){
      Jnt.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.J), nrow = r.J ))
    } else if (Loads == "Rademacher"){
      Jnt.Loads.All[[k]] = matrix(rsign(p[k]*(r.J)), nrow = r.J)
    } 
    D.J[[k]] = (equal.eig + 1)*diag(r.J:1) + equal.eig*diag(rep(1,r.J))
    
    Joint.Sigs[[k]] = JntScores%*%sqrt(D.J[[k]])%*%Jnt.Loads.All[[k]]
    
    if(SVD.plots){
      plot(svd(Joint.Sigs[[k]])$d, ylab = "Singular Values")
      title(paste0("SVD of Joint Signal from X", k))
    }
    # (k>1)*sum(P[1:(k-1)])+(1:P[k])
    # temp.IndScores = IndivScores[,grep(paste0("Indiv X", k), colnames(IndivScores))]
    temp.IndScores = IndivScores[,(k>1)*sum(r.I[1:(k-1)])+(1:r.I[k])]
    
    if(Loads == "Gaussian"){
      Indiv.Loads.All[[k]] = matrix(rnorm(n = p[k]*r.I[k]), nrow = r.I[k], ncol = p[k])
    } else if (Loads == "Fixed"){
      temp.fcn = function(x){x[sample(round(length(x)/4))] = 1; x}
      Indiv.Loads.All[[k]] = matrix(apply(matrix(-1, nrow = r.I[k], ncol = p[k]), 1, temp.fcn), nrow = r.I[k])
    } else if (Loads == "Double_Exp"){
      Indiv.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.I[k]), nrow = r.I[k] ))
    } else if (Loads == "Rademacher"){
      Indiv.Loads.All[[k]] = matrix(rsign(p[k]*(r.I[k])), nrow = r.I[k])
    } 
    
    D.I[[k]] = (equal.eig + 1)*diag(r.I[k]:1) + equal.eig*diag(rep(1,r.I[k]))
    
    Indiv.Sigs[[k]] = temp.IndScores%*%sqrt(D.I[[k]])%*%Indiv.Loads.All[[k]]
    
    if(SVD.plots){
      plot(svd(Indiv.Sigs[[k]])$d, ylab = "Singular Values")
      title(paste0("SVD of Individual Signal from X", k))
    }
    
    Sig.Mats[[k]] = Joint.Sigs[[k]] + Indiv.Sigs[[k]]
    
    Noise[[k]] = matrix(rnorm(n*p[k]), nrow  = n)
    
    Data.Mats[[k]] = AdjSigVarExp(Joint.Sigs[[k]], Indiv.Sigs[[k]], Noise[[k]],
                                  JntVarEx[k], IndVarEx[k])
    
    Joint.Sigs[[k]] = Data.Mats[[k]]$J
    Indiv.Sigs[[k]] = Data.Mats[[k]]$I
  }
  
  Blocks = lapply(Data.Mats, function(x) x[["Data"]])
  
  Dat.Comps = list(Joint.Sigs, Indiv.Sigs, Noise)
  names(Dat.Comps) = c("JointSignalMatrices", "IndivSignalMatrices", "NoiseMatrices")
  
  Scores = list(JntScores, IndivScores)
  names(Scores) = c("Joint", "Indiv")
  
  Loadings = list(mapply(function(x,y) x%*%y, D.J, Jnt.Loads.All), 
                  mapply(function(x,y) x%*%y, D.I, Indiv.Loads.All))
  names(Loadings) = c("Joint", "Indiv")
  
  out = list(Dat.Comps, Blocks, Scores, Loadings)
  names(out) = c("Data Components", "Data Blocks", "Scores", "Loadings")
  
  out
}

#############################################################################
###########          pJIVE ML estimation of JIVE model     ##################
###########          that uses EM algorithm                ##################
#############################################################################
ProJIVE_EM_Kge2=function(Y,P,Q,Max.iter=10000,diff.tol=1e-5,plots=TRUE,chord.tol=-1,sig_hat=NULL, init.loads = NULL){
  
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
  } else if(init.loads == "CJIVE"){
    dat.blocks=list(Y[,1:P[1]],Y[,P[1]+(1:P[2])])
    cjive.solution = cc.jive(dat.blocks = dat.blocks, signal.ranks = Q[1]+Q[-1], 
                             joint.rank = Q[1], perm.test = FALSE)
    
    WJ = lapply(cjive.solution$sJIVE$joint_matrices, function(x) x$v)
    WI = lapply(cjive.solution$sJIVE$indiv_matrices, function(x) x$v)
  }
  
  # Block specific loading matrices W_k
  wk_hat=wji_hat=list()
  for(k in 1:K){ # Selecting a sub matrx of the Cholesky Decomp solution L
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
  
  # Total noise matrices D
  d_hat=generate_d(sig_hat,P)
  
  # Initiate LogLik
  all_obs.LogLik=c(-Inf)
  all_complete.LogLik = c(-Inf)
  
  # Initiate Chordal Distance
  chord.dist = 1
  mu_hat=apply(as.matrix(Y),2,sum) / N   
  
  ## Store some values to save computation time
  Yc=sweep(Y, 2, mu_hat)
  S=t(Yc)%*%Yc
  
  Iq=diag(sum(Q))
  Ip=diag(sum(P))
  
  # Set initial iteration number:
  iter=0
  
  while (Max.iter>=iter & eval_converge(all_obs.LogLik,all_obs.LogLik,N*diff.tol) & (chord.dist[iter+1] >= chord.tol)) 
  {
    ################## START OF EM-ALGORITHM ######################
    w=w_hat
    
    c_solv=solve(Iq+t(w)%*%solve(d_hat)%*%w)
    
    U=S%*%solve(d_hat)%*%w%*%c_solv
    V=c_solv+c_solv%*%t(w)%*%solve(d_hat)%*%S%*%solve(d_hat)%*%w%*%c_solv
    
    ## Update d_tild
    
    d_tild=S-2*w%*%t(U)+w%*%V%*%t(w)   
    for(k in 1:K){
      
      ## Update wk_hat 
      wk_hat[[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))
      
      ## Update sigma_hat
      sig_hat[k]=mean(diag(A[[k]]%*%diag(diag(d_tild))%*%t(A[[k]])))
    }
    
    # w=w_hat
    # 
    # Cinv = solve(d_hat)%*%(Ip-w%*%solve(Iq+t(w)%*%solve(d_hat)%*%w)%*%t(w)%*%solve(d_hat))
    # 
    # # Compute subject scores
    # exp.theta = U = Y%*%Cinv%*%w
    # # Compute expected outer product
    # cov.theta = Iq - t(w)%*%Cinv%*%w
    # exp.theta2 = V = cov.theta + t(exp.theta)%*%exp.theta
    # 
    # ## Update d_tild
    # 
    # # d_tild=(S-2*w%*%t(U)%*%Y+w%*%V%*%t(w))/N   
    # d_tild=(S + 2*t(Y)%*%U%*%t(w) - w%*%V%*%t(w))/N   
    # for(k in 1:K){
    #   U = t(Y)%*%exp.theta
    #   
    #   ## Update wk_hat 
    #   wk_hat[[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))
    #   
    #   ## Update sigma_hat
    #   sig_hat[k]=ifelse(is.null(sig_hat),mean(diag(A[[k]]%*%diag(diag(d_tild))%*%t(A[[k]]))),sig_hat[k])
    # }
    # 
    ## Update w_hat 
    w_hat=wk_to_w(wk_hat, P, Q)
    chord.dist = c(chord.dist, chord.norm.diff(w, w_hat))
    
    ## Update d_hat 
    d_hat=generate_d(sig_hat,P)
    
    
    
    ################## End of EM-ALGORITHM ######################
    
    ## Update iter
    iter=iter + 1
    
    # Compute subject scores
    exp.theta = U = Y%*%solve(d_hat)%*%w%*%c_solv
    
    all_obs.LogLik=append(all_obs.LogLik,obs_LogLik(Y, mu_hat, w_hat, d_hat))  
    all_complete.LogLik=append(all_complete.LogLik, complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat))  
  }
  
  
  
  print(paste0("Total Iteration = ",toString(iter)))
  print(paste0("Observed Data Likelihood = ",toString(all_obs.LogLik[length(all_obs.LogLik)])))
  print(paste0("Complete Data Likelihood = ",toString(all_complete.LogLik[length(all_complete.LogLik)])))
  print(paste0("Chordal Norm = ",toString(chord.dist[length(chord.dist)])))
  
  if(plots){
    if(is.finite(min(all_obs.LogLik[-1]))) {
      layout(matrix(1:3, nrow = 1))
      plot(all_obs.LogLik[-1], ylab = "Log-Likelihood", main = "Observed Data Log-Likelihood")
    } else {
      layout(matrix(1:2, nrow = 1))
    }
    if(is.finite(min(all_complete.LogLik[-1]))) {plot(all_complete.LogLik[-1], ylab = "Log-Likelihood", 
                                                      main = "Complete Data Log-Likelihood")}
    plot(chord.dist, ylab = "Chordal Norm", main = "Distance between consecutive \n estimates of 'W'")
    layout(1)
  }
  theta.names = paste0("Joint_Score_", 1:Q[1])
  for(k in 1+1:K){theta.names = c(theta.names, paste0("Individual_Data", k-1, "_Score_", 1:Q[k]))}
  colnames(exp.theta) = theta.names
  names(Q) = c("Joint", "Indiv_1", "Indiv_2")
  out = list(exp.theta, w_hat, Q, sig_hat, chord.dist, all_complete.LogLik[-1], all_obs.LogLik[-1])
  names(out) = c("SubjectScoreMatrix", "LoadingMatrix", "Ranks", "ErrorVariances","ChordalDistances","Complete-Data-Log-Likelihood", "Observed-Data-Log-Likelihood")
  return(out)
}

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