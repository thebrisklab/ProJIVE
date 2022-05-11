###########################################################################################################################
#######                 Functions for JIVE-PNC Project(s)                              ####################################
#######                 Author: Raphiel J. Murden                                      ####################################
#######                 Supervised by Benjamin Risk                                    ####################################
###########################################################################################################################
require(rootSolve); require(Matrix); require(ggplot2); require(reshape2); require(fields); require(mvtnorm)
require(dplyr); require(xtable); require(optimx); require(gplots); require(MASS); require(r.jive); 
require(extraDistr)

ajive.dir = "r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

######################################################################################################################
###########   Generates 2 Simulated Datasets that follow JIVE Model using binary subject scores   ####################
######################################################################################################################
GenerateToyData <- function(n, p, JntVarEx, IndVarEx, jnt_rank = 1, equal.eig = F, ind_ranks, JntVarAdj = T, mix.probs = NULL,
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
      Jnt.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.J)), nrow = r.J)
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
      Indiv.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.I[k])), nrow = r.I[k])
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
  
  roots = multiroot(simul.quads, start, parms = parms)
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
  
  Method = substr(sim.names[grep("Scores",sim.names)], 1,7)
  Method = sub("E.O", "E", Method); Method = sub("E.", "E", Method)
  Method = rep(as.factor(Method), each = num.sims)
  Method = relevel(Method, "ProJIVE")
  
  Type = substring(sim.names[grep("Scores", sim.names)], 8)
  Type = sub("racle.", "", Type); Type = sub(".J", "J", Type); Type = sub(".I", "I", Type);
  Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type)
  Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type)
  Type = rep(as.factor(Type), each = num.sims)
  
  n.levs = nlevels(Type)*nlevels(Method)
  sim.score.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1.raw, each = n.levs), 
                                  JVE_2 = rep(JVE_2.raw, each = n.levs),
                                  IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  Norm = c(as.numeric(as.matrix(AllSims[,grep("Loads", sim.names)])))
  
  Method = rep(substr(sim.names[grep("Loads",sim.names)], 1,7), each = num.sims)
  Method = sub("E.O", "E", Method); Method = sub("E.", "E", Method)
  Method = as.factor(Method)
  Method = relevel(Method, "ProJIVE")
  
  Type = rep(substring(sim.names[grep("Loads", sim.names)], 7), num.sims)
  Type = sub("Oracle.", "", Type); Type = sub("E.", "", Type)
  Type = sub(".I", "I", Type);  Type = sub(".J", "J", Type);
  Type = sub("Loads", "Loadings", Type)
  Type = gsub("X", "X[1]", Type)
  Type = gsub("Y", "X[2]", Type)
  Type = as.factor(gsub("[.]", " ", Type))
  
  n.levs = nlevels(Type)*nlevels(Method)
  JVE_1 = rep(JVE_1.raw, each = n.levs)
  JVE_2 = rep(JVE_2.raw, each = n.levs)
  
  sim.load.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = JVE_1, JVE_2 = JVE_2, 
                                 IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))

  sim.all.norms.gg = rbind(sim.score.norms.gg, sim.load.norms.gg)
  out = sim.all.norms.gg 
  out
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
    w_k[[k]] = w[P[k]+(1:P[k+1]),c(1:Q[1],sum(Q[1:k])+1:Q[k+1])]
  }
  
  
  return(w_k)
}

## FUNC to generate list(G) of D matrix from sigma vectors
generate_d=function(sig_lst, p_vec){
  
  d=list()
  for(k in 1:length(p_vec)){
    d[[k]]=rep(sig_lst[k],p_vec[k])
  }
  D=diag(unlist(d))
  
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
      return((diff.ll>=(diff.tol)))
    }
  }
}
########################################################################
####################END OF PRE-DEFINED FUNCS############################
########################################################################

#############################################################################
###########          pJIVE ML estimation of JIVE model     ##################
###########          that uses EM algorithm                ##################
#############################################################################
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
    stop("Error: The number of feature blocks do not match the number of invididual scores")
  }
  
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
    WJ = list(matrix(rnorm(Q[1]*P[1]), nrow = P[1]),matrix(rnorm(Q[1]*P[2]), nrow = P[2]))
    WI = list(matrix(rnorm(Q[2]*P[1]), nrow = P[1]),matrix(rnorm(Q[3]*P[2]), nrow = P[2]))
  } else if (is.list(init.loads)){
    WJ = init.loads[[1]]
    WI = init.loads[[2]]
  } else if(init.loads == "AJIVE"){
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
    sig_hat=c(mean(svd(Y[,1:P[1]], nu = 0, nv = 0)$d[-(1:sum(Q[-3]))]),mean(svd(Y[,P[1]+(1:P[2])], nu = 0, nv = 0)$d[-(1:sum(Q[-2]))]))
  # } else if (is.numeric(sig_hat)){
  #   print("Using fixed values for noise variances:", round(sig_hat, 3))
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