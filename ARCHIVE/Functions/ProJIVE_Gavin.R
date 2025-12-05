ProJIVE_EM2=function(Y,P,Q,Max.iter=10000,diff.tol=1e-5,print=TRUE){
  

  ########################################################################
  ####################START OF PRE-DEFINED FUNCS##########################
  ########################################################################
  ## FUNC to Generate the corresponding Ak and Bk:
  generate_ab=function(m,n){
    return(matrix(rep(0,n*m),nrow=n,ncol=m))
  }
  
  
  ## FUNC to generate list(G) of W matrix with right dimension
  generate_w=function(k){
    r=lapply(Q,generate_ab,P[k])
    r[[1]]=as.matrix(wk_hat[[k]][1:P[k],1:Q[1]], nrow=P[k], ncol=Q[1])
    r[[k+1]]=as.matrix(wk_hat[[k]][1:P[k],Q[1]+1:Q[k+1]], nrow=P[k], ncol=Q[k+1])
    return(do.call('cbind',r))
  }
  
  ## FUNC to generate list(G) of W matrix from W_k matrices
  wk_to_w=function(wk){
    
    w=do.call('rbind',lapply(1:K,generate_w))
    
    return(w)
  }
  
  ## FUNC to generate list(G) of D matrix from sigma vectors
  generate_d=function(sig_lst, p_vec){
    
    d=list()
    for(k in 1:K){
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
    
    
    #         w=cbind(do.call(rbind,wj[[g]]),as.matrix(bdiag(wi[[g]])))
    s=w%*%t(w)+d
    
    lik=dmvnorm(Y,mu,s)
    
    
    LogLik=sum(log(lik))
    
    return(LogLik)
  }
  
  ## FUNC to evaluate convergence in the loop
  eval_converge=function(vals_vec){
    if(length(vals_vec)==1) {
      return(TRUE)
    }else{
      if(vals_vec[length(vals_vec)]==-Inf){
        return(TRUE)
      }else{
        return(abs(all_obs.LogLik[length(all_obs.LogLik)]-all_obs.LogLik[length(all_obs.LogLik)-1])>=(N*diff.tol))
      }
    }
  }
  ########################################################################
  ####################END OF PRE-DEFINED FUNCS############################
  ########################################################################
  
  
  
  
  # Total sample size
  N=dim(Y)[1]
  
  # Total number feature blocks
  if(length(P)==(length(Q)-1)){
    K=length(P)
  } else{
    stop("Error: The number of feature blocks does not match the number of invididual signals.")
  }
  
  # Selection matrices A_k and B_k
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
  
  # Block specific loading matrices W_k
  
  wji_hat=list()
  for(k in 1:K){ # Selecting a sub matrx of the Cholesky Decomp solution L
    if(k==1){
      L=t(chol(cov(Y[,1:P[k]])))
    }else{
      L=t(chol(cov(Y[,(sum(P[1:(k-1)])+1):sum(P[1:k])])))
    }
    wji_hat[[k]]=L[,1:(Q[1]+Q[k+1])]
  }
  wk_hat=wji_hat
  
  
  
  # wji_hat=list()
  #     for(k in 1:K){ # Selecting a sub matrx of the Cholesky Decomp solution L
  
  #         wji_hat[[k]]=matrix(rep(1,P[k]*(Q[1]+Q[k+1])),nrow=P[k],ncol=(Q[1]+Q[k+1]))
  
  #     }
  # wk_hat=wji_hat
  
  
  sig_hat=rep(1,K)                   # Initializationf of sig_hat
  
  # Total loading matrices W
  w_hat=wk_to_w(wk_hat)
  
  # Total noise matrices D
  d_hat=generate_d(sig_hat,P)
  
  # all_obs.LogLik=append(c(-Inf),obs_LogLik(Y, mu_hat, w_hat, d_hat, pie_hat))
  all_obs.LogLik=c(-Inf)
  
  # Set initial iteration number:
  iter=0
  
  
  
  while (Max.iter>=iter & eval_converge(all_obs.LogLik) ) 
  {
    
    
    
    ################## START OF M-STEP######################
    
    
    ## Update mu_hat
    mu_hat=list()
    
    
    ## Update mu_hat
    Y_star=Y # use Y_star in case we need to implement censoring later
    
    mu_hat=apply(as.matrix(Y_star),2,sum) / N   
    
    ## Store some values to save computation time
    Yc=sweep(Y, 2, mu_hat)
    S=t(Yc)%*%Yc/ N
    
    Iq=diag(sum(Q))
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
    
    
    ## Update w_hat 
    w_hat=wk_to_w(wk_hat)
    
    ## Update d_hat 
    d_hat=generate_d(sig_hat,P)
    
    ################## End of M-step######################
    
    ## Update iter
    iter=iter + 1
    
    all_obs.LogLik=append(all_obs.LogLik,obs_LogLik(Y, mu_hat, w_hat, d_hat))  
  }
  
  #     if(iter<Max.iter & diff<=diff.tol){
  if(iter<Max.iter){
    converge_=TRUE
  }else{
    converge_=FALSE
  }    
  
  if(print){
    print(paste0("Total Iteration = ",toString(iter)))
    print(paste0("Convergence = ",toString(converge_)))
    print(paste0("LogLik = ",toString(round(tail(all_obs.LogLik, n=1),3))))
    # print(pie_hat)
    print(mu_hat)
    print(wk_hat)
    print(sig_hat)
    #     plot(all_obs.LogLik[-1])
  }
  LogLik=tail(all_obs.LogLik, n=1)
  
  
  OUTPUT=list()
  
  OUTPUT[[1]]=iter
  OUTPUT[[2]]=(iter<Max.iter)
  OUTPUT[[3]]=LogLik
  OUTPUT[[4]]=mu_hat
  OUTPUT[[5]]=wk_hat
  OUTPUT[[6]]=sig_hat
  OUTPUT[[7]]=all_obs.LogLik
  
  names(OUTPUT)<-c("Iterations", "Converged", "LogLik", "Mu", "Wk", "Sig", "AllLogLik")
  return(OUTPUT)
}
