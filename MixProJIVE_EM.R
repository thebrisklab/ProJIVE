source("Util_Func.R")
source("ProJIVE_EM.R")
MixProJIVE_EM=function(Y,P,Q,G,Max.iter=10000,diff.tol=1e-5,pie_hat=NA,mu_hat=NA,wk_hat=NA,sig_hat=NA,print=TRUE){
  
  ## Primary Parameters    
  # Y=as.matrix(df[1:6]);G=3;P=c(3,3);Q=c(1,1,1);Max.iter=10000;diff.tol=1e-5
  # pie_hat=NA;mu_hat=NA;wk_hat=NA;sig_hat=NA;print=TRUE
  #p_1=3, p_2=3, q_J=1, q_1=1, q_2=1
  # P=c(3,3)   # Dimensions of feature blocks:p_1, p_2 
  # Q=c(1,1,1) # Dimensions of scores: q_J, q_1, q_2
  
  
  
  ## Secondary Parameters
  
  ## Total sample size
  N=dim(Y)[1]
  
  if(sum(P)!=dim(Y)[2]){
    stop("Dimension of Y NOT in accordance with P!")
  }
  
  if(is.null(colnames(Y))){
    colnames(Y)=paste(rep('Y',sum(P)), 1:sum(P), sep="_")
  }
  
  initial=FALSE
  
  if(all(is.na(pie_hat)==FALSE) & all(is.na(mu_hat)==FALSE) & all(is.na(wk_hat)==FALSE) & all(is.na(sig_hat)==FALSE)){
    
    initial=TRUE
    
    if(length(pie_hat)==length(mu_hat) & length(mu_hat)==length(wk_hat) & length(wk_hat)==length(sig_hat)){
      G=length(pie_hat)
    }else{
      stop("Intial pie, mu, w, sig lengths are different")
    }
    ## Initialization of W and D
    w_hat=wk_to_w(wk_hat,Q,P)
    
    d_hat=lapply(sig_hat,generate_d,P)
    
  }
  
  ## Total number feature blocks
  if(length(P)==(length(Q)-1)){
    K=length(P)
  } else{
    stop("Error: The number of feature blocks don not match the number of invididual scores")
  }
  
  ## Generate selection matrices A_knd B_k
  A=list()
  B=list()
  for(k in 1:K){
    up=lapply(Q,generate_nm,Q[1])
    up[[1]]=diag(Q[1])
    down=lapply(Q,generate_nm,Q[(k+1)])
    down[[(k+1)]]=diag(Q[(k+1)])
    B[[k]]=rbind(do.call("cbind",up),
                 do.call("cbind",down))
    
    a=lapply(P,generate_nm,P[k])
    a[[k]]=diag(P[k])
    A[[k]]=do.call("cbind",a)
  }
  
  
  all_obs.LogLik=c()
  
  ## Set initial iteration number:
  iter=0
  diff=Inf
  ## Set initial conditional probabilities:
  ind_density=matrix(NA,nrow=N,ncol=G) #individual contribution to likelihood (NxG)
  
  # while (Max.iter>=iter & eval_converge(all_obs.LogLik) ) 
  while (Max.iter>=iter & diff>diff.tol)     
  {
    ################## START OF E-STEP######################
    
    if(iter==0 & initial==FALSE){
      
      ## E-step: computing the conditional posterior probabilities tau_hat
      
      ### First assign every one to be 1/G or a very small value
      #             tau_hat=matrix(rep(1/G,N*G), nrow=N, ncol=G)
      tau_hat=matrix(rep(0,N*G), nrow=N, ncol=G)
      
      ### Then randomly assign every one to be 1/G
      subsample=sample(1:N,round(N/10),replace=FALSE)
      tau_hat[subsample,]=t(rmultinom(round(N/10), 1, rep(1,G)/G))
      
      pie_hat=apply(tau_hat,2,mean)
      
      mu_hat=list()
      wk_hat=list()
      sig_hat=list()
      for(g in 1:G){ # Start Loop for each g=1:G mixture components
        
        ## Update mu_hat
        Y_sub=Y[tau_hat[,g]==1,]
        
        subres=ProJIVE_EM(Y=Y_sub,P=P,Q=Q,Max.iter=1000,diff.tol=1e-3,print=FALSE)
        
        mu_hat[[g]]=subres$Mu
        wk_hat[[g]]=subres$Wk
        sig_hat[[g]]=subres$Sig
        
      }
      # Initialization of W and D
      w_hat=wk_to_w(wk_hat,Q,P)
      
      d_hat=lapply(sig_hat,generate_d,P)
      
    }else{
      for(g in 1:G){
        ## E-step: computing individual observation's contribution to the likelihood
        COV=w_hat[[g]]%*%t(w_hat[[g]])+d_hat[[g]]
        ind_density[,g]=pie_hat[g]*dmvnorm(Y,mu_hat[[g]],COV)
      }
      
      ## E-step: computing the conditional posterior probabilities tau_hat
      tau_hat=ind_density/apply(ind_density,1,sum) #(NxG)
      
      tau_hat[is.na(tau_hat)]=1/G # adjust the nans 
      tau_hat[tau_hat < 1e-10] = 1e-10
      ## Evaluate Likelihood using the ind_density for convenience
      #         obs.LogLik=sum(log(apply(ind_density,1,sum)))
      #         all_obs.LogLik=append(all_obs.LogLik,obs.LogLik)
      
      
      
      
      ################## END OF E-STEP######################
      
      
      ################## START OF M-STEP######################
      
      ## Update pie_hat
      #     pie_hat=apply(tau_hat,2,sum)/N
      
      new_pie_hat=apply(tau_hat,2,mean)
      diff=max(abs(pie_hat-new_pie_hat))
      pie_hat=new_pie_hat
      
      ## Update mu_hat
      new_mu_hat=list()
      for(g in 1:G){ # Start Loop for each g=1:G mixture components
        
        ## Update mu_hat
        Y_star=Y # use Y_star in case we need to implement censoring later
        
        new_mu_hat[[g]]=apply(tau_hat[,g] * as.matrix(Y_star),2,sum) / sum(tau_hat[,g])   
        
        ## Store some values to save computation time
        Yc=sweep(Y, 2, new_mu_hat[[g]])
        S=t(tau_hat[,g]*Yc)%*%Yc/ sum(tau_hat[,g])
        
        Iq=diag(sum(Q))
        w=w_hat[[g]]

        c_solv=solve(Iq+t(w)%*%solve(d_hat[[g]])%*%w)
        
        U=S%*%solve(d_hat[[g]])%*%w%*%c_solv
        V=c_solv+c_solv%*%t(w)%*%solve(d_hat[[g]])%*%S%*%solve(d_hat[[g]])%*%w%*%c_solv
        
        ## Update d_tild
        
        d_tild=S-2*w%*%t(U)+w%*%V%*%t(w)   
        for(k in 1:K){
          
          ## Update wk_hat 
          wk_hat[[g]][[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))
          
          ## Update sigma_hat
          sig_hat[[g]][k]=mean(diag(A[[k]]%*%diag(diag(d_tild))%*%t(A[[k]])))
        }
        
      }
      
      diff=max(diff,abs(unlist(new_mu_hat)-unlist(mu_hat)))
      mu_hat=new_mu_hat
      
      ## Update w_hat 
      new_w_hat=wk_to_w(wk_hat,Q,P)
      diff=max(diff,abs(unlist(new_w_hat)-unlist(w_hat)))
      w_hat=new_w_hat
      
      #     w_hat=lapply(wk_hat,combine_wk,Q,P)
      #     w_hat=wk_to_w()
      ## Update d_hat 
      new_d_hat=lapply(sig_hat,generate_d,P)
      diff=max(diff,abs(unlist(new_d_hat)-unlist(d_hat)))
      d_hat=new_d_hat
    } 
    ################## End of M-step######################
    
    ## Update iter
    iter=iter + 1
    
    all_obs.LogLik=append(all_obs.LogLik,obs_LogLik_MPJIVE(Y, pie_hat,mu_hat, w_hat, d_hat))  
  }
  
  if(iter<Max.iter & diff<=diff.tol){
    converge_=TRUE
  }else{
    converge_=FALSE
  }
  
  names(pie_hat)=paste0("pie",1:G)
  names(mu_hat)=paste0("mu",1:G)
  names(wk_hat)=paste0("wk",1:G)
  names(sig_hat)=paste0("sig",1:G)
  
  if(print){
    print(paste0("Total Iteration = ",toString(iter)))
    print(paste0("Convergence = ",toString(converge_)))
    print(paste0("LogLik = ",toString(all_obs.LogLik[length(all_obs.LogLik)])))
    print(pie_hat)
    print(mu_hat)
    print(wk_hat)
    print(sig_hat)
    plot(all_obs.LogLik[-1])
  }
  LogLik=tail(all_obs.LogLik, n=1)
  nparm=length(as.vector(unlist(mu_hat)))+length(as.vector(unlist(wk_hat)))+length(as.vector(unlist(sig_hat)))+length(pie_hat)-1
  
  POST_PROB=tau_hat
  POST_PROB[POST_PROB < 1e-8] = 1e-8
  entropy=-sum(POST_PROB*log(POST_PROB))
  entropy.rsquare=1 - entropy/(N*log(G))
  
  AIC=-2*LogLik+nparm*2
  BIC=-2*LogLik+nparm*log(N)
  ICL=BIC+2*entropy
  
  OUTPUT=list()
  
  OUTPUT[[1]]=iter
  OUTPUT[[2]]=converge_
  OUTPUT[[3]]=LogLik
  OUTPUT[[4]]=AIC
  OUTPUT[[5]]=BIC
  OUTPUT[[6]]=ICL
  
  OUTPUT[[7]]=pie_hat
  OUTPUT[[8]]=mu_hat
  OUTPUT[[9]]=wk_hat
  OUTPUT[[10]]=sig_hat
  OUTPUT[[11]]=tau_hat
  OUTPUT[[12]]=apply(tau_hat,1,eval_class)
  
  names(OUTPUT)<-c("Iterations", "Converged", "LogLik","AIC", "BIC", "ICL","Pie", "Mu", "Wk", "Sig","Posterior","Class")
  
  return(OUTPUT)
}