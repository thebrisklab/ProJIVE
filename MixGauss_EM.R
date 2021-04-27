MixGauss_EM=function(Y,G,Max.iter=10000,diff.tol=1e-5,pie_hat=NA,mu_hat=NA,sigma_hat=NA,print=TRUE){
  
  ## Primary Parameters    
  # Y=as.matrix(df[1:6]);G=3;P=c(3,3);Q=c(1,1,1);Max.iter=10000;diff.tol=1e-5
  # pie_hat=NA;mu_hat=NA;wk_hat=NA;sig_hat=NA;print=TRUE
  #p_1=3, p_2=3, q_J=1, q_1=1, q_2=1
  # P=c(3,3)   # Dimensions of feature blocks:p_1, p_2 
  # Q=c(1,1,1) # Dimensions of scores: q_J, q_1, q_2
  
  
  
  ## Secondary Parameters
  
  ## Total sample size
  N=dim(Y)[1]
  P=dim(Y)[2]

  
  if(is.null(colnames(Y))){
    colnames(Y)=paste(rep('Y',sum(P)), 1:sum(P), sep="_")
  }
  
  initial=FALSE
  
  if(all(is.na(pie_hat)==FALSE) & all(is.na(mu_hat)==FALSE) & all(is.na(sigma_hat)==FALSE)){
    
    initial=TRUE
    
    if(length(pie_hat)==length(mu_hat) & length(mu_hat)==length(sigma_hat) ){
      G=length(pie_hat)
    }else{
      stop("Intial pie, mu, w, sig lengths are different")
    }
    
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
      sigma_hat=list()
      for(g in 1:G){ # Start Loop for each g=1:G mixture components
        
        ## Update mu_hat
        Y_sub=Y[tau_hat[,g]==1,]
        
        mu_hat[[g]]=apply(Y_sub,2,mean)

        sigma_hat[[g]]=cov(Y_sub)
        
      }

    }else{
      for(g in 1:G){
        ## E-step: computing individual observation's contribution to the likelihood
        ind_density[,g]=pie_hat[g]*dmvnorm(Y,mu_hat[[g]],sigma_hat[[g]])
      }
      
      ## E-step: computing the conditional posterior probabilities tau_hat
      tau_hat=ind_density/apply(ind_density,1,sum) #(NxG)
      
      tau_hat[is.na(tau_hat)]=1/G # adjust the nans 
      tau_hat[tau_hat < 1e-10] = 1e-10
      # Evaluate Likelihood using the ind_density for convenience
      obs.LogLik=sum(log(apply(ind_density,1,sum)))
      all_obs.LogLik=append(all_obs.LogLik,obs.LogLik)
      
      
      
      
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
        
      }
      
      diff=max(diff,abs(unlist(new_mu_hat)-unlist(mu_hat)))
      mu_hat=new_mu_hat
      
      new_sigma_hat=list()
      for(g in 1:G){ # Start Loop for each g=1:G mixture components
        
        ## Update mu_hat
        Y_star=sweep(Y, 2, mu_hat[[g]]) # use Y_star in case we need to implement censoring later
        
        new_sigma_hat[[g]]=t(tau_hat[,g] * as.matrix(Y_star))%*%as.matrix(Y_star) / sum(tau_hat[,g])   
        
      }

      ## Update sigma_hat 
        
      diff=max(diff,abs(unlist(new_sigma_hat)-unlist(sigma_hat)))
      sigma_hat=new_sigma_hat
    } 
    ################## End of M-step######################
    
    ## Update iter
    iter=iter + 1
    
#     all_obs.LogLik=append(all_obs.LogLik,obs_LogLik(Y, pie_hat,mu_hat, w_hat, d_hat))  
  }
  
  if(iter<Max.iter & diff<=diff.tol){
    converge_=TRUE
  }else{
    converge_=FALSE
  }
  
  names(pie_hat)=paste0("pie",1:G)
  names(mu_hat)=paste0("mu",1:G)
  names(sigma_hat)=paste0("sigma",1:G)
  
  if(print){
    print(paste0("Total Iteration = ",toString(iter)))
    print(paste0("Convergence = ",toString(converge_)))
    print(paste0("LogLik = ",toString(all_obs.LogLik[length(all_obs.LogLik)])))
    print(pie_hat)
    print(mu_hat)
    print(sigma_hat)
    plot(all_obs.LogLik[-1])
  }
  LogLik=tail(all_obs.LogLik, n=1)
  nparm=length(as.vector(unlist(mu_hat)))+length(sigma_hat)*(dim(sigma_hat)[1]^2+dim(sigma_hat)[1])/2+length(pie_hat)-1
  
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
  OUTPUT[[9]]=sigma_hat
  OUTPUT[[10]]=tau_hat
  OUTPUT[[11]]=apply(tau_hat,1,eval_class)
  
  names(OUTPUT)<-c("Iterations", "Converged", "LogLik","AIC", "BIC", "ICL","Pie", "Mu", "Sigma","Posterior","Class")
  
  return(OUTPUT)
}