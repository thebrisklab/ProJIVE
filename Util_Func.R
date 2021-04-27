########################################################################
####################START OF PRE-DEFINED FUNCS##########################
########################################################################

## FUNC to Generate the corresponding Ak and Bk:
generate_nm=function(m,n){
  ############################################################################################
  # input:    -n      :an integer for the row number
  #           -m      :an integer for the col number
  #
  # output:   a nxm matrix of 0.
  ############################################################################################
  return(matrix(rep(0,n*m),nrow=n,ncol=m))
}




## FUNC to generate list(G) of W matrix from W_k matrices
wk_to_w=function(wk,q,p){
  
  rearrange_wk=function(k,wgk,q,p){
    c=lapply(q,generate_nm,p[k]) # Generate Col matrices accordingly
    c[[1]]=as.matrix(wgk[[k]][1:p[k],1:q[1]], nrow=p[k], ncol=q[1])
    c[[k+1]]=as.matrix(wgk[[k]][1:p[k],q[1]+1:q[k+1]], nrow=p[k], ncol=q[k+1])
    
    return(do.call('cbind',c))
  }
  
  combine_wk=function(wg,q,p){
    return(do.call('rbind',lapply(1:length(wg),rearrange_wk,wg,q,p)))
  }
  
  #     w=list()
  #     for(g in 1:G){
  # #         w[[g]]=do.call('rbind',lapply(1:K,generate_w,g))
  # #         w[[g]]=do.call('rbind',lapply(1:K,rearrange_wk,wk_hat[[g]]))
  #           w[[g]]=combine_wk(wk_hat[[g]],Q,P)
  #     }
  #     return(w)
  return(lapply(wk,combine_wk,q,p))   
}


## FUNC to generate list(G) of W_k matrix from W matrices
w_to_wk=function(w,p,q){
    extract_w=function(r,c,wg){
        return(as.matrix(wg[r,append((1:q[1]),c)]))
    }
    
    split_w=function(wg,p,q){
        rows=split(1:sum(p),rep(1:length(p),p))
        
        cols=split(1:sum(q[-1])+1,rep(1:length(q[-1]),q[-1]))
        
        return(mapply(extract_w,rows,cols,MoreArgs=list(wg=wg),SIMPLIFY=FALSE))
    }

    return(lapply(w,split_w,p,q))
}


## FUNC to generate diagonal D matrix from sig_hat
generate_d=function(sig, p){ 
  return(diag(rep(sig,p)))
}

## FUNC to generate diagonal sig  from D matrix
generate_sig=function(d, p){ 
  return(unname(unlist(lapply(split(diag(d),rep(1:length(p),p)),unique))))
}



obs_LogLik_MPJIVE<-function(y, pie, mu, w, d){
  ##############################################
  # input:    -Y     :a nxp data frame as the observations
  #           -pie   :a G length vector
  #           -mu    :a G list of p dimension vectors
  #           -w     :a G list of pxq matrices
  #           -d     :a G list of p length vector indicating the noise
  #           
  # output:   a real value of the log likelihood
  ##############################################
  
  N=dim(y)[1]
  
  if(length(mu)==length(d) & length(mu)==length(pie)){
    G=length(mu)
  } else {
    stop("Number of components G in 'mu', 'sig', 'pie' are different")
  }
  
  lik<-rep(0, N)
  
  for(g in 1:G){
    #         w=cbind(do.call(rbind,wj[[g]]),as.matrix(bdiag(wi[[g]])))
    s=w[[g]]%*%t(w[[g]])+d[[g]]
    
    lik=lik+pie[g]*dmvnorm(y,mu[[g]],s)
  }
  
  LogLik=sum(log(lik))
  
  return(LogLik)
}

obs_LogLik_GMM<-function(y, pie, mu, sigma){
  ##############################################
  # input:    -Y     :a nxp data frame as the observations
  #           -pie   :a G length vector
  #           -mu    :a G list of p dimension vectors
  #           -w     :a G list of pxq matrices
  #           -d     :a G list of p length vector indicating the noise
  #           
  # output:   a real value of the log likelihood
  ##############################################
  
  N=dim(y)[1]
  
  if(length(mu)==length(sigma) & length(mu)==length(pie)){
    G=length(mu)
  } else {
    stop("Number of components G in 'mu', 'sigma', 'pie' are different")
  }
  
  lik<-rep(0, N)
  
  for(g in 1:G){

    
    lik=lik+pie[g]*dmvnorm(y,mu[[g]],sigma[[g]])
  }
  
  LogLik=sum(log(lik))
  
  return(LogLik)
}


eval_class=function(p){
  return(which(p==max(p)))
}
########################################################################
####################END OF PRE-DEFINED FUNCS############################
########################################################################