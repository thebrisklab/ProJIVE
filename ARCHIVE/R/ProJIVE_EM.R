#' PJIVE ML estimation of JIVE model that uses EM algorithm
#'
#' @param Y Matrix. The observed data.
#' @param P Integer vector. The number of variables in each dataset.
#' @param Q Integer vector. The number of components in each dataset.
#' @param Max.iter The maximum number of iterations for the EM algorithm (default = 1000).
#' @param diff.tol Numeric. Specifies the tolerance for convergence based on log-likelihood (default = 1e-8).
#' @param plots Logical (TRUE/FALSE): Specifies if plots are generated (default = TRUE).
#' @param chord.tol Numeric. Tolerance for the chordal distance for convergence (default = -1).
#' @param sig_hat The error variances for each dataset.
#' @param init.loads Either a list of initial loadings or a string specifying the initialization method ("AJIVE" or "CJIVE").
#'
#' @return A list containing:
#'            - SubjectScoreMatrix: The matrix of subject scores.
#'            - LoadingMatrix: The final loading matrix.
#'            - Ranks: The ranks for joint and individual components.
#'            - VarianceExplained: The variance explained by the joint and individual components for each dataset.
#'            - ErrorVariances: The estimated error variances.
#'            - ChordalDistances: The distance between consecutive estimates of the loading matrix.
#'            - Complete-Data-Log-Likelihood: The log-likelihood for the complete data.
#'            - Observed-Data-Log-Likelihood: The log-likelihood for the observed data.
#'            - BIC: The Bayesian Information Criterion for the model.
#'            - AIC: The Akaike Information Criterion for the model.
#' @export
#'
#' @examples
#' Y <- matrix(rnorm(100), 20, 5)
#' P <- c(5, 5)
#' Q <- c(2, 3)
#' result <- ProJIVE_EM(Y, P, Q)
#' str(result)
ProJIVE_EM=function(Y,P,Q,Max.iter=10000,diff.tol=1e-8,plots=TRUE,chord.tol=-1,sig_hat=NULL, init.loads = NULL){

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


    #         w=cbind(do.call(rbind,wj[[g]]),as.matrix(bdiag(wi[[g]])))
    s=w%*%t(w)+d

    lik=mvtnorm::dmvnorm(Y,mu,s)


    LogLik=sum(log(lik))

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
  eval_converge=function(vals_vec, diff.tol){
    len = length(vals_vec)
    if(length(vals_vec)==1) {
      return(TRUE)
    }else{
      if(vals_vec[length(vals_vec)]==-Inf){
        return(TRUE)
      }else{
        diff.ll = vals_vec[len]-vals_vec[len-1]
        return((diff.ll>=(diff.tol)))
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
  Q.tot = Q[1] + Q[-1]

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

  # get initial estimates of loadings matrices via cc.jive
  wk_hat=wji_hat=list()
  if (is.null(init.loads)){
    #Default: Initialize loadings as a sub matrix of the Cholesky Decomposition solution L
    WJ = WI = list()
    for(k in 1:K){
      if((k==1)){
        L=t(chol(cov(Y[,1:P[k]])))
      } else if(k>1){
        L=t(chol(cov(Y[,(sum(P[1:(k-1)])+1):sum(P[1:k])])))
      }
      wji_hat[[k]]=L[,1:(Q[1]+Q[k+1])]
      WJ[[k]] = wji_hat[[k]][,1:Q[1]]
      WI[[k]] = wji_hat[[k]][,-(1:Q[1])]
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
    ajive.solution = CJIVE::cc.jive(dat.blocks, signal.ranks = Q[1]+Q[-1],joint.rank = Q[1], perm.test = FALSE)

    WJ = lapply(ajive.solution$sJIVE$joint_matrices, function(x) x$v)
    WI = lapply(ajive.solution$sJIVE$indiv_matrices, function(x) x$v)
  }

  # Block specific loading matrices W_k
  for(k in 1:K){
    wji_hat[[k]]=cbind(WJ[[k]], WI[[k]])
  }
  wk_hat=wji_hat

  # Total loading matrices W
  w_hat=wk_to_w(wk_hat,P,Q)

  # Initializationf of sig_hat
  # sig_hat=list(rep(1,K))
  if(is.null(sig_hat)){
    sig_hat = stats::rnorm(length(P))
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
  mu_hat=apply(as.matrix(Y),2,sum)/N

  Iq=diag(sum(Q))
  Ip=diag(sum(P))

  c_solv=solve(Iq+t(w_hat)%*%solve(d_hat)%*%w_hat)
  exp.theta =  Y%*%solve(d_hat)%*%w_hat%*%c_solv

  all_obs.LogLik=obs_LogLik(Y, mu_hat, w_hat, d_hat)
  all_complete.LogLik = complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat)

  # Set initial iteration number:
  iter=0

  #create flag to stop EM algorithm if/when necessary
  while (Max.iter>=iter
         & eval_converge(all_obs.LogLik,N*diff.tol)
         & eval_converge(all_complete.LogLik,N*diff.tol)
  )
  {
    ################## START OF EM-ALGORITHM ######################
    ## Store some values to save computation time
    Yc=sweep(Y, 2, mu_hat)
    S=t(Yc)%*%Yc/N


    w=w_hat
    d=d_hat

    c_solv=solve(Iq+t(w)%*%solve(d_hat)%*%w)

    U=S%*%solve(d_hat)%*%w%*%c_solv
    V=c_solv+c_solv%*%t(w)%*%solve(d_hat)%*%S%*%solve(d_hat)%*%w%*%c_solv

    ## Update d_tild

    d_tild=S-2*w%*%t(U)+w%*%V%*%t(w)
    wk_hat_old = wk_hat
    sig_hat_old = sig_hat
    for(k in 1:K){

      ## Update wk_hat
      wk_hat_temp = list()
      wk_hat[[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))

      ## Update sigma_hat
      sig_hat[k]=mean(diag(A[[k]]%*%diag(diag(d_tild))%*%t(A[[k]])))
    }


    w_hat=wk_to_w(wk_hat, P, Q)
    chord.dist = c(chord.dist, chord.norm.diff(w, w_hat))

    ## Update d_hat
    d_hat=generate_d(sig_hat,P)

    ################## End of EM-ALGORITHM ######################
    iter=iter + 1

    # Compute subject scores
    c_solv=chol2inv(chol(Iq+t(w_hat)%*%solve(d_hat)%*%w_hat))
    exp.theta = Y%*%solve(d_hat)%*%w_hat%*%c_solv

    all_obs.LogLik=append(all_obs.LogLik, obs_LogLik(Y, mu_hat, w_hat, d_hat))
    all_complete.LogLik=append(all_complete.LogLik, complete_LogLik(Y, exp.theta, mu_hat, w_hat, d_hat))

    #Take previous iteration of solution if current iteration decreases LogLik
    ## Update iter

  }

  # len = length(all_obs.LogLik)
  # if(len>1){
  #   temp.flag = !(all_obs.LogLik[len] < all_obs.LogLik[len-1] & all_complete.LogLik[len] < all_complete.LogLik[len-1])
  #
  #   if(temp.flag){
  #     w_hat=wk_to_w(wk_hat_old, P, Q)
  #     chord.dist = chord.dist[1:iter]
  #     all_obs.LogLik = all_obs.LogLik[1:(len-1)]
  #     all_complete.LogLik = all_complete.LogLik[1:(len-1)]
  #
  #     ## Update d_hat
  #     d_hat=generate_d(sig_hat_old,P)
  #
  #     # Compute subject scores
  #     c_solv=chol2inv(chol(Iq+t(w_hat)%*%solve(d_hat)%*%w_hat))
  #     exp.theta = Y%*%solve(d_hat)%*%w_hat%*%c_solv
  #   }
  # }
  ################## End of EM-ALGORITHM ######################

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
