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
#' @param isotropic.error Logical (TRUE/FALSE) Indicator of whether to assume isotropic error, i.e., one scalar parameter for error variance per data set, or non-isotropic error, i.e., P_k scalar parameters for error variance in data set k (default = TRUE)
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
ProJIVE_EM=function(Y, P, Q, Max.iter=10000, diff.tol=1e-8, plots=TRUE,
                    chord.tol=-1, sig_hat=NULL, init.loads = NULL,
                    isotropic.error = TRUE){

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

  ## FUNC to generate the length-p vector of D's diagonal entries from sigma
  ## vectors. D is always diagonal in this model (Sec 2.2), so it is kept as
  ## a plain vector rather than expanded into a dense pxp matrix: forming
  ## that matrix (and calling solve()/det() on it, as the previous
  ## implementation did every EM iteration) costs O(p^3) even though the
  ## true cost of inverting/taking the determinant of a diagonal matrix is
  ## O(p). For p in the thousands this difference is the dominant cost of
  ## the algorithm.
  generate_d=function(sig_lst, p_vec, iso.err){
    d=list()
    if(iso.err){
      for(k in 1:K){
        d[[k]]=rep(sig_lst[[k]],p_vec[k])
      }
    } else {
      for(k in 1:K){
        d[[k]] = sig_lst[[k]]
      }
    }

    return(unlist(d))
  }

  ## Efficient observed-data log-likelihood using the matrix determinant lemma
  ## and Woodbury identity. Avoids ever forming the dense pxp covariance
  ## C = W W^T + D (or its Cholesky/inverse), which is O(p^3) per call and
  ## infeasible for large p (e.g., p in the thousands). Instead all work is
  ## done with the low-rank Wp x q loadings and the diagonal noise vector,
  ## which is O(N*p*q + p*q^2).
  obs_LogLik<-function(Yc, w, d_vec){
    ##############################################
    # input:    -Yc    : an nxp matrix of centered observations
    #           -w     : the pxq stacked loading matrix
    #           -d_vec : length-p vector of (diagonal) noise variances
    #
    # output:   a real value of the log likelihood
    ##############################################
    N=nrow(Yc)
    p=ncol(Yc)
    q=ncol(w)

    A = w/d_vec                      # solve(D)%*%w, done row-wise (D diagonal)
    YcA = Yc%*%A                     # N x q
    SA = crossprod(Yc, YcA)/N        # S%*%A without ever forming S (p x p)
    M = diag(q) + crossprod(w, A)    # I_q + W^T D^-1 W
    c_solv = solve(M)

    logdetC = sum(log(d_vec)) + as.numeric(determinant(M, logarithm = TRUE)$modulus)
    trace_term = sum(colMeans(Yc^2)/d_vec) - sum(diag(c_solv%*%crossprod(A, SA)))

    LogLik = -(N/2)*(p*log(2*pi) + logdetC) - (N/2)*trace_term

    return(LogLik)
  }

  ## Efficient complete-data log-likelihood: vectorized over subjects (no
  ## per-observation loop), and log|D| computed as sum(log(diag)) instead of
  ## det() on a dense pxp matrix (also O(p^3) otherwise).
  complete_LogLik<-function(Y, theta, mu, w, d_vec){
    ##############################################
    # input:    -mu    : length-p mean vector
    #           -w     : the pxq stacked loading matrix
    #           -d_vec : length-p vector of (diagonal) noise variances
    #           -theta : nxq matrix of (expected) subject scores
    #           -Y     : an nxp data matrix as the observations
    #
    # output:   a real value of the log likelihood
    ##############################################

    N=dim(Y)[1]

    constant.term = -(N/2)*(sum(c(ncol(Y),ncol(w)))*log(2*pi) + sum(log(d_vec)))

    Resid = sweep(Y, 2, mu) - theta%*%t(w)
    kernel.term = sum(sweep(Resid^2, 2, d_vec, "/")) + sum(theta^2)

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
        L=t(chol(stats::cov(Y[,1:P[k]])))
      } else if(k>1){
        L=t(chol(stats::cov(Y[,(sum(P[1:(k-1)])+1):sum(P[1:k])])))
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
    if(K==2){
      ajive.solution = CJIVE::cc.jive(dat.blocks, signal.ranks = Q[1]+Q[-1],
                               joint.rank = Q[1], perm.test = FALSE)

      WJ = lapply(ajive.solution$sJIVE$joint_matrices, function(x) x$v)
      WI = lapply(ajive.solution$sJIVE$indiv_matrices, function(x) x$v)
    } else if(K>2){
      stop("Please provide inital loading estimates from the AJIVE package, available at https://github.com/idc9/r_jive.")
    }
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
    sig_hat = rep(NA, K)
    for(k in 1:K){
      temp.Y = Y[,(k>1)*sum(P[1:(k-1)])+(1:P[k])]
      sig_hat[k] = mean(svd(temp.Y)$d[-(1:Q.tot[k])]^2)
    }
    rm(temp.Y)
  }

  # Total noise matrices D
  if(isotropic.error){
    sig_hat = as.numeric(sig_hat)
    d_hat=generate_d(sig_hat,P,TRUE)
  } else{
    sig_lst = list()
    for(k in 1:K){
      sig_lst[[k]] = sig_hat[k] + stats::rnorm(P[k], sd = 0.1*sig_hat[k])
    }
    d_hat=generate_d(sig_lst,P,FALSE)
  }

  # Initiate Chordal Distance
  chord.dist = NULL
  mu_hat=apply(as.matrix(Y),2,sum)/N

  Iq=diag(sum(Q))

  # Data enter the algorithm only through Yc (centered data). mu_hat never
  # changes inside the EM loop, so Yc (and any block-index bookkeeping) is
  # computed once here rather than on every iteration.
  Yc=sweep(Y, 2, mu_hat)
  cumP = c(0, cumsum(P))
  block_idx = lapply(1:K, function(k) (cumP[k]+1):cumP[k+1])

  # d_hat is a length-p vector of the diagonal of D. Since D is diagonal,
  # solve(D)%*%M / D^{-1} is just M divided row-wise by d_hat -- no pxp
  # matrix (or its inverse/determinant) is ever formed.
  Dinv = function(M) M/d_hat

  c_solv=solve(Iq+crossprod(w_hat, Dinv(w_hat)))
  exp.theta =  Yc%*%Dinv(w_hat)%*%c_solv

  all_obs.LogLik=obs_LogLik(Yc, w_hat, d_hat)
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
    w=w_hat
    d=d_hat

    ## E-step conditional moments. S = t(Yc)%*%Yc/N (the pxp sample
    ## covariance) is never formed: any product S%*%X for a p x q matrix X
    ## is computed as crossprod(Yc, Yc%*%X)/N, which costs O(N*p*q) instead
    ## of O(p^2) to store S plus O(p^2*q) to multiply it -- this is what
    ## made the algorithm infeasible for datasets with thousands of
    ## features.
    A_mat = Dinv(w)                                  # solve(D)%*%w
    SA = crossprod(Yc, Yc%*%A_mat)/N                  # S%*%solve(D)%*%w
    c_solv=solve(Iq+crossprod(w, A_mat))
    U=SA%*%c_solv
    V=c_solv+c_solv%*%crossprod(A_mat, SA)%*%c_solv

    wk_hat_old = wk_hat

    ## Update wk_hat (M-step for W; does not depend on D, see manuscript Eq. 9)
    for(k in 1:K){
      wk_hat[[k]]=A[[k]]%*%U%*%t(B[[k]])%*%solve(B[[k]]%*%V%*%t(B[[k]]))
    }
    w_hat=wk_to_w(wk_hat, P, Q)
    chord.dist = c(chord.dist, CJIVE::chord.norm.diff(w, w_hat))

    ## Update sigma_hat using the *updated* w_hat, matching the closed-form
    ## solution in the manuscript (Eq. 9), which plugs the new W_k into the
    ## D_k formula. Only the diagonal of
    ##   d_tild = S - 2*w_hat%*%t(U) + w_hat%*%V%*%t(w_hat)
    ## is ever needed, so it is computed directly as a length-p vector
    ## instead of forming the full pxp matrix.
    diag_S = colMeans(Yc^2)
    diag_d_tild = diag_S - 2*rowSums(w_hat*U) + rowSums((w_hat%*%V)*w_hat)

    if(isotropic.error){
      for(k in 1:K){
        sig_hat[k]=mean(diag_d_tild[block_idx[[k]]])
      }
    } else{
      sig_hat = list()
      for(k in 1:K){
        sig_hat[[k]]=diag_d_tild[block_idx[[k]]]
      }
    }

    ## Update d_hat
    d_hat=generate_d(sig_hat,P,isotropic.error)

    ################## End of EM-ALGORITHM ######################
    iter=iter + 1

    # Compute subject scores
    c_solv=chol2inv(chol(Iq+crossprod(w_hat, Dinv(w_hat))))
    exp.theta = Yc%*%Dinv(w_hat)%*%c_solv

    all_obs.LogLik=append(all_obs.LogLik, obs_LogLik(Yc, w_hat, d_hat))
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
      graphics::layout(matrix(1:3, nrow = 1))
      plot(all_obs.LogLik, ylab = "Log-Likelihood", main = "Observed Data Log-Likelihood")
    } else {
      graphics::layout(matrix(1:2, nrow = 1))
    }
    if(is.finite(min(all_complete.LogLik))) {plot(all_complete.LogLik, ylab = "Log-Likelihood",
                                                  main = "Complete Data Log-Likelihood")}
    plot(chord.dist, ylab = "Chordal Norm", main = "Distance between consecutive \n estimates of 'W'")
    graphics::layout(1)
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
