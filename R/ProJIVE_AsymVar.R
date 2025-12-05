#' Calculate empirical observed information matrix
#'
#' @param W.mats A list of matrices that each represents the variable loadings of the dataset. Rows as variables and columns as latent components.
#' @param error.vars Variance of the noise for each dataset.
#' @param theta A matrix that represents the subject scores.
#' @param r.J Integer. The number of joint components in the PJIVE decomposition.
#' @param Y A matrix that represents the observed data. Rows as subjects and columns as features from all datasets combined.
#'
#' @return A list containing:
#'            - MeanScoreVector: The mean of the score vectors across subjects.
#'            - ObservedEmpiricalInformationMatrix: The empirical observed information matrix.
#'            - Inverse_ObservedEmpiricalInformationMatrix: The inverse of the observed information matrix.
#' @export
#'
#' @examples
#' W1 <- matrix(rnorm(50), 10, 5)
#' W2 <- matrix(rnorm(60), 12, 5)
#' W.mats <- list(W1, W2)
#' error.vars <- c(0.5, 0.8)
#' theta <- matrix(rnorm(100), 20, 5)
#' r.J <- 2
#' Y <- matrix(rnorm(400), 20, 20)
#' result <- ProJIVE_AsymVar(W.mats, error.vars, theta, r.J, Y)
#' str(result)
ProJIVE_AsymVar<-function(W.mats, error.vars, theta, r.J, Y){
  K1 = length(W.mats)
  K2 = length(error.vars)
  if(K1!=K2){stop("W matrices and error variances must have the same length.")}

  K = K1; rm(K1, K2)
  r = ncol(theta) ##r_J + r_I1 + ... + r_I2
  r_k = sapply(W.mats, ncol)
  p_k = sapply(W.mats, nrow)
  r_Ik = r_k - r.J
  n = nrow(Y)

  theta_k = list()
  theta_k[[1]] = theta[,1:r_k[1]]

  Y_k = list()
  Y_k[[1]] = Y[,1:p_k[1]]

  for(k in 2:K){
    theta_k[[k]] = theta[,c(1:r.J, r.J + cumsum(r_Ik)[k-1] + (1:r_Ik[k]))]
    Y_k[[k]] = Y[,cumsum(p_k)[k-1]+(1:p_k[k])]
  }

  score_i = list()
  for(i in 1:nrow(Y)){
    score_wk = score_sigmak = NULL
    for(k in 1:K){
      epsilon_ik = (Y_k[[k]][i,] - W.mats[[k]]%*%t(t(theta_k[[k]][i,])))
      score_wk = c(score_wk, c(epsilon_ik%*%t(theta_k[[k]][i,])/error.vars[k]))
      score_sigmak = c(score_sigmak, (t(epsilon_ik)%*%epsilon_ik/error.vars[k] - p_k[k])/(2*error.vars[k]))
    }
    score_i[[i]] = matrix(c(score_wk, score_sigmak),ncol = 1)
  }
  mean.score = info.mat = 0
  for(i in 1:n){
    mean.score = mean.score + score_i[[i]]/nrow(Y)
  }
  for(i in 1:n){
    info.mat = info.mat + (score_i[[i]]-mean.score)%*%t(score_i[[i]]-mean.score)
  }
  out = list(mean.score, info.mat)
  names(out) = c("MeanScoreVector", "ObservedEmpericalInformationMatrix")
  tryCatch(
    {
      inv.info = Matrix::chol2inv(Matrix::chol(info.mat))
      out = list(mean.score, info.mat, inv.info)
      names(out) = c("MeanScoreVector", "ObservedEmpericalInformationMatrix", "Inverse_ObservedEmpericalInformationMatrix")
    },
    error = function(cond){
      message("Warning: Information matrix is not positive definite.
              Generalized inverse used to estimate covariance matrix.")
      inv.info = MASS::ginv(info.mat)
      inv.info
    }
  )
  return(out)
}
