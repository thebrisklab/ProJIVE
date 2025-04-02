#' Bootstrap
#'
#' @param B Integer. The number of bootstrap samples desired to generate (default = 50).
#' @param P A integer vector that contains the number of variables in each dataset.
#' @param Q A integer vector that contains the number of components in each dataset.
#' @param theta.hat Matrix. The estimated subject scores from the ProJIVE model.
#' @param W.hat Matrix. The estimated loading matrix from the ProJIVE model.
#' @param error.vars A numeric vector containing the estimated noise variances for each dataset.
#'
#' @return A list containing:
#'            - Bootstrap_Covariance: The covariance matrix of the bootstrapped loadings.
#'            - Bootstrap_StdErrs: The standard errors of the loadings computed from bootstrap variance.
#' @export
#'
#' @examples
#' set.seed(123)
#' P <- c(10, 15)
#' Q <- c(3, 4)
#' theta.hat <- matrix(rnorm(200), 20, 10)
#' W.hat <- matrix(rnorm(100), 10, 10)
#' error.vars <- c(0.5, 0.7)
#' result <- ProJIVE_BootsratVar(B, P, Q, theta.hat, W.hat, error.vars)
#' str(result)
ProJIVE_BootsratVar = function(B = 50, P, Q, theta.hat, W.hat, error.vars){
  nobs = nrow(theta.hat)
  bstrap.mean = 0
  bstraps = list()
  for(b in 1:B){
    errors = NULL
    for(k in 1:length(P)){
      errors = cbind(errors, matrix(sqrt(error.vars[k])*rnorm(n=nobs*(P[k])),nrow=nobs))
    }

    tempdata = theta.hat%*%t(W.hat)+errors

    temp.res = ProJIVE_EM(Y=tempdata, P=P, Q=Q, Max.iter=10000, diff.tol=1e-7, sig_hat = "MLE",
                          init.loads = "CJIVE", plots = FALSE)#, verbose = FALSE)

    Procrustes.Loading = MCMCpack::procrustes(X = temp.res$LoadingMatrix, Xstar = W.hat)$X.new
    bstrap.mean = bstrap.mean + temp.res$LoadingMatrix/B
    bstraps[[b]] = temp.res$LoadingMatrix
  }

  bstrap.var = 0
  for(b in 1:B){
    bstrap.var = bstrap.var + matrix(bstraps[[b]] - bstrap.mean, ncol = 1)%*%matrix(bstraps[[b]] - bstrap.mean, nrow = 1)/(B-1)
  }
  round(bstrap.var)

  bstrap.se = sqrt(diag(bstrap.var))
  out = list(bstrap.var, bstrap.se);
  names(out) = c("Boostrap_Covariance", "Bootstrap_StdErrs")
  return(out)
}
