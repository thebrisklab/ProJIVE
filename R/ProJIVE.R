#' Wrapper function to conduct ProJIVE analyses (w option for multiple initial values for loadings), and calculate asymptotic variance
#'
#' @param Y Matrix. The observed data.
#' @param P Integer vector. The number of variables in each dataset.
#' @param Q Integer vector. The number of components in each dataset.
#' @param Max.iter The maximum number of iterations for the EM algorithm (default = 5000).
#' @param diff.tol Numeric. Specifies the tolerance for convergence based on log-likelihood (default = 1e-5).
#' @param plots Logical (TRUE/FALSE): Specifies if plots are generated (default = TRUE).
#' @param chord.tol Numeric. Tolerance for the chordal distance for convergence (default = -1).
#' @param sig_hat The error variances for each dataset.
#' @param init.loads Either a list of initial loadings or a string specifying the initialization method ("AJIVE" or "CJIVE").
#' @param center Logical (TRUE/FALSE): Specifies if the data matrix is centered before estimation (default = FALSE).
#' @param num.starts Integer. The number of different initializations (default = 5).
#' @param AsymVar Logical (TRUE/FALSE): If TRUE, computes the asymptotic variance using the observed empirical information matrix (default = FALSE).
#' @param return.all.starts Logical (TRUE/FALSE): If TRUE, returns results from all initializations; if FALSE, return only the best result (default: FALSE).
#'
#' @return A list containing:
#'            - ProJIVE_Results: The best ProJIVE estimation result (or all results if "return.all.starts = TRUE").
#'            - ObservedDataLogLikelihood: The observed log-likelihood for each initialization.
#'            - ObservedEmpiricalInfo (when "AsymVar = TRUE"): The observed empirical information matrix for asymptotic variance estimation.
#' @export
#'
#' @examples
#' Y <- matrix(rnorm(100), 20, 5)
#' P <- c(5, 5)
#' Q <- c(2, 3)
#' result <- ProJIVE(Y, P, Q)
#' str(result)
ProJIVE<-function(Y, P, Q, Max.iter=5000, diff.tol=1e-5, plots=TRUE,
                  chord.tol=-1, sig_hat=NULL, init.loads = NULL, center = FALSE,
                  num.starts = 5, AsymVar = FALSE, return.all.starts = FALSE){
  ProJIVE.res = list()
  obs.lik = NULL
  for(start in 1:num.starts){
    if(is.list(init.loads) | is.character(init.loads)){
      if(start==1){
        init.loads.in = init.loads
      } else if(start==2){
        init.loads.in = NULL
      } else if(start==3){
        init.loads.in = "CJIVE"
      } else if(start>=4){
        init.loads = list()
        WJ.init = WI.init = list()
        for(ind in 1:length(P)){
          WJ.init[[ind]] = matrix(rnorm(P[ind]*Q[1]), nrow = P[ind])
          WI.init[[ind]] = matrix(rnorm(P[ind]*Q[1+ind]), nrow = P[ind])
        }
        init.loads.in = list(WJ.init, WI.init)
      }
    } else if(!is.list(init.loads)){
      if(start==1){
        init.loads.in = NULL
      } else if(start==2){
        init.loads.in = "CJIVE"
      } else if(start>=3){
        init.loads = list()
        WJ.init = WI.init = list()
        for(ind in 1:length(P)){
          WJ.init[[ind]] = matrix(rnorm(P[ind]*Q[1]), nrow = P[ind])
          WI.init[[ind]] = matrix(rnorm(P[ind]*Q[1+ind]), nrow = P[ind])
        }
        init.loads.in = list(WJ.init, WI.init)
      }
    }

    ProJIVE.res[[start]] = ProJIVE_EM(Y, P, Q, Max.iter, diff.tol, plots,
                                      chord.tol, sig_hat, init.loads.in)
    cat(paste0("Start number ", start, " completed.\n",
               "************************************** \n"))
    obs.lik = c(obs.lik, tail(ProJIVE.res[[start]]$`Observed-Data-Log-Likelihood`, n=1))
  }
  out = list()

  if(return.all.starts){
    out[["ProJIVE_Results"]] = ProJIVE.res
    out[["ObservedDataLogLikelihood"]] = obs.lik
  } else if(!return.all.starts){
    out[["ProJIVE_Results"]] = ProJIVE.res[[which.max(obs.lik)]]
    out[["ObservedDataLogLikelihood"]] = obs.lik
  }

  if(AsymVar){
    PJIVE.res = ProJIVE.res[[which.max(obs.lik)]]
    PJIVE.scores = PJIVE.res$SubjectScoreMatrix
    PJIVE.loads.X = PJIVE.res$LoadingMatrix[1:p1,-(sum(Q):(sum(Q[-3])+1))]
    PJIVE.loads.Y = PJIVE.res$LoadingMatrix[-(1:p1),-(Q[1]+1:Q[2])]
    PJIVE.err.var = PJIVE.res$ErrorVariances

    out[["ObservedEmpireicalInfo"]] = ProJIVE_AsymVar(list(PJIVE.loads.X, PJIVE.loads.Y), PJIVE.err.var, PJIVE.scores, r.J, Y)
  }
  return(out)
}
