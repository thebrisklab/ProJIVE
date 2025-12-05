#' Generates K Simulated Datasets that follow JIVE Model using binary subject scores
#'
#' @param n Sample size.
#' @param p Number of features per block (input as a vector).
#' @param JntVarEx Proportions of joint variation explained (input as a vector).
#' @param IndVarEx Proportions of individual variation explained (input as a vector).
#' @param jnt_rank Joint Rank.
#' @param equal.eig Logical (TRUE/FALSE), which allows the user to specify whether components within a data-blocks joint (or individual) signal should be equally waited. The default is FALSE.
#' @param ind_ranks Individual rank for each block (input as a vector).
#' @param JntVarAdj Logical (TRUE/FALSE): Specify whether signal matrices should be weighted to achieve the desired proportions of variation attributable to the joint signal.
#' @param mix.probs A numerical vector that specifies the number of individual components (ranks) for each data block. Each entry describes the probability in the Gaussian Mixture model for the score. Entries need to add to 1.
#' @param SVD.plots Logical (TRUE/FALSE): Should plots of signal matrix singular values be produced to verify ranks?
#' @param Error Logical (TRUE/FALSE): Should the data be noise-contaminated?
#' @param print.cor logical (TRUE/FALSE), Print the correlation matrix for the scores? (Allows one to assess orthoganility between scores/compnents)
#' @param Loads char: Toy data can be generated with loadings from 'Gaussian', 'Rademacher', or 'Double_Exp'  (double exponential) distributions. Loadings can also be fixed at binary (0/1) values assigned to half of the variables with keyword 'Fixed'.
#' @param Scores char: Joint subject scores can be randomly generated from 'Gaussian', 'Binomial', or 'Gaussian_Mixture' distributions. The last refers to a  mixture of Gaussians with unit variance, where 20% have mean -4, 50% have mean 0, and 30% have mean 4. In all cases individual scores are standard Gaussian.
#' @param error.variances A numeric vector specifying the variance of the noise for each dataset.
#'
#' @return K simulated datasets that follow JIVE Model.
#' Contains four main elements:
#' - Data Components: A list containing three sub-lists (each of length K), JointSignalMatrices, IndivSignalMatrices, and NoiseMatrices.
#' - Data Blocks: A list of length K, where each element is the final observed dataset (joint + individual + noise) of dimension n x p(k).
#' - Scores: A list with two matrices, Joint (n x jnt_rank) and Indiv (n x sum(ind_ranks)).
#' - Loadings: A list containing Joint and Indiv loadings.
#' @export
#'
#' @examples
#'rep_number = 1
#'r.J = 3
#'r.I1 = 2
#'r.I2 = 2
#'n = 1000
#'p1 = 20
#'p2 = 200
#'JntVarEx1 = 0.1
#'JntVarEx2 = 0.1
#'IndVarEx1 = 0.25
#'IndVarEx2 = 0.25
#'nparams = p1*(r.J+r.I1)+p2*(r.J+r.I2)+2
#'prop = n/nparams
#'JntScores = matrix(rdunif(n*r.J, 0, 5), nrow = n, ncol = r.J)
#'IndivScore.X = matrix(rweibull(n*r.I1, shape = 1), nrow = n)
#'IndivScore.Y = matrix(rhnorm(n*r.I2), nrow = n)
#'Scores = cbind(JntScores, IndivScore.X, IndivScore.Y)
#'true_signal_ranks = r.J + c(r.I1,r.I2)
#'ToyDat = GenerateToyData(n = n, p = c(p1, p2), JntVarEx = c(JntVarEx1, JntVarEx2), equal.eig = F, IndVarEx = c(IndVarEx1, IndVarEx2), jnt_rank = r.J, ind_ranks = c(r.I1, r.I2), JntVarAdj = F, SVD.plots = F, Error = T, print.cor = TRUE, Loads = "Rademacher", Scores = "Gaussian_Mixture", error.variances = c(1,1))
GenerateToyData <- function(n, p, JntVarEx, IndVarEx, jnt_rank = 1, equal.eig = FALSE, ind_ranks, JntVarAdj = TRUE, mix.probs = NULL,
                            SVD.plots = TRUE, Error = TRUE, print.cor = TRUE, Loads = "Rademacher", Scores = "Gaussian_Mixture",
                            error.variances = NULL){

  r.J = jnt_rank
  r.I = ind_ranks
  K = length(p)
  Scores.text = ifelse(is.numeric(Scores), "as defined by user", paste("from", Scores, "distributions"))
  Loads.text = ifelse(is.list(Loads), "as defined by user.", paste("from", Loads, "distributions."))
  cat(paste0("Generating Scores ", Scores.text, " and Loadings ", Loads.text, ". \n"))

  if(is.numeric(Scores)){
    JntScores = Scores[,1:r.J, drop = FALSE]
    IndivScores = Scores[,-(1:r.J), drop = FALSE]
  } else if(Scores=="Binomial"){
    JntScores = matrix(MCMCpack::rbinom(n*r.J, size=1, prob=0.2), nrow = n, ncol = r.J)

    b = MCMCpack::rbinom(n*sum(r.I), size=1, prob=0.4)
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
  } else if (is.numeric(Scores)){
    if(nrow(Scores == n) & ncol(Scores = r.J+sum(r.I))){
      JntScores = Scores[,1:r.J]
      IndivScores = Scores[,-(1:r.J)]
    } else{
      cat("Matrix of subject scores must have n rows and jnt_rank+sum(ind_ranks) columns.")
    }
  } else {
    message("Please use one of the following three options to generate data with corresponding subject scores: 1) 'Gaussian' 2) 'Gaussian_Mixture' 3) Binomial.\n")
    message("Alternatively, the user can enter a matrix of subject scores with n rows and jnt_rank+sum(ind_ranks) columns.\n")
  }

  colnames(JntScores) = paste("Jnt Score", 1:r.J)
  IndivScores.names = NULL

  for(k in 1:K){
    IndivScores.names = c(IndivScores.names, paste0("Indiv X", k, " Score ", 1:r.I[k]))
  }
  colnames(IndivScores) = IndivScores.names

  if(print.cor){cat("The correlation between subject scores is given by"); print(round(cor(cbind(JntScores, IndivScores)),4))}

  Jnt.Loads.All = list()
  Indiv.Loads.All = list()
  D.I = list()
  Noise = Joint.Sigs = Indiv.Sigs = list()
  Sig.Mats =  Data.Mats = list()
  temp.fcn = function(x){x[sample(round(length(x)/2))] = 1; x}

  if(is.null(error.variances)){
    error.variances = rep(1,K)
  } else if(is.numeric(error.variances) & length(error.variances) != K){
    cat("Number of entries in 'error.variances' must equal number of data sets and thus have the same length as 'p'.\n")
  }

  for(k in 1:K){
    if(length(Loads)==1){
      if(Loads == "Gaussian"){
        Jnt.Loads.All[[k]] = matrix(rnorm(r.J*p[k]), nrow = r.J, ncol = p[k])
      } else if (Loads == "Fixed"){
        Jnt.Loads.All[[k]] = matrix(apply(matrix(0, nrow = r.J, ncol = p[k]), 1, temp.fcn), nrow = r.J)
      } else if (Loads == "Double_Exp"){
        Jnt.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.J)), nrow = r.J)
      } else if (Loads == "Rademacher"){
        Jnt.Loads.All[[k]] = matrix(rsign(p[k]*(r.J)), nrow = r.J)
      }
    } else if (is.list(Loads)){
      Jnt.Loads.All[[k]] = t(Loads[[1]][[k]])
    }
    D.J = (1 - equal.eig)*diag(r.J:1) + equal.eig*diag(rep(1,r.J))

    # Joint.Sigs[[k]] = JntScores%*%sqrt(D.J[[k]])%*%Jnt.Loads.All[[k]]
    Joint.Sigs[[k]] = JntScores%*%D.J%*%Jnt.Loads.All[[k]]

    if(SVD.plots){
      plot(svd(Joint.Sigs[[k]])$d, ylab = "Singular Values")
      title(paste0("SVD of Joint Signal from X", k))
    }
    temp.IndScores = IndivScores[,(k>1)*sum(r.I[1:(k-1)])+(1:r.I[k])]

    if(length(Loads)==1){
      if(Loads == "Gaussian"){
        Indiv.Loads.All[[k]] = matrix(rnorm(n = p[k]*r.I[k]), nrow = r.I[k], ncol = p[k])
      } else if (Loads == "Fixed"){
        temp.fcn = function(x){x[sample(round(length(x)/4))] = 1; x}
        Indiv.Loads.All[[k]] = matrix(apply(matrix(-1, nrow = r.I[k], ncol = p[k]), 1, temp.fcn), nrow = r.I[k])
      } else if (Loads == "Double_Exp"){
        Indiv.Loads.All[[k]] = matrix(rlaplace(p[k]*(r.I[k])), nrow = r.I[k])
      } else if (Loads == "Rademacher"){
        Indiv.Loads.All[[k]] = matrix(extraDistr::rsign(p[k]*(r.I[k])), nrow = r.I[k])
      }
    } else if (is.list(Loads)){
      Indiv.Loads.All[[k]] = t(Loads[[2]][[k]])
    }

    D.I[[k]] = (1 - equal.eig)*diag(r.I[k]:1) + equal.eig*diag(rep(1,r.I[k]))

    Indiv.Sigs[[k]] = temp.IndScores%*%D.I[[k]]%*%Indiv.Loads.All[[k]]

    if(SVD.plots){
      plot(svd(Indiv.Sigs[[k]])$d, ylab = "Singular Values")
      title(paste0("SVD of Individual Signal from X", k))
    }

    Sig.Mats[[k]] = Joint.Sigs[[k]] + Indiv.Sigs[[k]]

    Noise[[k]] = matrix(rnorm(n*p[k], sd = sqrt(error.variances[k])), nrow  = n)

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

  Loadings = list(list(),list())
  names(Loadings) = c("Joint", "Indiv")
  for(k in 1:K){
    Loadings[["Joint"]][[k]] = D.J%*%Jnt.Loads.All[[k]]
    Loadings[["Indiv"]][[k]] = D.I[[k]]%*%Indiv.Loads.All[[k]]
  }

  out = list(Dat.Comps, Blocks, Scores, Loadings)
  names(out) = c("Data Components", "Data Blocks", "Scores", "Loadings")

  return(out)
}
