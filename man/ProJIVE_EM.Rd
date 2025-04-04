% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProJIVE_EM.R
\name{ProJIVE_EM}
\alias{ProJIVE_EM}
\title{PJIVE ML estimation of JIVE model that uses EM algorithm}
\usage{
ProJIVE_EM(
  Y,
  P,
  Q,
  Max.iter = 10000,
  diff.tol = 1e-08,
  plots = TRUE,
  chord.tol = -1,
  sig_hat = NULL,
  init.loads = NULL
)
}
\arguments{
\item{Y}{Matrix. The observed data.}

\item{P}{Integer vector. The number of variables in each dataset.}

\item{Q}{Integer vector. The number of components in each dataset.}

\item{Max.iter}{The maximum number of iterations for the EM algorithm (default = 1000).}

\item{diff.tol}{Numeric. Specifies the tolerance for convergence based on log-likelihood (default = 1e-8).}

\item{plots}{Logical (TRUE/FALSE): Specifies if plots are generated (default = TRUE).}

\item{chord.tol}{Numeric. Tolerance for the chordal distance for convergence (default = -1).}

\item{sig_hat}{The error variances for each dataset.}

\item{init.loads}{Either a list of initial loadings or a string specifying the initialization method ("AJIVE" or "CJIVE").}
}
\value{
A list containing:
- SubjectScoreMatrix: The matrix of subject scores.
- LoadingMatrix: The final loading matrix.
- Ranks: The ranks for joint and individual components.
- VarianceExplained: The variance explained by the joint and individual components for each dataset.
- ErrorVariances: The estimated error variances.
- ChordalDistances: The distance between consecutive estimates of the loading matrix.
- Complete-Data-Log-Likelihood: The log-likelihood for the complete data.
- Observed-Data-Log-Likelihood: The log-likelihood for the observed data.
- BIC: The Bayesian Information Criterion for the model.
- AIC: The Akaike Information Criterion for the model.
}
\description{
PJIVE ML estimation of JIVE model that uses EM algorithm
}
\examples{
Y <- matrix(rnorm(100), 20, 5)
P <- c(5, 5)
Q <- c(2, 3)
result <- ProJIVE_EM(Y, P, Q)
str(result)
}
