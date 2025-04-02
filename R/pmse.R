#' Adaptation of PMSE function from singR Package
#'
#' @param M1 A numeric matrix (d x p) representing an estimated structure.
#' @param M2 A numeric matrix (q x p) representing a reference structure.
#' @param S1 A numeric matrix (n x p) representing estimated subject scores.
#' @param S2 A numeric matrix (n x q) representing reference subject scores.
#' @param standardize Logical (TRUE/FALSE): If TRUE, standardizes matrices to unit variance before computing PMSE (default = FALSE).
#'
#' @return A numeric value. The Procrustes mean squared error.
#' @export
#'
#' @examples
#' M1 <- matrix(rnorm(20), 5, 4)
#' M2 <- matrix(rnorm(20), 5, 4)
#' pmse(M1 = M1, M2 = M2)
#'
#' S1 <- matrix(rnorm(30), 10, 3)
#' S2 <- matrix(rnorm(30), 10, 3)
#' pmse(S1 = S1, S2 = S2)
pmse<-function (M1 = NULL, M2 = NULL, S1 = NULL, S2 = NULL, standardize = FALSE){
  tfun = function(x) all(x == 0)
  if (is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2))
    stop("need to supply either M1 and M2 or S1 and S2")
  if (!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if (!is.null(M1) && nrow(M1) > ncol(M1)){
    stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")
  }
  if (is.null(M1)) {
    if(is.vector(S1)) {S1 = matrix(S1, ncol = 1)}
    if(is.vector(S2)) {S2 = matrix(S2, ncol = 1)}

    nS = nrow(S1)
    if (nS != nrow(S2))
      stop("S1 and S2 must have the same number of rows")
    if (sum(apply(S1, 2, tfun)) + sum(apply(S2, 2, tfun)))
      stop("pmse not defined when S1 or S2 has a column of all zeros")
    if (standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if (p < q) {
      S1 = cbind(S1, matrix(0, nS, (q - p)))
    }
    if (q < p) {
      S2 = cbind(S2, matrix(0, nS, (p - q)))
    }
    Stemp = matchICA.2(S = S1, template = S2)
    n.comp = max(q, p)
    indices = c(1:n.comp)[!(apply(Stemp, 2, tfun) | apply(S2,
                                                          2, tfun))]
    return(sqrt(sum((Stemp[, indices] - S2[, indices])^2))/sqrt(nS *
                                                                  min(p, q)))
  }
  else {
    if (sum(apply(M1, 1, tfun)) + sum(apply(M2, 1, tfun)))
      stop("pmse not defined when M1 or M2 has a row of all zeros")
    if (standardize) {
      temp = diag((diag(M1 %*% t(M1)))^(-1/2))
      M1 = temp %*% M1
      temp = diag((diag(M2 %*% t(M2)))^(-1/2))
      M2 = temp %*% M2
    }
    p = ncol(M1)
    if (p != ncol(M2))
      stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp = max(d, q)
    if (n.comp > p)
      warning("M should be d x p")
    if (d < q) {
      M1 = rbind(M1, matrix(0, (q - d), p))
    }
    if (q < d) {
      M2 = rbind(M2, matrix(0, (d - q), p))
    }
    l2.mat1 = l2.mat2 = matrix(NA, nrow = n.comp, ncol = n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        l2.mat1[i, j] = sum((M2[i, ] - M1[j, ])^2)
        l2.mat2[i, j] = sum((M2[i, ] + M1[j, ])^2)
      }
    }
    l2.mat1 = sqrt(l2.mat1)
    l2.mat2 = sqrt(l2.mat2)
    l2.mat = l2.mat1 * (l2.mat1 <= l2.mat2) + l2.mat2 * (l2.mat2 <
                                                           l2.mat1)
    map = as.vector(solve_LSAP(l2.mat))
    l2.1 = ifelse(nrow(l2.mat1)==1 && ncol(l2.mat1)==1 && l2.mat1[1,1] == 0, 0, diag(l2.mat1[, map]))
    l2.2 = ifelse(nrow(l2.mat2)==1 && ncol(l2.mat2)==1 && l2.mat2[1,1] == 0, 0, diag(l2.mat2[, map]))
    sign.change = -1 * (l2.2 < l2.1) + 1 * (l2.1 <= l2.2)
    perm = diag(n.comp)[, map] %*% diag(sign.change)
    M.perm = t(perm) %*% M1
    indices = c(1:n.comp)[!(apply(M.perm, 1, tfun) | apply(M2,
                                                           1, tfun))]
    return(sqrt(sum((M.perm[indices, ] - M2[indices, ])^2))/sqrt(p *
                                                                   min(d, q)))
  }
}
