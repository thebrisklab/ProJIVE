#' Adjust Dataset Components to get Desired R^2 Values
#'
#' @param J The joint component of the dataset (a matrix).
#' @param I The individual component of the dataset (a matrix).
#' @param E The noise component of the dataset (a matrix).
#' @param JntVarEx Proportions of joint variation explained (input as a vector).
#' @param IndVarEx Proportions of individual variation explained (input as a vector).
#'
#' @return A list containing:
#'            - J: The adjusted joint component matrix.
#'            - I: The adjusted individual component matrix.
#'            - Data: The new dataset matrix after adjustment.
#' @export
#'
#' @examples
#' J <- matrix(rnorm(100), 10, 10)
#' I <- matrix(rnorm(100), 10, 10)
#' E <- matrix(rnorm(100), 10, 10)
#' result <- AdjSigVarExp(J, I, E, JntVarEx = 0.3, IndVarEx = 0.5)
#' str(result)
AdjSigVarExp <-function(J, I, E, JntVarEx, IndVarEx){
  simul.quads = function(x, parms){
    JJ = parms[1]
    II = parms[2]
    EE = parms[3]
    JE = parms[4]
    IE = parms[5]
    R_J = parms[6]
    R_I = parms[7]

    y1 = x[1]^2*II*(1 - R_I) - 2*x[1]*IE*R_I - R_I*(x[2]^2*JJ + 2*x[2]*JE + EE)
    y2 = x[2]^2*JJ*(1 - R_J) - 2*x[2]*JE*R_J - R_J*(x[1]^2*II + 2*x[1]*IE + EE)

    y = c(y1,y2)
    return(y)
  }

  JJ = CJIVE::MatVar(J)
  II = CJIVE::MatVar(I)
  EE = CJIVE::MatVar(E)
  JE = sum(diag(J%*%t(E)))
  IE = sum(diag(I%*%t(E)))
  R_J = JntVarEx
  R_I = IndVarEx

  parms = c(JJ, II, EE, JE, IE, R_J, R_I)

  A = J + I
  AA = CJIVE::MatVar(A)
  EE = CJIVE::MatVar(E)
  AE = sum(diag(A%*%t(E)))

  d0 = IndVarEx + JntVarEx
  a = AA*(1 - d0)
  b = -2*AE
  c = -d0*EE

  d.A = (-b+sqrt(b^2 - 4*a*c))/(2*a)

  start = c(0.5, 0.5)*d.A

  roots = rootSolve::multiroot(simul.quads, start, parms = parms)
  c = roots$root[1]
  d = roots$root[2]

  J.d = d*J; I.c = c*I;
  Dat = J.d + I.c + E

  res = list(J.d, I.c, Dat)
  names(res) = c("J", "I", "Data")
  return(res)
}
