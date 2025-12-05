###########################################################################################################################
#######   Conduct D-CCCA analysis on author's (Shu) Github page using the R code for   ####################################
#######           D-CCCA that I adapted from Shu's Python version on their github      ####################################
#######                 Author: Raphiel J. Murden                                      ####################################
#######                 Supervised by Benjamin Risk                                    ####################################
###########################################################################################################################
############  Notes: 
############ 05DEC2022: D-CCA was conducted using the Python version on the RSPH cluster. Results were saved as .csv files
############            This script imports those results, conducts D-CCA on the example data provided therein, and then 
############            compares the imported results to those generated via R
###########################################################################################################################

require(dplyr); require(mgcv); require(r.jive); require(ggplot2); require(psych); require(r.jive); require(stringr)
set.seed(0)
prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
dat.dir = "C:/Users/rmurden/Downloads/D-CCA-master"
proc.dat.dir = "H:/My Documents/P-JIVE/Data/Processed"
py.res.dir = "C:/Users/rmurden/Downloads/D-CCA-master/Output/"

source(file.path(prog.dir, "Functions_for_PJIVE.R"))
source(file.path(prog.dir, "Functions_for_CJIVE.R"))

################ Import data and cluster results
dat.1 = as.matrix(read.delim(file = file.path(dat.dir, "Y_1.txt"), header = FALSE, sep = " "), nrow = 900, ncol = 300)
dat.2 = as.matrix(read.delim(file = file.path(dat.dir, "Y_2.txt"), header = FALSE, sep = " "), nrow = 300, ncol = 300)

X1.py = as.matrix(read.csv(file.path(py.res.dir, "X1_hat.csv"), header = FALSE))
X2.py = as.matrix(read.csv(file.path(py.res.dir, "X2_hat.csv"), header = FALSE))

C1.py = as.matrix(read.csv(file.path(py.res.dir, "C1_hat.csv"), header = FALSE))
C2.py = as.matrix(read.csv(file.path(py.res.dir, "C2_hat.csv"), header = FALSE))

D1.py = as.matrix(read.csv(file.path(py.res.dir, "D1_hat.csv"), header = FALSE))
D2.py = as.matrix(read.csv(file.path(py.res.dir, "D2_hat.csv"), header = FALSE))


################ Conduct D-CCA in R
r_1 = 3; r_2 = 5; r_12 = 1 
res.R = dCCA(dat.1, dat.2, r_1 = r_1,  r_2 = r_2, r_12 = r_12)

###Signal matrices
X1.R = res.R$SignalMatrix_X1
norm(X1.R - X1.py, "F")
hist(c(abs(X1.R - X1.py)), breaks = 1000, xlab = "Abs. Diff. Between X1", main = paste0("Frobenius Norm = ", round(norm(X1.R - X1.py, "F"), 2)))

X2.R = res.R$SignalMatrix_X2
norm(X2.R - X2.py, "F")
hist(c(abs(X2.R - X2.py)), breaks = 1000, xlab = "Abs. Diff. Between X2", main = paste0("Frobenius Norm = ", round(norm(X2.R - X2.py, "F"), 2)))

## Common Signal matrices
C1.R = res.R$Common_SignalMatrix_X1
norm(C1.R - C1.py, "F")
hist(c(abs(C1.R - C1.py)), breaks = 1000, xlab = "Abs. Diff. Between C1", main = paste0("Frobenius Norm = ", round(norm(C1.R - C1.py, "F"), 2)))

C2.R = res.R$Common_SignalMatrix_X2
norm(C2.R - C2.py, "F")
hist(c(abs(C2.R - C2.py)), breaks = 1000, xlab = "Abs. Diff. Between C2", main = paste0("Frobenius Norm = ", round(norm(C2.R - C2.py, "F"), 2)))

## DistinDt Signal matriDes
D1.R = res.R$Distinct_SignalMatrix_X1
norm(D1.R - D1.py, "F")
hist(c(abs(D1.R - D1.py)), breaks = 1000, xlab = "Abs. Diff. Between D1", main = paste0("Frobenius Norm = ", round(norm(D1.R - D1.py, "F"), 2)))

D2.R = res.R$Distinct_SignalMatrix_X2
norm(D2.R - D2.py, "F")
hist(c(abs(D2.R - D2.py)), breaks = 1000, xlab = "Abs. Diff. Between D2", main = paste0("Frobenius Norm = ", round(norm(D2.R - D2.py, "F"), 2)))

