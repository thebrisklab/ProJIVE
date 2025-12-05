#############################################################################################################################
### Examine the null distribution of the chordal norm between random pairs of subject scores ###############################
### Author: Raphiel J. Murden                                                               #################################
### Notes:
#############################################################################################################################

rep_number = 0
# NOTE: Change the location of 'prog.dir' to the location where you have saved the file 'Functions_to_SimulateData.R'
prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
prog.gipca.dir = "H:/My Documents/P-JIVE/Programs/GeneralizedIntegrativePCA-master/Functions"
source(file.path(prog.dir, "Functions_for_CJIVE.R"))
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
gipca.files = list.files(prog.gipca.dir,full.names = TRUE)
lapply(gipca.files, source)

r.J = 3
r.I1 = 2
r.I2 = 2
#outdir = args[2]
n = 200
p1 = 30
p2 = 50 ####Note that we have p2 = 1000 here as opposed to p2 = 10,000 in simulations
JntVarEx1 = 0.5
JntVarEx2 = 0.5
#files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25
N = 200
#######JIVE Implementations for Toy Data No 1################################################################################
###Setup input parameters for JIVE implementations
true_signal_ranks = r.J + c(r.I1,r.I2)                          ##true jranks of overall signals

norms = NULL
for(i in 1:N){
  set.seed(i)
  ToyDat = GenerateToyData(n = n, p1 = p1, p2 = p2, JntVarEx1 = JntVarEx1, JntVarEx2 = JntVarEx2, 
                           IndVarEx1 = IndVarEx1, IndVarEx2 =  IndVarEx2, jnt_rank = r.J,
                           equal.eig = F,ind_rank1 = r.I1, ind_rank2 = r.I2, JntVarAdj = T, SVD.plots = F,
                           Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian")
  
  JntScores = ToyDat[['Scores']][['Joint']]
  
  set.seed(2*i)
  ToyDat = GenerateToyData(n = n, p1 = p1, p2 = p2, JntVarEx1 = JntVarEx1, JntVarEx2 = JntVarEx2, 
                           IndVarEx1 = IndVarEx1, IndVarEx2 =  IndVarEx2, jnt_rank = r.J,
                           equal.eig = F,ind_rank1 = r.I1, ind_rank2 = r.I2, JntVarAdj = T, SVD.plots = F,
                           Error = T, print.cor = F, Loads = "Gaussian", Scores = "Gaussian_Mixture")
  
  JntScores2 = ToyDat[['Scores']][['Joint']]
  norms = c(norms, chord.norm.diff(JntScores,JntScores2))
}

hist(norms, breaks = 20, main = paste0("JVE_1=",round(JntVarEx1,2),"; JVE_2=",round(JntVarEx2,2)))
summary(norms)
sd(norms)
