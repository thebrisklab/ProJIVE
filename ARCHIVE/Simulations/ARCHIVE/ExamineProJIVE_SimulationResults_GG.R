#############################################################################################################################
### compare methods of implementing JIVE analysis on simulated datasets                  ####################################
### Author: Raphiel J. Murden                                                            ####################################
### Supervised by Benjamin Risk                                                          ####################################
#############################################################################################################################
require(r.jive); require(ggplot2); require(xtable); require(reshape)
require(tidyr); require(arsenal); require(dplyr); require(stringr); require(gridExtra)

results.dir_rJ1 = "H:/My Documents/P-JIVE/Results/Simulation_Results/JointRank1_GaussGauss_n1000"
results.dir_rJ3 = "H:/My Documents/P-JIVE/Results/Simulation_Results/JointRank3_GaussGauss_n1000"
imgs.fldr_rJ1 = results.dir_rJ1
imgs.fldr_rJ3 = results.dir_rJ3
# imgs.fldr = "H:/My Documents/Apps/Overleaf/CJIVE Manuscript/Images"

prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
source(file.path(prog.dir, "Functions_for_CJIVE.R"))

ajive.dir = "H:/My Documents/Applications2/r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb.cols = cbPalette[c(2:5,7)]

########################################################################################################################
########################################################################################################################
#########       Rank 1                     #############################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 20, p2 = 20

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_rJ1, "SimBin_P120_P220.05.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_rJ1,"SimBin_P120_P220.005.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_rJ1,"SimBin_P120_P220.05.005"), 20, 20)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_rJ1,"SimBin_P120_P220.005.005"), 20, 20)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

JntVar200.Table = aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[3:4] = c("Mean X_1","S.D. X_1")
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)[3])
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[5:6] = c("Mean X_2","S.D. X_2")
JntVar200.Table = round(JntVar200.Table,3)

Time.Table.Mean = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                            ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) mean(x/60))
Time.Table.SD = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                          ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) sd(x/60))

sim.gg.data = ConvSims_gg_ProJIVE(AllSims)
# temp = ConvSims_gg(AllSims)
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cb.cols, text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cb.cols, text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cb.cols, text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_rJ1_P120_P220_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ1_P120_P220_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_rJ1_P120_P220_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ1_P120_P220_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
print(load.plot)
dev.off()
########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 20, p2 = 200

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_rJ1, "SimBin_P120_P2200.05.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_rJ1,"SimBin_P120_P2200.005.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_rJ1,"SimBin_P120_P2200.05.005"), 20, 200)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_rJ1,"SimBin_P120_P2200.005.005"), 20, 200)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

JntVar200.Table = aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[3:4] = c("Mean X_1","S.D. X_1")
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)[3])
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[5:6] = c("Mean X_2","S.D. X_2")
JntVar200.Table = round(JntVar200.Table,3)

Time.Table.Mean = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                            ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) mean(x/60))
Time.Table.SD = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                          ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) sd(x/60))

sim.gg.data = ConvSims_gg_ProJIVE(AllSims)
# temp = ConvSims_gg(AllSims)
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cb.cols, text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cb.cols, text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cb.cols, text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_rJ1_P120_P2200_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ1_P120_P2200_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()


plot.name = paste("GG_SimBin_rJ1_P120_P2200_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
print(scor.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ1_P120_P2020_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ1, plot.name))
print(norm.plot)
dev.off()
########################################################################################################################
########################################################################################################################
#########       Rank 3                     #############################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
######Rank 3: p1 = 20, p2 = 20

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_rJ3, "SimBin_P120_P220.05.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_rJ3,"SimBin_P120_P220.005.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_rJ3,"SimBin_P120_P220.05.005"), 20, 20)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_rJ3,"SimBin_P120_P220.005.005"), 20, 20)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

JntVar200.Table = aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[3:4] = c("Mean X_1","S.D. X_1")
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)[3])
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[5:6] = c("Mean X_2","S.D. X_2")
JntVar200.Table = round(JntVar200.Table,3)

Time.Table.Mean = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                            ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) mean(x/60))
Time.Table.SD = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                          ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) sd(x/60))

sim.gg.data = ConvSims_gg_ProJIVE(AllSims)
# temp = ConvSims_gg(AllSims)
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cb.cols, text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cb.cols, text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cb.cols, text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_rJ3_P120_P220_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ3_P120_P220_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_rJ3_P120_P220_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ3_P120_P220_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
print(load.plot)
dev.off()
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######Rank 3: p1 = 20, p2 = 200

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_rJ3, "SimBin_P120_P2200.05.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_rJ3,"SimBin_P120_P2200.005.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_rJ3,"SimBin_P120_P2200.05.005"), 20, 200)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_rJ3,"SimBin_P120_P2200.005.005"), 20, 200)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

JntVar200.Table = aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.X ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[3:4] = c("Mean X_1","S.D. X_1")
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = mean)[3])
JntVar200.Table = cbind(JntVar200.Table, aggregate(formula = Joint.Var.Exp.Y ~ JntVarEx1 + JntVarEx2, data = AllSims.rows, FUN = sd)[3])
colnames(JntVar200.Table)[5:6] = c("Mean X_2","S.D. X_2")
JntVar200.Table = round(JntVar200.Table,3)

Time.Table.Mean = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                            ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) mean(x/60))
Time.Table.SD = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time) 
                          ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) sd(x/60))

sim.gg.data = ConvSims_gg_ProJIVE(AllSims)
# temp = ConvSims_gg(AllSims)
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cb.cols, text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cb.cols, text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cb.cols, text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_rJ3_P120_P2200_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ3_P120_P2200_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_rJ3_P120_P2200_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_rJ3_P120_P2200_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_rJ3, plot.name))
print(load.plot)
dev.off()
################################################################################################################################
######################################## Investigate associations between proportion of correct rank and #######################
############################################## Indiv Var Ex values #############################################################
# ################################################################################################################################
# JVE.Ranks = rbind(sim.results.0505[,c(grep("Var", sim.names), grep("Rank", sim.names))],
#                   sim.results.05005[,c(grep("Var", sim.names), grep("Rank", sim.names))],
#                   sim.results.00505[,c(grep("Var", sim.names), grep("Rank", sim.names))],
#                   sim.results.005005[,c(grep("Var", sim.names), grep("Rank", sim.names))])
# colnames(JVE.Ranks)
# 
# ###No clear association from scatter plots
# plot(JVE.Ranks[JVE.Ranks$JntVarEx2==0.5, c("Indiv.Var.Exp.Y", "Joint.Var.Exp.Y")])
# plot(JVE.Ranks$Indiv.Var.Exp.Y, JVE.Ranks$CC.95..Joint.Rank)
# ###THE R_{I,Y}^2 value is 0.3 
# 
# 
# plot(JVE.Ranks$Indiv.Var.Exp.X, JVE.Ranks$Joint.Var.Exp.X)
# plot(round(JVE.Ranks$Indiv.Var.Exp.Y, 1), JVE.Ranks$Joint.Var.Exp.Y)
# 
# plot(JVE.Ranks$Indiv.Var.Exp.X, JVE.Ranks$iJIVE.Joint.Rank)
# plot(JVE.Ranks$Indiv.Var.Exp.X, JVE.Ranks$CC.Elbow.Joint.Rank)
# plot(JVE.Ranks$Indiv.Var.Exp.X, JVE.Ranks$aJIVE.Correct.Joint.Rank)
# plot(JVE.Ranks$Indiv.Var.Exp.X, JVE.Ranks$aJIVE.Wrong.Joint.Rank)
# 
# plot(JVE.Ranks$Indiv.Var.Exp.X, JVE.Ranks$CC.95..Joint.Rank)
# 
# 
# 

