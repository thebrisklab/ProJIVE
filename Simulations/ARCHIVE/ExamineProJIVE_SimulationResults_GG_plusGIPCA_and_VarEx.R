#############################################################################################################################
### Compare methods of implementing JIVE analysis on simulated datasets                  ####################################
### Author: Raphiel J. Murden                                                            ####################################
### Supervised by Benjamin Risk                                                          ####################################
#############################################################################################################################
require(r.jive); require(ggplot2); require(xtable); require(reshape)
require(tidyr); require(arsenal); require(dplyr); require(stringr); require(gridExtra)

results.dir_n500_rJ1 = "H:/My Documents/P-JIVE/Results/Simulation_Results/VarEx_Sims/JointRank1_GaussGauss_n500"
results.dir_n500_rJ3 = "H:/My Documents/P-JIVE/Results/Simulation_Results/VarEx_Sims/JointRank3_GaussGauss_n500"
results.dir_n1000_rJ1 = "H:/My Documents/P-JIVE/Results/Simulation_Results/VarEx_Sims/JointRank1_GaussGauss_n1000"
results.dir_n1000_rJ3 = "H:/My Documents/P-JIVE/Results/Simulation_Results/VarEx_Sims/JointRank3_GaussGauss_n1000"
imgs.fldr_n500_rJ1 = results.dir_n500_rJ1
imgs.fldr_n500_rJ3 = results.dir_n500_rJ3
imgs.fldr_n1000_rJ1 = results.dir_n1000_rJ1
imgs.fldr_n1000_rJ3 = results.dir_n1000_rJ3

prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
source(file.path(prog.dir, "Functions_for_CJIVE.R"))

ajive.dir = "H:/My Documents/Applications2/r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb.cols.vx = cbPalette[c(8,2,3,5)]
cb.cols = cbPalette[c(2:5,7)]

########################################################################################################################
########################################################################################################################
#########       Rank 1 - n=500             #############################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 20, p2 = 20

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n500_rJ1, "SimBin_P120_P220.05.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n500_rJ1,"SimBin_P120_P220.005.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n500_rJ1,"SimBin_P120_P220.05.005"), 20, 20)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n500_rJ1,"SimBin_P120_P220.005.005"), 20, 20)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 500; AllSims.rows[,"r.J"] = 1;
Time.Table_p220_n500_rJ1 = aggregate(cbind(AJIVE_Time, R.JIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                     ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                     FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))

sim.gg = ConvSims_gg_ProJIVE2(AllSims, 500)
sim.gg.data = sim.gg$Norms

# temp = ConvSims_gg(AllSims)
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n500_rJ1_P120_P220_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ1_P120_P220_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ1_P120_P220_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ1_P120_P220_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
print(load.plot)
dev.off()


VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 20, 500)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n500_rJ1, paste0("GG_VarEx_n500_rJ1_p2_20_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")

########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 20, p2 = 200

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n500_rJ1, "SimBin_P120_P2200.05.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n500_rJ1,"SimBin_P120_P2200.005.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n500_rJ1,"SimBin_P120_P2200.05.005"), 20, 200)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n500_rJ1,"SimBin_P120_P2200.005.005"), 20, 200)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 500; AllSims.rows[,"r.J"] = 1;
Time.Table_p2200_n500_rJ1 = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                      ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                      FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))
Time.Table_n500_rJ1 = rbind(Time.Table_p220_n500_rJ1, Time.Table_p2200_n500_rJ1)

sim.gg = ConvSims_gg_ProJIVE2(AllSims, 500)
sim.gg.data = sim.gg$Norms

# temp = ConvSims_gg(AllSims)
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n500_rJ1_P120_P2200_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ1_P120_P2200_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ1_P120_P2200_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ1_P120_P2200_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ1, plot.name))
print(load.plot)
dev.off()

# VarEx.dat.gg = sim.gg$VarEx
# varex.plot = gg.varex.plot(VarEx.dat.gg, cb.cols, show.legend = T)
# varex.plot


VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 200, 500)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n500_rJ1, paste0("GG_VarEx_n500_rJ1_p2_200_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")

########################################################################################################################
########################################################################################################################
#########       Rank 1 - n=1000             #############################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 20, p2 = 20

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1, "SimBin_P120_P220.05.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1,"SimBin_P120_P220.005.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1,"SimBin_P120_P220.05.005"), 20, 20)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1,"SimBin_P120_P220.005.005"), 20, 20)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 1000; AllSims.rows[,"r.J"] = 1;
Time.Table_p220_n1000_rJ1 = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                      ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                      FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))

sim.gg.data = ConvSims_gg_ProJIVE2(AllSims, 1000)$Norms
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P220_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P220_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P220_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P220_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
print(load.plot)
dev.off()

Time.Table_p220_rJ1.SD.Mean = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                        ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) mean(x/60))
Time.Table_p220_rJ1.SD.SD = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                      ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) sd(x/60))

VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 20, 1000)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n1000_rJ1, paste0("GG_VarEx_n1000_rJ1_p2_20_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")

########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 20, p2 = 200

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1, "SimBin_P120_P2200.05.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1,"SimBin_P120_P2200.005.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1,"SimBin_P120_P2200.05.005"), 20, 200)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ1,"SimBin_P120_P2200.005.005"), 20, 200)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 1000; AllSims.rows[,"r.J"] = 1;
Time.Table_p2200_n1000_rJ1 = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                       ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                       FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))
Time.Table_n1000_rJ1 = rbind(Time.Table_p220_n1000_rJ1, Time.Table_p2200_n1000_rJ1)
Time.Table_rJ1 = rbind(Time.Table_n500_rJ1, Time.Table_n1000_rJ1)

sim.gg.data = ConvSims_gg_ProJIVE2(AllSims, 1000)$Norms
sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P2200_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P2200_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P2200_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ1_P120_P2200_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ1, plot.name))
print(load.plot)
dev.off()

VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 200, 1000)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n1000_rJ1, paste0("GG_VarEx_n1000_rJ1_p2_200_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")
########################################################################################################################
########################################################################################################################
#########       Rank 3 - n=500             #############################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
######Rank 3: p1 = 20, p2 = 20

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n500_rJ3, "SimBin_P120_P220.05.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n500_rJ3,"SimBin_P120_P220.005.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n500_rJ3,"SimBin_P120_P220.05.005"), 20, 20)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n500_rJ3,"SimBin_P120_P220.005.005"), 20, 20)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 500; AllSims.rows[,"r.J"] = 3;
Time.Table_p220_n500_rJ3 = aggregate(cbind(AJIVE_Time, R.JIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                     ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                     FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))

sim.gg = ConvSims_gg_ProJIVE2(AllSims, 500)
sim.gg.data = sim.gg$Norms

sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n500_rJ3_P120_P220_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ3_P120_P220_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ3_P120_P220_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ3_P120_P220_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
print(load.plot)
dev.off()


VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 20, 500)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n500_rJ3, paste0("GG_VarEx_n500_rJ3_p2_20_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######Rank 3: p1 = 20, p2 = 200

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n500_rJ3, "SimBin_P120_P2200.05.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n500_rJ3,"SimBin_P120_P2200.005.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n500_rJ3,"SimBin_P120_P2200.05.005"), 20, 200)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n500_rJ3,"SimBin_P120_P2200.005.005"), 20, 200)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 500; AllSims.rows[,"r.J"] = 3;
Time.Table_p2200_n500_rJ3 = aggregate(cbind(AJIVE_Time, R.JIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                      ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                      FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))
Time.Table_n500_rJ3 = rbind(Time.Table_p220_n500_rJ3, Time.Table_p2200_n500_rJ3)

sim.gg = ConvSims_gg_ProJIVE2(AllSims, 500)
sim.gg.data = sim.gg$Norms

sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n500_rJ3_P120_P2200_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ3_P120_P2200_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ3_P120_P2200_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n500_rJ3_P120_P2200_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n500_rJ3, plot.name))
print(load.plot)
dev.off()


Time.Table_p220_rJ3.SD.Mean = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                        ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) mean(x/60))
Time.Table_p220_rJ3.SD.SD = aggregate(cbind(R.JIVE_Time, AJIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                      ~ JntVarEx1 + JntVarEx2 + p2, data = AllSims.rows, FUN = function(x) sd(x/60))


VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 200, 500)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n500_rJ3, paste0("GG_VarEx_n500_rJ3_p2_200_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")

########################################################################################################################
########################################################################################################################
#########       Rank 3 - n=1000             #############################################################################
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
######Rank 3: p1 = 20, p2 = 20

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3, "SimBin_P120_P220.05.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3,"SimBin_P120_P220.005.05"), 20, 20)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3,"SimBin_P120_P220.05.005"), 20, 20)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3,"SimBin_P120_P220.005.005"), 20, 20)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 1000; AllSims.rows[,"r.J"] = 3;
Time.Table_p220_n1000_rJ3 = aggregate(cbind(AJIVE_Time, R.JIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                      ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                      FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))
sim.gg = ConvSims_gg_ProJIVE2(AllSims, 1000)
sim.gg.data = sim.gg$Norms

sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n1000_rJ3_P120_P220_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ3_P120_P220_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ3_P120_P220_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ3_P120_P220_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
print(load.plot)
dev.off()


VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 20, 1000)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25)  
ggsave(file.path(imgs.fldr_n1000_rJ3, paste0("GG_VarEx_n1000_rJ3_p2_20_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")

########################################################################################################################
########################################################################################################################
######Rank 3: p1 = 20, p2 = 200

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results.0505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3, "SimBin_P120_P2200.05.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in data set 1 and 0.5 in dataset 2 #######
##############################################################################################################################
sim.results.00505 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3,"SimBin_P120_P2200.005.05"), 20, 200)

##############################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets                  #######
##############################################################################################################################
sim.results.05005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3,"SimBin_P120_P2200.05.005"), 20, 200)

####################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.05 in both data sets (Total R^2 = 0.75)
####################################################################################################################
sim.results.005005 = GetSimResults_Dir(file.path(results.dir_n1000_rJ3,"SimBin_P120_P2200.005.005"), 20, 200)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
AllSims = cbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)
AllSims[,"Indiv.Var.Exp.X"] = 0.25
AllSims[,"Indiv.Var.Exp.Y"] = 0.25
AllSims.rows = rbind(sim.results.005005, sim.results.05005, sim.results.00505, sim.results.0505)

AllSims.rows[,"n"] = 1000; AllSims.rows[,"r.J"] = 3;
Time.Table_p2200_n1000_rJ3 = aggregate(cbind(AJIVE_Time, R.JIVE_Time, ProJIVE_Time, GIPCA_Time) 
                                       ~ JntVarEx1 + JntVarEx2 + p2 + n + r.J, data = AllSims.rows, 
                                       FUN = function(x) paste0(round(mean(x/60),1), " (", round(sd(x/60),3), ")"))
Time.Table_n1000_rJ3 = rbind(Time.Table_p220_n1000_rJ3, Time.Table_p2200_n1000_rJ3)
Time.Table_rJ3 = rbind(Time.Table_n500_rJ3, Time.Table_n1000_rJ3)

sim.gg = ConvSims_gg_ProJIVE2(AllSims, 1000)
sim.gg.data = sim.gg$Norms

sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loadings", sim.gg.data[,"Type"]),]

norm.plot = gg.norm.plot(sim.gg.data, cbPalette[2:5], text.size = 15)
score.plot = gg.score.norm.plot(sim.scores, cbPalette[2:5], text.size = 15, show.legend = TRUE)
load.plot = gg.load.norm.plot(sim.loads, cbPalette[2:5], text.size = 15, show.legend = TRUE)

plot.name = paste("GG_SimBin_n1000_rJ3_P120_p2200_ScoreAndLoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
print(norm.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ3_P120_p2200_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
gg.norm.plot(sim.gg.data, cb.cols, text.size = 15, show.legend = TRUE)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ3_P120_p2200_ScorePlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
print(score.plot)
dev.off()

plot.name = paste("GG_SimBin_n1000_rJ3_P120_p2200_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr_n1000_rJ3, plot.name))
print(load.plot)
dev.off()

VarEx.dat.gg = MakeVarEx.data.gg(AllSims.rows, 20, 200, 1000)
x.labs = c(expression("Indiv"*"X"[1]), expression("Indiv"*"X"[2]), 
           expression("Joint"*"X"[1]), expression("Joint"*"X"[2]))

ggplot(data = VarEx.dat.gg, aes(y = Mean_EmpJntVarEx, x = Type, fill = Method)) + 
  geom_bar(position = "dodge", stat = "identity") + ylab("Mean Variance Explained") +  
  scale_fill_manual(name = "JIVE Method", values = cb.cols.vx) +
  scale_x_discrete(labels = x.labs) + 
  facet_grid(JntVarEx1.fac~JntVarEx2.fac, scales = "free", labeller = label_parsed) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 01, angle = 70, size = 9),
        text = element_text(size = 11)) + 
  geom_errorbar(aes(ymin=Mean_EmpJntVarEx-SD_EmpJntVarEx, ymax=Mean_EmpJntVarEx+SD_EmpJntVarEx), position=position_dodge(), size = 0.25) 
ggsave(file.path(imgs.fldr_n1000_rJ3, paste0("GG_VarEx_n1000_rJ3_p2_200_", lubridate::today(), ".pdf")),
       width = 5, height = 4, units = "in")


GG.Time.Table = rbind(Time.Table_rJ1, Time.Table_rJ3)[,c(5:1,6:9)]
write.csv(GG.Time.Table, file = "H:/My Documents/P-JIVE/Results/Simulation_Results/VarEx_Sims/GG_Time_Table.csv", row.names = FALSE)

# apply(GG.Time.Table[,6:9], 2, function(x) str_extract(x, "(?<=.)[:space:]"))

GG.Mean.Time.Table = cbind(GG.Time.Table[,1:5], str_extract(GG.Time.Table[,6:9], "^[:space:]"))
lm(cbind(AJIVE_Time, R.JIVE_Time, ProJIVE_Time, GIPCA_Time)~n , data = as.data.frame(GG.Time.Table))
