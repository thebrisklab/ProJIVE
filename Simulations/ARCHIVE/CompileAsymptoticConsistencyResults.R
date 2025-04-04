#############################################################################################################################
### Compile and evaluate ProJIVE analyses of simulated data to check statistical consistency  ###############################
### Author: Raphiel J. Murden                                                                 ###############################
### Supervised by Benjamin Risk                                                               ###############################
#############################################################################################################################
prog.dir = "H:/My Documents/P-JIVE/Programs/Functions"
require(CJIVE); require(singR); require(ggplot2); require(lubridate)
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
res.rJ1.rI1.dir = "H:/My Documents/P-JIVE/Results/Simulation_Results/ConsistencyCheck/ConsistencyCheckJointRank_rJ1_rI1_Results"
res.rJ1.rI2.dir = "H:/My Documents/P-JIVE/Results/Simulation_Results/ConsistencyCheck/ConsistencyCheckJointRank_rJ1_rI2_Results"
res.rJ3.rI1.dir = "H:/My Documents/P-JIVE/Results/Simulation_Results/ConsistencyCheck/ConsistencyCheckJointRank_rJ3_rI1_Results"
res.rJ3.rI2.dir = "H:/My Documents/P-JIVE/Results/Simulation_Results/ConsistencyCheck/ConsistencyCheckJointRank_rJ3_rI2_Results"

#############################################################################################
########### Results: Joint Rank = 1; Indiv Ranks = 1; p2 = 20 ###############################
########### n = 200
count = 1
res.rJ1.rI1.p20.n200 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p220"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p20.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n200 = sapply(res.rJ1.rI1.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n200 = sapply(res.rJ1.rI1.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])


########### n = 1000
count = 1
res.rJ1.rI1.p20.n1000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p220"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p20.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n1000 = sapply(res.rJ1.rI1.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n1000 = sapply(res.rJ1.rI1.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ1.rI1.p20.n2000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p220"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p20.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n2000 = sapply(res.rJ1.rI1.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n2000 = sapply(res.rJ1.rI1.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ1.rI1.p20.n5000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p220"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p20.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n5000 = sapply(res.rJ1.rI1.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n5000 = sapply(res.rJ1.rI1.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ1.rI1.p20.n10000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p220"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p20.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n10000 = sapply(res.rJ1.rI1.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n10000 = sapply(res.rJ1.rI1.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ1.rI1.p20 = data.frame(PMSE = c(pmse.X1.p20.n200, pmse.X2.p20.n200, pmse.X1.p20.n1000, pmse.X2.p20.n1000,
                                                 pmse.X1.p20.n2000, pmse.X2.p20.n2000, pmse.X1.p20.n5000, pmse.X2.p20.n5000,
                                                 pmse.X1.p20.n10000, pmse.X2.p20.n10000
                                           ),
                                  SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p20.n200), length(pmse.X1.p20.n1000), 
                                                                             length(pmse.X1.p20.n2000), length(pmse.X1.p20.n5000), 
                                                                             length(pmse.X1.p20.n10000))),
                                                        labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                  DataSet = factor(rep(c("X1", "X2"), each = 10)))

pmse.rJ1.rI1.p20.plot = ggplot(data = dat.pmse.rJ1.rI1.p20) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
ggsave(plot = pmse.rJ1.rI1.p20.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ1.rI1.dir, paste0("Asmy.PMSE_rJ1_rI1_p20_", today(), ".pdf"))))

##############################################################################################
########### Results: Joint Rank = 1; Indiv Ranks = 1; p2 = 200 ###############################
########### n = 200
count = 1
res.rJ1.rI1.p200.n200 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p2200"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p200.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n200 = sapply(res.rJ1.rI1.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n200 = sapply(res.rJ1.rI1.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ1.rI1.p200.n1000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p2200"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p200.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n1000 = sapply(res.rJ1.rI1.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n1000 = sapply(res.rJ1.rI1.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ1.rI1.p200.n2000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p2200"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p200.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n2000 = sapply(res.rJ1.rI1.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n2000 = sapply(res.rJ1.rI1.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ1.rI1.p200.n5000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p2200"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p200.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n5000 = sapply(res.rJ1.rI1.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n5000 = sapply(res.rJ1.rI1.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ1.rI1.p200.n10000 = list()
for(nm in list.files(file.path(res.rJ1.rI1.dir, "SimBin_p2200"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ1.rI1.p200.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n10000 = sapply(res.rJ1.rI1.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n10000 = sapply(res.rJ1.rI1.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ1.rI1.p200 = data.frame(PMSE = c(pmse.X1.p200.n200, pmse.X2.p200.n200, pmse.X1.p200.n1000, pmse.X2.p200.n1000,
                                                 pmse.X1.p200.n2000, pmse.X2.p200.n2000, pmse.X1.p200.n5000, pmse.X2.p200.n5000,
                                                 pmse.X1.p200.n10000, pmse.X2.p200.n10000),
                                  SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p200.n200), length(pmse.X1.p200.n1000), 
                                                                           length(pmse.X1.p200.n2000), length(pmse.X1.p200.n5000), 
                                                                           length(pmse.X1.p200.n10000))), 
                                                          labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                  DataSet = factor(rep(c("X1", "X2"), each = 10)))
pmse.rJ1.rI1.p200.plot = ggplot(data = dat.pmse.rJ1.rI1.p200) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
ggsave(plot = pmse.rJ1.rI1.p200.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ1.rI1.dir, paste0("Asmy.PMSE_rJ1_rI1_p200_", today(), ".pdf"))))

# load(file.path(res.rJ1.rI1.dir, ))
#############################################################################################
########### Results: Joint Rank = 1; Indiv Ranks = 2; p2 = 20 ###############################
########### n = 200
count = 1
res.rJ1.rI2.p20.n200 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p220"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p20.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n200 = sapply(res.rJ1.rI2.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n200 = sapply(res.rJ1.rI2.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ1.rI2.p20.n1000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p220"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p20.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n1000 = sapply(res.rJ1.rI2.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n1000 = sapply(res.rJ1.rI2.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ1.rI2.p20.n2000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p220"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p20.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n2000 = sapply(res.rJ1.rI2.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n2000 = sapply(res.rJ1.rI2.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ1.rI2.p20.n5000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p220"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p20.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n5000 = sapply(res.rJ1.rI2.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n5000 = sapply(res.rJ1.rI2.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ1.rI2.p20.n10000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p220"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p20.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n10000 = sapply(res.rJ1.rI2.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n10000 = sapply(res.rJ1.rI2.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ1.rI2.p20 = data.frame(PMSE = c(pmse.X1.p20.n200, pmse.X2.p20.n200, pmse.X1.p20.n1000, pmse.X2.p20.n1000,
                                                 pmse.X1.p20.n2000, pmse.X2.p20.n2000, pmse.X1.p20.n5000, pmse.X2.p20.n5000,
                                                 pmse.X1.p20.n10000, pmse.X2.p20.n10000),
                                  SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p20.n200), length(pmse.X1.p20.n1000), 
                                                                           length(pmse.X1.p20.n2000), length(pmse.X1.p20.n5000), 
                                                                           length(pmse.X1.p20.n10000))),
                                                      labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                  DataSet = factor(rep(c("X1", "X2"), each = 10)))
pmse.rJ1.rI2.p20.plot = ggplot(data = dat.pmse.rJ1.rI2.p20) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
ggsave(plot = pmse.rJ1.rI2.p20.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ1.rI2.dir, paste0("Asmy.PMSE_rJ1_rI2_p20_", today(), ".pdf"))))

#############################################################################################
########### Results: Joint Rank = 1; Indiv Ranks = 2; p2 = 200 ###############################
########### n = 200
count = 1
res.rJ1.rI2.p200.n200 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p2200"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p200.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n200 = sapply(res.rJ1.rI2.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n200 = sapply(res.rJ1.rI2.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ1.rI2.p200.n1000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p2200"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p200.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n1000 = sapply(res.rJ1.rI2.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n1000 = sapply(res.rJ1.rI2.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ1.rI2.p200.n2000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p2200"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p200.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n2000 = sapply(res.rJ1.rI2.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n2000 = sapply(res.rJ1.rI2.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ1.rI2.p200.n5000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p2200"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p200.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n5000 = sapply(res.rJ1.rI2.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n5000 = sapply(res.rJ1.rI2.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ1.rI2.p200.n10000 = list()
for(nm in list.files(file.path(res.rJ1.rI2.dir, "SimBin_p2200"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ1.rI2.p200.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n10000 = sapply(res.rJ1.rI2.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n10000 = sapply(res.rJ1.rI2.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ1.rI2.p200 = data.frame(PMSE = c(pmse.X1.p200.n200, pmse.X2.p200.n200, pmse.X1.p200.n1000, pmse.X2.p200.n1000,
                                                  pmse.X1.p200.n2000, pmse.X2.p200.n2000, pmse.X1.p200.n5000, pmse.X2.p200.n5000,
                                                  pmse.X1.p200.n10000, pmse.X2.p200.n10000),
                                   SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p200.n200), length(pmse.X1.p200.n1000), 
                                                                            length(pmse.X1.p200.n2000), length(pmse.X1.p200.n5000), 
                                                                            length(pmse.X1.p200.n10000))),
                                                       labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                   DataSet = factor(rep(c("X1", "X2"), each = length(pmse.X1.p200.n200))))
pmse.rJ1.rI2.p200.plot = ggplot(data = dat.pmse.rJ1.rI2.p200) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
ggsave(plot = pmse.rJ1.rI2.p200.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ1.rI2.dir, paste0("Asmy.PMSE_rJ1_rI2_p200_", today(), ".pdf"))))

#############################################################################################
########### Results: Joint Rank = 3; Indiv Ranks = 1; p2 = 20 ###############################
########### n = 200
count = 1
res.rJ3.rI1.p20.n200 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p220"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p20.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n200 = sapply(res.rJ3.rI1.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n200 = sapply(res.rJ3.rI1.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ3.rI1.p20.n1000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p220"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p20.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n1000 = sapply(res.rJ3.rI1.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n1000 = sapply(res.rJ3.rI1.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ3.rI1.p20.n2000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p220"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p20.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n2000 = sapply(res.rJ3.rI1.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n2000 = sapply(res.rJ3.rI1.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ3.rI1.p20.n5000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p220"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p20.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n5000 = sapply(res.rJ3.rI1.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n5000 = sapply(res.rJ3.rI1.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ3.rI1.p20.n10000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p220"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p20.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n10000 = sapply(res.rJ3.rI1.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n10000 = sapply(res.rJ3.rI1.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ3.rI1.p20 = data.frame(PMSE = c(pmse.X1.p20.n200, pmse.X2.p20.n200, pmse.X1.p20.n1000, pmse.X2.p20.n1000,
                                            pmse.X1.p20.n2000, pmse.X2.p20.n2000, pmse.X1.p20.n5000, pmse.X2.p20.n5000,
                                            pmse.X1.p20.n10000, pmse.X2.p20.n10000),
                                   SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p20.n200), length(pmse.X1.p20.n1000), 
                                                                            length(pmse.X1.p20.n2000), length(pmse.X1.p20.n5000), 
                                                                            length(pmse.X1.p20.n10000))),
                                                       labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                   DataSet = factor(rep(c("X1", "X2"), each = length(pmse.X1.p20.n200))))
pmse.rJ3.rI1.p20.plot = ggplot(data = dat.pmse.rJ3.rI1.p20) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
ggsave(plot = pmse.rJ3.rI1.p20.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ3.rI1.dir, paste0("Asmy.PMSE_rJ1_rI2_p20_", today(), ".pdf"))))

#############################################################################################
########### Results: Joint Rank = 3; Indiv Ranks = 1; p2 = 200 ###############################
########### n = 200
count = 1
res.rJ3.rI1.p200.n200 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p2200"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p200.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n200 = sapply(res.rJ3.rI1.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n200 = sapply(res.rJ3.rI1.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ3.rI1.p200.n1000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p2200"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p200.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n1000 = sapply(res.rJ3.rI1.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n1000 = sapply(res.rJ3.rI1.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ3.rI1.p200.n2000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p2200"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p200.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n2000 = sapply(res.rJ3.rI1.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n2000 = sapply(res.rJ3.rI1.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ3.rI1.p200.n5000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p2200"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p200.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n5000 = sapply(res.rJ3.rI1.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n5000 = sapply(res.rJ3.rI1.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ3.rI1.p200.n10000 = list()
for(nm in list.files(file.path(res.rJ3.rI1.dir, "SimBin_p2200"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ3.rI1.p200.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n10000 = sapply(res.rJ3.rI1.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n10000 = sapply(res.rJ3.rI1.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ3.rI1.p200 = data.frame(PMSE = c(pmse.X1.p200.n200, pmse.X2.p200.n200, pmse.X1.p200.n1000, pmse.X2.p200.n1000,
                                           pmse.X1.p200.n2000, pmse.X2.p200.n2000, pmse.X1.p200.n5000, pmse.X2.p200.n5000,
                                           pmse.X1.p200.n10000, pmse.X2.p200.n10000),
                                  SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p200.n200), length(pmse.X1.p200.n1000), 
                                                                           length(pmse.X1.p200.n2000), length(pmse.X1.p200.n5000), 
                                                                           length(pmse.X1.p200.n10000))),
                                                      labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                  DataSet = factor(rep(c("X1", "X2"), each = length(pmse.X1.p200.n200))))
pmse.rJ3.rI1.p200.plot = ggplot(data = dat.pmse.rJ3.rI1.p200) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
ggsave(plot = pmse.rJ3.rI1.p200.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ3.rI1.dir, paste0("Asmy.PMSE_rJ1_rI2_p200_", today(), ".pdf"))))

#############################################################################################
########### Results: Joint Rank = 3; Indiv Ranks = 2; p2 = 20 ###############################
########### n = 200
count = 1
res.rJ3.rI2.p20.n200 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p220"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p20.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n200 = sapply(res.rJ3.rI2.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n200 = sapply(res.rJ3.rI2.p20.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ3.rI2.p20.n1000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p220"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p20.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n1000 = sapply(res.rJ3.rI2.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n1000 = sapply(res.rJ3.rI2.p20.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ3.rI2.p20.n2000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p220"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p20.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n2000 = sapply(res.rJ3.rI2.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n2000 = sapply(res.rJ3.rI2.p20.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ3.rI2.p20.n5000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p220"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p20.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n5000 = sapply(res.rJ3.rI2.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n5000 = sapply(res.rJ3.rI2.p20.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ3.rI2.p20.n10000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p220"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p20.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p20.n10000 = sapply(res.rJ3.rI2.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p20.n10000 = sapply(res.rJ3.rI2.p20.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ3.rI2.p20 = data.frame(PMSE = c(pmse.X1.p20.n200, pmse.X2.p20.n200, pmse.X1.p20.n1000, pmse.X2.p20.n1000,
                                           pmse.X1.p20.n2000, pmse.X2.p20.n2000, pmse.X1.p20.n5000, pmse.X2.p20.n5000,
                                           pmse.X1.p20.n10000, pmse.X2.p20.n10000),
                                  SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p20.n200), length(pmse.X1.p20.n1000), 
                                                                           length(pmse.X1.p20.n2000), length(pmse.X1.p20.n5000), 
                                                                           length(pmse.X1.p20.n10000))),
                                                      labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                  DataSet = factor(rep(c("X1", "X2"), each = length(pmse.X1.p20.n200))))
pmse.rJ3.rI2.p20.plot = ggplot(data = dat.pmse.rJ3.rI2.p20) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
pmse.rJ3.rI2.p20.plot + labs(title = "p2 = 20")
ggsave(plot = pmse.rJ3.rI2.p20.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ3.rI2.dir, paste0("Asmy.PMSE_rJ1_rI2_p20_", today(), ".pdf"))))

#############################################################################################
########### Results: Joint Rank = 3; Indiv Ranks = 2; p2 = 200 ###############################
########### n = 200
count = 1
res.rJ3.rI2.p200.n200 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p2200"), pattern = "n200_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p200.n200[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n200 = sapply(res.rJ3.rI2.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n200 = sapply(res.rJ3.rI2.p200.n200, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 1000
count = 1
res.rJ3.rI2.p200.n1000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p2200"), pattern = "n1000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p200.n1000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n1000 = sapply(res.rJ3.rI2.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n1000 = sapply(res.rJ3.rI2.p200.n1000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 2000
count = 1
res.rJ3.rI2.p200.n2000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p2200"), pattern = "n2000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p200.n2000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n2000 = sapply(res.rJ3.rI2.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n2000 = sapply(res.rJ3.rI2.p200.n2000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 5000
count = 1
res.rJ3.rI2.p200.n5000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p2200"), pattern = "n5000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p200.n5000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n5000 = sapply(res.rJ3.rI2.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n5000 = sapply(res.rJ3.rI2.p200.n5000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

########### n = 10000
count = 1
res.rJ3.rI2.p200.n10000 = list()
for(nm in list.files(file.path(res.rJ3.rI2.dir, "SimBin_p2200"), pattern = "n10000_", full.names = TRUE)){
  load(nm); res.rJ3.rI2.p200.n10000[[count]] = out; rm(out)
  count = count + 1
}
rm(count)

pmse.X1.p200.n10000 = sapply(res.rJ3.rI2.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X1"]])
pmse.X2.p200.n10000 = sapply(res.rJ3.rI2.p200.n10000, FUN = function(x) x[["PMSE.PJIVE.Loads.X2"]])

dat.pmse.rJ3.rI2.p200 = data.frame(PMSE = c(pmse.X1.p200.n200, pmse.X2.p200.n200, pmse.X1.p200.n1000, pmse.X2.p200.n1000,
                                            pmse.X1.p200.n2000, pmse.X2.p200.n2000, pmse.X1.p200.n5000, pmse.X2.p200.n5000,
                                            pmse.X1.p200.n10000, pmse.X2.p200.n10000),
                                   SampleSize = factor(rep(1:5, times = 2*c(length(pmse.X1.p200.n200), length(pmse.X1.p200.n1000), 
                                                                            length(pmse.X1.p200.n2000), length(pmse.X1.p200.n5000), 
                                                                            length(pmse.X1.p200.n10000))),
                                                       labels = paste("n =", c(200, 1000, 2000, 5000, 10000))),
                                   DataSet = factor(rep(c("X1", "X2"), each = length(pmse.X1.p200.n200))))
pmse.rJ3.rI2.p200.plot = ggplot(data = dat.pmse.rJ3.rI2.p200) + geom_boxplot(aes(SampleSize, PMSE, fill = DataSet))  
pmse.rJ3.rI2.p200.plot + labs(title = "p2 = 200")
ggsave(plot = pmse.rJ3.rI2.p200.plot,  width = 4, height = 4,
       filename = file.path(file.path(res.rJ3.rI2.dir, paste0("Asmy.PMSE_rJ1_rI2_p200_", today(), ".pdf"))))








