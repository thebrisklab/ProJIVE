#############################################################################################################################
### compare methods of implementing JIVE analysis on simulated datasets                  ####################################
### Author: Raphiel J. Murden                                                            ####################################
### Supervised by Benjamin Risk                                                          ####################################
#############################################################################################################################
require(r.jive); require(ggplot2); require(xtable); require(reshape)
require(tidyr); require(arsenal); require(dplyr); require(stringr); require(gridExtra)

results.dir = "H:/My Documents/ProJIVE/Results/Simulation_Results/SimFeng_P1100_P21000"
imgs.fldr = results.dir

prog.dir = "H:/My Documents/ProJIVE/Programs/Functions"
source(file.path(prog.dir, "Functions_for_PJIVE.R"))
source(file.path(prog.dir, "Functions_for_CJIVE.R"))

ajive.dir = "H:/My Documents/Applications2/r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb.cols = cbPalette[c(2:5,7)]

norm.print = function(x) paste0(format(round(mean(x),2), nsmall = 2),
                                " (", format(round(sd(x),2), nsmall = 2), ")")

########################################################################################################################
########################################################################################################################
######Rank 1: p1 = 100, p2 = 1000

########################################################################################################################
#############################Simulation Results for Joint Variance Explained = 0.5 in both data sets (Total R^2 = 0.75)#
########################################################################################################################
sim.results = GetSimResults_Dir(file.path(results.dir), 100, 1000)

##############################################################################################################################
#############################Make Plots and Images    ########################################################################
##############################################################################################################################
Time.Summary.Table = summary(tableby(formula =  ~ R.JIVE_Time + AJIVE_Time +
                               ProJIVE_Time + GIPCA_Time + dCCA_Time,
                             data = sim.results))

sim.gg.data = ConvSims_gg_ProJIVE(sim.results)

sim.scores = sim.gg.data[grep("Score", sim.gg.data[,"Type"]),]
sim.loads = sim.gg.data[grep("Loads", sim.gg.data[,"Type"]),]

labs.nm.s = levels(sim.gg.data$Type)[grep("Score",levels(sim.gg.data$Type))][c(3,1,2)]
labs.ex.s = c("Joint Subj Scores", expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]))
names(labs.ex.s) = labs.nm.s

names(cb.cols) = levels(sim.scores$Method)[-1]
score.plot =
  sim.scores %>%
  filter(Method != levels(sim.scores$Method)[1]) %>%
  transform(Method  = factor(as.numeric(Method), levels = 2:6, labels = levels(Method)[-1])) %>%
  ggplot(aes(x = Type, y = Norm)) +
  geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = TRUE, 
               fatten = 0.5) +
  labs(y = "Chordal Norm", x = "Type") +
  scale_x_discrete(limits = labs.nm.s, labels = labs.ex.s) +
  scale_fill_manual(values=cb.cols) +
  theme_bw() + coord_cartesian(ylim = c(0, 1)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = 12),
        text = element_text(size = 15))
plot.name = paste("SimFeng_ScorePlot", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr, plot.name))
print(score.plot)
dev.off()

labs.nm.l = levels(sim.gg.data$Type)[grep("Load",levels(sim.gg.data$Type))][c(3,4,1,2)]
labs.ex.l = c(expression("Joint Variable Loadings"*"X"[1]),expression("Joint Variable Loadings"*"X"[2]),
            expression("Indiv Variable Loadings"*"X"[1]), expression("Indiv Variable Loadings"*"X"[2]))
names(labs.ex.l) = labs.nm.l

load.plot =
  sim.loads %>%
  filter(Method != levels(Method)[1]) %>% 
  transform(Method  = factor(as.numeric(Method), levels = 2:6, labels = levels(Method)[-1])) %>%
  ggplot(aes(x = Type, y = Norm)) +
  geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = TRUE, 
               fatten = 0.5) +
  labs(y = "Chordal Norm", x = "Type") +
  scale_x_discrete(limits = labs.nm.l, labels = labs.ex.l) +
  scale_fill_manual(values=cb.cols) +
  theme_bw() + coord_cartesian(ylim = c(0, 1)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = 12),
        text = element_text(size = 15))

plot.name = paste("SimFeng_LoadNormPlots", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr, plot.name))
print(load.plot)
dev.off()

labs.nm = levels(sim.gg.data$Type)[c(3,6,7,1,2,4,5)]
labs.ex = c("Joint Subj Scores", expression("Joint Loadings"*"X"[1]), expression("Joint Loadings"*"X"[2]), 
            expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]),
            expression("Indiv Loadings"*"X"[1]), expression("Indiv Loadings"*"X"[2]))
names(labs.ex) = labs.nm

norm.plot =
  sim.gg.data %>%
  filter(Method != levels(Method)[1]) %>% 
  transform(Method  = factor(as.numeric(Method), levels = 2:6, labels = levels(Method)[-1])) %>%
  ggplot(aes(x = Type, y = Norm)) +
  geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = TRUE, 
               fatten = 0.5) +
  labs(y = "Chordal Norm", x = "Type") +
  scale_x_discrete(limits = labs.nm, labels = labs.ex) +
  scale_y_log10(breaks = c(0.2*(0:5), 0.05)) +
  scale_fill_manual(values=cb.cols) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = 12),
        text = element_text(size = 15))


plot.name = paste("SimFeng_ScoreAndLoadNormPlots_wLegend", Sys.Date(), ".pdf", sep = "")
pdf(file = file.path(imgs.fldr, plot.name), width = 5, height = 6)
print(norm.plot)
dev.off()

tab.out =
  sim.gg.data %>% 
  filter(Method != "ProJIVE Oracle") %>%
  aggregate(Norm ~ Method + Type, FUN = norm.print) %>%
  mutate(rJ = 1, p2 = 1000) %>%
  pivot_wider(names_from = Method, values_from = Norm) %>% 
  relocate(rJ, p2)

require(writexl)
write_xlsx(tab.out,
           path = paste0("H:/My Documents/ProJIVE/Results/Simulation_Results/Feng1k_Table_",today(),".xlsx"))


