###########################################################################################################################
#######                 Functions for JIVE-PNC Project(s)                              ####################################
#######                 Author: Raphiel J. Murden                                      ####################################
#######                 Supervised by Benjamin Risk                                    ####################################
###########################################################################################################################
require(rootSolve); require(Matrix); require(ggplot2); require(reshape2); require(fields); require(mvtnorm)
require(dplyr); require(xtable); require(optimx); require(gplots); require(MASS); require(r.jive); 
require(extraDistr)

######################################################################################################################
###########   Generates 2 Simulated Datasets that follow JIVE Model using binary subject scores   ####################
######################################################################################################################
GenToyDatBinRank <- function(n, p1, p2, JntVarEx1, JntVarEx2, IndVarEx1, IndVarEx2, jnt_rank = 1, equal.eig = F,
                             ind_rank1 = 2, ind_rank2 = 2, JntVarAdj = T, SVD.plots = T, Error = T, print.cor = T){
  
  #Write out both joint and indiv subject scores for both data sets first
  r.J = jnt_rank
  JntScores = matrix(rbinom(n*r.J, size=1, prob=0.2), nrow = n, ncol = r.J)
  colnames(JntScores) = paste("Jnt Score", 1:r.J)
  
  r.I1 = ind_rank1
  r.I2 = ind_rank2
  
  b = rbinom(n*(r.I1 + r.I2), size=1, prob=0.4)
  b = 1 - 2*b
  IndivScores = matrix(b, nrow = n, ncol = (r.I1 + r.I2))
  colnames(IndivScores) = c(paste("Ind X Score", 1:r.I1), paste("Ind Y Score", 1:r.I2))
  
  if(print.cor){print("The correlation between subject scores is given by"); print(round(cor(cbind(JntScores, IndivScores)),4))}
  
  ##############################Define X Dataset##############################
  ##Then write each of the 3 variable loading vectors for the first dataset
  ##Note that the joint signal has rank=1
  AdjJntLoad.X = matrix(rnorm(r.J*p1), nrow = r.J, ncol = p1)
  
  #change relavent scaling of joint components
  D.J = (equal.eig + 1)*diag(r.J:1) + equal.eig*diag(rep(1,r.J))
  JX = JntScores%*%sqrt(D.J)%*%AdjJntLoad.X
  
  if(SVD.plots){
    plot(svd(JX)$d, ylab = "Singular Values")
    title("SVD of Joint Signal from X")
  }
  
  ##Note that the individual signal has rank = 2 as well
  IndScores.X = IndivScores[,1:r.I1]
  IndLoad.X = matrix(rnorm(n = p1*r.I1), nrow = r.I1, ncol = p1)
  D.IX = (equal.eig + 1)*diag((r.I1:1)) + equal.eig*diag(rep(1,r.I1))
  
  IX = IndScores.X%*%sqrt(D.IX)%*%IndLoad.X
  
  if(SVD.plots){
    plot(svd(IX)$d, ylab = "Singular Values")
    title("SVD of Individual Signal from X")
  }
  
  AX = JX + IX
  
  ##############################Define Y Dataset##############################
  AdjJntLoad.Y = matrix(rnorm(r.J*p2), nrow = r.J, ncol = p2)
  
  ##Note that the joint signal has rank = 3
  JY = JntScores%*%sqrt(D.J)%*%AdjJntLoad.Y
  
  if(SVD.plots){
    plot(svd(JY)$d, ylab = "Singular Values")
    title("SVD of Joint Signal from Y")
  }
  
  IndScores.Y = IndivScores[,(r.I1 + 1:r.I2)]
  IndLoad.Y = matrix(rnorm(r.I2*p2), nrow = r.I2, ncol = p2)
  D.IY = (equal.eig + 1)*diag((r.I2:1)) + equal.eig*diag(rep(1,r.I2)) 
  
  ##Note that the individual signal has rank=2
  IY = IndScores.Y%*%sqrt(D.IY)%*%IndLoad.Y
  
  if(SVD.plots){
    plot(svd(IY)$d, ylab = "Singular Values")
    title("SVD of Individual Signal from Y")
  }
  
  ##Error matrix
  EX = matrix(rnorm(n*p1), nrow=n, ncol=p1)*Error
  
  ##Error matrix
  EY = matrix(rnorm(n*p2), nrow=n, ncol=p2)*Error
  
  Dat.X = AdjSigVarExp(JX, IX, EX, JntVarEx1, IndVarEx1)
  JX = Dat.X$J
  IX = Dat.X$I
  
  Dat.Y = AdjSigVarExp(JY, IY, EY, JntVarEx2, IndVarEx2)
  JY = Dat.Y$J
  IY = Dat.Y$I
  
  Blocks = list(Dat.X[["Data"]], Dat.Y[["Data"]])
  
  Dat.Comps = list(JX, JY, IX, IY, EX, EY)
  names(Dat.Comps) = c("J1", "J2", "I1", "I2", "E1", "E2" )
  
  Scores = list(JntScores, IndScores.X, IndScores.Y)
  names(Scores) = c("Joint", "Indiv_1", "Indiv_2")
  
  Loadings = list(AdjJntLoad.X, IndLoad.X, AdjJntLoad.Y, IndLoad.Y)
  names(Loadings) = c("Joint_1", "Indiv_1", "Joint_2", "Indiv_2")
  
  out = list(Dat.Comps, Blocks, Scores, Loadings)
  names(out) = c("Data Components", "Data Blocks", "Scores", "Loadings")
  
  out
}


####################################################################################
###########   Adjust Dataset Components to get Desired R^2 Values  #################
####################################################################################
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
  
  JJ = MatVar(J)
  II = MatVar(I)
  EE = MatVar(E)
  JE = sum(diag(J%*%t(E)))
  IE = sum(diag(I%*%t(E)))
  R_J = JntVarEx
  R_I = IndVarEx
  
  parms = c(JJ, II, EE, JE, IE, R_J, R_I)
  
  A = J + I
  AA = MatVar(A)
  EE = MatVar(E)
  AE = sum(diag(A%*%t(E)))
  
  ##Desired Total Variance Explained
  d0 = IndVarEx + JntVarEx
  a = AA*(1 - d0)
  b = -2*AE
  c = -d0*EE
  
  d.A = (-b+sqrt(b^2 - 4*a*c))/(2*a)
  
  start = c(0.5, 0.5)*d.A
  
  roots = multiroot(simul.quads, start, parms = parms)
  c = roots$root[1]
  d = roots$root[2]
  
  J.d = d*J; I.c = c*I;
  Dat = J.d + I.c + E
  
  res = list(J.d, I.c, Dat)
  names(res) = c("J", "I", "Data")
  return(res)
}


####################################################################################
###########           Frobenius Norm of Data Matrix Values         #################
####################################################################################
MatVar2 = function(X){
sum(X^2)
}


##################################################################
###########           Sum of Eigenvalues         #################
##################################################################
MatVar = function(X){
  sum(diag(X%*%t(X)))
}


#####################################################################################################
###########          Use PC Scores to calculate joint rank based on cancor         ##################
#####################################################################################################
perm.jntrank <- function(dat.blocks, signal.ranks = NULL, nperms = 500, perc.var = 0.95, alpha = 0.05, center = T){
  n.r.1 = nrow(dat.blocks[[1]])
  n.r.2 = nrow(dat.blocks[[2]])
  if(n.r.1 != n.r.2){stop("The number of rows in each data matrix must match")}
  n = n.r.1
  
  K = length(dat.blocks)
  ##Column center data blocks
  if(center){
    cent.blocks = lapply(dat.blocks, function(D){scale(D, scale = F)})
  }  else{
    cent.blocks = dat.blocks
  }
  
  if(is.null(signal.ranks)){
    all.singvals = sapply(cent.blocks, function(x) svd(x$d))
    for(k in 1:K){
      d = all.singvals
      signal.ranks = c(signal.ranks, which.min(cumsum(d^2) <= perc.var*sum(d^2)) )
    }
  }
  
  x.all.svd = list()
  for(k in 1:K){
    x.all.svd[[k]] = svd(cent.blocks[[k]], nu = signal.ranks[k], nv = signal.ranks[k])
  }

  U.all = lapply(x.all.svd, function(x) scale(x$u, scale = FALSE))

  U1tU2 = t(U.all[[1]])%*%U.all[[2]]
  orig.corrs = svd(U1tU2)$d
  
  perm.corrs = NULL
  
  for (i in 1:nperms){
    U1tU2.perm = t(U.all[[1]])%*%U.all[[2]][sample(n),]
    perm.svd = svd(U1tU2.perm)
    perm.corrs = rbind(perm.corrs, perm.svd$d)
  }
  
  test.val = quantile(perm.corrs[,1],probs =  1-alpha)
  r.j = which.min(orig.corrs >= test.val) - 1
  ranks = c(r.j, signal.ranks)
  p.vals = NULL
  for(i in 1:r.j){
    p.vals = c(p.vals, sum(orig.corrs[i] <= perm.corrs[,1])/nperms)
    }
  names(ranks) = c("Joint", "Total Signal 1", "Total Signal 2")
  res = list(ranks, orig.corrs[1:r.j], p.vals)
  names(res) = c("Ranks", "Canonical Corrs", "P-values")
  return(res)
}

######################################################################################
###########       Chordal norm for matrices with diff ranks         ##################
######################################################################################
chord.norm.diff = function(X,Y, tol = 1E-8){
  svd.X = svd(X)
  svd.Y = svd(Y)
  if(svd.Y$d[1]>0){
    Uy = svd.Y$u
  } else {
      Uy = matrix(0, dim(svd.Y$u))
  }
  if(svd.X$d[1]>0){
    Ux = svd.X$u
  } else {
    Ux = matrix(0, dim(svd.X$u))
  }
  k = max(min(sum(svd.X$d > tol), sum(svd.Y$d > tol)),1)
  inner.prod = t(Ux)%*%Uy
  sig = round(svd(inner.prod)$d[1:k], 8)
  sqrt(sum((1 - sig^2))/k)
}

#############################################################################################################
###########   CC.JIVE uses permutation test based on PC scores to find joint rank          ##################
###########   and estimates joint subject scores as scaled average of canonical variables   #################
#############################################################################################################
cc.jive<-function(dat.blocks, signal.ranks = NULL, joint.rank = 1, perc.var = 0.95, perm.test = T, center = F, nperms = 1000){

  n.1 = nrow(dat.blocks[[1]])
  n.2 = nrow(dat.blocks[[2]])
  if(n.1 != n.2){stop("The number of rows in each data matrix must match")}
  n = n.1
  
  if(center){
    cent.blocks = lapply(dat.blocks, scale)
  } else {
    cent.blocks = dat.blocks
  }
    
  if(is.null(signal.ranks)){
    d1 = svd(cent.blocks[[1]])$d
    d2 = svd(cent.blocks[[2]])$d
    
    r.1 = which.min(cumsum(d1^2) <= perc.var*sum(d1^2)) 
    r.2 = which.min(cumsum(d2^2) <= perc.var*sum(d2^2))
    
    signal.ranks = c(r.1,r.2)
  }
  
  r.1 = signal.ranks[1]
  r.2 = signal.ranks[2]
  
  r.J = joint.rank
  
  if(perm.test){
    cca.perm.res = perm.jntrank(cent.blocks, signal.ranks = c(r.1,r.2), nperms = nperms)
    r.J = cca.perm.res$Ranks[1]
  }
  
  x1.svd = svd(cent.blocks[[1]], nu = signal.ranks[1], nv = signal.ranks[1])
  x2.svd = svd(cent.blocks[[2]], nu = signal.ranks[2], nv = signal.ranks[2])
  
  U.1 = x1.svd$u
  U.2 = x2.svd$u
  
  if(r.J > 0){
    U1tU2 = t(U.1)%*%U.2
    U1tU2.svd = svd(U1tU2, nu = r.J, nv = r.J)
    orig.corrs = U1tU2.svd$d
    
    if(r.J == 1){
      jnt.scores = (U.1%*%U1tU2.svd$u +U.2%*%U1tU2.svd$v)*(2*(1 + orig.corrs[r.J]))^-.5
    } else {
      jnt.scores = (U.1%*%U1tU2.svd$u +U.2%*%U1tU2.svd$v)%*%diag((2*(1 + orig.corrs[1:r.J]))^-.5)
    }
  } else {
    jnt.scores = rep(0, nrow(cent.blocks[[1]]))
    U1tU2 = t(U.1)%*%U.2
    U1tU2.svd = svd(U1tU2)
    orig.corrs = U1tU2.svd$d
  }
  sjive.out = sjive(cent.blocks, c(r.1, r.2), r.J, jnt.scores)
    
  if(perm.test){
    res = list(r.J, jnt.scores, list(orig.corrs,cca.perm.res$'P-values'), list(U1tU2.svd$u, U1tU2.svd$v), signal.ranks)
    names(res) = c("Jnt_Rank","Jnt_Scores", "Canonical_Correlations", "Loadings", "Signal Ranks")
  } else {
    res = list(r.J, jnt.scores, orig.corrs, list(U1tU2.svd$u, U1tU2.svd$v), signal.ranks)
    names(res) = c("Jnt_Rank","Jnt_Scores", "Canonical_Correlations", "Loadings", "Signal Ranks")
  } 
  
  res.out = list(res, sjive.out)
  names(res.out) = c("CanCorRes","sJIVE")
  return(res.out)
}

#############################################################################################################
#########   Compute predicted joint scores for new subjects using current CCJIVE Joint scores  s #############
#############################################################################################################
cc.jive.pred<-function(orig.dat.blocks, new.subjs, joint.rank = 1, tot.signal.ranks, cc.jive.loadings){
  
  r1 = tot.signal.ranks[1]
  r2 = tot.signal.ranks[2]
  X1.svd = svd(orig.dat.blocks[[1]], nu = r1, nv = r1)
  X2.svd = svd(orig.dat.blocks[[2]], nu = r2, nv = r2)
  
  U1.Ldngs = cc.jive.loadings[[1]]
  U2.Ldngs = cc.jive.loadings[[2]]
  
  Pred.CanVar.1 = new.subjs[[1]]%*%X1.svd$v%*%diag(X1.svd$d[1:r1]^-1)%*%U1.Ldngs
  Pred.CanVar.2 = new.subjs[[2]]%*%X2.svd$v%*%diag(X2.svd$d[1:r2]^-1)%*%U2.Ldngs
  pred.jnt.scores = sqrt(1/2)*(Pred.CanVar.1 + Pred.CanVar.2)
  
  return(pred.jnt.scores)
}

#####################################################################################p########################
##############            Retrieve Simulation Results stored in a directory                ##################
#############################################################################################################
GetSimResults_Dir = function(sim.dir, p1, p2, Preds=FALSE){
#  setwd(sim.dir)
  files = list.files(sim.dir, pattern = ".csv")
  num_sims = length(files)
  
  JVEs = as.numeric(paste("0",strsplit(sim.dir, "[.]0")[[1]][-1], sep = "."))
  JntVarEx1 = JVEs[1]
  JntVarEx2 = JVEs[2]
  
  results = NULL
  for (nm in files){
    temp = read.csv(file = file.path(sim.dir, nm), header = T)
    varnames = temp[,1]; values = temp[,2]
    to.keep = !grepl("Pred",varnames)
    if(Preds==F & length(to.keep)>0){
      values = values[to.keep]
      varnames = varnames[to.keep]
      }
    results = rbind(results,c(values, p1, p2))
  }
  
  colnames(results) = c(varnames, "p1", "p2")
  
  results = data.frame(cbind(results, JntVarEx1, JntVarEx2))
}


###############################################################
############    Convert Vector to Network       ###############
###############################################################
vec2net = function(invector) {
  #invector: p x 1, where p is the number of edges
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[upper.tri(outNet,diag = F)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  #  diag(outNet) = 1
  outNet
}

vec2net.2 = function(invector) {
  #invector: p x 1, where p is the number of edges
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet,diag = F)] = invector
  outNet = t(outNet)
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  #  diag(outNet) = 1
  outNet
}

#############################################################################################################
############    Visually display heatmap of matrix: adapted from Erick Lock's show.image  ###################
############    function employed in r.jive package                                       ###################
#############################################################################################################
show.image.2<-function (Image, ylab = "", xlab="",net=F,main="", sub="", colorbar = T) 
{
  #if(net){Image=Image-diag(Image)}
  
  #Image = Image*upper.tri(Image)
  lower = mean(Image) - 5 * sd(Image)
  upper = mean(Image) + 5 * sd(Image)
  Image[Image < lower] = lower
  Image[Image > upper] = upper
  if(colorbar){
    image.plot(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
               zlim = c(lower, upper), axes = FALSE, col = bluered(100), 
               xlab = xlab, ylab = ylab)
    title(main=main, sub = sub)
  }
  
  else if(!colorbar){
    image(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
               zlim = c(lower, upper), axes = FALSE, col = bluered(100), 
               xlab = xlab, ylab = ylab)
    title(main=main, sub = sub)
  }
  
  if(net){
    mod = read.table("C:/Users/Raphiel\'s PC/Dropbox/JIVE-PNC/PNC-Data/roi_id.txt")
    mod = mod[,1]
    
    Count.m = table(sort(mod))
    k = 0
    Grid.m = vector()
    for(i in 1:(length(Count.m)-1))
    {   
      k = k + Count.m[i]
      Grid.m[i] = k
    }
    abline(v=c(0,Grid.m[-1]-32,232)+.5,
           h=c(0,232-Grid.m[-1]+32,232)+.5)
  }
}

#####################################################
###########          sJIVE         ##################
#####################################################
sjive<-function (blocks, initial_signal_ranks, joint_rank, joint_scores = NULL) 
{
  
  j=joint_rank
  K <- length(blocks)
  if (K < 2) {
    stop("ajive expects at least two data matrices.")
  }
  if (sum(sapply(blocks, function(X) any(is.na(X)))) > 0) {
    stop("Some of the blocks has missing data -- ajive expects full data matrices.")
  }
  
  block_svd <- list()
  scaled.centered <- list()
  
  for (k in 1:K) {
    
    #Scale and center each block
    temp<-blocks[[k]]
    
    for (i in 1:dim(blocks[[k]])[2]) {
      temp[,i]<-blocks[[k]][,i]-mean(blocks[[k]][,i]) ##subtract the mean from each column
    }
    
    scaled.centered[[k]] <- temp
    block_svd[[k]] <- svd(scaled.centered[[k]],nu=initial_signal_ranks[k]) ##Take SVD of each data block
    
    if (k==1) {
      stacked_mat<-block_svd[[k]]$u[,1:initial_signal_ranks[k]] ##initialize matrix that results from stacking SVD matrices: contains first k vectors of U
    }
    else if (k>1){
      stacked_mat <- cbind(stacked_mat,block_svd[[k]]$u[,1:initial_signal_ranks[k]]) ##stack the rest of the U matrices together with the first
    } 
  }
  if(is.null(joint_scores)){
    if(j>0){
      joint_bases<-svd(stacked_mat,nu=j)$u ##take SVD of stacked matrix: nu=j indicates that the joint structure has rank=j
      joint_proj<-(joint_bases)%*%t(joint_bases) ##orthogonal projection matrix onto joint space based on stacked bases
    } else {
      joint_bases = rep(0, nrow(scaled.centered[[1]]))
      joint_proj<-(joint_bases)%*%t(joint_bases) 
    }
  } else {
    joint_proj<-(joint_scores)%*%t(joint_scores) 
  }
  indiv_structure <- list()
  joint_structure <- list()
  
  for (k in K:1) {
    joint_structure[[k]] = list()
    joint_structure[[k]][['full']]<-joint_proj%*%as.matrix(scaled.centered[[k]])
    temp.svd = svd(joint_structure[[k]][['full']], nu = joint_rank, nv = joint_rank)
    
    joint_structure[[k]][['u']]<-temp.svd[['u']]
    joint_structure[[k]][['v']]<-temp.svd[['v']]
    joint_structure[[k]][['d']]<-temp.svd[['d']][1:joint_rank]
    
    dat_proj<-block_svd[[k]]$u%*%t(block_svd[[k]]$u)
    # indiv_structure[[k]]<-(dat_proj-joint_proj)%*%as.matrix(scaled.centered[[k]])
    temp = (diag(nrow(scaled.centered[[k]]))-joint_proj)%*%as.matrix(scaled.centered[[k]])
    indiv.rank = initial_signal_ranks[k]-joint_rank
    temp.svd = svd(temp)
    
    indiv_structure[[k]] = list() 
    
    indiv_structure[[k]][['u']] = temp.svd[['u']][,1:indiv.rank, drop = FALSE]
    indiv_structure[[k]][['v']] = temp.svd[['v']][,1:indiv.rank, drop = FALSE]
    indiv_structure[[k]][['d']] = temp.svd[['d']][1:indiv.rank]
    indiv_structure[[k]][['full']] = indiv_structure[[k]][['u']]%*%
      diag(indiv_structure[[k]][['d']], nrow = indiv.rank, ncol = indiv.rank)%*%
      t(indiv_structure[[k]][['v']])
  }
  
  out = list(joint_structure,indiv_structure,stacked_mat,joint_proj)
  names(out) = c("joint_matrices", "indiv_matrices", "stacked_mat", "joint_projection")
  return(out)
}


#######################################################################################################################
#####################         Convert Simulation Results to a form that will allow ggplot       #######################
#######################################################################################################################
ConvSims_gg<-function(AllSims){
  sim.names = colnames(AllSims)
  #CompEvals.names = c(sim.names[grep("Scores", sim.names)], sim.names[grep("Loads", sim.names)])
  num.sims = dim(AllSims)[1]
  
  p1 = unique(AllSims$p1)
  p2 = unique(AllSims$p2)
  
  JVE_1 =  c(as.matrix(AllSims[,grep("JntVarEx1",sim.names)]))
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))
  JVE_1 = factor(JVE_1, labels = JVE1.labs, levels = c(0.05, 0.5))
  
  JVE_2 = c(as.matrix(AllSims[,grep("JntVarEx2",sim.names)]))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)), bquote("R"[J2]^2*"=0.5, p"[2]*"="*.(p2)) )
  JVE_2 = factor(JVE_2, labels = JVE2.labs, levels = c(0.05, 0.5))
  
  IVE_1 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.X", sim.names)]))
  IVE_1 = round(as.numeric(IVE_1), 2)
  
  IVE_2 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.Y", sim.names)]))
  IVE_2 = round(as.numeric(IVE_2), 2)
  
  if("True.Joint.Rank" %in% sim.names){
    sim.names.2 = sim.names[-grep("True.Joint.Rank", sim.names)]
    AllSims.2 = AllSims[,which(sim.names %in% sim.names.2)]
    Method = sim.names.2[grep("Joint.Rank", sim.names.2)]
    temp = colnames(AllSims.2)
    Rank = factor(c(as.matrix(AllSims.2[,temp[grep("Joint.Rank", temp)]])))
  } else {
    Method = sim.names[grep("Joint.Rank", sim.names)]
    Rank = factor(c(as.matrix(AllSims[,grep("Joint.Rank", sim.names)])))
  }
  
  Method = sub(".Joint.Rank", "", Method, fixed = T);  Method = sub("aJI", "AJI", Method, fixed = T)
  Method = sub("Correct", "", Method, fixed = T); Method = sub(".Wrong", "-Over", Method, fixed = T)
  Method = sub(".oracle", "-Oracle", Method, fixed = T); Method = sub(".over", "-Over", Method, fixed = T)
  Method = sub("Elbow", "Oracle", Method, fixed = T); Method = sub(".95.", "-Over", Method, fixed = T)
  Method = sub("iJI", "R.JI", Method, fixed = T); Method = sub("AJIVE.", "AJIVE-", Method, fixed = T)
  Method = sub("CC", "CJIVE", Method, fixed = T); Method = sub("cJIVE", "CJIVE", Method, fixed = T)
  Method[which(Method == "CJIVE.Oracle")] = "CJIVE-Oracle"
  Method[which(Method %in% c("CJIVE.Over", "cJIVE.Over"))] = "CJIVE-Over"
  Method = as.factor(Method)
  
  sim.ranks.gg = data.frame(Rank = Rank, Method = rep(Method, each = num.sims), JVE_1 = rep(JVE_1, each = length(levels(Method))),
                            JVE_2 = rep(JVE_2, each = length(levels(Method))), 
                            IVE_1 = rep(IVE_1, each = length(levels(Method))), 
                            IVE_2 = rep(IVE_2, each = length(levels(Method))))
  if("True.Joint.Rank" %in% sim.names){
    True_Rank = c(as.numeric(as.matrix(AllSims[,grep("True.Joint.Rank", sim.names)])))
    True_Rank = factor(True_Rank)
    sim.ranks.gg$True_Rank = rep(True_Rank, each = length(levels(Method)))
  }
  
  Norm = c(as.numeric(as.matrix(AllSims[,grep("Scores", sim.names)])))
  
  Method = rep(substr(sim.names[grep("Scores",sim.names)], 1,12), each = num.sims)
  Method = sub(".Joint.", "", Method, fixed = T); Method = sub(".Indiv.", "", Method, fixed = T);
  # Method = sub(".Joi", "", Method); Method = sub("E.w.J", "E.w", Method)
  # Method = sub(".Ind", "", Method); Method = sub("E.w.I", "E.w", Method)
  Method = sub("aJI", "AJI", Method, fixed = T); Method = sub("iJI", "R.JI", Method, fixed = T);
  Method = sub("r.J", "r", Method, fixed = T); Method = sub("r.I", "r", Method, fixed = T)
  Method[which(Method == "AJIVE")] = "AJIVE-Oracle"; Method = sub(".w", "-Over", Method, fixed = T); 
  Method[which(Method == "CC.Oracle")] = "CJIVE-Oracle"; Method[which(Method == "CC-Over")] = "CJIVE-Over"
  Method = sub("JIVE.", "JIVE-", Method, fixed = T); Method = sub("[[:punct:]]$", "", Method, fixed = T)
  Method = as.factor(Method)
  
  Type = rep(substring(sim.names[grep("Scores", sim.names)], 7), each = num.sims)
  Type = sub("w.", "", Type)
  Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type)
  Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type)
  Type = gsub("Over ", "", Type)
  Type = gsub("Ora", "", Type)
  Type = factor(Type)
  
  n.levs = nlevels(Type)*nlevels(Method)
  sim.score.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1, each = n.levs), 
                                  JVE_2 = rep(JVE_2, each = n.levs),
                                  IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  if("True.Joint.Rank" %in% sim.names){
    True_Rank = c(as.numeric(as.matrix(AllSims[,grep("True.Joint.Rank", sim.names)])))
    True_Rank = factor(True_Rank)
    sim.score.norms.gg$True_Rank = rep(True_Rank, each = length(levels(Method)))
  }
  
  Norm = c(as.numeric(as.matrix(AllSims[,grep("Loads", sim.names)])))
  
  Method = rep(substr(sim.names[grep("Loads",sim.names)], 1,12), each = num.sims)
  Method = sub(".Joint.", "", Method, fixed = T); Method = sub(".Indiv.", "", Method, fixed = T);
  # Method = sub(".Joi", "", Method); Method = sub("E.w.J", "E.w", Method)
  # Method = sub(".Ind", "", Method); Method = sub("E.w.I", "E.w", Method)
  Method = sub("aJI", "AJI", Method, fixed = T); Method = sub("iJI", "R.JI", Method, fixed = T);
  Method = sub("r.J", "r", Method, fixed = T); Method = sub("r.I", "r", Method, fixed = T)
  Method[which(Method == "AJIVE")] = "AJIVE-Oracle"; Method = sub(".w", "-Over", Method, fixed = T); 
  Method[which(Method == "CC.Oracle")] = "CJIVE-Oracle"; Method[which(Method == "CC-Over")] = "CJIVE-Over"
  Method = sub("JIVE.", "JIVE-", Method, fixed = T); Method = sub("[[:punct:]]$", "", Method, fixed = T)
  Method = as.factor(Method)
  
  Type = factor(rep(substring(sim.names[grep("Loads", sim.names)], 7), each = num.sims))
  Type = sub("w.", "", Type)
  Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type)
  Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type)
  Type = gsub("Over ", "", Type)
  Type = gsub("Ora", "", Type)
  Type = factor(Type)
  
  n.levs = nlevels(Type)*nlevels(Method)
  sim.load.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1, each = n.levs), 
                                 JVE_2 = rep(JVE_2, each = n.levs), IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))
  
  if("True.Joint.Rank" %in% sim.names){
    True_Rank = c(as.numeric(as.matrix(AllSims[,grep("True.Joint.Rank", sim.names)])))
    True_Rank = factor(True_Rank)
    sim.load.norms.gg$True_Rank = rep(True_Rank, each = length(levels(Method)))
  }
  
  sim.all.norms.gg = rbind(sim.score.norms.gg, sim.load.norms.gg)
  
  out = list(sim.ranks.gg, sim.all.norms.gg)
  names(out) = c("Ranks", "Subj and Ldg Norms")
  out
}


##############Author: Ben Risk, PhD
# Function for plotting networks with ggplot
create.graph.long = function(gmatrix,sort_indices=NULL) {
  nnode = nrow(gmatrix)
  X1 = c(1:nnode)%x%rep(1,nnode)
  X2 =  rep(1,nnode)%x%c(1:nnode)
  X1 = factor(X1, )
  if (!is.null(sort_indices)) {
    gmatrix = gmatrix[sort_indices,sort_indices]
  }
  value = as.vector(as.matrix(gmatrix))
  data.frame(X1,X2,value)
}

#################################################################################################
#####################         Plot CJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.norm.plot<-function(norm.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)]
  labs.ex = c("Joint Subj Scores", expression("Joint Loadings"*"X"[1]), expression("Joint Loadings"*"X"[2]), 
              expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]),
              expression("Indiv Loadings"*"X"[1]), expression("Indiv Loadings"*"X"[2]))
  ggplot(data = norm.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)], labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}

#################################################################################################
#####################         Plot CJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.score.norm.plot<-function(norm.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[grep("Score",levels(norm.dat$Type))][c(3,1,2)]
  labs.ex = c("Joint Subj Scores", expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]))
  ggplot(data = norm.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs, labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}

#################################################################################################
#####################         Plot CJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.load.norm.plot<-function(norm.dat, cols, show.legend = F, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[grep("Load",levels(norm.dat$Type))][c(3,4,1,2)]
  labs.ex = c(expression("Joint Variable Loadings"*"X"[1]),expression("Joint Variable Loadings"*"X"[2]),
              expression("Indiv Variable Loadings"*"X"[1]), expression("Indiv Variable Loadings"*"X"[2]))
  ggplot(data = norm.dat, aes(x = Type, y = Norm)) +
    geom_boxplot(aes(fill = Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    # geom_boxplot(aes(color = Method),
    #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
    #              show.legend = F) +
    labs(y = "Chordal Norm", x = "Type") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_x_discrete(limits = labs, labels = labs.ex) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + coord_cartesian(ylim = c(0, y.max)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = element_text(size = text.size))
}

#################################################################################################
#####################         Plot CJIVE Norms from Simulation Study      #######################
#################################################################################################
gg.corr.plot<-function(cor.dat, cols, show.legend = F, text.size){
  ggplot(data = cor.dat, aes(x = Component, y = Prediction_Correlation)) +
    geom_boxplot(aes(fill = Type), position = "dodge", outlier.alpha = 0, show.legend = show.legend, fatten = 0.5) +
    # geom_boxplot(aes(color = Type, fill = Type),
    #              fatten = NULL, coef = 0, outlier.alpha = 0,
    #              show.legend = T) +
    labs(y = "Absolute Pearson Correlation", x = "Component") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_fill_manual(values=cols) +
    scale_colour_manual(values=cols) + 
    theme_bw() + 
    theme(text = element_text(size = text.size))
}

#################################################################################################
#####################         Plot Ranks CJIVE from Simulation Study      #######################
#################################################################################################
gg.rank.plot<-function(rank.dat, cols, show.legend = F, text.size){
  ggplot(data = rank.dat, aes(x = Rank, y = n/nrow(AllSims), fill=Method)) +
    geom_bar(position = position_dodge(), stat = "identity", width = 0.6, show.legend = show.legend) +
    labs(y = "Proportion", x = "Selected Rank") +
    facet_grid(JVE_2 ~ JVE_1, labeller = label_parsed) +
    scale_fill_manual(values=cols) + 
    theme_bw() + 
    theme(text = element_text(size = text.size))
}

############################################################c#####################################
##############   Box-plots of chordal norms from CJIVE from Simulation Study      ################
##################################################################################################
Melt.Sim.Cors<-function(sim.dat,r.J,p1,p2){
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)), bquote("R"[J2]^2*"=0.5, p"[1]*"="*.(p2)))
  
  A.PredCors.melt = melt(sim.dat, id.vars = c("JntVarEx1", "JntVarEx2"), 
                         measure.vars = paste("AJIVE_Pred_Corr", 1:r.J, sep = ""), variable_name = "Component")
  colnames(A.PredCors.melt)[which(colnames(A.PredCors.melt)=="variable")]="Component"
  colnames(A.PredCors.melt)[which(colnames(A.PredCors.melt)=="value")]="Prediction_Correlation"
  
  levels(A.PredCors.melt$Component) = paste("No.", 1:r.J)
  A.PredCors.melt$Type = "AJIVE"
  
  CC.PredCors.melt = melt(sim.dat, id.vars = c("JntVarEx1", "JntVarEx2"), 
                          measure.vars = paste("CC_Pred_Corr", 1:r.J, sep = ""), variable_name = "Component", 
                          value.name = "Prediction_Correlation")
  colnames(CC.PredCors.melt)[which(colnames(CC.PredCors.melt)=="variable")]="Component"
  colnames(CC.PredCors.melt)[which(colnames(CC.PredCors.melt)=="value")]="Prediction_Correlation"
  
  levels(CC.PredCors.melt$Component) = paste("No.", 1:r.J)
  CC.PredCors.melt$Type = "CJIVE"
  
  Pred.Cors.melt = rbind(A.PredCors.melt, CC.PredCors.melt)
  Pred.Cors.melt$Prediction_Correlation = abs(Pred.Cors.melt$Prediction_Correlation)
  Pred.Cors.melt$JVE_1 = factor(Pred.Cors.melt$JntVarEx1, labels = JVE1.labs)
  Pred.Cors.melt$JVE_2 = factor(Pred.Cors.melt$JntVarEx2, labels = JVE2.labs)
  Pred.Cors.melt
}

############################################################c#####################################
##############         Function to sign-correct and scale variable loadings       ################
##################################################################################################
scale.loadings = function(loading.comp){
  x = loading.comp
  x1 = x/max(abs(x), na.rm = TRUE)
  # x2 = x/max(x, na.rm = TRUE)
  # pos = max(x1, na.rm = TRUE)==max(x2, na.rm = TRUE)
  # ((-1)^(pos))*x1
  pos = sign(skew(x1, na.rm = TRUE))
  ((-1)^(pos))*x1
}
