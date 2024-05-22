library(matrixStats)

DCCA=NULL
JIVE=NULL
R.JIVE=NULL
DISCO=NULL
AJIVE=NULL
OnPLS=NULL
COBE=NULL
for(seed in 1:1000){
  DCCA=rbind(as.numeric(as.matrix(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1.txt',sep='')))[3,]) ,DCCA)
  JIVE=rbind(t(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_JIVE.txt',sep=''))) ,JIVE)
  R.JIVE=rbind(t(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_RJIVE.txt',sep=''))) ,R.JIVE)
  AJIVE=rbind(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_AJIVE.txt',sep='')), AJIVE)
  OnPLS=rbind(as.numeric(t(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_OnPLS.txt',sep='')))),OnPLS)
  DISCO=rbind(t(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_DISCO.txt',sep=''))) , DISCO)
  COBE=rbind(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_COBE.txt',sep='')),COBE)
}


#c(errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2)
round(colMeans(DCCA)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(DCCA)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

round(colMeans(JIVE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(JIVE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

round(colMeans(R.JIVE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(R.JIVE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

AJIVE=as.matrix(AJIVE)
round(colMeans(AJIVE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(AJIVE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

round(colMeans(OnPLS)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(OnPLS)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

round(colMeans(DISCO)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(DISCO)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

COBE=as.matrix(COBE)
round(colMeans(COBE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)
round(colSds(COBE)[c(2,4,1,3,6,8,5,7,10,12,9,11)],3)

## test zero correlations
max(abs(AJIVE[,13:14])) #0

max(abs(R.JIVE[,13]))#1.021627e-12
max(abs(DCCA[,13]))

method='JIVE'
psc=numeric(1000) # proportion of significant correlations between d1 and d2
for(seed in 1:1000){
  pvalue=exp(-scan(paste('siml2_seed',seed,'_p900_angle45_nstd1_',method,'_nlogP.txt',sep='')))
  psc[seed]=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)
}
mean(psc)*100 #in %
sd(psc)*100 #in %



##canonical angle/correlation estimated by D-CCA
p=c(100,600,900,1500)
theta= c(0,45,60,75)
sigma= c(0.1,1,3,4)

rlt1.mean.theta=matrix(0,4,4)
rlt1.sd.theta=matrix(0,4,4)
rlt1.mean.ccor=matrix(0,4,4)
rlt1.sd.ccor=matrix(0,4,4)
for(i in 1:4){
  for(j in 1:4){
    data=matrix(0,1000,2)
    for(k in 1:1000){
      data[k,]=as.matrix( read.table(paste('siml2_theta_ED_seed',k,'_p',p[i],'_angle',theta[j],'_nstd1.txt',sep='') ))[,1]
    }
    rlt1.mean.ccor[i,j]=mean(data[,1])
    rlt1.sd.ccor[i,j]=sd(data[,1])
    rlt1.mean.theta[i,j]=mean(data[,2])
    rlt1.sd.theta[i,j]=sd(data[,2])
  }
}


rlt2.mean.theta=matrix(0,4,4)
rlt2.sd.theta=matrix(0,4,4)
rlt2.mean.ccor=matrix(0,4,4)
rlt2.sd.ccor=matrix(0,4,4)
for(i in 1:4){
  for(j in 1:4){
    data=matrix(0,1000,2)
    for(k in 1:1000){
      data[k,]=as.matrix( read.table(paste('siml2_theta_ED_seed',k,'_p900_angle',theta[j],'_nstd',sigma[i],'.txt',sep='') ))[,1]
    }
    rlt2.mean.ccor[i,j]=mean(data[,1])
    rlt2.sd.ccor[i,j]=sd(data[,1])
    rlt2.mean.theta[i,j]=mean(data[,2])
    rlt2.sd.theta[i,j]=sd(data[,2])
  }
}

cos(c(0,45,60,75)/180*pi)

signif(rlt1.mean.theta,3)
round(rlt1.sd.theta,2)
round(rlt1.mean.ccor,3)
round(rlt1.sd.ccor,3)

signif(rlt2.mean.theta,3)
round(rlt2.sd.theta,2)
round(rlt2.mean.ccor,3)
round(rlt2.sd.ccor,3)


####Tables for Block-GDFM
GDFM.45=NULL
for(seed in 1:1000){
 GDFM.45=rbind(read.table(paste('siml2_seed',seed,'_p900_angle45_nstd1_gdfm.txt',sep='')),GDFM.45)

#  (output=[max(max(abs(phi_y1))), max(max(abs(phi_y2))), ...
#          errF_chi_xy1, errF_chi_xy2, err2_chi_xy1, err2_chi_xy2, ...
#          errF_chi_y1, errF_chi_y2, err2_chi_y1, err2_chi_y2, ...
#          normF_nu_y1_X, normF_nu_y2_X, norm2_nu_y1_X, norm2_nu_y2_X, ...
#          normF_nu_y1_chi, normF_nu_y2_chi, norm2_nu_y1_chi, norm2_nu_y2_chi];  

}

GDFM.45=as.matrix(GDFM.45)

round(colMeans(GDFM.45)[c(7:10,3:6,15:18)],3)
round(colSds(GDFM.45)[c(7:10,3:6,15:18)],3)

## test zero correlations
psc.45.1=numeric(1000)# proportion of significant correlations between d1 and d2
psc.45.2=numeric(1000)
for(seed in 1:1000){
  #weekly common
  pvalue=exp(-scan(paste('siml2_seed',seed,'_p900_angle45_nstd1_gdfm_nlogP_1.txt',sep='')))
  psc.45.1[seed]=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)
  
  #weekly common + weekly idiosyncratic
  pvalue=exp(-scan(paste('siml2_seed',seed,'_p900_angle45_nstd1_gdfm_nlogP_2.txt',sep='')))
  psc.45.2[seed]=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)
}


mean(psc.45.1)*100 #in %
sd(psc.45.1)*100 #in %

mean(psc.45.2)*100 #in %
sd(psc.45.2)*100 #in %








