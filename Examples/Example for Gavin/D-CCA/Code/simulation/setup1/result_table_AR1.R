library(matrixStats)
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
      data[k,]=as.matrix( read.table(paste('siml1_AR1_theta_ED_seed',k,'_p',p[i],'_angle',theta[j],'_nstd1.txt',sep='') ))[,1]
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
      data[k,]=as.matrix( read.table(paste('siml1_AR1_theta_ED_seed',k,'_p900_angle',theta[j],'_nstd',sigma[i],'.txt',sep='') ))[,1]
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

