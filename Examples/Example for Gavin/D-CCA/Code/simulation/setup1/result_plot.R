library("matrixStats")

setwd('/Volumes/My Passport/Hai_Windows/work/Hai_paper1/simulation/setup1_paper_revision1/result')
###############################################
#####plot1
subtitle_1=c("(1a)","(1b)","(1c)")
#subtitle_2=c("(2a)","(2b)","(2c)")
#subtitle_3=c("(3a)","(3b)","(3c)")
x.lab=c(expression('p'[1]),expression(theta[1]),expression(sigma[e]^2))
est.matrix=c(expression("Averge relative norm error of  "*widehat(bold('X'))[1]),expression("Averge relative norm error of  "*widehat(bold('X'))[2]),expression("Averge relative norm error of  "*widehat(bold('C'))[1]),expression("Averge relative norm error of  "*widehat(bold('C'))[2]),expression("Averge relative norm error of  "*widehat(bold('D'))[1]),expression("Averge relative norm error of  "*widehat(bold('D'))[2]),'Proportion')
x.tick=list(c('100','300','600','900','1200','1500'),
            c('0','15','30','45','60','75'),
            c('0.01','0.25','1','4','9','16'))


#pdf(file="plot1.pdf" ,width=9, height=9)
#16*8 inch plot
par(mfrow=c(3,6),mar=c(5, 5, 1, 1))
#par(mfrow=c(2,2))
x6=seq(1,6)

##########################################
######## set1: sigma_sq=1, theta=1, vary p
rlt1=matrix(0,1000,12)
rlt2=matrix(0,1000,12)
rlt3=matrix(0,1000,12)
rlt1.mean=matrix(0,6,12)     
rlt2.mean=matrix(0,6,12)   
rlt3.mean=matrix(0,6,12) 
p=c(100,300,600,900,1200,1500)
for(j in 1:6){
  for(i in 1:1000){
    data=as.matrix( read.table(paste('siml1_seed',i,'_p',p[j],'_angle45_nstd1.txt',sep='') ))[,1:12]
    rlt1[i,]=data[1,]
    rlt2[i,]=data[3,]
    rlt3[i,]=data[7,]
  }
  rlt1.mean[j,]=colMeans(rlt1)
  rlt2.mean[j,]=colMeans(rlt2)
  rlt3.mean[j,]=colMeans(rlt3)
}


for(j in 0:5)
{
  plot(x=x6,y=rlt1.mean[,j*2+1],ylim=c(0.8*min(   min(rlt1.mean[,c(j*2+1,j*2+2)]), min(rlt2.mean[,c(j*2+1,j*2+2)])  ), 1.1*max(   max(rlt1.mean[,c(j*2+1,j*2+2)]), max(rlt2.mean[,c(j*2+1,j*2+2)])  ) ),xlab=x.lab[1],
       ylab=est.matrix[j+1],pch=2,type="o",lty=4,lwd=1,xaxt="n",cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1.5)
  axis(side=1,at=x6,labels=x.tick[[1]])
  #title(sub=subtitle_1[i+1])
  lines(x=x6,y=rlt1.mean[,j*2+2],pch=1,type="o",lty=2,lwd=1)
  lines(x=x6,y=rlt2.mean[,j*2+1],pch=17,type="o",lty=1,lwd=1)
  lines(x=x6,y=rlt2.mean[,j*2+2],pch=16,type="o",lty=3,lwd=1)
}


##########################################
######## set2: sigma_sq=1, p=900, vary theta
rlt1=matrix(0,1000,12)
rlt2=matrix(0,1000,12)
rlt3=matrix(0,1000,12)
rlt1.mean=matrix(0,6,12)     
rlt2.mean=matrix(0,6,12)   
rlt3.mean=matrix(0,6,12)  
theta= c(0,15,30,45,60,75)
for(j in 1:6){
  for(i in 1:1000){
    data=as.matrix( read.table(paste('siml1_seed',i,'_p900_angle',theta[j],'_nstd1.txt',sep='') ))[,1:12]
    rlt1[i,]=data[1,]
    rlt2[i,]=data[3,]
    rlt3[i,]=data[7,]
  }
  rlt1.mean[j,]=colMeans(rlt1)
  rlt2.mean[j,]=colMeans(rlt2)
  rlt3.mean[j,]=colMeans(rlt3)
}


for(j in 0:5)
{
  plot(x=x6,y=rlt1.mean[,j*2+1],ylim=c(0.8*min(   min(rlt1.mean[,c(j*2+1,j*2+2)]), min(rlt2.mean[,c(j*2+1,j*2+2)])  ), 1.1*max(   max(rlt1.mean[,c(j*2+1,j*2+2)]), max(rlt2.mean[,c(j*2+1,j*2+2)])  ) ),xlab=x.lab[2],
       ylab=est.matrix[j+1],pch=2,type="o",lty=4,lwd=1,xaxt="n",cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1.5)
  axis(side=1,at=x6,labels=x.tick[[2]])
  #title(sub=subtitle_1[i+1])
  lines(x=x6,y=rlt1.mean[,j*2+2],pch=1,type="o",lty=2,lwd=1)
  lines(x=x6,y=rlt2.mean[,j*2+1],pch=17,type="o",lty=1,lwd=1)
  lines(x=x6,y=rlt2.mean[,j*2+2],pch=16,type="o",lty=3,lwd=1)
}


##########################################
######## set3: theta=45, p=900, vary sigma_2
rlt1=matrix(0,1000,12)
rlt2=matrix(0,1000,12)
rlt3=matrix(0,1000,12)
rlt1.mean=matrix(0,6,12)     
rlt2.mean=matrix(0,6,12)   
sigma= c(0.1,0.5,1,2,3,4)
for(j in 1:6){
  for(i in 1:1000){
    data=as.matrix( read.table(paste('siml1_seed',i,'_p900_angle45_nstd',sigma[j],'.txt',sep='') ))[,1:12]
    rlt1[i,]=data[1,]
    rlt2[i,]=data[3,]
    rlt3[i,]=data[7,]
  }
  rlt1.mean[j,]=colMeans(rlt1)
  rlt2.mean[j,]=colMeans(rlt2)
  rlt3.mean[j,]=colMeans(rlt3)
}


for(j in 0:5)
{
  plot(x=x6,y=rlt1.mean[,j*2+1],ylim=c(0.8*min(   min(rlt1.mean[,c(j*2+1,j*2+2)]), min(rlt2.mean[,c(j*2+1,j*2+2)])  ), 1.1*max(   max(rlt1.mean[,c(j*2+1,j*2+2)]), max(rlt2.mean[,c(j*2+1,j*2+2)])  ) ),xlab=x.lab[3],
       ylab=est.matrix[j+1],pch=2,type="o",lty=4,lwd=1,xaxt="n",cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1.5)
  axis(side=1,at=x6,labels=x.tick[[3]])
  #title(sub=subtitle_1[i+1])
  lines(x=x6,y=rlt1.mean[,j*2+2],pch=1,type="o",lty=2,lwd=1)
  lines(x=x6,y=rlt2.mean[,j*2+1],pch=17,type="o",lty=1,lwd=1)
  lines(x=x6,y=rlt2.mean[,j*2+2],pch=16,type="o",lty=3,lwd=1)
}



#dev.off()

######################################
####set 2: proportion of matrix norm


p=c(100,300,600,900,1200,1500)
theta= c(0,15,30,45,60,75)
sigma= c(0.1,0.5,1,2,3,4)

rlt1=matrix(0,1000,12)
rlt1.mean=matrix(0,6,12)     

for(j in 1:6){
  for(i in 1:1000){
    data=as.matrix( read.table(paste('siml1_seed',i,'_p',p[j],'_angle45_nstd1.txt',sep='') ))[,1:12]
    rlt1[i,]=data[7,]
    
  }
  rlt1.mean[j,]=colMeans(rlt1)
  
}
rlt1.mean





rlt2=matrix(0,1000,12)
rlt2.mean=matrix(0,6,12)     

for(j in 1:6){
  for(i in 1:1000){
    data=as.matrix( read.table(paste('siml1_seed',i,'_p900_angle',theta[j],'_nstd1.txt',sep='') ))[,1:12]
    rlt2[i,]=data[7,]
    
  }
  rlt2.mean[j,]=colMeans(rlt2)
  
}
rlt2.mean


rlt3=matrix(0,1000,12)
rlt3.mean=matrix(0,6,12)     

for(j in 1:6){
  for(i in 1:1000){
    data=as.matrix( read.table(paste('siml1_seed',i,'_p900_angle45_nstd',sigma[j],'.txt',sep='') ))[,1:12]
    rlt3[i,]=data[7,]
    
  }
  rlt3.mean[j,]=colMeans(rlt3)
  
}
rlt3.mean



## The accuracy of rank selection
rlt1=matrix(0,1000*6,12)
rlt2=matrix(0,1000*6,12)
rlt3=matrix(0,1000*6,12)
p=c(100,300,600,900,1200,1500)
theta= c(0,15,30,45,60,75)
sigma= c(0.1,0.5,1,2,3,4)
for(j in 1:6){
  for(i in 1:1000){
    data1=as.matrix( read.table(paste('siml1_seed',i,'_p',p[j],'_angle45_nstd1.txt',sep='') ))[,1:12]
    data2=as.matrix( read.table(paste('siml1_seed',i,'_p900_angle',theta[j],'_nstd1.txt',sep='') ))[,1:12]
    data3=as.matrix( read.table(paste('siml1_seed',i,'_p900_angle45_nstd',sigma[j],'.txt',sep='') ))[,1:12]
    
    rlt1[i+1000*(j-1),]=data1[4,]
    rlt2[i+1000*(j-1),]=data2[4,]
    rlt3[i+1000*(j-1),]=data3[4,]
  }
}

rlt=rbind(rlt1,rlt2,rlt3)

sum(rlt[,1]==3)/18000  #17921,  0.9956111
sum(rlt[,2]==3)/18000  #17941,  0.9967222
sum(rlt[,3]==1)/18000  #17975,  0.9986111


