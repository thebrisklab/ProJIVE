library('RegularizedSCA')#v0.5.3
source('test_zero_corr.R')
args <- commandArgs(TRUE)
seed=args[1]

Y_1=t(as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_Y_1.txt',sep=''))))
Y_2=t(as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_Y_2.txt',sep=''))))

p_1=dim(Y_1)[2]
p_2=dim(Y_2)[2]

norm.Y_1=norm(Y_1,'f')
norm.Y_2=norm(Y_2,'f')

Y_1.scaled=Y_1/norm.Y_1
Y_2.scaled=Y_2/norm.Y_2

Y_all.scaled=cbind(Y_1.scaled,Y_2.scaled)

R=6

Jk=c(p_1,p_2)

set.seed(0)
invisible(capture.output(  rlt<-DISCOsca(Y_all.scaled, R, Jk)   ))

X_1.DISCO= t(rlt$Trot_best[[1]]%*%t(rlt$Prot_best[[1]][1:p_1,]))  * norm.Y_1
X_2.DISCO= t(rlt$Trot_best[[1]]%*%t(rlt$Prot_best[[1]][(p_1+1):(p_1+p_2),]))  * norm.Y_2


col.match=NULL
for(i in 1:R){
  if(sum(rlt$comdist[[1]][,i]==c(1,1))==2)
  {
    col.match=c(col.match,i)
  }
}


C_1.DISCO= t(rlt$Trot_best[[1]][,col.match]%*%t(rlt$Prot_best[[1]][1:p_1,col.match]))   * norm.Y_1
C_2.DISCO= t(rlt$Trot_best[[1]][,col.match]%*%t(rlt$Prot_best[[1]][(p_1+1):(p_1+p_2),col.match]))  * norm.Y_2


D_1.DISCO=X_1.DISCO-C_1.DISCO
D_2.DISCO=X_2.DISCO-C_2.DISCO


X_1=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_X_1.txt',sep='')))
X_2=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_X_2.txt',sep='')))

C_1=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_C_1.txt',sep='')))
C_2=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_C_2.txt',sep='')))

D_1=X_1-C_1
D_2=X_2-C_2

errF_X_1=base::norm(X_1-X_1.DISCO,"F")/base::norm(X_1,"F")
errF_X_2=base::norm(X_2-X_2.DISCO,"F")/base::norm(X_2,"F")

errF_C_1=base::norm(C_1-C_1.DISCO,"F")/base::norm(C_1,"F")
errF_C_2=base::norm(C_2-C_2.DISCO,"F")/base::norm(C_2,"F")

errF_D_1=base::norm(D_1-D_1.DISCO,"F")/base::norm(D_1,"F")
errF_D_2=base::norm(D_2-D_2.DISCO,"F")/base::norm(D_2,"F")

err2_X_1=base::norm(X_1-X_1.DISCO,"2")/base::norm(X_1,"2")
err2_X_2=base::norm(X_2-X_2.DISCO,"2")/base::norm(X_2,"2")

err2_C_1=base::norm(C_1-C_1.DISCO,"2")/base::norm(C_1,"2")
err2_C_2=base::norm(C_2-C_2.DISCO,"2")/base::norm(C_2,"2")

err2_D_1=base::norm(D_1-D_1.DISCO,"2")/base::norm(D_1,"2")
err2_D_2=base::norm(D_2-D_2.DISCO,"2")/base::norm(D_2,"2")

output=c(errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2,  max(abs(D_1.DISCO%*%t(D_2.DISCO))) )

write.table(output,file=paste('siml1_seed',seed,'_p900_angle45_nstd1_DISCO.txt',sep=''),row.names = FALSE,col.names = FALSE)

#test correlations between D_1.DISCO and D_2.DISCO
p1=nrow(D_1.DISCO)
p2=nrow(D_2.DISCO)
nlogP=numeric(p1*p2)
t=1
for(i in 1:p1){
  for(j in 1:p2){
    nlogP[t]=-log(test_zero_corr(D_1.DISCO[i,],D_2.DISCO[j,]))
    t=t+1
  }
}

write.table(nlogP,file=paste('siml1_seed',seed,'_p900_angle45_nstd1_DISCO_nlogP.txt',sep=''),row.names = FALSE,col.names = FALSE)



