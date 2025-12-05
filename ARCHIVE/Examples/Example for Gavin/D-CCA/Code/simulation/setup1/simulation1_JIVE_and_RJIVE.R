library(r.jive) #v2.1
source('test_zero_corr.R')
args <- commandArgs(TRUE)
seed=args[1]

Y_1=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_Y_1.txt',sep='')))
Y_2=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_Y_2.txt',sep='')))

X_1=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_X_1.txt',sep='')))
X_2=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_X_2.txt',sep='')))

C_1=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_C_1.txt',sep='')))
C_2=as.matrix(read.table(paste('/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision2/data/datafile/siml1_seed',seed,'_p900_angle45_nstd1_C_2.txt',sep='')))

D_1=X_1-C_1
D_2=X_2-C_2

data=list(Y_1,Y_2)
###JIVE without the orthongonal constraint on distinct matrices (JIVE)
set.seed(0)
Results = jive(data,center=F,orthIndiv =F)

Results$converged

C_1.jive=Results$joint[[1]]*Results$scale$`Scale Values`[1]
C_2.jive=Results$joint[[2]]*Results$scale$`Scale Values`[2]

D_1.jive=Results$individual[[1]]*Results$scale$`Scale Values`[1]
D_2.jive=Results$individual[[2]]*Results$scale$`Scale Values`[2]

X_1.jive=C_1.jive+D_1.jive
X_2.jive=C_2.jive+D_2.jive

errF_X_1=base::norm(X_1-X_1.jive,"F")/base::norm(X_1,"F")
errF_X_2=base::norm(X_2-X_2.jive,"F")/base::norm(X_2,"F")

errF_C_1=base::norm(C_1-C_1.jive,"F")/base::norm(C_1,"F")
errF_C_2=base::norm(C_2-C_2.jive,"F")/base::norm(C_2,"F")

errF_D_1=base::norm(D_1-D_1.jive,"F")/base::norm(D_1,"F")
errF_D_2=base::norm(D_2-D_2.jive,"F")/base::norm(D_2,"F")

err2_X_1=base::norm(X_1-X_1.jive,"2")/base::norm(X_1,"2")
err2_X_2=base::norm(X_2-X_2.jive,"2")/base::norm(X_2,"2")

err2_C_1=base::norm(C_1-C_1.jive,"2")/base::norm(C_1,"2")
err2_C_2=base::norm(C_2-C_2.jive,"2")/base::norm(C_2,"2")

err2_D_1=base::norm(D_1-D_1.jive,"2")/base::norm(D_1,"2")
err2_D_2=base::norm(D_2-D_2.jive,"2")/base::norm(D_2,"2")

output=c(errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2)

write.table(output,file=paste('siml1_seed',seed,'_p900_angle45_nstd1_JIVE.txt',sep=''),row.names = FALSE,col.names = FALSE)

#test correlations between D_1.jive and D_2.jive
p1=nrow(D_1.jive)
p2=nrow(D_2.jive)
nlogP=numeric(p1*p2)
t=1
for(i in 1:p1){
  for(j in 1:p2){
    nlogP[t]=-log(test_zero_corr(D_1.jive[i,],D_2.jive[j,]))
    t=t+1
  }
}

write.table(nlogP,file=paste('siml1_seed',seed,'_p900_angle45_nstd1_JIVE_nlogP.txt',sep=''),row.names = FALSE,col.names = FALSE)


###JIVE without the orthongonal constraint on distinct matrices (rJIVE)
set.seed(0)
Results = jive(data,center=F,orthIndiv =T)

Results$converged

C_1.rjive=Results$joint[[1]]*Results$scale$`Scale Values`[1]
C_2.rjive=Results$joint[[2]]*Results$scale$`Scale Values`[2]

D_1.rjive=Results$individual[[1]]*Results$scale$`Scale Values`[1]
D_2.rjive=Results$individual[[2]]*Results$scale$`Scale Values`[2]


X_1.rjive=C_1.rjive+D_1.rjive
X_2.rjive=C_2.rjive+D_2.rjive

errF_X_1=base::norm(X_1-X_1.rjive,"F")/base::norm(X_1,"F")
errF_X_2=base::norm(X_2-X_2.rjive,"F")/base::norm(X_2,"F")

errF_C_1=base::norm(C_1-C_1.rjive,"F")/base::norm(C_1,"F")
errF_C_2=base::norm(C_2-C_2.rjive,"F")/base::norm(C_2,"F")

errF_D_1=base::norm(D_1-D_1.rjive,"F")/base::norm(D_1,"F")
errF_D_2=base::norm(D_2-D_2.rjive,"F")/base::norm(D_2,"F")

err2_X_1=base::norm(X_1-X_1.rjive,"2")/base::norm(X_1,"2")
err2_X_2=base::norm(X_2-X_2.rjive,"2")/base::norm(X_2,"2")

err2_C_1=base::norm(C_1-C_1.rjive,"2")/base::norm(C_1,"2")
err2_C_2=base::norm(C_2-C_2.rjive,"2")/base::norm(C_2,"2")

err2_D_1=base::norm(D_1-D_1.rjive,"2")/base::norm(D_1,"2")
err2_D_2=base::norm(D_2-D_2.rjive,"2")/base::norm(D_2,"2")

output=c(errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2, max(abs(D_1.rjive%*%t(D_2.rjive))) )

write.table(output,file=paste('siml1_seed',seed,'_p900_angle45_nstd1_RJIVE.txt',sep=''),row.names = FALSE,col.names = FALSE)

#test correlations between D_1.rjive and D_2.rjive
p1=nrow(D_1.rjive)
p2=nrow(D_2.rjive)
nlogP=numeric(p1*p2)
t=1
for(i in 1:p1){
  for(j in 1:p2){
    nlogP[t]=-log(test_zero_corr(D_1.rjive[i,],D_2.rjive[j,]))
    t=t+1
  }
}

write.table(nlogP,file=paste('siml1_seed',seed,'_p900_angle45_nstd1_RJIVE_nlogP.txt',sep=''),row.names = FALSE,col.names = FALSE)





