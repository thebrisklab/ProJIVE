library(r.jive) #v2.1
source('test_zero_corr.R')

Y_1=read.table('Y_1.txt')
Y_2=read.table('Y_2.txt')
  
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

write.table(C_1.jive,file='C_1_jive.txt',row.names = FALSE,col.names = FALSE)
write.table(C_2.jive,file='C_2_jive.txt',row.names = FALSE,col.names = FALSE)


write.table(D_1.jive,file='D_1_jive.txt',row.names = FALSE,col.names = FALSE)
write.table(D_2.jive,file='D_2_jive.txt',row.names = FALSE,col.names = FALSE)


write.table(X_1.jive,file='X_1_jive.txt',row.names = FALSE,col.names = FALSE)
write.table(X_2.jive,file='X_2_jive.txt',row.names = FALSE,col.names = FALSE)

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

write.table(nlogP,file=paste('siml3_JIVE_nlogP.txt',sep=''),row.names = FALSE,col.names = FALSE)


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

write.table(C_1.rjive,file='C_1_rjive.txt',row.names = FALSE,col.names = FALSE)
write.table(C_2.rjive,file='C_2_rjive.txt',row.names = FALSE,col.names = FALSE)


write.table(D_1.rjive,file='D_1_rjive.txt',row.names = FALSE,col.names = FALSE)
write.table(D_2.rjive,file='D_2_rjive.txt',row.names = FALSE,col.names = FALSE)


write.table(X_1.rjive,file='X_1_rjive.txt',row.names = FALSE,col.names = FALSE)
write.table(X_2.rjive,file='X_2_rjive.txt',row.names = FALSE,col.names = FALSE)

max(abs(D_1.rjive%*%t(D_2.rjive)))

