library(qvalue)#‘2.10.1’
source('test_zero_corr.R')
method='OnPLS'
D_1=as.matrix(read.table(paste('D_1_hat_',method,'_above90.txt',sep='')))
D_2=as.matrix(read.table(paste('D_2_hat_',method,'_above90.txt',sep='')))

p1=nrow(D_1)
p2=nrow(D_2)
nlogP=numeric(p1*p2)
t=1
for(i in 1:p1){
  for(j in 1:p2){
    nlogP[t]=-log(test_zero_corr(D_1[i,],D_2[j,]))
    t=t+1
  }
}

write.table(nlogP,file=paste(method,'_nlogP_above90.txt',sep=''),row.names = FALSE,col.names = FALSE)


method='OnPLS'
pvalue=exp(-scan(paste(method,'_nlogP_above90.txt',sep='')))
psc=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)# proportion of significant correlations between d1 and d2
psc


#Block-GDFM1
D_1=t(as.matrix(read.table('psi_y1_above90.txt')))
D_2=t(as.matrix(read.table('psi_y2_above90.txt')))

p1=nrow(D_1)
p2=nrow(D_2)
nlogP=numeric(p1*p2)
t=1
for(i in 1:p1){
  for(j in 1:p2){
    nlogP[t]=-log(test_zero_corr(D_1[i,],D_2[j,]))
    t=t+1
  }
}

write.table(nlogP,file='Block-GDFM1_nlogP_above90.txt',row.names = FALSE,col.names = FALSE)

pvalue=exp(-scan('Block-GDFM1_nlogP_above90.txt'))
psc=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)# proportion of significant correlations between d1 and d2
psc



#Block-GDFM2
D_1=t(as.matrix(read.table('psi_y1_above90.txt')  +    read.table('nu_y1_above90.txt')  ) )
D_2=t(as.matrix(read.table('psi_y2_above90.txt')  +    read.table('nu_y2_above90.txt')  ) )

p1=nrow(D_1)
p2=nrow(D_2)
nlogP=numeric(p1*p2)
t=1
for(i in 1:p1){
  for(j in 1:p2){
    nlogP[t]=-log(test_zero_corr(D_1[i,],D_2[j,]))
    t=t+1
  }
}

write.table(nlogP,file='Block-GDFM2_nlogP_above90.txt',row.names = FALSE,col.names = FALSE)

pvalue=exp(-scan('Block-GDFM2_nlogP_above90.txt'))
psc=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)# proportion of significant correlations between d1 and d2
psc


