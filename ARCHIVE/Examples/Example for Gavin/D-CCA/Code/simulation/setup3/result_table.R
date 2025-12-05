source('test_zero_corr.R')
D_1.DISCO=as.matrix(read.table('D_1_DISCO.txt'))
D_2.DISCO=as.matrix(read.table('D_2_DISCO.txt'))

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

write.table(nlogP,file=paste('siml3_DISCO_nlogP.txt',sep=''),row.names = FALSE,col.names = FALSE)


method='DISCO'
pvalue=exp(-scan(paste('siml3_',method,'_nlogP.txt',sep='')))
psc=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)# proportion of significant correlations between d1 and d2
psc


###gdfm
pvalue=exp(-scan('siml3_gdfm_nlogP_1.txt'))
psc=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)# proportion of significant correlations between d1 and d2
psc

pvalue=exp(-scan('siml3_gdfm_nlogP_2.txt'))
psc=sum(p.adjust(pvalue,method='BH')<0.05)/length(pvalue)# proportion of significant correlations between d1 and d2
psc

