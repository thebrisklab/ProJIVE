# R version 3.3.3
setwd('E:/Hai Windows/work/Hai_paper1/realdata')
#setwd('/Volumes/My Passport/Hai Windows/work/Hai_paper1/realdata')

num.IQR.outlier<-function(data){
  data.IQR=IQR(data)
  data.Q1=quantile(data,0.25)
  data.Q3=quantile(data,0.75)
  data.min=data.Q1-1.5*data.IQR
  data.max=data.Q3+1.5*data.IQR
  return(sum(data>data.max | data<data.min))
}

SWISS_1d<-function(x,label){
  category=unique(label)
  n.class=length(category)
  WISS=0
  for(i in 1:n.class)
  {
    x.class.i=x[ which(label==category[i]) ]
    WISS=WISS+sum( (x.class.i-mean(x.class.i))^2 )
  }
  SST=sum((x-mean(x))^2)
  return(WISS/SST)
}

#### find differentially variables between normals and tumors
#mRNA Expression is obtained from https://tcga-data.nci.nih.gov/docs/publications/brca_2015/
Exp=read.delim('BRCA817_20140528_log_medcntr.txt',header=T)
rownames(Exp)=Exp[,1]
Exp=Exp[-c(1,2),-1]
#DNA methylation is obtained from https://tcga-data.nci.nih.gov/docs/publications/brca_2012/
Meth=read.table('BRCA.methylation.27k.450k.txt',header=T)

Exp.gene=rownames(Exp)
Meth.probe=rownames(Meth)

################ Match columns (samples) between sources##############################
Exp.id = colnames(Exp)
Meth.id = colnames(Meth)

Exp.id = substr(Exp.id,1,16)

MatchMeth = match(Exp.id,Meth.id,nomatch=0)

#same subject id in the two datasets
index.Exp=which(MatchMeth!=0)
index.Meth=setdiff(MatchMeth,0)

length(Meth.id[index.Meth])
length(unique(Meth.id[index.Meth]))

sum(Meth.id[index.Meth] == Exp.id[index.Exp]) ##all match

Exp.id=Exp.id[index.Exp]

Exp.mat =as.matrix(Exp[,index.Exp])
storage.mode(Exp.mat)<-"numeric"
Meth.mat =as.matrix(Meth[,index.Meth])

################### Subtype labels#############################
library('xlsx')
#obtained from https://tcga-data.nci.nih.gov/docs/publications/brca_2015/
subtype= read.xlsx("BRCA_freeze_3.26.2014_ver06102014.xlsx", sheetIndex=1,header = TRUE) 
subtype=subtype[,c(1,4)]
dim(subtype)# 

length(subtype[,1])#817
length(unique(subtype[,1]))#817


subtype.id=as.character(subtype[,1])
subtype.id=gsub("[[:punct:]]", ".", subtype.id)

Exp.id=substr(Exp.id,1,12)

Matchsubtype = match(Exp.id,subtype.id,nomatch=0)

index.Exp=which(Matchsubtype!=0)
index.subtype=setdiff(Matchsubtype,0)

length(index.Exp)
length(unique(index.Exp))
sum(Exp.id[index.Exp]==subtype.id[index.subtype])#675

Exp.mat = Exp.mat[,index.Exp]
Meth.mat = Meth.mat[,index.Exp]
subtype=subtype[index.subtype,]

dim(subtype)

label=numeric(dim(subtype)[1])
label[which(subtype[,2]=='Basal')]=1
label[which(subtype[,2]=='LumA')]=2
label[which(subtype[,2]=='LumB')]=3
label[which(subtype[,2]=='Her2')]=4
label[which(subtype[,2]=='Normal')]=5

sum(label==1)#112
sum(label==2)#331
sum(label==3)#162
sum(label==4)#55
sum(label==5)#15

index=c(which(label==1),which(label==2),which(label==3),which(label==4))

Exp.mat = Exp.mat[,index]
Meth.mat =Meth.mat[,index]
subtype=subtype[index,]
label=label[index]

length(label)#660

################### Data preprocessing#############################
#check missing values
#remove genes with over 20% missing values
Exp.na=apply(Exp.mat,1,function(x) sum(is.na(x))/dim(Exp.mat)[2])
sum(Exp.na>0.2)#3944

Exp.mat=Exp.mat[-which(Exp.na>0.2),]
Exp.gene=Exp.gene[-which(Exp.na>0.2)]

Exp.NA=list()
Meth.NA=list()

Exp.i=list()
Meth.i=list()

for(i in 1:4)
{
  num.i=sum(label==i)
  Exp.NA[[i]]=numeric(num.i)
  Meth.NA[[i]]=numeric(num.i)
  
  Exp.i[[i]]=Exp.mat[,label==i]
  Meth.i[[i]]=Meth.mat[,label==i]
  
  for(j in 1:num.i)
  {
    Exp.NA[[i]][j]=sum(is.na(Exp.i[[i]][j,]))
    Meth.NA[[i]][j]=sum(is.na( Meth.i[[i]][j,] ))
  }
}

max.Exp.NA.i=numeric(4)
max.Meth.NA.i=numeric(4)
for(i in 1:4)
{
  max.Exp.NA.i[i]=max(Exp.NA[[i]])
  max.Meth.NA.i[i]=max(Meth.NA[[i]])
}
max.Exp.NA.i #39 72 32 13
max.Meth.NA.i # 2 14  5  1


#impute missing values
#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute)#version 1.46.0
Exp.imputed=impute.knn(Exp.mat,rng.seed=0)$data 
Meth.imputed=impute.knn(Meth.mat,rng.seed=0)$data 


Exp.num.out=apply(Exp.imputed,1,num.IQR.outlier)
Meth.num.out=apply(Meth.imputed,1,num.IQR.outlier)

quantile(Exp.num.out,seq(0.1,1,by=0.01))
quantile(Meth.num.out,seq(0.1,1,by=0.01))

#Exp.prep = Exp.imputed[Exp.num.out<25,] 
#Meth.prep = Meth.imputed[Meth.num.out<25,]

Exp.SWISS=apply(Exp.imputed,1,function(x) SWISS_1d(x,label))
quantile(Exp.SWISS,seq(0,1,by=0.01))

Meth.SWISS=apply(Meth.imputed,1,function(x) SWISS_1d(x,label))
quantile(Meth.SWISS,seq(0,1,by=0.01))

Exp.sd=apply(Exp.imputed,1,sd)
Meth.sd=apply(Meth.imputed,1,sd)
quantile(Exp.sd,seq(0,1,by=0.05))
quantile(Meth.sd,seq(0,1,by=0.05))


Exp.mad=apply(Exp.imputed,1,mad)
Meth.mad=apply(Meth.imputed,1,mad)

quantile(Exp.mad,seq(0,1,by=0.05))
quantile(Meth.mad,seq(0,1,by=0.05))


Exp.median=apply(Exp.imputed,1,median)
quantile(Exp.median,seq(0,1,by=0.05))

Meth.median=apply(Meth.imputed,1,median)
quantile(Meth.median,seq(0,1,by=0.05))


library(fBasics)
Exp.skew=apply(Exp.imputed,1,skewness)
Exp.kurtosis=apply(Exp.imputed,1,kurtosis)


Meth.skew=apply(Meth.imputed,1,skewness)
Meth.kurtosis=apply(Meth.imputed,1,kurtosis)


  
#Exp.prep= Exp.imputed[Exp.sd>quantile(Exp.sd,0.75) &Exp.SWISS<0.8 & Exp.num.out<=33  & abs(Exp.skew)<2 & abs(Exp.kurtosis)<5, ]
#Meth.prep = Meth.imputed[Meth.sd>quantile(Meth.sd,0.9) & Meth.SWISS>0.9 & Meth.num.out<=33 & abs(Meth.skew)<2 & abs(Meth.kurtosis)<5,]

  
#Exp.prep= Exp.imputed[Exp.sd>=quantile(Exp.sd,0.75) & Exp.SWISS<=0.8 & Exp.num.out<=33  & abs(Exp.skew)<=2 & abs(Exp.kurtosis)<=5, ]

Exp.prep= Exp.imputed[Exp.sd>=quantile(Exp.sd,0.85) & Exp.SWISS<=0.9 & Exp.num.out<=33  & abs(Exp.skew)<=2 & abs(Exp.kurtosis)<=5, ]

Meth.prep.above90 = Meth.imputed[Meth.sd>=quantile(Meth.sd,0.9) & Meth.SWISS>0.9 & Meth.num.out<=33 & abs(Meth.skew)<=2 & abs(Meth.kurtosis)<=5,]
Meth.prep.below90 = Meth.imputed[Meth.sd>=quantile(Meth.sd,0.9) & Meth.SWISS<=0.9 & Meth.num.out<=33 & abs(Meth.skew)<=2 & abs(Meth.kurtosis)<=5,]

Exp.gene=Exp.gene[Exp.sd>=quantile(Exp.sd,0.85) & Exp.SWISS<=0.9 & Exp.num.out<=33  & abs(Exp.skew)<=2 & abs(Exp.kurtosis)<=5]
Meth.probe.above90=Meth.probe[Meth.sd>=quantile(Meth.sd,0.9) & Meth.SWISS>0.9 & Meth.num.out<=33 & abs(Meth.skew)<=2 & abs(Meth.kurtosis)<=5]
Meth.probe.below90=Meth.probe[Meth.sd>=quantile(Meth.sd,0.9) & Meth.SWISS<=0.9 & Meth.num.out<=33 & abs(Meth.skew)<=2 & abs(Meth.kurtosis)<=5]

dim(Exp.prep) #1195  660
dim(Meth.prep.above90)#1202  660
dim(Meth.prep.below90)#881 660

length(Exp.gene)
length(Meth.probe.above90)
length(Meth.probe.below90)

#centering
library(matrixStats)
Expression=t(scale(t(Exp.prep), center = T, scale =apply(Exp.prep, 1, sd, na.rm = TRUE)))
max(abs(rowMeans(Expression)))
rowSds(Expression)

Exp.svd=svd(Expression%*%t(Expression)/dim(Expression)[2])
#plot(Exp.svd$d)

Methylation.above90=t(scale(t(Meth.prep.above90), center = T, scale =apply(Meth.prep.above90, 1, sd, na.rm = TRUE)))
max(abs(rowMeans(Methylation.above90)))
rowSds(Methylation.above90)

Meth.svd.above90=svd(Methylation.above90%*%t(Methylation.above90)/dim(Methylation.above90)[2])
#plot(Meth.svd.above90$d)

Methylation.below90=t(scale(t(Meth.prep.below90), center = T, scale =apply(Meth.prep.below90, 1, sd, na.rm = TRUE)))
max(abs(rowMeans(Methylation.below90)))
rowSds(Methylation.below90)

Meth.svd.below90=svd(Methylation.below90%*%t(Methylation.below90)/dim(Methylation.below90)[2])
#plot(Meth.svd.below90$d)

par(mfrow=c(1,3),mar=c(5, 5, 3, 1))#10*3.5 inch^2
plot(Exp.svd$d[1:100],ylab='Eigenvalue',main='EXP90')
plot(Meth.svd.below90$d[1:100],ylab='Eigenvalue',main='METH90b')
plot(Meth.svd.above90$d[1:100],ylab='Eigenvalue',main='METH90a')


write.table(Expression,file='Expression.txt',row.names = FALSE,col.names = FALSE)
write.table(Methylation.above90,file='Methylation_above90.txt',row.names = FALSE,col.names = FALSE)
write.table(Methylation.below90,file='Methylation_below90.txt',row.names = FALSE,col.names = FALSE)

write.table(label,file='PAM50_label.txt',row.names = FALSE,col.names = FALSE)

write.table(Exp.gene,file='Exp_gene.txt',row.names = FALSE,col.names = FALSE)
write.table(Meth.probe.above90,file='Meth_probe_above90.txt',row.names = FALSE,col.names = FALSE)
write.table(Meth.probe.below90,file='Meth_probe_below90.txt',row.names = FALSE,col.names = FALSE)

