########################################################
##Parc differential correlation analysis
##Lara Mangravite and Barbara Engelhardt, Oct 2012
########################################################
setwd("/gluster/home/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/DiffCorNetworks")
library(matrixStats)
library(lattice)
library(synapseClient)
synapseCacheDir('~/synapseCache/')

##Grab normalized data from Synapse
PARCnorm<-loadEntity("syn1417424")
names(PARCnorm$objects)
cov<-PARCnorm$objects$A932cov
dat<-PARCnorm$objects$A932norm
bdat<-dat[,cov$Treated=="Control"]
tdat<-dat[,cov$Treated=="Treated"]

################################################
#Start analysis with a subset of candidate genes
#################################################
#Chol biosynth genes, genes diff correl with HMGCR, genes not diff correl with HMGCR
#chol bio genes
genes<-c("HMGCR","LDLR","HMGCS1","MVK","MVD","LSS","DHCR7")

#load gene-gene correlation analysis output performed using cor.test(dat)$p.value
# in order to id genes for inifical analysis.

#load from server on belltown
load(file="/gluster/home/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/GeneGeneCorrelation/PvalResponseCorrelations.RData")
#load from server using local computer
#diffcor<-load(file="/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/GeneGeneCorrelation/PvalResponseCorrelations.RData")

g<-grep("GATM",names(combop)) 
hist(combop[[g]])
sum(combop[[g]]==0) #18
h<-grep("HMGCR",names(combop)) 
hist(combop[[h]])
sum(combop[[h]]==0) #137
l<-grep("LDLR",names(combop)) 
l<-l[1]
hist(combop[[l]],breaks=100)
sum(combop[[l]]==0) #713
genes<-unique(c(genes,
         names(combop[[h]][order(combop[[h]])][c(1:137,((10195-137):10195))]),
         names(combop[[l]][order(combop[[l]])][c(1:713,((10195-713):10195))]),
         names(combop[[g]][order(combop[[g]])][c(1:100,((10195-100):10195))])))

bdata<-bdat[which(rownames(bdat) %in% genes),]
tdata<-tdat[which(rownames(tdat) %in% genes),]

source("/gluster/home/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")

#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#seen in the control-treated samples but not in the simvastatin-treated samples
d1not2<-build_correlation_graph_differential_1not2(gex1=bdata,gex2=tdata,gexnames=genes)
save(d1not2,file="AbridgedControlspecificcorrelation.Rdat")
cor<-as.numeric(d1not2$corg)[as.numeric(d1not2$corg) != 1]
e<-as.numeric(d1not2$edges)[as.numeric(d1not2$corg) != 1]
summary(abs(cor[e == 1]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1041  0.1471  0.2188  0.2653  0.3551  0.9285  
summary(abs(cor[e == 0]))
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#5.000e-08 2.156e-02 4.456e-02 4.693e-02 7.096e-02 1.041e-01  


#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#seen in the control-treated samples but not in the simvastatin-treated samples
d2not1<-build_correlation_graph_differential_2not1(gex1=bdata,gex2=tdata,gexnames=genes)
save(d2not1,file="AbridgedTreatedspecificcorrelatin.Rdat")
cor<-as.numeric(d2not1$corg)[as.numeric(d2not1$corg) != 1]
e<-as.numeric(d2not1$edges)[as.numeric(d2not1$corg) != 1]
summary(abs(cor[e == 1]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1042  0.1463  0.2134  0.2550  0.3353  0.9249 
summary(abs(cor[e == 0]))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#1.000e-08 2.178e-02 4.486e-02 4.715e-02 7.124e-02 1.042e-01 
#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#seen in the control-treated samples but not in the simvastatin-treated samples
dop<-build_correlation_graph_differential_opposite(gex1=bdata,gex2=tdata,gexnames=genes)
save(dop,file="AbridgedDifferentialcorrelation.Rdat")
cor<-as.numeric(dop$corg)[as.numeric(dop$corg) != 1]
e<-as.numeric(dop$edges)[as.numeric(dop$corg) != 1]
summary(abs(cor[e == 1]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1141  0.1411  0.1861  0.2284  0.2788  0.9086 
summary(abs(cor[e == 0]))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#7.000e-08 2.045e-02 4.276e-02 4.686e-02 7.052e-02 1.141e-01 

cont<-build_correlation_graph(gex1=bdata,gexnames=genes)
save(cont,file="AbridgedControlCorrelation.Rdat")
treat<-build_correlation_graph(gex1=tdata,gexnames=genes)
save(treat,file="AbridgedTreatmentCorrelation.Rdat")

############################
#Run full dataset
##############################
genes<-rownames(bdat)
source("/gluster/home/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")

#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#seen in the control-treated samples but not in the simvastatin-treated samples
d1not2<-build_correlation_graph_differential_1not2(gex1=bdat,gex2=tdat,gexnames=genes)
save(d1not2,file="controlspecificcorrelationFull.Rdat")
cor<-as.numeric(d1not2$corg)[as.numeric(d1not2$corg) != 1]
e<-as.numeric(d1not2$edges)[as.numeric(d1not2$corg) != 1]
summary(abs(cor[e == 1]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1056  0.1378  0.1819  0.2105  0.2546  0.9731 
summary(abs(cor[e == 0]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02263 0.04651 0.04852 0.07330 0.10560 
sum(e == 1)/(10195*10195)
#0.4378163 percent have significant edges with this model.


#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#seen in the control-treated samples but not in the simvastatin-treated samples
d2not1<-build_correlation_graph_differential_2not1(gex1=bdat,gex2=tdat,gexnames=genes)
save(d2not1,file="treatedspecificcorrelationFull.Rdat")
cor<-as.numeric(d2not1$corg)[as.numeric(d2not1$corg) != 1]
e<-as.numeric(d2not1$edges)[as.numeric(d2not1$corg) != 1]
summary(abs(cor[e == 1]))
summary(abs(cor[e == 0]))
sum(e == 1)/(10195*10195)
rm(d2not1,cor,e)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1059  0.1373  0.1802  0.2074  0.2503  0.9714 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02264 0.04656 0.04860 0.07342 0.10590 
#[1] 0.4301667

#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#seen in the control-treated samples but not in the simvastatin-treated samples
dop<-build_correlation_graph_differential_opposite(gex1=bdat,gex2=tdat,gexnames=genes)
save(dop,file="differentialcorrelationFull.Rdat")
cor<-as.numeric(dop$corg)[as.numeric(dop$corg) != 1]
e<-as.numeric(dop$edges)[as.numeric(dop$corg) != 1]
summary(abs(cor[e == 1]))
summary(abs(cor[e == 0]))
sum(e == 1)/(10195*10195)
rm(dop,cor,e)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1174  0.1383  0.1694  0.1959  0.2264  0.9147 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02086 0.04370 0.04795 0.07209 0.11740 
#[1] 0.2143046

#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#in control samples
control<-build_correlation_graph(gex1=bdat,gexnames=genes)
save(control,file="ControlCorrelationFull.Rdat")
cor<-as.numeric(control$corg)[as.numeric(control$corg) != 1]
e<-as.numeric(control$edges)[as.numeric(control$corg) != 1]
summary(abs(cor[e == 1]))
summary(abs(cor[e == 0]))
sum(e == 1)/(10195*10195)
rm(control,cor,e)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1057  0.1371  0.1799  0.2073  0.2498  0.9789 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02278 0.04677 0.04870 0.07355 0.10570 
#[1] 0.4349275

#contcyto<-build_cytoscape_input()


#Correlation matrix and edges (FDR<0.05) for gene-gene correlation structure
#in treated samples
treat<-build_correlation_graph(gex1=tdat,gexnames=genes)
save(treat,file="ControlCorrelationFull.Rdat")
cor<-as.numeric(treat$corg)[as.numeric(treat$corg) != 1]
e<-as.numeric(treat$edges)[as.numeric(treat$corg) != 1]
summary(abs(cor[e == 1]))
summary(abs(cor[e == 0]))
sum(e == 1)/(10195*10195)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1060  0.1368  0.1786  0.2050  0.2465  0.9800 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02277 0.04678 0.04875 0.07363 0.10600 
#[1] 0.4281149

i<-grep("HMGCR",rownames(bdat))
hist(treat$corg[,i],breaks=40,col="grey")
hist(cor,breaks=40,col="grey")

##how many significant edges per gene.
avee<-rowSums(treat$edges)
png("HistTotalGGCsPerGeneTreated.png")
hist(avee,breaks=40,col="grey",main="Treated",xlim=c(0,10000),xlab="N, significant correlations per gene")
dev.off()
summary(avee)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#155    3474    4553    4365    5464    7112 
#st dev: 1448.295

#Candidate genes, N corr> 0.2
test<-c("HMGCR","LDLR","GATM","RHOA","SQLE","MVK","MVD")
avee[rownames(tdat) %in% test]
#GATM HMGCR  LDLR   MVD   MVK  RHOA  SQLE 
#3763 3350 5588 5693 4537 6190 5263

##how many sigificant edges with cor > 0.2 per gene?
temp<-abs(treat$corg)>=0.2
ave<-rowSums(temp)
hist(ave,col="grey")
png("HistCorrGreater2Treated.png")
hist(ave,breaks=40,col="grey",main="Treated",xlim=c(0,10000),xlab="N, correlations > 0.2 per gene")
dev.off()
summary(ave)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0   793.5  1639.0  1763.0  2635.0  4856.0 
#st dev: 1175.03
#Candidate genes, N corr> 0.2
test<-c("HMGCR","LDLR","GATM","RHOA","SQLE","MVK","MVD")
ave[rownames(tdat) %in% test]
#GATM HMGCR  LDLR   MVD   MVK  RHOA  SQLE 
#955   810  2709  2866  1682  3628  2389 

grep("LDLR",rownames(tdat))
#4046, 4776
temp1<-treat$corg[rownames(tdat) == "HMGCR",]
png("HMGCRcorVsLDLRcor.png")
plot(treat$corg[rownames(tdat) == "HMGCR",c(1:4045,4047:4775,4777:10195)]
     ~treat$corg[rownames(tdat) == "LDLR",c(1:4045,4047:4775,4777:10195)],xlab="LDLR ggc",ylab="HMGCR ggc")
abline(lm(treat$corg[rownames(tdat) == "HMGCR",c(1:4045,4047:4775,4777:10195)]~treat$corg[rownames(tdat) == "LDLR",c(1:4045,4047:4775,4777:10195)]),col="red")
dev.off()
cor.test(treat$corg[rownames(tdat) == "HMGCR",c(1:4045,4047:4775,4777:10195)],treat$corg[rownames(tdat) == "LDLR",c(1:4045,4047:4775,4777:10195)])
#r2=0.58, p=2.2x10-16

i<-c(grep("HMGCR",rownames(tdat)),grep("GATM",rownames(tdat)))
png("HMGCRcorVsGAATMcor.png")
plot(treat$corg[rownames(tdat) == "HMGCR",-c(i)]
     ~treat$corg[rownames(tdat) == "GATM",-c(i)],xlab="GATM ggc", ylab="HMGCR ggc")
abline(lm(treat$corg[rownames(tdat) == "HMGCR",-c(i)]~treat$corg[rownames(tdat) == "GATM",-c(i)]),col="red")
dev.off()
cor.test(treat$corg[rownames(tdat) == "HMGCR",-c(i)],treat$corg[rownames(tdat) == "GATM",-c(i)])
#r2=0.22, p=2.2x10-16

#####Formally, assess the correlation of HMGCR ggc to all other genes on the array.
corofcor<-c()
pofcor<-c()
  for (i in 1:10195) {
  j<-c(4046,i)
  temp<-cor.test(treat$corg[4046,-c(j)],treat$corg[i,-c(j)])
  corofcor<-c(corofcor,temp$estimate)
  pofcor<-c(pofcor,temp$p.value)
}
HMGCRcor<-cbind(corofcor,pofcor)
rownames(HMGCRcor)<-rownames(tdat)
colnames(HMGCRcor)<-c("cor","pval")
library(qvalue)
HMGCRcor<-cbind(HMGCRcor,qvalue(HMGCRcor[,2])$q)
HMGCRcor<-cbind(HMGCRcor,(HMGCRcor[,1]*HMGCRcor[,1]))
colnames(HMGCRcor)[3:4]<-c("qval","r2")
HMGCRcor<-HMGCRcor[order(HMGCRcor[,1],decreasing=T),]
HMGCRcor<-cbind(HMGCRcor,c(1:10195))
HMGCRcor[which(rownames(HMGCRcor) %in% test),]
#cor         pval         qval           r2     
#HMGCR  1.00000000 0.000000e+00 0.000000e+00 1.0000000000    1
#LDLR   0.57705510 0.000000e+00 0.000000e+00 0.3329925894   42
#MVK    0.40314237 0.000000e+00 0.000000e+00 0.1625237731  550
#MVD    0.28486747 0.000000e+00 0.000000e+00 0.0811494762 1643
#GATM   0.22459503 0.000000e+00 0.000000e+00 0.0504429297 2365
#SQLE   0.02164992 2.883229e-02 8.888218e-04 0.0004687189 5012
#RHOA  -0.16847861 8.840625e-66 3.727478e-67 0.0283850433 7376

#####Formally, assess the correlation of LDLR ggc to all other genes on the array.
corofcor<-c()
pofcor<-c()
for (i in 1:10195) {
  j<-c(4776,i)
  temp<-cor.test(treat$corg[4776,-c(j)],treat$corg[i,-c(j)])
  corofcor<-c(corofcor,temp$estimate)
  pofcor<-c(pofcor,temp$p.value)
}
LDLRcor<-cbind(corofcor,pofcor)
rownames(LDLRcor)<-rownames(tdat)
colnames(LDLRcor)<-c("cor","pval")
library(qvalue)
LDLRcor<-cbind(LDLRcor,qvalue(LDLRcor[,2])$q)
LDLRcor<-cbind(LDLRcor,(LDLRcor[,1]*LDLRcor[,1]))
colnames(LDLRcor)[3:4]<-c("qval","r2")
LDLRcor<-LDLRcor[order(LDLRcor[,1],decreasing=T),]
LDLRcor<-cbind(LDLRcor,c(1:10195))
LDLRcor[which(rownames(LDLRcor) %in% test),]
#cor          pval          qval         r2     
#LDLR   1.0000000  0.000000e+00  0.000000e+00 1.00000000    1
#HMGCR  0.5770551  0.000000e+00  0.000000e+00 0.33299259 1486
#GATM   0.5515804  0.000000e+00  0.000000e+00 0.30424088 1676
#MVK   -0.2800522 5.081697e-183 6.790689e-185 0.07842924 7016
#MVD   -0.3122058 2.829708e-229 3.867101e-231 0.09747247 7193
#RHOA  -0.3565077 3.433827e-303 4.870543e-305 0.12709774 7478
#SQLE  -0.3608318  0.000000e+00  0.000000e+00 0.13019955 7504

#####Formally, assess the correlation of GATM ggc to all other genes on the array.
corofcor<-c()
pofcor<-c()
for (i in 1:10195) {
  j<-c(4776,i)
  temp<-cor.test(treat$corg[4776,-c(j)],treat$corg[i,-c(j)])
  corofcor<-c(corofcor,temp$estimate)
  pofcor<-c(pofcor,temp$p.value)
}
GATMcor<-cbind(corofcor,pofcor)
rownames(GATMcor)<-rownames(tdat)
colnames(GATMcor)<-c("cor","pval")
library(qvalue)
GATMcor<-cbind(GATMcor,qvalue(GATMcor[,2])$q)
GATMcor<-cbind(GATMcor,(GATMcor[,1]*GATMcor[,1]))
colnames(GATMcor)[3:4]<-c("qval","r2")
GATMcor<-GATMcor[order(GATMcor[,1],decreasing=T),]
GATMcor<-cbind(GATMcor,c(1:10195))
GATMcor[which(rownames(GATMcor) %in% test),]
#cor          pval          qval         r2     
#LDLR   1.0000000  0.000000e+00  0.000000e+00 1.00000000    1
#HMGCR  0.5770551  0.000000e+00  0.000000e+00 0.33299259 1486
#GATM   0.5515804  0.000000e+00  0.000000e+00 0.30424088 1676
#MVK   -0.2800522 5.081697e-183 6.790689e-185 0.07842924 7016
#MVD   -0.3122058 2.829708e-229 3.867101e-231 0.09747247 7193
#RHOA  -0.3565077 3.433827e-303 4.870543e-305 0.12709774 7478
#SQLE  -0.3608318  0.000000e+00  0.000000e+00 0.13019955 7504

#############Formally do this analyssi for all genes.  But how to assess significance?

##################################################################
##Cluster genes based on their gene-gene correlation structures.##
##################################################################

####Pull this information out so that we can use it in cytoscape.
####Will be much faster to do so with a more stringent cutt off.
temp<-abs(treat$corg)>=0.2 

  modcor<-c()
  for (i in 11:2000){
    temp<-c()
    for (j in ((i+1):10195)) {
      if (edges[i,j] == 1 & abs(cor[i,j]) >= 0.2) {
        temp<-data.frame(Node1=genenames[i],Node2=genenames[j],CorCoeff=as.numeric(cor[i,j]),AbsCorCoeff=as.numeric(abs(cor[i,j])))
        modcor<-rbind(modcor,temp)
      } else  {
        next
      }
    }
  }
  #  colnames(modcor)<-c("Node1","Node2","CorCoeff","AbsCorCoeff")
  modcor<-cbind(modcor,modcor[,4]>=0.5,modcor[,4]>=0.4,modcor[,3]>=0.3)
  colnames(modcor)[5:7]<-c("CorCoeffGreater5","CorCoeffGreater4","CorCoeffGreater3")
  return(modcor)
}

modcor10<-modcor

modcor<-c()
for (i in 11:2000){
  temp<-c()
  for (j in ((i+1):10195)) {
    if (edges[i,j] == 1 & abs(cor[i,j]) >= 0.2) {
      temp<-data.frame(Node1=genenames[i],Node2=genenames[j],CorCoeff=as.numeric(cor[i,j]),AbsCorCoeff=as.numeric(abs(cor[i,j])))
      modcor<-rbind(modcor,temp)
    } else  {
      next
    }
  }
}
modcorend<-rbind(modcor10,modcor)
write.table(modcorend,file="treatForCyto.txt",header=T)

modcor<-c()
for (i in 2001:5000){
  temp<-c()
  for (j in ((i+1):10195)) {
    if (edges[i,j] == 1 & abs(cor[i,j]) >= 0.2) {
      temp<-data.frame(Node1=genenames[i],Node2=genenames[j],CorCoeff=as.numeric(cor[i,j]),AbsCorCoeff=as.numeric(abs(cor[i,j])))
      modcor<-rbind(modcor,temp)
    } else  {
      next
    }
  }
}
modcorend<-rbind(modcorend,modcor)
write.table(modcorend,file="treatForCyto.txt",header=T)

modcor<-c()
for (i in 5001:8000){
  temp<-c()
  for (j in ((i+1):10195)) {
    if (edges[i,j] == 1 & abs(cor[i,j]) >= 0.2) {
      temp<-data.frame(Node1=genenames[i],Node2=genenames[j],CorCoeff=as.numeric(cor[i,j]),AbsCorCoeff=as.numeric(abs(cor[i,j])))
      modcor<-rbind(modcor,temp)
    } else  {
      next
    }
  }
}
modcorend<-rbind(modcorend,modcor)
write.table(modcorend,file="treatForCyto.txt",header=T)

modcor<-c()
for (i in 8001:10195){
  temp<-c()
  for (j in ((i+1):10195)) {
    if (edges[i,j] == 1 & abs(cor[i,j]) >= 0.2) {
      temp<-data.frame(Node1=genenames[i],Node2=genenames[j],CorCoeff=as.numeric(cor[i,j]),AbsCorCoeff=as.numeric(abs(cor[i,j])))
      modcor<-rbind(modcor,temp)
    } else  {
      next
    }
  }
}
modcorend<-rbind(modcorend,modcor)
write.table(modcorend,file="treatForCyto.txt",header=T)

#treatCyto<-build_cytoscape_input(cor=treat$corg,edges=treat$edges,genenames<-colnames(treat$corg))
#write.table(treatCyto,file="TreatCorrCytoInput.txt",header=T,row.names=F)
#rm(treat,cor,e)















