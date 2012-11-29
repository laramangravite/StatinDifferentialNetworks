################################################################
##Renormalization CAP data for network analysis.
#Aug-Sept 2012
################################################################
####load packages and obtain raw PARC expression and covariate data from Synapse.
setwd("/gluster/home/lmangrav/PARC/CoexpNetworks/CoexpnJuly2012")
library(sva)
library(gcrma)
library(matrixStats)
library(WGCNA)
library(lattice)
library(synapseClient)
synapseCacheDir('~/synapseCache/')
synapseLogin(user="lara.mangravite@sagebase.org")
#Load data
PARCraw<-loadEntity("syn1417162")
names(PARCraw$objects)
Araw<-PARCraw$objects$A960raw
Acov<-PARCraw$objects$All960cov
Araw<-log2(Araw)
rownames(Araw)<-Acov$SampleID

################################################################
##Prepare files for normalization.
################################################################

################################################################
##REmove samples represented as the only "pass" sample on an array. Remove 10 samples.
################################################################
####Using the below defined variables as covariates, reduce samples so that model matrix is not singular and there are no unpaired samples.
covariates<-Acov
remove=model.matrix(~covariates$ExposureBatch+covariates$RNAbatch+covariates$Array+covariates$Cell.Count
                    +covariates$Gender+covariates$Age+covariates$BMI+covariates$Smoking.Status)
X<-remove
b1 <- solve(t(X) %*% X) %*% t(X) %*% (data)
##Model Matrix is Singular.
colnames(remove)[colSums(remove)==1]
#[1] "covariates$Array4482545049" "covariates$Array4256290076" "covariates$Array4256290098"
#[4] "covariates$Array4256425106" "covariates$Array4669609056"
Acov[which(Acov$Array %in% c("4256290098","4256290076","4482545049","4256425106","4669609056")),]
#CAP03.1361P mangravite_12.23.09.csv  4482545049 44
#CAP02.1786S mangravite_12.23.09.csv 4669609056  44
#CAP03.1135S            DoseResponse 4256290076  1
#CAP04.0727S mangravite_11.25.08.csv 4256425106  5
#CAP04.1958S            DoseResponse  4256290098 1

##Remove these from the data.
order<-data.frame(unique(Acov$Array),unique(as.numeric(Acov$Array)))
covariates$Array<-as.numeric(covariates$Array)
c<-c("4256290098","4256290076","4482545049","4256425106","4669609056")
covariates<-covariates[-(which(covariates$Array %in% order[which(order[,1] %in% c),2])),]
temp<-(as.character((Acov$SampleID[which(Acov$Array %in% c)])))
data<-Araw[-c(which(rownames(Araw) %in% (as.character(Acov$SampleID[which(Acov$Array %in% c)])))),]
covariates$Array<-as.factor(covariates$Array)


##Now remove the paired sample for each of the five that were just removed.
covP<-covariates[covariates$Treated==0,1:2]
covT<-covariates[covariates$Treated==1,1:2]
covP<-covP[order(covP[,2]),]
covP<-cbind(covP,rep(1,479))
colnames(covP)<-c("Psamples","PID","HasPsample")
covT<-covT[order(covT[,2]),]
covT<-cbind(covT,rep(1,476))
colnames(covT)<-c("Tsamples","PID","HasTsample")
covPT<-merge(covP,covT,by.x=2,by.y=2,all.x=T,all.y=T)
covPT<-cbind(covPT,(as.numeric(covPT[,3])+as.numeric(covPT[,5])))
#           PID    Psamples HasPsample    Tsamples HasTsample
#75  CAP02.1786 CAP02.1786P          1        <NA>         NA
#166 CAP03.1135 CAP03.1135P          1        <NA>         NA
#181 CAP03.1361        <NA>         NA CAP03.1361S          1
#400 CAP04.0727 CAP04.0727P          1        <NA>         NA
#472 CAP04.1958 CAP04.1958P          1        <NA>         NA

c<-c("CAP02.1786P","CAP03.1135P","CAP03.1361S","CAP04.0727P","CAP04.1958P")
A950norm<-A955Norm[-(which(rownames(A955Norm) %in% c)),]
A950cov<-covariates[-(which(covariates$SampleID %in% c)),]
A950raw<-raw[-(which(rownames(raw) %in% c)),]
##Does this make the model matrix singular?
remove=model.matrix(~A950cov$ExposureBatch+A950cov$RNAbatch+A950cov$Array+A950cov$Cell.Count
                    +A950cov$Gender+A950cov$Age+A950cov$BMI+A950cov$Smoking.Status)
X<-remove
b1 <- solve(t(X) %*% X) %*% t(X) %*% (A950raw)
##Not singular!

################################################################
##Assess if we should remove the few samples left fromt the 0point RNA plate/batch . Remove 0 samples.
################################################################
###Now remove the few samples left fromt the 0point RNA plate/batch to see if this is the source of the RNA batch variance
###############This  not the source so leave these samples in.########################
#A950cov[A950cov$RNA.Plate==0,]$SampleID
#CAP04.1433P CAP04.1433S
#A950cov[A950cov$RNAbatch==1,]$SampleID
#same samples.
#str(A950cov)
#order<-data.frame(unique(A950cov$RNA.Plate),unique(as.numeric(A950cov$RNA.Plate)))
#A950cov$RNA.Plate<-as.numeric(A950cov$RNA.Plate)
#temp<-A950cov[-c(which(A950cov$RNA.Plate==order[order[,1]==0,2])),]
#temp$RNA.Plate<-as.factor(temp$RNA.Plate)
#temp$Array<-as.numeric(temp$Array)
#temp$Array<-as.factor(temp$Array)
#temp$RNAbatch<-as.factor(temp$RNAbatch)
#tempraw<-Araw[c(which(rownames(Araw) %in% as.character(temp$SampleID))),]
#source("/gluster/home/lmangrav/PARC/CoexpNetworks/CoexpnJuly2012/normParcData.R")
#A948Norm<-normParcData(tempraw,temp)
#A948cov<-temp

################################################################
##Remove samples that looked like the entire batch was reversed. Remove 8 samples.
################################################################

out<-as.character(A950cov[A950cov$Date.Exposed=="90122",]$SampleID)
#[1] CAP04.0441P CAP04.0451P CAP04.0508P CAP04.0672P CAP04.0441S CAP04.0451S CAP04.0508S
#[8] CAP04.0672S
unique(A950cov[A950cov$Date.Exposed=="90122",]$ExposureBatch)  #batch 27
order<-data.frame(unique(A950cov$ExposureBatch),unique(as.numeric(A950cov$ExposureBatch)))
A950cov$ExposureBatch<-as.numeric(A950cov$ExposureBatch)
A942cov<-A950cov[-c(which(A950cov$ExposureBatch==order[order[,1]==27,2])),]
A942cov$ExposureBatch<-as.factor(A942cov$ExposureBatch)
hist(as.numeric(A942cov$Array),breaks=length(unique(as.numeric(A942cov$Array))))
sum(as.numeric(A942cov$Array)==4)
A942cov[as.numeric(A942cov$Array)==4,1:4]
#temp$Array<-as.numeric(temp$Array)
#temp$Array<-as.factor(temp$Array)
A942raw<-Araw[c(which(rownames(Araw) %in% as.character(A942cov$SampleID))),]

################################################################
##Reassess gender mismatch using a different gene than the original. Remove 10 samples.
################################################################
gender<-c("ASMT","IL3RA","SLC25A6","CSF2RA","CD99","ASMTL","VAMP7","RBMY3AP","PCDH11Y","JARID1D")
rownames(A942norm)[which(rownames(A942norm) %in% gender)]
gender<-c("ASMTL","CD99","IL3RA","JARID1D","VAMP7")
i<-4
hist(A942norm[rownames(A942norm)==gender[i],])
hist(A942raw[,colnames(A942raw)==gender[i]])
##Of these, only gender[4] has a bimodal distribution and then only in the raw data.
xyplot(A942raw[,colnames(A942raw)==gender[i]]~A942cov$Gender)
summary(A942raw[A942cov$Gender==0,colnames(A942raw)==gender[i]])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.418  10.870  11.110  11.110  11.340  12.030 
tempm<-A942raw[A942cov$Gender==0,]
tempm<-tempm[,colnames(A942raw)==gender[i]]
(tempm[order(tempm,decreasing=F)])[1:10]
#CAP03.0925P CAP04.0509P CAP04.0509S CAP03.0190P CAP02.1051P CAP04.0483P CAP02.1051S CAP02.1743P 
#8.417733    9.216815    9.510668    9.564954    9.932274   10.112705   10.137939   10.338387 
#CAP04.1116P CAP04.1032P 
#10.368821   10.408825

summary(A942raw[A942cov$Gender==1,colnames(A942raw)==gender[i]])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#7.956   8.185   8.251   8.290   8.323  12.090 
tempf<-A942raw[A942cov$Gender==1,]
tempf<-tempf[,colnames(A942raw)==gender[i]]
(tempf[order(tempf,decreasing=T)])[1:10]
#CAP03.2486P CAP03.2174S CAP04.0523S CAP03.2174P CAP03.1901S CAP02.1783S CAP03.1113S CAP03.0229S 
#12.087603   11.771995   11.475950   11.295348    8.739233    8.727898    8.717620    8.653318 
#CAP03.0778S CAP03.2245S CAP03.2056S CAP02.1206S CAP04.0160P CAP03.1714S CAP03.1749S CAP03.0583S 
#8.632616    8.629725    8.614396    8.613495    8.605760    8.569264    8.560948    8.548729

A942raw[which(rownames(A942raw) %in% c("CAP03.0925P","CAP03.0925S")), 
        colnames(A942raw)==gender[i]]
#CAP03.0925P CAP03.0925S 
#8.417733   11.212314 
A942raw[which(rownames(A942raw) %in% c("CAP04.0509P","CAP04.0509S")), 
        colnames(A942raw)==gender[i]]
#CAP04.0509P CAP04.0509S 
#9.216815    9.510668 
A942raw[which(rownames(A942raw) %in% c("CAP03.0190P","CAP03.0190S")), 
        colnames(A942raw)==gender[i]]
#CAP03.0190P CAP03.0190S 
#9.564954   11.162569 
A942raw[which(rownames(A942raw) %in% c("CAP02.1051P","CAP02.1051S")), 
        colnames(A942raw)==gender[i]]
#CAP02.1051P CAP02.1051S 
#9.932274   10.137939 
A942raw[which(rownames(A942raw) %in% c("CAP03.2486P","CAP03.2486S")), 
        colnames(A942raw)==gender[i]]
#CAP03.2486P CAP03.2486S 
#12.087603    8.326346 
A942raw[which(rownames(A942raw) %in% c("CAP04.0523P","CAP04.0523S")), 
        colnames(A942raw)==gender[i]]
#CAP04.0523P CAP04.0523S 
#8.143851   11.475950 
A942raw[which(rownames(A942raw) %in% c("CAP03.2174P","CAP03.2174S")), 
        colnames(A942raw)==gender[i]]
#CAP03.2174P CAP03.2174S 
#11.29535    11.77199 

###According to this, there are a few gender mismatches.
##We should remove:
out<-c("CAP03.2174","CAP04.0523","CAP03.2486","CAP03.0190","CAP03.0925")

A942cov$PID<-as.character(A942cov$PID)
A942cov$SampleID<-as.character(A942cov$SampleID)
A932cov<-A942cov[-c(which(A942cov$PID %in% out)),]
A932cov$ExposureBatch<-as.numeric(A932cov$ExposureBatch)
A932cov$Array<-as.numeric(A932cov$Array)

A932cov$PID<-as.factor(A932cov$PID)
A932cov$SampleID<-as.factor(A932cov$SampleID)
A932cov$ExposureBatch<-as.factor(A932cov$ExposureBatch)
A932cov$Array<-as.factor(A932cov$Array)

hist(as.numeric(A932cov$Array),breaks=length(unique(as.numeric(A932cov$Array))))

A932raw<-Araw[c(which(rownames(Araw) %in% as.character(A932cov$SampleID))),]

################################################################
##Normalize 932 samples.
################################################################
source("/gluster/home/lmangrav/PARC/CoexpNetworks/CoexpnJuly2012/normParcData.R")
A932norm<-normParcData(A932raw,A932cov)
temp<-A932norm


################################################################
#Assessing quality of normalization.
################################################################

pcs2<-eigen(cov(A932norm))
plot(pcs2$vectors[,1]~A932cov$RNA.Plate)
plot(pcs2$vectors[,1]~A932cov$RNAbatch)
plot(pcs2$vectors[,1]~A932cov$Treated)
plot(pcs2$vectors[,1]~A932cov$Array)
plot(pcs2$vectors[,1]~A932cov$ExposureBatch)
xyplot(pcs2$vectors[,1]~A932cov$ExposureBatch, groups=A932cov$Treated)
xyplot(pcs2$vectors[,1]~A932cov$Array, groups=A932cov$Treated)
xyplot(pcs2$vectors[,1]~A932cov$RNA.Plate, groups=A932cov$Treated)
temp<-pcs2
(temp$values/sum(temp$values))[1:10]
# [1] 0.9827473121 0.0014544109 0.0011629168 0.0010953139 0.0007060413 0.0005980230 0.0004918448
# [8] 0.0004394740 0.0003614691 0.0003020580
p<-NULL
coef<-NULL
for (i in 2:14) {
  p<-c(p,cor.test(temp$vectors[,1],as.numeric(A932cov[,i]))$p.value)
  coef<-c(coef,cor.test(temp$vectors[,1],as.numeric(A932cov[,i]))$estimate)
}
eig1cor<-cbind(p,coef)
rownames(eig1cor)<-colnames(A932cov)[2:14]
eig1cor
#                          p          coef
#PID            0.0009614574 -1.079751e-01
#ExposureBatch  0.9999999998  6.464781e-12
#RNAbatch       0.2969478141  3.420054e-02
#HybBatch       0.9967349500  1.342228e-04
#Treated        0.0000000000  3.253716e-01
#Cell.Count     1.0000000000 -2.556283e-13
#Gender         0.9999999999 -4.391169e-12
#Age            0.9999999999 -2.114850e-12
#BMI            0.9999999999  4.562396e-12
#Smoking.Status 0.9999999999  2.116950e-12
#Date.Exposed   0.9999999998 -8.846086e-12
#RNA.Plate      0.5319245837  2.050024e-02
#Array          1.0000000000  1.493602e-12

cov<-NULL
for (i in 2:14){
  cov<-c(cov,anova(lm(temp$vectors[,1]~A932cov[,i]))$P[1])
}
names(cov)<-colnames(A950cov)[2:14]
cov
#           PID  ExposureBatch       RNAbatch       HybBatch        Treated     Cell.Count 
#9.067997e-10   1.000000e+00   2.969478e-01   1.000000e+00   2.011000e-24   1.000000e+00 
#Gender            Age            BMI Smoking.Status   Date.Exposed      RNA.Plate 
#1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00   5.319246e-01 
#Array 
#1.000000e+00  
A932cov$Treated<-factor(as.numeric(A932cov$Treated),labels=c("Control","Treated"))
png("PARCnormEig1ByExposureBatch.png")
xyplot(pcs2$vectors[,1]~A932cov$ExposureBatch, groups=A932cov$Treated,
       auto.key=list(columns=2),xlab="ExposureBatch",ylab="PC1 from Normed Data", )
dev.off()
png("PARCnormEig1ByRNABatch.png")
xyplot(pcs2$vectors[,1]~A932cov$RNAbatch, groups=A932cov$Treated,
       auto.key=list(columns=2),xlab="RNAbatch",ylab="PC1 from Normed Data", )
dev.off()
png("PARCnormEig1ByArray.png")
xyplot(pcs2$vectors[,1]~A932cov$Array, groups=A932cov$Treated,
       auto.key=list(columns=2),xlab="Array",ylab="PC1 from Normed Data", )
dev.off()
png("PARCnormEig1ByRNAPlate.png")
xyplot(pcs2$vectors[,1]~A932cov$RNA.Plate, groups=A932cov$Treated,
       auto.key=list(columns=2),xlab="RNAPlate",ylab="PC1 from Normed Data", )
dev.off()
png("PARCnormEig1ByHybBatch.png")
xyplot(pcs2$vectors[,1]~A932cov$HybBatch, groups=A932cov$Treated,
       auto.key=list(columns=2),xlab="Hybridization Batch",ylab="PC1 from Normed Data", )
dev.off()
png("PARCnormEig1ByPID.png")
xyplot(pcs2$vectors[,1]~A932cov$PID, groups=A932cov$Treated,
       auto.key=list(columns=2),xlab="Individual",ylab="PC1 from Normed Data", )
dev.off()


xyplot(pcs2$vectors[,2]~A932cov$ExposureBatch, groups=A932cov$Treated)
xyplot(pcs2$vectors[,2]~A932cov$Array, groups=A932cov$Treated)
xyplot(pcs2$vectors[,2]~A932cov$RNA.Plate, groups=A932cov$Treated)
p<-NULL
coef<-NULL
for (i in 2:14) {
  p<-c(p,cor.test(pcs2$vectors[,2],as.numeric(A932cov[,i]))$p.value)
  coef<-c(coef,cor.test(pcs2$vectors[,2],as.numeric(A932cov[,i]))$estimate)
}
eig2cor<-cbind(p,coef)
rownames(eig2cor)<-colnames(A932cov)[2:14]
eig2cor
#                          p          coef
#PID            9.717329e-01 -1.162270e-03
#ExposureBatch  1.000000e+00 -5.827167e-15
#RNAbatch       2.787194e-04  1.187837e-01
#HybBatch       9.969408e-01 -1.257612e-04
#Treated        2.168823e-23 -3.183138e-01
#Cell.Count     1.000000e+00  6.088162e-16
#Gender         1.000000e+00  1.676282e-16
#Age            1.000000e+00  4.810339e-15
#BMI            1.000000e+00 -1.065175e-14
#Smoking.Status 1.000000e+00  1.215846e-15
#Date.Exposed   1.000000e+00  1.925090e-14
#RNA.Plate      2.666966e-02  7.259956e-02
#Array          1.000000e+00  4.761706e-15

################################################################
#Place normalized data in Synapse.
################################################################
norm<-Data(list(name="CAPexpression-NormalizedForNetworksSept272012",parentId="syn301168"))
norm<-addObject(norm,A932norm)
norm<-addObject(norm,A932cov)
norm<-storeEntity(norm)
###Parent Project:  Network Analysis of Statin Response
##https://synapse.sagebase.org/#Synapse:syn301168
propertyValue(norm,"id")
onWeb(propertyValue(norm,"id"))
##Synapse ID for the normalized data:  "syn1417424"
##To view in Synapse:
##onWeb(propertyValue(norm,"id"))

