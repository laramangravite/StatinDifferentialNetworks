setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012")
library(matrixStats)
library(WGCNA)
library(lattice)
library(synapseClient)
synapseCacheDir('~/synapseCache/')
#synapseLogin(user="lara.mangravite@sagebase.org")
#Load data
PARCnorm<-loadEntity("syn1417424")
names(PARCnorm$objects)
cov<-PARCnorm$objects$A932cov
dat<-PARCnorm$objects$A932norm
bdat<-dat[,cov$Treated=="Control"]
tdat<-dat[,cov$Treated=="Treated"]

##########################
##Assess correlation of modules with clinical traits.
##############################
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules")
tnet<-read.table("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/TreatmentExpnModules/TreatmentWgcnaOct31Modules.csv",sep=",")
cnet<-read.table("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules/ControlWgcnaOct31Modules.csv",sep=",")
load(file="/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/TreatmentExpnModules/TreatmentModuleEigens.Rdat")
load(file="/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules/ControlModuleEigens.Rdat")
subjects<-colnames(bdat)
subjects<-gsub("P","",subjects)
subjects<-gsub("CA","CAP",subjects)

pheno<-read.table("/Volumes/lmangrav/PARC/CoexpNetworks/Data/D480phenotypes.TXT",sep=",",header=T)
pheno<-pheno[which(pheno[,1] %in% subjects),]
phenos<-pheno[,13:28]
cov<-pheno[,9:12]
mod<-model.matrix()
rownames(cov)<-pheno[,1]
cov[is.na(cov[,3]),3]
#Gender   Age BMI
#CAP04.1480      1 43.61  NA
grep("CAP04.1480",rownames(cov))
cov[448,3]<-(-3.674639 + 0.3357048*84) #impute missing BMI from waist.
X<-model.matrix(~cov$Gender+cov$Age+cov$BMI+cov$Smoking.Status)
data<-phenos
data[is.na(data[,2]),]
#LDLR.B LDLC.B ApoB.B TC.B TG.B HDLC.B ApoAI.B CRP.B   LDLR.T LDLC.T ApoB.T  TC.T  TG.T
#CAP02.1818 2.282508     NA     NA   NA   NA     NA      NA    NA 4.393827     82     67 173.5 144.5
#HDLC.T ApoAI.T CRP.T
#CAP02.1818   62.8     128  0.25
rownames(data)<-sub("CAP02.1818","CAP02.1808",rownames(data))
rownames(bdat)<-sub("CAP02.1818","CAP02.1808",rownames(bdat))
data[rownames(data)=="CAP02.1808",c(2:8)]<- c(132,90,222,90.5,71.7,149,0.3)
clindat<-data[,c(2:7,10:15)]
clindat<-log2(clindat)
clindat<-cbind(clindat,(clindat[,7]/clindat[,1]),(clindat[,8]/clindat[,2]),
               (clindat[,9]/clindat[,3]),(clindat[,10]/clindat[,4]),(clindat[,11]/clindat[,5]),
               (clindat[,12]/clindat[,6]))
clindat<-as.matrix(clindat)
b1 <- solve(t(X) %*% X) %*% t(X) %*% (clindat)
res1 <- t(clindat)-t(X %*% b1)
res2=res1+t(b1)[,1] #add intercept back
clindatr<-res2 #phenos by samples
rownames(clindatr)[13:18]<-c("LDLC.D","ApoB.D","TC.D","TG.D","HDLC.D","ApoAI.D")
modassocr<-c()
modassocp<-c()
library(qvalue)
for (i in 1:ncol(ceig)) {
  temp<-apply(clindatr,1,function(x) {cor.test(x,ceig[,i])$estimate})
  modassocr<-cbind(modassocr,temp)
  temp<-apply(clindatr,1,function(x) {cor.test(x,ceig[,i])$p.value})
  modassocp<-cbind(modassocp,temp)
}
apply(modassocp,1,function(x){sum(x<=0.05)})
modassocq<-t(apply(modassocp,1,function(x) {qvalue(x)$q}))
rownames(modassocr)<-rownames(clindatr)
colnames(modassocr)<-c(paste("mod",c(1:24),sep=""))
rownames(modassocp)<-rownames(clindatr)
colnames(modassocp)<-c(paste("mod",c(1:24),sep=""))
i<-16
hist(modassocp[i,],breaks=20,col="grey")
modassocp[i,]
#ApoB.B with mod21. p:0.05568555 salmon module
#LDLC.B with mod12. p: 0.04026196 lightcyan
#TG.B with mod13.
#LDLC.T with mod12. p: 0.02347116
#ApoB.T with mod2,mod21. p:0.02571987 and 0.01857225
#TG.T with mod2,mod13,mod19,mod21.
#ApoB.D with mod8. p:0.02181154, r:0.106235, q:0.523477 green module
#TG.D with mod3.
rownames(modassocp)
modassocq[1,4]

##LDLR.
sum(is.na(data[,1]))
#[1] 43
LDLR<-data[,c(1,9)]
LDLR<-LDLR[!is.na(LDLR[,1]),]
LDLR<-LDLR[!is.na(LDLR[,2]),]
LDLR<-log2(LDLR)
LDLR<-cbind(LDLR,(LDLR[,2]/LDLR[,1]))
colnames(LDLR)[3]<-"LDLR.D"
LDLR<-as.matrix(LDLR)
rownames(X)<-rownames(data)
X<-X[which(rownames(X) %in% colnames(LDLR)),]
rownames(ceig)<-rownames(data)
ceigl<-ceig
ceigl<-ceigl[which(rownames(ceigl) %in% colnames(LDLR)),]
b1 <- solve(t(X) %*% X) %*% t(X) %*% t(LDLR)
res1 <- (LDLR)-t(X %*% b1)
res2=res1+t(b1)[,1] #add intercept back
LDLR<-t(res2) #phenos by samples
ldlrr<-c()
ldlrp<-c()
for (i in 1:ncol(ceigl)) {
  temp<-apply(LDLR,2,function(x) {cor.test(x,ceigl[,i])$estimate})
  ldlrr<-cbind(ldlrr,temp)
  temp<-apply(LDLR,2,function(x) {cor.test(x,ceigl[,i])$p.value})
  ldlrp<-cbind(ldlrp,temp)
}
rownames(ldlrp)[3]<-"LDLR.D"
apply(ldlrp,1,function(x){sum(x<=0.05)})
#LDLR.B LDLR.T LDLR.D 
#4      2      0 
#LDLR.B is correlated with mod12,mod15,mod18,mod23.
#lightcyan, magenta, purple, turquoise.
#which contains GATM, none, LDLR, HMGCR.
ldlrr[1,c(12,15,18,23)]
ldlrp[1,c(12,15,18,23)]
modassocq[19,c(12,15,18,23)]
#cor:   -0.1002295  0.1239932  0.1291734 -0.1036812 
#pvals: 0.044335111 0.012737178 0.009432482 0.037476811 
#qvals: 0.15206305 0.08737337 0.08737337 0.15206305 

modassocp<-rbind(modassocp,ldlrp)
modassocr<-rbind(modassocr,ldlrr)
ldlrq<-t(apply(ldlrp,1,function(x){qvalue(x)$q}))
modassocq<-rbind(modassocq,ldlrq)
modassocp[19,]
cmodnames<-c("black","blue","brown","cyan","darkgreen,","darkred","darkturquoise","green",
            "greenyellow","grey","grey60","lightcyan","lightgreen","lightyellow","magenta",
            "midnightblue","pink","purple","red","royalblue","salmon","tan","turquoise","yellow")
colnames(modassocq)<-colnames(modassocp)<-colnames(modassocr)<-cmodnames
baselineModuleCorrelations<-list(pval=modassocp,qval=modassocq,cor=modassocr)
save(baselineModuleCorrelations,file="baselineModuleCorrelations.Rdat")

################
##Traetment module associations
####################
teig<-tMEList$eigengenes
modassocr<-c()
modassocp<-c()
#library(qvalue)
for (i in 1:ncol(teig)) {
  temp<-apply(clindatr,1,function(x) {cor.test(x,teig[,i])$estimate})
  modassocr<-cbind(modassocr,temp)
  temp<-apply(clindatr,1,function(x) {cor.test(x,teig[,i])$p.value})
  modassocp<-cbind(modassocp,temp)
}
modassocq<-t(apply(modassocp,1,function(x) {qvalue(x)$q}))
rownames(modassocr)<-rownames(clindatr)
rownames(modassocp)<-rownames(clindatr)
rownames(modassocq)<-rownames(clindatr)
tmodnames<-c("black","blue","brown","cyan","darkgreen","darkgrey","darkorange","darkred",
             "darkturquoise","green","greenyellow","grey","grey60","lightcyan","lightgreen",
             "lightyellow","magenta","midnightblue","orange","pink","purple","red","royalblue",
             "salmon","tan","turquoise","yellow")
colnames(modassocr)<-tmodnames
colnames(modassocp)<-tmodnames
colnames(modassocq)<-tmodnames

teigl<-teig
rownames(teigl)<-rownames(data)
teigl<-teigl[which(rownames(teigl) %in% rownames(LDLR)),]
ldlrr<-c()
ldlrp<-c()
for (i in 1:ncol(teigl)) {
  temp<-apply(LDLR,2,function(x) {cor.test(x,teigl[,i])$estimate})
  ldlrr<-cbind(ldlrr,temp)
  temp<-apply(LDLR,2,function(x) {cor.test(x,teigl[,i])$p.value})
  ldlrp<-cbind(ldlrp,temp)
}
ldlrq<-t(apply(ldlrp,1,function(x) {qvalue(x)$q}))
modassocp<-rbind(modassocp,ldlrp)
modassocq<-rbind(modassocq,ldlrq)
modassocr<-rbind(modassocr,ldlrr)


###
apply(modassocp,1,function(x){sum(x<=0.05)})
#LDLC.B  ApoB.B    TC.B    TG.B  HDLC.B ApoAI.B  LDLC.T  ApoB.T    TC.T    TG.T  HDLC.T ApoAI.T 
#2       0       2       0       0       0       0       0       0       0       0       1 
#LDLC.D  ApoB.D    TC.D    TG.D  HDLC.D ApoAI.D   LDLR.B  LDLR.T  LDLR.D 
#2       0       0       0       0       0        7       2       0 
apply(modassocq,1,function(x){sum(x<=0.05)})  
#LDLC.B  ApoB.B    TC.B    TG.B  HDLC.B ApoAI.B  LDLC.T  ApoB.T    TC.T    TG.T  HDLC.T ApoAI.T 
#0       0       0       0       0       0       0       0       0       0       0       0 
#LDLC.D  ApoB.D    TC.D    TG.D  HDLC.D ApoAI.D  LDLR.B  LDLR.T  LDLR.D 
#0       0       0       0       0       0       9       0       0 
modassocp[19,modassocp[19,]<=0.05]
#Baseline LDLR
#darkgreen    darkgrey     darkred       green   lightcyan         red         tan 
#0.011651296 0.037491082 0.017256499 0.035436838 0.006163726 0.019620659 0.025226880 
modassocp[20,modassocp[20,]<=0.05]
#Treated LDLR
#darkgreen      green 
#0.02211558 0.04465244 
modassocp[1,modassocp[1,]<=0.05]
#BAseline LDLC
#darkred  turquoise 
#0.03870458 0.04318315 
modassocp[13,modassocp[13,]<=0.05]
#Delta LDLC
#darkturquoise         green 
#0.04973938    0.04620837 

treatmentModuleCorrelations<-list(pval=modassocp,qval=modassocq,cor=modassocr)
save(treatmentModuleCorrelations,file="TreatmentModuleCorrelations.Rdat")

##########################################
##Comparison of modules across networks###
##########################################