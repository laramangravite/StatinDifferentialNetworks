#######August 30 2012
#######Normalizing with less-granular Hyb and RNA labeling batches.
setwd("/gluster/home/lmangrav/PARC/CoexpNetworks/CoexpnJuly2012")
library(sva)
library(gcrma)
library(matrixStats)
library(WGCNA)
library(lattice)
library(synapseClient)
synapseCacheDir('~/synapseCache/')
synapseLogin(user="lara.mangravite@sagebase.org")
#################################August 16th: Resumming Normalization Analysis.
PARCraw <- loadEntity("syn299515")
Braw<-PARCraw$objects$Br
Bcov<-PARCraw$objects$B480Cov
Traw<-PARCraw$objects$Tr
Tcov<-PARCraw$objects$T480Cov
Araw<-rbind(Braw,Traw)
Araw<-log2(Araw)
Bcov<-cbind(Bcov,as.factor(c(rep(0,480))))
Tcov<-cbind(Tcov,as.factor(c(rep(1,480))))
colnames(Bcov)[9]<-"Treated"
colnames(Tcov)[9]<-"Treated"
Bcov<-cbind(Bcov,rownames(Bcov))
Tcov<-cbind(Tcov,rownames(Tcov))
colnames(Bcov)[10]<-"PID"
colnames(Tcov)[10]<-"PID"
rownames(Bcov)<-paste(rownames(Bcov),"P",sep="")
rownames(Tcov)<-paste(rownames(Tcov),"S",sep="")
Acov<-rbind(Bcov,Tcov)
Acov<-cbind(Acov,rownames(Acov))
colnames(Acov)[11]<-"SampleID"
###############Prep the hybridization batch data.
files<-c(
  #"mangravite_11.12.08.csv",
  "mangravite_11.18.08.csv","mangravite_11.20.08.csv","mangravite_11.24.08.csv","mangravite_11.25.08.csv","mangravite_12.01.08.csv","mangravite_12.02.08.csv","mangravite_12.03.08.csv","mangravite_12.04.08.csv",
         "mangravite_3.09.09.csv","mangravite_3.10.09.csv","mangravite_3.11.09.csv","mangravite_3.12.09.csv","mangravite_3.23.09.csv","mangravite_3.24.09.csv",
         "mangravite_3.25.09.csv","mangravite_3.26.09.csv","mangravite_3.30.09.csv","mangravite_3.31.09.csv","mangravite_4.01.09.csv","mangravite_4.02.09.csv","mangravite_4.13.09.csv","mangravite_8.24.09.csv",
         "mangravite_8.25.09.csv","mangravite_8.26.09.csv","mangravite_8.27.09.csv","mangravite_8.31.09.csv","mangravite_9.01.09.csv","mangravite_9.02.09.csv","mangravite_9.03.09.csv","mangravite_9.14.09.csv",
         "mangravite_9.15.09.csv","mangravite_9.16.09.csv","mangravite_9.17.09.csv","mangravite_9.21.09.csv","mangravite_9.22.09.csv","mangravite_9.23.09.csv","mangravite_9.24.09.csv","mangravite_9.28.09.csv",
         "mangravite_9.29.09.csv","mangravite_9.30.09.csv","mangravite_11.23.09.csv","mangravite_12.22.09.csv","mangravite_12.23.09.csv")
temp<-NULL
arrays<-NULL
for (i in 1:length(files)){
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
colnames(temp)<-c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position")
temp[,8]<-files[i]
temp[,9]<-i
arrays<-rbind(arrays,temp)
}
arrays<-arrays[,-c(1:5)]
colnames(arrays)[3:4]<-c("HybDate","HybBatch")
array<-unique(arrays[,c(1,3:4)])

Acov<-merge(Acov,array,by.x=4,by.y=1,all.x=T,all.y=F)            
Acov[is.na(Acov[,12]),12]<-"DoseResponse"
Acov[is.na(Acov[,13]),13]<-0
Acov[,13]<-Acov[,13]+1
Acov[,13]<-as.factor(Acov[,13])
rownames(Acov)<-Acov$SampleID
summary(Acov$HybBatch)
##Hybridization batches.
#Batch    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 
#Narrays 17 14 23 24 11 25 21 27 13 27 15 23 22 16 30 30 25 28 27 30 23 13 24 27 28 28 23 24 26 23 23 17 

#Batch   33 34 35 36 37 38 39 40 41 42 43 44 
#Narrays 14 21 24 25 23 27 24 27  6 14 26  2 

dbatch<-as.matrix(unique(Acov$Date.Exposed[order(as.numeric(Acov$Date.Exposed))]))
dbatch<-cbind(dbatch,c(1:46))
colnames(dbatch)<-c("date","ExposureBatch")
Acov<-merge(Acov,dbatch,by.x=2,by.y=1,all.x=T,all.y=F)
summary(Acov$ExposureBatch)
#Ebatch    1 10 11 12 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 35 36 37 38 
#Nsamples 10 18 24 24 28 18 34 26 30 25 35 12 34 40 28 29 24 28 28  8 32 20 14 28 18 34 20 26 34 14 32 15 

#Ebatch   39  4 40 41 42 43 44 45 46  5  6  7  8  9 
#Nsamples 22  8 24  2 14  2 13  2  2 16 18 12 20 15 

##Convert RNA plate to RNA batch.
rbatch<-(unique(Acov$RNA.Plate))
rbatch<-rbatch[order(rbatch)]
rbatch<-cbind(rbatch,c(1,rep(2,4),rep(3,4),rep(4,4),rep(5,3)))
rbatch[,1]<-rbatch[,1]-1
colnames(rbatch)<-c("RNAplate","RNAbatch")
Acov<-merge(Acov,rbatch,by.x=4,by.y=1,all.x=T,all.y=F)
rownames(Acov)<-Acov$SampleID
Acov<-Acov[,c(11,10,14,15,13,9,4,5,6,7,8,2,1,3,12)]
write.table(Acov,file="AllSampleCovariatesAug302012.txt",row.names=T,col.names=T)
rm(dbatch,rbatch,array,arrays,files,i)

####################Analysis of Normalization with E,R,H batches.
png("ExposureBatchVsRNABatch.png")
plot(as.numeric(Acov[,4]),as.numeric(Acov[,3]),xlab=colnames(Acov)[4],ylab=colnames(Acov)[3])
abline(lm(as.numeric(Acov[,3])~as.numeric(Acov[,4])),col="red")
dev.off()
png("ExposureBatchVsArray.png")
plot(as.numeric(Acov[,4]),as.numeric(Acov[,5]),xlab=colnames(Acov)[4],ylab=colnames(Acov)[5])
abline(lm(as.numeric(Acov[,4])~as.numeric(Acov[,1])))
dev.off()
png("ArrayVsRNABatch.png")
plot(as.numeric(Acov[,3]),as.numeric(Acov[,5]),xlab=colnames(Acov)[3],ylab=colnames(Acov)[5])
abline(lm(as.numeric(Acov[,4])~as.numeric(Acov[,3])))
dev.off()

#################Check batch data and update for accuracy.  See "UpdatingBatchInfo.R".
##There are some samples for which the array batch predates the RNA labeling batch, which is not possible.
##These were probably samples that were repeated, with the RNA labeling but not array batch info having been updated.
##Load data from Synapse
#outliers from RNA batch vs. array batch plot.
#CAP03.1371S #This is a dose response sample, so is exp batch 1 but a late RNA and array batch.
#CAP04.1201P #This is a dose response sample, so is exp batch 1 but a late RNA and array batch.
#CAP03.1371P #This is a dose response sample, so is exp batch 1 but a late RNA and array batch.
#CAP04.1201S #This is a dose response sample, so is exp batch 1 but a late RNA and array batch.
#CAP03.1918S #Is correct
#CAP03.1358P #Updated
#CAP03.0704P #Updated
#CAP03.2615P #Updated
#CAP04.0243P #Updated

##########Assess the correct batch info for each of these samples using the array hyb run sheets.
files<-c(
  #"mangravite_11.12.08.csv",
  "mangravite_11.18.08.csv","mangravite_11.20.08.csv","mangravite_11.24.08.csv","mangravite_11.25.08.csv","mangravite_12.01.08.csv","mangravite_12.02.08.csv","mangravite_12.03.08.csv","mangravite_12.04.08.csv",
  "mangravite_3.09.09.csv","mangravite_3.10.09.csv","mangravite_3.11.09.csv","mangravite_3.12.09.csv","mangravite_3.23.09.csv","mangravite_3.24.09.csv",
  "mangravite_3.25.09.csv","mangravite_3.26.09.csv","mangravite_3.30.09.csv","mangravite_3.31.09.csv","mangravite_4.01.09.csv","mangravite_4.02.09.csv","mangravite_4.13.09.csv","mangravite_8.24.09.csv",
  "mangravite_8.25.09.csv","mangravite_8.26.09.csv","mangravite_8.27.09.csv","mangravite_8.31.09.csv","mangravite_9.01.09.csv","mangravite_9.02.09.csv","mangravite_9.03.09.csv","mangravite_9.14.09.csv",
  "mangravite_9.15.09.csv","mangravite_9.16.09.csv","mangravite_9.17.09.csv","mangravite_9.21.09.csv","mangravite_9.22.09.csv","mangravite_9.23.09.csv","mangravite_9.24.09.csv","mangravite_9.28.09.csv",
  "mangravite_9.29.09.csv","mangravite_9.30.09.csv","mangravite_11.23.09.csv","mangravite_12.22.09.csv","mangravite_12.23.09.csv")
temp<-NULL
arrays<-NULL

##search for outlier sample 1.
k<-"04-0243"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  cat(files[i],"   ") 
  hold<-
    j<-(grep("04-0243",temp[,1]))
  if(length(j>0)){
    cat(length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}
#k<-"04-0243"
#"mangravite_11.18.08.csv"  32 04-0243_P NA NA All NA 4256290107  H
#"mangravite_11.24.08.csv"  32 04-0243_S NA NA All NA 4256425105  H
#"mangravite_9.30.09.csv"   1 04-0243P D3 NA All NA 4482545215  A
##Placebo sample was run in array batch 2 and 41.  RNA batch was 5 (plate 14), so it must be array batch 41.
#Change array batch to 41 within Acov file.
l<-"CAP04.0243P"
B<-41
f<-"mangravite_9.30.09.csv"
a<-"4482545215"
Acov[Acov$SampleID==l,]
Acov[Acov$SampleID==l,]$HybBatch<-B
Acov[Acov$SampleID==l,]$HybDate<-f
Acov[Acov$SampleID==l,]$Array<-a
Acov[Acov$SampleID==l,]

k<-"03-0704"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  cat(files[i],"   ") 
  hold<-
    j<-(grep(k,temp[,1]))
  if(length(j>0)){
    cat(i, length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}

temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]

#k<-"03-0704"
#mangravite_3.09.09.csv    30 03-0704_P NA NA All NA 4482545422  D 
#mangravite_3.25.09.csv    14 03-0704_S NA NA All NA 4521621077  D
#mangravite_8.31.09.csv    20 03-0704P NA NA All NA 4482545382  G

##Placebo sample was run in array batch 10 and 27.  RNA batch was 5 (plate 13), so it must be array batch 27.
#Change array batch to 41 within Acov file.
l<-"CAP03.0704P"
B<-27
f<-"mangravite_8.31.09.csv"
a<-"4482545382"
Acov[Acov$SampleID==l,]
Acov[Acov$SampleID==l,]$HybBatch<-B
Acov[Acov$SampleID==l,]$HybDate<-f
Acov[Acov$SampleID==l,]$Array<-a
Acov[Acov$SampleID==l,]

k<-"03-1358"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  cat(files[i],"   ") 
  hold<-
    j<-(grep(k,temp[,1]))
  if(length(j>0)){
    cat(i, length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}

temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]

#k<-"03-1358"
#mangravite_3.09.09.csv    15 03-1358_S NA NA All NA 4482545424  F 
#mangravite_3.11.09.csv    3 03-1358_P NA NA All NA 4457315028  E 
#mangravite_8.27.09.csv    28 03-1358P NA NA All NA 4482545298  G 
##Placebo sample was run twice.  RNA batch was 5 (plate 13), so it must be the later array batch.
#Change array batch to 41 within Acov file.
l<-"CAP03.1358P"
Acov[Acov$SampleID==l,]
f<-"mangravite_8.27.09.csv"
B<-as.numeric(unique(Acov$HybBatch[Acov$HybDate==f]))
a<-"4482545382"
Acov[Acov$SampleID==l,]$HybBatch<-B
Acov[Acov$SampleID==l,]$HybDate<-f
Acov[Acov$SampleID==l,]$Array<-a
Acov[Acov$SampleID==l,]

k<-"03-2615"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  j<-(grep(k,temp[,1]))
  if(length(j>0)){
    cat(files[i],"   ") 
    cat(i, length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}

temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
#mangravite_12.03.08.csv    23 03-2615_S NA NA All NA 4482545233  F
#mangravite_12.03.08.csv    31 03-2615_S1 NA NA All NA 4482545261  F 
#mangravite_12.04.08.csv     7 03-2615_P NA NA All NA 4482545160  F 
#mangravite_9.17.09.csv     22 03-2615S NA NA All NA 4521621003  D
##These are both wrong.  
#Placebo:  says RNA batch 13 but array batch 9.  Doesn't look to have been re-run. Change RNA batch back to plate 2.
#SVA:  we reran the sample but it appears that the data wasn't used. Leave as is. (RNA batch 2, HybBatch 8)
l<-"CAP03.2615P"
Acov[Acov$SampleID==l,]
p<-2
r<-as.numeric(unique(Acov$RNAbatch[Acov$RNA.Plate==p]))
Acov[Acov$SampleID==l,]$RNAbatch<-r
Acov[Acov$SampleID==l,]$RNA.Plate<-p
Acov[Acov$SampleID==l,]

k<-"04-1201"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  j<-(grep(k,temp[,1]))
  if(length(j>0)){
    cat(files[i],"   ") 
    cat(i, length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}
i<-32
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
i<-35
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
#04-1201S NA NA All NA 4521621035  D
#04-1201P NA NA All NA 4521621091  G 
##This is one of the original dose response samples so has early exposure batch.
##We ran these on the later arrays.  I'm sure we relabeled them.  Change label to last batch/plate.
l<-"CAP04.1201S"
Acov[Acov$SampleID==l,]
p<-15
r<-as.numeric(unique(Acov$RNAbatch[Acov$RNA.Plate==p]))
Acov[Acov$SampleID==l,]$RNAbatch<-r
Acov[Acov$SampleID==l,]$RNA.Plate<-p
Acov[Acov$SampleID==l,]
l<-"CAP04.1201P"
Acov[Acov$SampleID==l,]
Acov[Acov$SampleID==l,]$RNAbatch<-r
Acov[Acov$SampleID==l,]$RNA.Plate<-p
Acov[Acov$SampleID==l,]

##CAP03.1918S
k<-"03-1918"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  j<-(grep(k,temp[,1]))
  if(length(j>0)){
    cat(files[i],"   ") 
    cat(i, length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}
i<-33
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
i<-7
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
#20  03-1918_P NA NA All NA 4482545233  G
#28 03-1918_P1 NA NA All NA 4482545261  G
#30 03-1918S_1 NA NA All NA 4521621017  D
#Look at both.
l<-"CAP03.1918"
Acov[Acov$PID==l,]
##SVA sample failed RNA labeling or array first time.  Was rerun at a later date. Leave data as is.

k<-"03-1371"
for (i in 1:length(files)){
  temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
  j<-(grep(k,temp[,1]))
  if(length(j>0)){
    cat(files[i],"   ") 
    cat(i, length(grep(k,temp[,1])),"\n")
  }else{
    cat("none", "\n")
  }
}
i<-38
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
files[i]
i<-30
temp<-read.table(file=paste("~/PARC/CoexpNetworks/Data/PARCexpnBatches/",files[i],sep=""),sep=",",na.strings="",skip=9)
temp[grep(k,temp[,1]),]
files[i]
#"mangravite_9.28.09.csv" 03-1371P NA NA All NA 4482545411  H
# "mangravite_9.14.09.csv" 16 03-1371S NA NA All NA 4513250088  H
l<-"CAP03.1371"
Acov[Acov$PID==l,]
##This is also a dose response sample.  Exposure is early but RNA and array batch are late.
##Update RNA batch to be 15/5.
p<-15
r<-as.numeric(unique(Acov$RNAbatch[Acov$RNA.Plate==p]))
Acov[Acov$PID==l,]$RNAbatch<-r
Acov[Acov$PID==l,]$RNA.Plate<-p
Acov[Acov$PID==l,]


plot(as.numeric(Acov[,4]),as.numeric(Acov[,3]),xlab=colnames(Acov)[4],ylab=colnames(Acov)[3])
abline(lm(as.numeric(Acov[,3])~as.numeric(Acov[,4])),col="red")

plot(as.numeric(Acov[,4]),as.numeric(Acov[,5]),xlab=colnames(Acov)[4],ylab=colnames(Acov)[5])
abline(lm(as.numeric(Acov[,4])~as.numeric(Acov[,1])))

plot(as.numeric(Acov[,3]),as.numeric(Acov[,5]),xlab=colnames(Acov)[3],ylab=colnames(Acov)[5])
abline(lm(as.numeric(Acov[,4])~as.numeric(Acov[,3])))


Bcov<-Acov[Acov$Treated==0,]
Tcov<-Acov[Acov$Treated==1,]
Bcov<-Bcov[order(Bcov$PID),]
Tcov<-Tcov[order(Tcov$PID),]
Acov<-rbind(Bcov,Tcov)
PARCraw<-deleteObject(PARCraw, "B480Cov")
PARCraw<-deleteObject(PARCraw, "T480Cov")
PARCraw<-storeEntity(PARCraw)

B480cov<-Bcov
T480cov<-Tcov
All960cov<-Acov
A960rawLogged<-Araw
PARCraw<-addObject(PARCraw,"B480cov")
PARCraw<-addObject(PARCraw,"T480cov")
PARCraw<-addObject(PARCraw,"All960cov")
PARCraw<-addObject(PARCraw,"A960rawLogged")
PARCraw<-storeEntity(PARCraw)
PARC<-loadEntity("syn299515")
PARC<-deleteObject(PARC, "B480cov")
PARC<-deleteObject(PARC, "T480cov")
PARC<-deleteObject(PARC, "A960rawLogged")
PARC<-deleteObject(PARC, "All960cov")
PARC<-deleteObject(PARC, "B480Cov")
PARC<-deleteObject(PARC, "T480Cov")
PARC<-storeEntity(PARC)
synapseCacheDir('~/synapseCache/try4')
try4<-loadEntity("syn299515")
try4<-addObject(try4,"All960cov")
try4<-updateEntity(try4)
try4<-storeEntity(try4)
B480raw<-PARCraw$objects$Br
T480raw<-PARCraw$objects$Tr
dat<-Data(list(name="RawExpressionCovariateData",parentId="syn299511"))
dat<-createEntity(dat)
dat<-addObject(dat,Araw)
dat<-deleteObject(dat,"Araw")
A960raw<-rbind(B480raw,T480raw)
dat<-addObject(dat,A960raw)
dat<-addObject(dat,All960cov)
dat<-addObject(dat,Bcov)
dat<-addObject(dat,Tcov)
dat<-addObject(dat,B480raw)
dat<-addObject(dat,T480raw)
dat<-updateEntity(dat)
dat<-storeEntity(dat)
synapseCacheDir('~/synapseCache/try5')
reload<-loadEntity("syn1417157")
##try to update previous entity objects
synapseCacheDir('~/synapseCache/try5')
PARCdat<-loadEntity("syn299515")
names(PARCdat$objects)
dim(PARCdat$objects$A960rawLogged)
dim(PARCdat$objects$Br)
PARCdat<-deleteObject(PARCdat, "A960rawLogged")
PARCdat<-deleteObject(PARCdat, "B480Cov")
PARCdat<-deleteObject(PARCdat, "Br")
PARCdat<-deleteObject(PARCdat, "T480Cov")
PARCdat<-deleteObject(PARCdat, "Tr")
PARCdat<-updateEntity(PARCdat)

A960raw<-rbind(B480raw,T480raw)
PARCdat<-addObject(PARCdat,A960raw)
PARCdat<-addObject(PARCdat,All960cov)
PARCdat<-addObject(PARCdat,Bcov)
PARCdat<-addObject(PARCdat,Tcov)
PARCdat<-addObject(PARCdat,B480raw)
PARCdat<-addObject(PARCdat,T480raw)
PARCdat<-storeEntity(PARCdat)

synapseCacheDir('~/synapseCache/try6')
PARCdat2<-loadEntity("syn299515")
PARCdat2<-deleteEntity(PARCdat2,"B480Cov")
##create a new entity and add objects.

dat<-Data(list(name="RawExpressionDataWithCovariates",parentId="syn299511"))
dat<-createEntity(dat)
A960raw<-rbind(Braw,Traw)
dat<-addObject(dat,A960raw)
dat<-addObject(dat,All960cov)
dat<-addObject(dat,B480cov)
dat<-addObject(dat,T480cov)
dat<-addObject(dat,B480raw)
dat<-addObject(dat,T480raw)
dat<-storeEntity(dat)
i<-propertyValue(dat,"id")
synapseCacheDir('~/synapseCache/try7')
dati<-loadEntity(i) ###"syn1417162"

##Remove old entities.
PARCraw<-loadEntity("syn299515")
PARCraw<-deleteEntity(PARCraw)
dat<-loadEntity("syn1417157")
dat<-deleteEntity(dat)


#####This is the PARC-CAP data with updated batch information.
synapseCacheDir('~/synapseCache')
PARCraw<-loadEntity("syn1417162")

