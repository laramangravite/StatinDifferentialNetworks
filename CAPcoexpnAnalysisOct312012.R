######################################################
#Running Coexpression Networks on CAP treatment data
##Sept 2012, Lara Mangravite, Sage Bionetworks
######################################################
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

######################################################
#WGCNA on Baseline data.
######################################################

library(cluster)
softPower = 4;
bdat<-t(bdat)
badjacency = adjacency(bdat, power = softPower);

TOM = TOMsimilarity(badjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
#Detecting clusters...
#..cutHeight not given, setting it to 0.995  ===>  99% of the (truncated) height range in dendro.
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#black          blue         brown          cyan     darkgreen 
#293          1606          1051           146            48 
#darkred darkturquoise         green   greenyellow          grey 
#72            45           423           187            20 
#grey60     lightcyan    lightgreen   lightyellow       magenta 
#111           128           107            93           213 
#midnightblue          pink        purple           red     royalblue 
#136           262           204           372            78 
#salmon           tan     turquoise        yellow 
#159           159          3356           926 

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
png("controlWGCNATree.png")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
cgeneMod<-cbind(rep(1,10195),colnames((bdat)),dynamicColors)
write.table(cgeneMod,file="ControlWgcnaOct31Modules.csv", sep=",",col.names=F,row.names=F)
source("/gluster/home/lmangrav/PARC/CoexpNetworks/Scripts/FastGO.R")
source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/FastGO.R")
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules")
FastModuleOntologyAnalysis ("ControlWgcnaOct31Modules.csv","/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/TGI-GO-HumanV3-16773gsymbols_trimed.xls",identifier_col=1,gene_symbolcol=2, signifLevel=0.001, ctopNrows=20, background=0, ReportGenes=T)

cgeneMod[grep("HMGCR",cgeneMod[,2]),]
#turquoise
cgeneMod[grep("LDLR",cgeneMod[,2]),]
#purple
cgeneMod[grep("GATM",cgeneMod[,2]),]
#lightcyan
cgeneMod[grep("SQLE",cgeneMod[,2]),]
#brown

MEList = moduleEigengenes(bdat, colors = dynamicColors)
cMEs = MEList$eigengenes
write.table(cMEs,"ControlWGCNAModulePCs.txt")
save(MEList,file="ControlModuleEigens.Rdat")
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
png("ModuleEigenCluster-controldata.png")
plot(METree, main = "Clustering of module eigengenes: cutHeight 0.98, Power 4, B466",
     xlab = "", sub = "")
dev.off()

######################################################
#WGCNA on Treatment data.
######################################################

tdat<-t(tdat)
tadjacency = adjacency(tdat, power = softPower);
TOM = TOMsimilarity(tadjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
#Detecting clusters...
#..cutHeight not given, setting it to 0.995  ===>  99% of the (truncated) height range in dendro.
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#black          blue         brown          cyan     darkgreen      darkgrey 
#336          1444          1105           263            82            45 
#darkorange       darkred darkturquoise         green   greenyellow          grey 
#31           108            77           558           277            14 
#grey60     lightcyan    lightgreen   lightyellow       magenta  midnightblue 
#160           160           154           153           322           175 
#orange          pink        purple           red     royalblue        salmon 
#41           326           285           526           115           268 
#tan     turquoise        yellow 
#276          1914           980 

# Plot the dendrogram and colors underneath
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/TreatmentExpnModules")
sizeGrWindow(8,6)
png("treatmentWGCNATree.png")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
tgeneMod<-cbind(rep(1,10195),colnames((tdat)),dynamicColors)
write.table(tgeneMod,file="TreatmentWgcnaOct31Modules.csv", sep=",",col.names=F,row.names=F)
source("/gluster/home/lmangrav/PARC/CoexpNetworks/Scripts/FastGO.R")
source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/FastGO.R")

FastModuleOntologyAnalysis ("TreatmentWgcnaOct31Modules.csv","/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/TGI-GO-HumanV3-16773gsymbols_trimed.xls",identifier_col=1,gene_symbolcol=2, signifLevel=0.001, ctopNrows=20, background=0, ReportGenes=T)

tgeneMod[grep("HMGCR",tgeneMod[,2]),]
#red
tgeneMod[grep("LDLR",tgeneMod[,2]),]
#darkred
tgeneMod[grep("GATM",tgeneMod[,2]),]
#yellow
tgeneMod[grep("SQLE",tgeneMod[,2]),]
#brown

tMEList = moduleEigengenes(tdat, colors = dynamicColors)
tMEs = tMEList$eigengenes
write.table(tMEs,"TreatmentWGCNAModulePCs.txt")
save(tMEList,file="TreatmentModuleEigens.Rdat")
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(tMEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
png("ModuleEigenCluster-Treateddata.png")
plot(METree, main = "Clustering of module eigengenes: cutHeight 0.98, Power 4, B466",
     xlab = "", sub = "")
dev.off()

#######################
##Which modules are significantly correlated with each other?
######################

##Control data
dim(cMEs)
pccor<-cor(cMEs)
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules")
#source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")
#cmodcor<-build_correlation_graph_fdr(gex=cMEs,gexnames=colnames(cMEs))
qvalue_approx <- function(pval)
{
  m <- length(pval)
  
  u <- sort(pval,index.return=TRUE)
  qval <- lapply(1:m, function(x) { return((u$x[x]*m)/x) } )
  qval <- unlist(qval)
  names(qval)<-u$ix
  qval<- qval[order(as.integer(names(qval)))]
  names(qval)<-NULL
  return(qval)
}

# input is an n x p matrix where n is samples and p is gene-residuals.
gex<-cMEs
  n = nrow(gex) #number of samples
  p = ncol(gex) #number of genes
  betahat = cor((gex)) #gene by gene correlation matrix
  se.betahat = sqrt((1-betahat^2)/n)  
  pval = 2*pnorm(-abs((betahat/se.betahat)[se.betahat>6.902830e-10])) #vector of length p*(p-1) with pvalues for all correlations except self-correlations (gene1 x gene1)
  qval <- qvalue_approx(as.vector(pval)) #vector of length p*(p-1) of qvals
  qval <- cbind(0, matrix(qval<0.05,p,p-1,byrow=T)) 
  #builds matrix with 1st col is zeros representing genex x genex correlations.
  #then qvalues start filling in by column with column 1 representing gene 1 correlations, etc
  #there are p (N gene) rows per column and p (N gene) columns.
  print(dim(qval))
  for (i in 2:p){
    if (i < p){
      qval[i,]<-c(qval[i,2:i],qval[i,1],qval[i,(i+1):p]) 
    } 
    else if (i == p) {
      qval[i,]<-c(qval[i,2:i],qval[i,1]) }
  }
cModCor<-list(edges=qval,corm=betahat)
sum(cModCor$edges > 0)
#[1] 350
length(cModCor$edges)
#[1] 576
summary(abs(as.numeric(cModCor$corm)))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000749 0.0622700 0.1909000 0.2763000 0.4230000 1.0000000 
summary(abs(as.numeric(cModCor$corm)[as.vector(cModCor$edges == 1)]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.09922 0.19380 0.34140 0.36120 0.48890 0.89170
c<-as.vector(cModCor$corm)[as.vector(cModCor$corm) <= 0.99]

e<-as.vector(cModCor$edges)[as.vector(cModCor$corm) <= 0.99]
summary(abs(c[e == 0]))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#7.494e-05 2.213e-02 3.985e-02 4.320e-02 6.653e-02 9.655e-02 
sum(abs(c) >= 0.2)
#260, meaning 130 connections.  what does that look like?
save(cModCor,file="BetweenModuleCorControl.Rdat")

##Treatment data
dim(tMEs)
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/TreatmentExpnModules")
#source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")
#cmodcor<-build_correlation_graph_fdr(gex=cMEs,gexnames=colnames(cMEs))
qvalue_approx <- function(pval)
{
  m <- length(pval)
  
  u <- sort(pval,index.return=TRUE)
  qval <- lapply(1:m, function(x) { return((u$x[x]*m)/x) } )
  qval <- unlist(qval)
  names(qval)<-u$ix
  qval<- qval[order(as.integer(names(qval)))]
  names(qval)<-NULL
  return(qval)
}

# input is an n x p matrix where n is samples and p is gene-residuals.
gex<-tMEs
n = nrow(gex) #number of samples
p = ncol(gex) #number of genes
betahat = cor((gex)) #gene by gene correlation matrix
se.betahat = sqrt((1-betahat^2)/n)  
pval = 2*pnorm(-abs((betahat/se.betahat)[se.betahat>1.380566e-09])) #vector of length p*(p-1) with pvalues for all correlations except self-correlations (gene1 x gene1)
qval <- qvalue_approx(as.vector(pval)) #vector of length p*(p-1) of qvals
qval <- cbind(0, matrix(qval<0.05,p,p-1,byrow=T)) 
#builds matrix with 1st col is zeros representing genex x genex correlations.
#then qvalues start filling in by column with column 1 representing gene 1 correlations, etc
#there are p (N gene) rows per column and p (N gene) columns.
print(dim(qval))
for (i in 2:p){
  if (i < p){
    qval[i,]<-c(qval[i,2:i],qval[i,1],qval[i,(i+1):p]) 
  } 
  else if (i == p) {
    qval[i,]<-c(qval[i,2:i],qval[i,1]) }
}
tModCor<-list(edges=qval,corm=betahat)
sum(tModCor$edges > 0)
#[1] 484
length(tModCor$edges)
#[1] 729
summary(abs(as.numeric(tModCor$corm)[as.vector(tModCor$edges == 1)]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.09843 0.17420 0.31180 0.35700 0.49490 0.94220
t<-as.vector(tModCor$corm)[as.vector(tModCor$corm) <= 0.99]
e<-as.vector(tModCor$edges)[as.vector(tModCor$corm) <= 0.99]
summary(abs(t[e == 0]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00207 0.01930 0.04445 0.04411 0.06191 0.09635 
length(t[e == 0])
#218
sum(abs(t[e == 1]) >= 0.2)
#342
save(tModCor,file="BetweenModuleCorTreated.Rdat")

##############################
####Prepare for cytoscape.
###############################

##Control data.
modcon<-c()
for (i in 1:p){
  temp<-c()
  for (j in 1:p){
    if (cModCor$edges[(i+1),j] == 1) {
      temp<-c(rownames(cModCor$corm)[(i+1)],colnames(cModCor$corm)[j],cModCor$corm[(i+1),j])
    } else {
      temp<-NULL
    }
    modcon<-rbind(modcon,temp)
  }
}

###########################
##Prepare HMGCR-containing modules for cytoscape analysis.
############################
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules")

cgeneMod<-read.table(file="/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules/ControlWgcnaOct31Modules.csv",sep=",",header=F)
turq<-cgeneMod[cgeneMod[,3]=="turquoise",]
load("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/DiffCorNetworks/ControlCorrelationFull.Rdat")

library(synapseClient)
synapseCacheDir('~/synapseCache/')
#synapseLogin(user="lara.mangravite@sagebase.org")
#Load data
PARCnorm<-loadEntity("syn1417424")
names(PARCnorm$objects)
cov<-PARCnorm$objects$A932cov
dat<-PARCnorm$objects$A932norm
bdat<-dat[,cov$Treated=="Control"]
bdat<-bdat[which(rownames(bdat) %in% turq[,2]),]
genes<-rownames(bdat)
source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")
control<-build_correlation_graph(gex1=bdat,gexnames=genes)
save(control,file="ControlCorrelationTurqModule.Rdat")
ccor<-control$corg
cedges<-control$edges
ccor<-ccor[which(rownames(dat) %in% turq[,2]),which(rownames(dat) %in% turq[,2])]
cedges<-cedges[which(rownames(dat) %in% turq[,2]),which(rownames(dat) %in% turq[,2])]
genenames<-red
#prep for cytoscape
modcor<-c()
for (i in 1:3356){
  temp<-c()
  for (j in ((i+1):3356)) {
    if (cedges[i,j] == 1 & abs(ccor[i,j]) >= 0.2) {
      temp<-data.frame(Node1=genenames[i],Node2=genenames[j])
      modcor<-rbind(modcor,temp)
    } else  {
      next
    }
  }
}
write.table(modcor,file="controlTurquoiseModuleForCyto.txt")


###treated dataset.
tgeneMod<-read.table(file="/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/TreatmentExpnModules/TreatmentWgcnaOct31Modules.csv",sep=",",header=F)
red<-tgeneMod[tgeneMod[,3]=="red",]
length(turq[which(turq[,2] %in% red[,2]),2])
#427 are shared between.
#turq is 3356
#red is 526
findt<-tgeneMod[which(tgeneMod[,2] %in% turq[,2]),]
summary(findt[,3])
#black          blue         brown          cyan     darkgreen      darkgrey    darkorange       darkred darkturquoise 
#48           235           505            20             3             9            25            35             7 
#green   greenyellow          grey        grey60     lightcyan    lightgreen   lightyellow       magenta  midnightblue 
#113            33             7            22            36            25            18            52            61 
#orange          pink        purple           red     royalblue        salmon           tan     turquoise        yellow 
#3           114            32           427            12            30            93          1050           341 
load("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/DiffCorNetworks/TreatmentCorrelationFull.Rdat")
names(treat)
edges<-treat$edges
red<-unique(as.character(red[,2]))
edges<-edges[which(rownames(dat) %in% red),]
edges<-edges[,which(rownames(dat) %in% red)]
cor<-treat$corg
cor<-cor[which(rownames(dat) %in% red),which(rownames(dat) %in% red)]
genenames<-red
#prep for cytoscape
modcor<-c()
for (i in 1:526){
  temp<-c()
  for (j in ((i+1):526)) {
    if (edges[i,j] == 1 & abs(cor[i,j]) >= 0.2) {
      temp<-data.frame(Node1=genenames[i],Node2=genenames[j],CorCoeff=as.numeric(cor[i,j]),AbsCorCoeff=as.numeric(abs(cor[i,j])))
      modcor<-rbind(modcor,temp)
    } else  {
      next
    }
  }
}

write.table(modcor,file="treatRedModuleForCyto.txt")
