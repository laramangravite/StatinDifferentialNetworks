##############################
##WGCNA Modules:  Connectivity
##############################
setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/ControlExpnModules")
source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")

cpc<-read.table(file="ControlWGCNAModulePCs.txt",header=T)
cpc<-t(cpc)
cm<-build_correlation_graph(gex1=cpc,gexnames=colnames(cpc))
cor<-cm$corg
edges<-cm$edges
modcor<-build_cytoscape_input(cor,edges,genenames=colnames(cor))
save(cm,file="ControlModuleConnectivity.Rdat")
write.table(modcor,file="ControlModConnectCytoscape.txt",row.names=F)

setwd("/Volumes/lmangrav/PARC/CoexpNetworks/CoexpnSept2012/TreatmentExpnModules")
source("/Volumes/lmangrav/PARC/CoexpNetworks/Scripts/DifferentialCorrelationAnalysisFunctions.R")

tpc<-read.table(file="TreatmentWGCNAModulePCs.txt",header=T)
tpc<-t(tpc)
tm<-build_correlation_graph(gex1=tpc,gexnames=colnames(tpc))
cor<-tm$corg
edges<-tm$edges
tmodcor<-build_cytoscape_input(cor,edges,genenames=colnames(cor))
tmodcor[,1]<-as.character(tmodcor[,1])
tmodcor[,2]<-as.character(tmodcor[,2])
#tmodcor<-tmodcor[-(tmodcor[,2]=="MEgrey"),]
tmodcor<-tmodcor[-132,]
save(tm,file="TreatmentModuleConnectivity.Rdat")
write.table(tmodcor,file="TreatmentModConnectCytoscape.txt",row.names=F)

