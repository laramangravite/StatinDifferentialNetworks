# assumes pi0=1
# does not take into account ties. use unique or p<=x for that (but
#will take longer).
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
build_correlation_graph_fdr <- function(gex, gexnames)
{
  n = nrow(gex) #number of samples
  p = ncol(gex) #number of genes
  betahat = cor((gex)) #gene by gene correlation matrix
  se.betahat = sqrt((1-betahat^2)/n)  
  pval = 2*pnorm(-abs((betahat/se.betahat)[se.betahat>0])) #vector of length p*(p-1) with pvalues for all correlations except self-correlations (gene1 x gene1)
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
  return(list(edges=qval, corg=betahat))
}

build_correlation_graph_differential_2not1 <- function(gex1, gex2, gexnames)
{
  n = nrow(gex1)
  p = ncol(gex1)
  gex = c()
  #Build quantile normalized residuals of control gene expression data with treated expression data removed on a gene-by-gene basis
  for(i in 1:n) {
    gex = cbind(gex, qqnorm(lm(gex2[i,]~gex1[i,])$residual,plot.it=FALSE)$x)
  }
  return(build_correlation_graph_fdr(gex, gexnames))
}

build_correlation_graph_differential_1not2 <- function(gex1, gex2, gexnames)
{
  n = nrow(gex1) #number of genes
  p = ncol(gex1) #number of samples
  gex = c()
  #Build quantile normalized residuals of treated gene expression data with control expression data removed on a gene-by-gene basis
  for(i in 1:n) {
    gex = cbind(gex, qqnorm(lm(gex1[i,]~gex2[i,])$residual,plot.it=FALSE)$x)
  }#returns a p x n (sample x gene) matrix.
  return(build_correlation_graph_fdr(gex, gexnames))
}


build_correlation_graph_differential_opposite <- function(gex1, gex2, gexnames)
{
  n = nrow(gex1) #number of genes
  p = ncol(gex1) #number of samples
  gex = c()
  #Build quantile normalized residuals of differential gene expression data with averaged expression data removed on a gene-by-gene basis
  for(i in 1:n) {
    gex = cbind(gex,
                qqnorm(lm((gex2[i,]-gex1[i,])~(gex2[i,]+gex1[i,]))$residual,plot.it=FALSE)$x)
  } #returns a p x n (sample x gene) matrix.
  return(build_correlation_graph_fdr(gex, gexnames))
}


build_correlation_graph <- function(gex1, gexnames)
{
  n = nrow(gex1) #number of genes
  p = ncol(gex1) #number of samples
  #Quantile normalize  
  gex<-apply(gex1,1,function(x) {return(qqnorm(x,plot.it=F)$x)})
  return(build_correlation_graph_fdr(gex, gexnames))
}

###########Transform output for cytoscape
build_cytoscape_input<-function(cor, edges,genenames) {
  p<-nrow(cor)
  modcor<-c()
  for (i in 1:(p-1)){
    temp<-c()
    for (j in ((i+1):p)) {
      if (edges[i,j] == 1) {
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
