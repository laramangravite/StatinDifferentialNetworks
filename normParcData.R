
normParcData<-function(data,covariates){
############CAP normalization
  remove=model.matrix(~covariates$ExposureBatch+covariates$Array+covariates$Cell.Count
                      +covariates$Gender+covariates$Age+covariates$BMI+covariates$Smoking.Status)
  X<-remove
  b1 <- solve(t(X) %*% X) %*% t(X) %*% (data)
  res1 <- t(data)-t(X %*% b1)
  res2=res1+t(b1)[,1] #add intercept back
return(res2)
}