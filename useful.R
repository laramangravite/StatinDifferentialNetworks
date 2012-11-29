library(corpcor)
fs <- function(dat) {
u <- fast.svd(sweep(dat,1,rowMeans(dat)))
u$d <- u$d^2 / sum(u$d^2)
u
}
fast.svd <- function( x, nu = min(n, p), nv = min(n, p)) {
    x <- as.matrix(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
 
    if(  p <= n )
      return( svd( x, nu, nv ) )
    else
      {
      s <- svd( t(x), nu=nv, nv=nu)
      retval <- list()
      retval$d <- s$d
      retval$u <- s$v
      retval$v <- s$u
      return(retval)
    }
  }




cor.pvalue <- function(X,method="pearson", use="complete" ) { 
dfr = nrow(X) - 2 
R <- cor(X,method=method, use= use) 
above <- row(R) < col(R) 
r2 <- R[above]^2 
Fstat <- r2 * dfr / (1 - r2) 
R[above] <- 1 - pf(Fstat, 1, dfr) 
R 
} 

ls.objects=function(pos = 1, pattern, order.by) 
{ 
    napply <- function(names, fn) sapply(names, function(x) 
fn(get(x, 
        pos = pos))) 
    names <- ls(pos = pos, pattern = pattern) 
    obj.class <- napply(names, function(x) as.character(class(x))[1]) 
    obj.mode <- napply(names, mode) 
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class) 
    obj.size <- napply(names, object.size) 
    obj.dim <- t(napply(names, function(x) 
as.numeric(dim(x))[1:2])) 
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function") 
    obj.dim[vec, 1] <- napply(names, length)[vec] 
    out <- data.frame(obj.type, obj.size, obj.dim) 
    names(out) <- c("Type", "Size", "Rows", "Columns") 
    if (!missing(order.by)) 
        out <- out[order(out[[order.by]]), ] 
    out 
} 

pcor.mat=function(x,y,z,method="p",na.rm=T){

        x <- c(x)
        y <- c(y)
        z <- as.data.frame(z)

        if(dim(z)[2] == 0){
                stop("There should be given data\n")
        }

        data <- data.frame(x,y,z)

        if(na.rm == T){
                data = na.omit(data)
        }

        xdata <- na.omit(data.frame(data[,c(1,2)]))
        Sxx <- cov(xdata,xdata,m=method)

        xzdata <- na.omit(data)
        xdata <- data.frame(xzdata[,c(1,2)])
        zdata <- data.frame(xzdata[,-c(1,2)])
        Sxz <- cov(xdata,zdata,m=method)

        zdata <- na.omit(data.frame(data[,-c(1,2)]))
        Szz <- cov(zdata,zdata,m=method)

        # is Szz positive definite?
        zz.ev <- eigen(Szz)$values
        if(min(zz.ev)[1]<0){
                stop("\'Szz\' is not positive definite!\n")
        }

        # partial correlation
        Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
        
        rxx.z <- cov2cor(Sxx.z)[1,2]

        rxx.z
}

pcor.test=function(x,y,z,use="mat",method="p",na.rm=T){
	# The partial correlation coefficient between x and y given z
	#
	# pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
	#
	# x and y should be vectors
	#
	# z can be either a vector or a matrix
	#
	# use: There are two methods to calculate the partial correlation coefficient.
	#	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
	#	 Default is "mat".
	#
	# method: There are three ways to calculate the correlation coefficient, 
	#	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
	# 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
	#	    Default is "p".
	#
	# na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
	#        If not, the missing samples will be removed just when the correlation coefficient is calculated.
	#	   However, the number of samples for the p-value is the number of samples after removing 
	#	   all the missing samples from the whole dataset.
	#	   Default is "T".

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(use == "mat"){
		p.use <- "Var-Cov matrix"
		pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
	}else if(use == "rec"){
		p.use <- "Recursive formula"
		pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
	}else{
		stop("\'use\' should be either \"rec\" or \"mat\"!\n")
	}

	# print the method
	if(gregexpr("p",method)[[1]][1] == 1){
		p.method <- "Pearson"
	}else if(gregexpr("s",method)[[1]][1] == 1){
		p.method <- "Spearman"
	}else if(gregexpr("k",method)[[1]][1] == 1){
		p.method <- "Kendall"
	}else{
		stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
	}

	# sample number
	n <- dim(na.omit(data.frame(x,y,z)))[1]
	
	# given variables' number
	gn <- dim(z)[2]

	# p-value
	if(p.method == "Kendall"){
		statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  		p.value <- 2*pnorm(-abs(statistic))
	}

	data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}

pcor.rec=function(x,y,z,method="p",na.rm=T){
	# 

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	# recursive formula
	if(dim(z)[2] == 1){
		tdata <- na.omit(data.frame(data[,1],data[,2]))
		rxy <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]))
		rxz <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]))
		ryz <- cor(tdata[,1],tdata[,2],m=method)

		rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
		
		return(rxy.z)
	}else{
		x <- c(data[,1])
		y <- c(data[,2])
		z0 <- c(data[,3])
		zc <- as.data.frame(data[,-c(1,2,3)])

		rxy.zc <- pcor.rec(x,y,zc,method=method,na.rm=na.rm)
		rxz0.zc <- pcor.rec(x,z0,zc,method=method,na.rm=na.rm)
		ryz0.zc <- pcor.rec(y,z0,zc,method=method,na.rm=na.rm)
		
		rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
		return(rxy.z)
	}			
}
cor.prob <- function(X, dfr = nrow(X) - 2) {
	 R <- cor(X)
	 above <- row(R) < col(R)
	 r2 <- R[above]^2
	 Fstat <- r2 * dfr / (1 - r2)
	 R[above] <- 1 - pf(Fstat, 1, dfr)
	 R
}