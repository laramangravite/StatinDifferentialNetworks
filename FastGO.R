##############################################################################################
##             Fast Gene Ontology enrichment analysis for modules using Fisher's exact test
##        
##     Author         :  Zhi Wang
##     Contact        :  Sage Bionetworks
##     Email          :  zhi.wang@sagebase.org
##     Date           :  03/03/2011

##------------------------------------------------------------------------------------------------------------------
## inputfname contains: gene information, expression data, module information (last column) 
## identifier_col     : unique gene identifier column
##                      =1 for  normal case (probeset)
##                      =3 for bxh (locus link number)
##                      actually, it will be ignored
##
## gene_symbolcol     : we assume the first col as probset and this column gives the real gene names
## ontologyfnlist     : GeneFisher ontology file, ONLY ACCEPT ONE FILE
## maxSignifLevel     : report only the categories with FisherExactTest Pvalue < maxSignifLevel
## outdir             : =="", put the results under the same directory as the input, otherwise under the new directory
## useAllModules	  : whether use all genes as one module
## background         : set background size (if set, it should be larger than number of genes in the input file unions ontology file)
## ReportGenes        : Whether report overlap genes. It is ~7 times slower for this calculation.

FastModuleOntologyAnalysis <- function(inputfname, ontologyfnlist,
									   identifier_col=1, gene_symbolcol=2,
									   outdir="",
									   useAllModules=F, background=0, minPopHits=5,
									   signifLevel=0.001, ctopNrows=20, ReportGenes=T
									  )
{
	## read input module file
	allMatrix <- read.delim(inputfname, sep=",", header=F)
	if(outdir=="") {
		mfname        = getFileName(inputfname)
    } else{
		mfname        = getFileNameNopath(inputfname)
		mfname        = paste(outdir, "/", mfname, sep="")
    }
	
	if(!useAllModules){ # consider modules
		geneSet <- allMatrix[, c(gene_symbolcol, ncol(allMatrix))]

    }  else{ # consider all gene in the list
		geneSet <- cbind(allMatrix[, gene_symbolcol], "all")
    }
	
	geneSet[,1] <- toupper(as.character(geneSet[,1]));	geneSet[,2] <- as.character(geneSet[,2]);
	
	## remove NA, "NA" and ""
	is.rm <- is.na(geneSet[,1]) | is.na(geneSet[,2]) | (geneSet[,1] == "") | (geneSet[,2] == "") | (geneSet[,1] == "NA") | (geneSet[,2] == "NA")
	geneSet <- geneSet[!is.rm, ]
	
	## remove redundant rows
	geneSet <- unique(geneSet)
	
	## split into list
	geneList <- split(geneSet[,1], geneSet[,2])
	
	## read ontology file
	ontoTab <- read.delim(ontologyfnlist, sep="\t", header=T, as.is=T)[,c(1,2,5)]
	
	## change ontology system to shorter names
	ontoTab[ontoTab[,1] == "GO Biological Process",1] <- "BP"
	ontoTab[ontoTab[,1] == "GO Molecular Function",1] <- "MF"
	ontoTab[ontoTab[,1] == "GO Cellular Component",1] <- "CC"
	ontoTab[ontoTab[,1] == "Panther Biological Process",1] <- "Panther BP"
	ontoTab[ontoTab[,1] == "Panther Molecular Function",1] <- "Panther MF"
	ontoTab[ontoTab[,1] == "Panther Cellular Component",1] <- "Panther CC"
	
	## combine ontology system and category
	ontoTab <- cbind(paste(ontoTab[,1], ontoTab[,2], sep='\t'), ontoTab[,3])
	
	ontoList <- list()
	
	for(i in 1:nrow(ontoTab))
	{
		genes <- toupper(unlist(strsplit(ontoTab[i,2], "; ")))
		#if(length(genes) > minPopHits)
		#{
			ontoList[[i]] <- genes
		#}
	}
	names(ontoList) <- ontoTab[,1]

	overlap <- GeneSetListOverlap(ontoList, geneList,
								  Background='set1', BackgroundSize=0,
								  ReportGenes=ReportGenes
								 )
	
	
	modNames <- names(overlap$Set2Size)
	ontNames <- names(overlap$Set1Size)
	enrichTab <- overlap$EnrichTable
	
			# tabnames <- c("Set1", "Set2", 
					  # "Overlap size", "Sampling size", "Positive size", "Background size", 
					  # "Fold enrichment", "P value",
					  # "Overlap genes"
					 # )
	if(ReportGenes)
	{
		colnames(enrichTab) <- c("System\tOntology category", "Module",
								 "Overlap size", "Module size", "Category size", "Background size", 
								 "Fold enrichment", "P value",
								 "Overlap genes"
								)
	}   else
	{
		colnames(enrichTab) <- c("System\tOntology category", "Module",
								 "Overlap size", "Module size", "Category size", "Background size", 
								 "Fold enrichment", "P value"
								)
	}
	enrichTab <- enrichTab[,c(2,4,1,3,5:ncol(enrichTab))]
	
	## filter unsignificant rows
	is.sig <- (as.numeric(enrichTab[,"P value"]) <= signifLevel & as.numeric(enrichTab[,"Overlap size"]) >= minPopHits)
	enrichTab <- enrichTab[is.sig,]

	
	ontTab <- NULL
	## output each module results to file
	for(mod in modNames)
	{
		outFile <- paste(mfname, "_Ontology_",  mod, ".xls",  sep='')
		modTab <- enrichTab[enrichTab[,"Module"] == mod, ]
		
		
		write.table(modTab, file=outFile, sep='\t', quote=F, col.names=T, row.names=F)
		
		if(nrow(modTab) > ctopNrows)
		{
			ontTab <- rbind(ontTab, modTab[1:ctopNrows,])
		}   else
		{
			ontTab <- rbind(ontTab, modTab)
		}
		
	}
	
	ontTab <- ontTab[order(-1* as.numeric(ontTab[,"Module size"]), ontTab[,"Module"], as.numeric(ontTab[,"P value"])), ]
	
	## output total result
	totalOutFile <- paste(mfname, "_Ontology.xls",  sep='')
	write.table(ontTab, file=totalOutFile, sep='\t', quote=F, col.names=T, row.names=F)	
}


##-----------------------------------------------------------------------------------------------
## Internal functions
##-----------------------------------------------------------------------------------------------

## New implementation of the test enrichment of overlaps between gene set1 and gene set2
## Input:
## 		GeneSet*:		    2 column data frame, 1st col Gene Identifiers, 2nd col Category Identifiers
##                          NOTICE: NA, "NA", "" are removed for Gene and Category Identifiers
##      Background:			which dataset used as background, usually set1 is background, and set2 is query set
##		BackgroundSize:		Specified background gene set size, only useful if Background="large",
##                     		if specify BackgroundSize, it should be larger than union(GeneSet1, GeneSet2)
##
## Return:
##		EnrichTable:		Table of enrichments, ranked by p-value
##		PVal:				P value matrix
##		Fold:               Fold enrichment matrix

GeneSetOverlap <- function(GeneSet1, GeneSet2,
                           Background="set1", BackgroundSize=0,
						   ReportGenes=T
						  )
{
	## determine input parameter 
	Background <- char.expand(Background, c("set1", "union", "intersect", "large"))
	if (length(Background) > 1 || is.na(Background)) 			stop("Background must be \"set1\", \"union\", \"intersect\", \"large\"")
	
	
	####################################
	## clean up data
	
	## convert datasets to characters
	GeneSet1[,1] <- toupper(as.character(GeneSet1[,1]));	GeneSet1[,2] <- as.character(GeneSet1[,2]);
	GeneSet2[,1] <- toupper(as.character(GeneSet2[,1]));	GeneSet2[,2] <- as.character(GeneSet2[,2]);
	
	## remove NA, "NA" and ""
	is.rm <- is.na(GeneSet1[,1]) | is.na(GeneSet1[,2]) | (GeneSet1[,1] == "") | (GeneSet1[,2] == "") | (GeneSet1[,1] == "NA") | (GeneSet1[,2] == "NA")
	GeneSet1 <- GeneSet1[!is.rm, ]
	is.rm <- is.na(GeneSet2[,1]) | is.na(GeneSet2[,2]) | (GeneSet2[,1] == "") | (GeneSet2[,2] == "") | (GeneSet2[,1] == "NA") | (GeneSet2[,2] == "NA")
	GeneSet2 <- GeneSet2[!is.rm, ]

	## remove redundant rows
	GeneSet1 <- unique(GeneSet1);			GeneSet2 <- unique(GeneSet2)
	
	###################################
	## determine background set, and do adjustments
	CatGene1 <- split(GeneSet1[,1], GeneSet1[,2]);	CatGene2 <- split(GeneSet2[,1], GeneSet2[,2])
	
	overlap <- GeneSetListOverlap(CatGene1, CatGene2,
								  Background=Background, BackgroundSize=BackgroundSize,
								  ReportGenes=ReportGenes
								 )
	return(overlap)
}

## New implementation of the test enrichment of overlaps between gene set1 and gene set2
## Input:
## 		CatGene*:		    list format of genes within each category
##                          names(list) are category names, list[[i]] are genes within category i
##                          NOTICE: NA, "NA", "" should have already been removed for Gene and Category Identifiers, upper/lower cases are unified
##      Background:			which dataset used as background, usually set1 is background, and set2 is query set
##		BackgroundSize:		Specified background gene set size, only useful if Background="large",
##                     		if specify BackgroundSize, it should be larger than union(GeneSet1, GeneSet2)
##      ReportGenes:        Whether report overlap genes. It is ~7 times slower for this calculation.  It will be improved later.
##
## Return:
##		EnrichTable:		Table of enrichments, ranked by p-value
##		EnrichMatrices:		Module-by-module enrichment matrices

GeneSetListOverlap <- function(CatGene1, CatGene2,
							   Background="set1", BackgroundSize=0,
							   ReportGenes=T
							  )
{
    ############################
	## debug use
	if(F)
	{
		g1 <- paste('g', c(1:10, 1), sep='')
		g2 <- paste('g', c(1:11), sep='')
		c1 <- paste('c', c(rep(1,3), rep(2,5), rep(3,2), 2), sep='')
		c2 <- paste('c', c(rep(2,3), rep(3,4), rep(1,3), rep(4,1)), sep='')
		GeneModuleSet1 <- GeneSet1 <- cbind(g1, c1)
		GeneModuleSet2 <- GeneSet2 <- cbind(g2, c2)
		
		Background="set1"
		BackgroundSize=0
		ReportGenes=T
		
		CatGene1 <- split(GeneSet1[,1], GeneSet1[,2]);	CatGene2 <- split(GeneSet2[,1], GeneSet2[,2])
	}
	## end of debug
	############################

	## determine input parameter 
	#st <- proc.time()

	Background <- char.expand(Background, c("set1", "union", "intersect", "large"))
	if (length(Background) > 1 || is.na(Background)) 			stop("Background must be \"set1\", \"union\", \"intersect\", \"large\"")
	
	## Unique genes & sets
	UnqGene1 <- unique(unlist(CatGene1));			UnqGene2 <- unique(unlist(CatGene2))
	CatList1 <- names(CatGene1);					CatList2 <- names(CatGene2)
	
	## create relationship matrices
	## rows are genes, columns are categories
	gsMat1 <- matrix(F, nrow=length(UnqGene1), ncol=length(CatList1), dimnames = list(UnqGene1, CatList1))
	for(cat1 in 1:length(CatGene1))
	{
		gsMat1[CatGene1[[cat1]], cat1] <- T
	}
	
	gsMat2 <- matrix(F, nrow=length(UnqGene2), ncol=length(CatList2), dimnames = list(UnqGene2, CatList2))
	for(cat2 in 1:length(CatGene2))
	{
		gsMat2[CatGene2[[cat2]], cat2] <- T
	}
	
	#print(proc.time()-st)
	#st <- proc.time()

	## determine background set
	if(Background=="set1")
	{
		totalGenes <- UnqGene1
		totalBalls <- length(totalGenes)
		
		#UnqGene1   <- UnqGene1[UnqGene1 %in% totalGenes]
		UnqGene2   <- UnqGene2[UnqGene2 %in% totalGenes]
		
		#tempMat1    <- gsMat1
		#gsMat1     <- matrix(F, nrow=totalBalls, ncol=length(CatList1), dimnames = list(totalGenes, CatList1))
		#gsMat1[UnqGene1,] <- tempMat1[UnqGene1,]
		
		tempMat2    <- gsMat2
		gsMat2     <- matrix(F, nrow=totalBalls, ncol=length(CatList2), dimnames = list(totalGenes, CatList2))
		gsMat2[UnqGene2,] <- tempMat2[UnqGene2,]
		
	} else if(Background=="union")
	{
		totalGenes <- union(UnqGene1, UnqGene2)
		totalBalls <- length(totalGenes)

		#UnqGene1   <- UnqGene1[UnqGene1 %in% totalGenes]
		#UnqGene2   <- UnqGene2[UnqGene2 %in% totalGenes]
		
		tempMat1    <- gsMat1
		gsMat1     <- matrix(F, nrow=totalBalls, ncol=length(CatList1), dimnames = list(totalGenes, CatList1))
		gsMat1[UnqGene1,] <- tempMat1[UnqGene1,]
		
		tempMat2    <- gsMat2
		gsMat2     <- matrix(F, nrow=totalBalls, ncol=length(CatList2), dimnames = list(totalGenes, CatList2))
		gsMat2[UnqGene2,] <- tempMat2[UnqGene2,]
		
	} else if(Background=="intersect")
	{
		totalGenes <- intersect(UnqGene1, UnqGene2)
		totalBalls <- length(totalGenes)

		UnqGene1   <- UnqGene1[UnqGene1 %in% totalGenes]
		UnqGene2   <- UnqGene2[UnqGene2 %in% totalGenes]
		
		tempMat1    <- gsMat1
		gsMat1     <- matrix(F, nrow=totalBalls, ncol=length(CatList1), dimnames = list(totalGenes, CatList1))
		gsMat1[UnqGene1,] <- tempMat1[UnqGene1,]
		
		tempMat2    <- gsMat2
		gsMat2     <- matrix(F, nrow=totalBalls, ncol=length(CatList2), dimnames = list(totalGenes, CatList2))
		gsMat2[UnqGene2,] <- tempMat2[UnqGene2,]
		
	} else if(Background=="large")
	{
		totalGenes <- union(UnqGene1, UnqGene2)
		totalBalls <- length(totalGenes)
		
		#UnqGene1   <- UnqGene1[UnqGene1 %in% totalGenes]
		#UnqGene2   <- UnqGene2[UnqGene2 %in% totalGenes]
		
		tempMat1    <- gsMat1
		gsMat1     <- matrix(F, nrow=totalBalls, ncol=length(CatList1), dimnames = list(totalGenes, CatList1))
		gsMat1[UnqGene1,] <- tempMat1[UnqGene1,]
		
		tempMat2    <- gsMat2
		gsMat2     <- matrix(F, nrow=totalBalls, ncol=length(CatList2), dimnames = list(totalGenes, CatList2))
		gsMat2[UnqGene2,] <- tempMat2[UnqGene2,]
		
		if(totalBalls < length(totalGenes))
		{
			warning("Defined BackgroundSize is too small, use union(GeneSet1, GeneSet2) instead!")
		}   else
		{
			totalBalls <- BackgroundSize
		}
		
	}
	
	if(totalBalls == 0)	stop("Your background set has problems!")

	#print(proc.time()-st)
	#st <- proc.time()
	
	## category size	
	CatSize1 <- colSums(gsMat1);					CatSize2 <- colSums(gsMat2)
	
	## overlap matrix
	CountMat <- t(gsMat1) %*% gsMat2
	
	## p-value matrix
	totalWhite <- matrix(CatSize1, nrow=length(CatSize1), ncol=length(CatSize2), byrow=F)
	totalBlack <- totalBalls - totalWhite
	totalDrawn <- matrix(CatSize2, nrow=length(CatSize1), ncol=length(CatSize2), byrow=T)
	PMat       <- phyper(CountMat-1, totalWhite, totalBlack, totalDrawn, lower.tail=F)
	FoldMat    <- (CountMat/totalDrawn)/(totalWhite/totalBalls)
	
	set1Mat    <- matrix(CatList1, nrow=length(CatSize1), ncol=length(CatSize2), byrow=F)
	set2Mat    <- matrix(CatList2, nrow=length(CatSize1), ncol=length(CatSize2), byrow=T)
		
	#print(proc.time()-st)
	#st <- proc.time()
	
	if(ReportGenes)
	{
		overGeneMat <- matrix("", nrow=length(CatSize1), ncol=length(CatSize2))
		for(i in 1:length(CatList1))
		{
			for(j in 1:length(CatList2))
			{
				overlapGenes <- totalGenes[gsMat1[,i] & gsMat2[,j]]
				overGeneMat[i,j] <- paste(overlapGenes, collapse=";")
			}
		}
		
		tabnames <- c("Set1", "Set2", 
					  "Overlap size", "Sampling size", "Positive size", "Background size", 
					  "Fold enrichment", "P value",
					  "Overlap genes"
					 )
		EnrichTable <- data.frame(as.character(set1Mat), as.character(set2Mat),
								  as.character(CountMat), as.character(totalDrawn), as.character(totalWhite), as.character(totalBalls),
								  as.character(FoldMat), as.character(PMat),
								  as.character(overGeneMat),
								  stringsAsFactors=F
								 )
		colnames(EnrichTable) <- tabnames
		
	}   else
	{
		tabnames <- c("Set1", "Set2", 
					  "Overlap size", "Sampling size", "Positive size", "Background size", 
					  "Fold enrichment", "P value"
					 )
		EnrichTable <- data.frame(as.character(set1Mat), as.character(set2Mat),
								  as.character(CountMat), as.character(totalDrawn), as.character(totalWhite), as.character(totalBalls),
								  as.character(FoldMat), as.character(PMat),
								  stringsAsFactors=F
								 )
		colnames(EnrichTable) <- tabnames
	}
	
	EnrichTable[is.na(EnrichTable)] <- 0
	EnrichTable <- EnrichTable[order(as.numeric(EnrichTable[,"P value"])) ,]

	#print(proc.time()-st)
		
	## return value
	return(list(Set1Size=CatSize1, Set2Size=CatSize2,
				Count=CountMat, Fold=FoldMat, Pval=PMat,
				EnrichTable=EnrichTable
			   )
		  )
}


#get the filename without extension
#
getFileExtension=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    
    if( length(splitted) >1){
        return (splitted[length(splitted)])
    } else{
        return ("")
    }
}

#get the filename without extension
getFileName=function(fullfname){
    ext=getFileExtension(fullfname)
    if(ext ==""){
        return (fullfname)
    }
    extd = paste(".", ext, sep="")
    splitted=splitString(fullfname, extd)
    
    splitted[1]
}

#get the filename without extension and path information
getFileNameNopath=function(fullfname){
    myfilename = getFileName(fullfname)
    splitted=unlist( strsplit(myfilename, "/") )
    splitted[length(splitted) ]
}

#get the filename without path information
getFileFullNameNopath=function(fullfnames){
    res = NULL
    for(each in fullfnames) {
        splitted=unlist( strsplit(each, "/") )
        res= c(res, splitted[length(splitted) ])
    }
    return (res)
}

# to split "abc|123", use sep="\\|", "abc.123" use "\\."
splitString =function(mystring, separator="; "){
    splitted = NULL
    for (each in mystring){
        if (is.na(each) | is.null(each)){
            next
        }
        a=unlist( strsplit(each, separator) )
        splitted =c(splitted, a)
    }
    #a=unlist( strsplit(mystring, separator) )
    return(splitted )
}
