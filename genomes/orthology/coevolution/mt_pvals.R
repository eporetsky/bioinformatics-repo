#Loading arguments
args <- commandArgs(TRUE)
if(length(args) != 6){
	cat("Usage: Rscript pvalMatrix.R <taxidsFile> <distancesDir> <groups> <probability> <iterations> <mtResults>\n")
	cat("\t<taxidsFile>: Text file containing a list of taxids of interest separated by \\n\n")
	cat("\t<distancesDir>: Directory containing a file of distances for each protein. The files are in format taxid1\\ttaxid2\\tdistance\\n\n")
	cat("\t<groups>: Number of discretized groups that the algorithm will create in order to calculate the p-values\n")
	cat("\t<probability>: Probability to switch a pair of branches between 2 trees\n")
	cat("\t<iterations>: Number of switch iterations for each group\n")
	cat("\t<mtResults>: path to output file\n")
	q()
}

#Saving arguments
taxidListFile <- args[1]
dir <- args[2]
groups <- as.numeric(args[3])
probability <- as.numeric(args[4])
iterations <- as.numeric(args[5])
outputFile <- args[6]

#Constants
debug <- TRUE	#prints help on progress
printDistances <- TRUE	#prints compacted distances in file
decimals <- 6	#number of decimal points to be printed in p-value matrix
minorgs <- 2	#minorgs in the group creation

#settings to avoid scientific notations 
options("scipen"=100, "digits"=4)


#####################
#	FUNCTIONS
#####################

# reads a file containing distances between leaves in a tree (taxid1\ttaxid2\tdistance\n)
# taxidList can be provided to subset the distances to a set of taxids
# It returns the distances in a matrix
readDistances <- function(file, taxidList=NA){
	colClasses <- c("character", "character", "numeric")
	colnames <- c("tax1", "tax2", "distance")
	#Read distances table
	tabAll <- read.table(file, header = FALSE, colClasses = colClasses, comment.char="", col.names=colnames)
	if(! is.na(taxidList[1])){
		tabAll <- tabAll[tabAll$tax1 %in% taxidList & tabAll$tax2 %in% taxidList, ]
	}
	if(nrow(tabAll) > 0){
		#Create mirror values
		inverted <- tabAll[ ,c(2,1,3)]
		names(inverted) <- colnames
		tabAll <- rbind(tabAll, inverted)
		matrix <- xtabs(distance~tax1+tax2, data=tabAll)
		if(debug){
			cat(file, "readed\n",sep=" ")
		}
		return(matrix)
	}else{
		return(FALSE)
	}
}

# Reads a full directory to batch readDistances
# taxidList can be provided to subset by organisms
# It returns a list of matrixes
readFolder <- function(folder, taxidList){
	files <- list.files(folder, pattern="*.dist")
	filenames <- unlist(strsplit(files, '.dist'))
	matrixes <- list()
	for(i in 1:length(files)){
		if(is.matrix(thismatrix <- readDistances(paste(folder,files[i],sep=""), taxidList))){
			matrixes[[filenames[i]]] <- thismatrix
		}
	}
	return(matrixes)
}

#Get vector of common organisms between 2 matrixes 
getCommonOrgs <- function(matrix1, matrix2){
	common <- intersect(colnames(matrix1), colnames(matrix2))
	return(common)
}

#Gets mirrortree correlation between 2 matrixes
getMtCorr <- function(matrix1, matrix2, common){
	if(length(common) > 2){
		m1 <- matrix1[common, common]
		m2 <- matrix2[common, common]
		v1 <- m1[upper.tri(m1)]
		v2 <- m2[upper.tri(m2)]
		if((sd(v1) > 0) & (sd(v2) >0)){
			c <- cor(v1,v2)
		}else{
			c <- NA
		}
	}else{
		c <- NA
	}
	return(c)
}

# This function will divide all the possible pairs of matrixes
# based on a discretize number of organisms in common
# It should take some time with many organisms
getPairsByCommon <- function(matrixes, taxidList, groups, minorgs){
	if(debug){
		cat("Getting pairs discretized by common organisms...\n")
	}
	#phylogenetic profiles of this sequences-organisms
	phyprofiles <- sapply(matrixes, function(x) taxidList %in% rownames(x))
	row.names(phyprofiles) <- taxidList
	# pairwise orgs in commoon
	orgsInCommon <- apply(phyprofiles, 2, function(x) apply(phyprofiles, 2, function(y) length(which(x & y))))
	#discretize levels
	if(groups > (max(orgsInCommon))){
		groups <- max(orgsInCommon)
	}
	levels <- cut(1:max(orgsInCommon), breaks=unique(c(0,unique(sapply(seq(log10(minorgs), log10(max(orgsInCommon)), (log10(max(orgsInCommon)) / groups)), function(x) round(10^x))), max(orgsInCommon))))
	orgsInCommonNoZero <- orgsInCommon
	orgsInCommonNoZero[orgsInCommonNoZero == 0] <- 1
	discretized <- apply(orgsInCommonNoZero, 1, function(x) levels[x])
	row.names(discretized) <- row.names(orgsInCommonNoZero)
	#only upper triangle without diagonal
	discretized[lower.tri(discretized)] <- NA
	diag(discretized) <- NA
	pairsByCommon <- sapply(unique(as.character(levels)), function(x) which(discretized == x, arr.ind=TRUE, useNames=FALSE))
	# sapply(unique(as.character(levels)), function(x) length(which(discretized == x)))
	return(pairsByCommon)
}

# This function switch distances in a couple of matrixes with a specified "probability"
# To switch the distances the matrixes are standarized before the switch and de-standarized after the switch
# using the previous mean and sd.
# WARNING: All the changed matrixes are saved in a list outside the function
mix <- function(namex, namey, listofMatrixes, probability){
	matx <- listofMatrixes[[namex]]
	maty <- listofMatrixes[[namey]]
	#gets orgs in common
	common <- getCommonOrgs(matx, maty)
	#get common submatrixes
	submatx <- matx[common, common]
	submaty <- maty[common, common]
	#means and sds
	submatx.mean <- mean(submatx[upper.tri(submatx)])
	submaty.mean <- mean(submaty[upper.tri(submaty)])
	submatx.sd <- sd(submatx[upper.tri(submatx)])
	submaty.sd <- sd(submaty[upper.tri(submaty)])
	#standarized matrixes
	if(submatx.sd != 0){
		submatx.std <- as.matrix((as.dist(submatx) - submatx.mean) / submatx.sd)
	}else{
		submatx.std <- submatx
	}
	if(submaty.sd != 0){
		submaty.std <- as.matrix((as.dist(submaty) - submaty.mean) / submaty.sd)
	}else{
		submaty.std <- submaty
	}
	#get organisms col/row to switch and Switch
	toSwitch <- replicate(length(common), sample(c(TRUE, FALSE),1,replace=TRUE, prob=c(probability, 1-probability)))
	cols.x <- submatx.std[ ,toSwitch]
	rows.x <- submatx.std[toSwitch, ]
	submatx.std[ ,toSwitch] <- submaty.std[ ,toSwitch]
	submatx.std[toSwitch,] <- submaty.std[toSwitch,]
	submaty.std[ ,toSwitch] <- cols.x
	submaty.std[toSwitch,] <- rows.x
	#de-standarize
	if(submatx.sd != 0){
		submatx.destd <- as.matrix((as.dist(submatx.std) * submatx.sd) + submatx.mean)
	}else{
		submatx.destd <- submatx.std
	}
	if(submaty.sd != 0){
		submaty.destd <- as.matrix((as.dist(submaty.std) * submaty.sd) + submaty.mean)
	}else{
		submaty.destd <- submaty.std			
	}
	#reconstruct
	matx[common, common] <- submatx.destd
	maty[common, common] <- submaty.destd
	#return
	thisMatrixes[[namex]] <<- matx
	thisMatrixes[[namey]] <<- maty
	return()
}

# It calculates the distribution of correlations based on permuted distances matrixes
# creates a list with correlations per categorie
getPvalList <- function(pairsByCommon,iterations, probability){
	results <- list()
	if(debug){
		cat("Getting p-values...\n")
	}
	for(i in c(2:length(names(pairsByCommon)))){
		if(debug){
			cat("p-values - ",(i-1)," out of ",length(names(pairsByCommon))," - ",names(pairsByCommon)[i],"\n")
		}
		categorie <- names(pairsByCommon)[i]
		#pairs between matrixes in this categorie
		thisPairs <- pairsByCommon[[categorie]]
		#list of submatrixes
		thisMatrixes <<- matrixes[unique(as.vector(thisPairs))]
		#list of names of submatrixes
		allNames <- names(matrixes)
	
		#random iterations over our set of matrixes in this categorie
		randomOrder <- sample(nrow(pairsByCommon[[categorie]]), iterations, replace=TRUE)

                if(nrow(thisPairs)>0){
                    #we create n iterations of switched matrixes
                    out <- apply(thisPairs[randomOrder, ],1,function(x) mix(allNames[x[1]],allNames[x[2]], thisMatrixes, probability))
                    if(debug){
                            cat(i-1,"of",length(names(pairsByCommon))-1,"-",iterations, " permutations in class ", categorie,"\n")
                    }
                    #we calculate correlations of mixed
                    corrs <- apply(thisPairs[randomOrder, ],1,function(x) {getMtCorr(thisMatrixes[[allNames[x[1]]]],thisMatrixes[[allNames[x[2]]]],getCommonOrgs(thisMatrixes[[allNames[x[1]]]], thisMatrixes[[allNames[x[2]]]]))})
                    if(debug){
                            cat(i-1,"of",length(names(pairsByCommon))-1,"- correlations calculated in class ", categorie,"\n")
                    }
                  }else{
                    if(debug){
                            cat("Not enough data to run permutations in class ", categorie,"\n")
                    }
                    corrs <- NA
                }
		
		results[[categorie]] <- corrs		
		gc()
	}
	return(results)
}

# It gets the distance between a pair of matrixes
getDistance <- function(matrix, combination){
	if(sum(combination %in% row.names(matrix)) == 2){
		result <- matrix[combination[1], combination[2]]
	}else{
		result <- NA
	}
	return(result)
}

mirrorTree <- function(matrixes, pairsByCommon, pvalDistribution){
	if(debug){
		cat("Running mirrortree...\n")
	}
	result <- data.frame(p1 = character(0), p2 = character(0), n=numeric(0), r=numeric(0), p.value=numeric(0))
	allNames <- names(matrixes)
	for(i in 1:length(names(pairsByCommon))){
		categorie <- names(pairsByCommon)[i]
		thisPairs <- pairsByCommon[[categorie]]
		if(length(allNames[thisPairs[ ,1]]) > 0){
			if(debug){
				cat(i, "out of ",length(names(pairsByCommon)),"...\n")
			}			
			thisResult <- data.frame(p1=as.character(allNames[thisPairs[ ,1]]), p2=as.character(allNames[thisPairs[ ,2]]))
			thisResult$commons <- apply(thisResult, 1, function(x) as.list(getCommonOrgs(matrixes[[as.character(x[[1]])]], matrixes[[as.character(x[[2]])]])))
			thisResult$n <- sapply(thisResult$commons, length)
			#cases [0-2) in common
			if(i == 1){
				thisResult$r <- rep(NA, length(thisResult$n))
				thisResult$p.value <- rep(1, length(thisResult$n))
			}else{
				thisResult$r <- apply(thisResult, 1, function(x) getMtCorr(matrixes[[as.character(x[[1]])]],matrixes[[as.character(x[[2]])]],as.vector(as.character(x[[3]]))))
				thisPvals <- pvalDistribution[[categorie]]
				thisResult$p.value <- sapply(thisResult$r, function(x) if(!is.na(x)){return(length(which(thisPvals >= x)) / length(thisPvals))}else{return(1)})
			}
			result <- rbind(result, thisResult[ ,c("p1","p2","n","r","p.value")])
		}
	}
	return(result)
}


#####################
# MAIN
#####################

#List of taxids we are interested in
taxidList <- unique(as.character(read.table(taxidListFile)$V1))

#list of matrixes in the folder
matrixes <- readFolder(dir, taxidList)

#pairs of organisms
pairsByCommon <- getPairsByCommon(matrixes, taxidList, groups, minorgs)

#Calculating p-value matrix
thisMatrixes <- list()
pvalDistribution <- getPvalList(pairsByCommon,iterations, probability)
rm(thisMatrixes)

mtResults <- mirrorTree(matrixes, pairsByCommon, pvalDistribution)
write.table(mtResults, file=outputFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

q()
