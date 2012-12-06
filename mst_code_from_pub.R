findRuns <- function (mst,mst2,nv,nv1){
	mstWM <- get.adjacency(mst, type = "lower", attr="weight")
	edgeind <- which(mstWM != 0, arr.ind = TRUE, useNames = FALSE) 
	run_count <- 1 + sum(V(mst2)[edgeind[,1]-1]$color != V(mst2)[edgeind[,2]-1]$color) # run_count = no. of deleted edges + 1
	run_count
}

#########################################################################
#									#
#	This function performs the multivariate Wald-Wolfowitz test	#
#									#
#########################################################################
# 'data' is a p-by-nv matrix that supply the gene expression profiles to the function
# 'nv' is the total number of samples in both groups
# 'nv1' is the number of samples of the first group. Obviously (nv-nv1) is the number of samples of the second group
# 'p' is the number of genes in the gene set

# The syntax for calling this function is:
# p_value <- MVWWtest(data,nv1,print_decision=TRUE)   (print_decision can be either TRUE or FALSE)

MVWWtest <- function (data,nv1,print_decision = FALSE) {
	library(igraph)
	library(combinat)
	dimensions <- dim(data)
	p <- dimensions[1]
	nv <- dimensions[2]
	number_perm <- 1000			# number of random permutations
	combinations <- factorial(nv) / (factorial(nv-nv1)*factorial(nv1))
	gt <- aperm(data, c(2,1))
	Wmat <- as.matrix(dist(gt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
	gr <- graph.adjacency(Wmat, weighted = TRUE, mode = "undirected")
	V(gr)[c(0:(nv1-1))]$color <- "green"
	V(gr)[c(nv1:(nv-1))]$color <- "red"
	mst <- minimum.spanning.tree(gr)
	permutations <- matrix(V(mst)$color, nv, 1)
	domain <- V(mst)$color
	
	#####################################################################################
	# if combinations > 3000, do an approximate permutation test (randomly picked       #
	# combinations withno replica). If nv1=nv/2, the condition is satisfied for nv>=14. #
	# if combinations < 3000, do full permutation test (all possible combinations with  #
	# no replica). Random search for 1000 distictive samples out of a total number of   #
	# samples which is little bit > 1000 will result in a prolonged loop execution	    #
	#####################################################################################
	if (combinations > 3000) {
		runs <- array(0,c(1,number_perm))	
		for (itr in 1:number_perm) 	{ 
			randperm <- sample(domain, replace = FALSE)
			mst2 <- mst
			while (sum(colSums(permutations == randperm) == nv) > 0)
			{
				randperm <- sample(domain, replace = FALSE)
			}
			permutations <- cbind(permutations,randperm)					
			V(mst2)$color <- randperm
			runs[itr] <- findRuns(mst,mst2,nv,nv1) 
		}
	} else {
		runs <- array(0,c(1,combinations))	# a vector to save results
		co <- combn(c(1:nv), nv1)
		for (itr in 1:combinations) 	{ 
			mst2 <- mst
			V(mst2)$color <- "red"
			V(mst2)[co[,itr]-1]$color <- "green"
			runs[itr] <- findRuns(mst,mst2,nv,nv1) 
		}
	}
	
	alpha <- 0.05
	runs_a <- findRuns(mst,mst,nv,nv1)
	p_value <- (sum(runs < runs_a) + 1) / (length(runs) + 1)
	if (print_decision == TRUE) {
		if (p_value < alpha) print(paste("p_value = ", p_value, "      Reject Ho"), quote=FALSE) else print(paste("p_value = ", p_value, "      Fail to reject Ho"), quote=FALSE)
	}
	p_value
}



findTestStatRKS <- function (radial_ranking,randperm,nv,nv1){
	ri <- 0
	si <- 0
	di <- array(0, c(1,nv))
	#radial_ranking <- radial_ranking + 1
	for (i in 1:nv) {
		ri <- sum(randperm[radial_ranking[1:i]] == "green")
		si <- sum(randperm[radial_ranking[1:i]] == "red")
		di[i] <- (ri/nv1) - (si/(nv-nv1))
	}
	D <- sqrt((nv1 * (nv-nv1)) / (nv1 + (nv-nv1))) * max(abs(di))
	D
}

#################################################################################
#										#
#	This function performs the multivariate radial Kolmogorov-Smirnov test	#
#										#
#################################################################################
# 'data' is a p-by-nv matrix that supply the gene expression profiles to the function
# 'nv' is the total number of samples in both groups
# 'nv1' is the number of samples of the first group. Obviously (nv-nv1) is the number of samples of the second group
# 'p' is the number of genes in the gene set

# The syntax for calling this function is:
# p_value <- MVKStest(data,nv1,print_decision=TRUE)	(print_decision can be either TRUE or FALSE)

MVradialKStest <- function (data,nv1,print_decision = FALSE) {
	library(igraph)
	library(combinat)
	dimensions <- dim(data)
	p <- dimensions[1]
	nv <- dimensions[2]
	number_perm <- 1000
	combinations <- factorial(nv) / (factorial(nv-nv1)*factorial(nv1))
	
	gt <- aperm(data, c(2,1))
	Wmat <- as.matrix(dist(gt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
	gr <- graph.adjacency(Wmat, weighted = TRUE, mode = "undirected")
	V(gr)[c(0:(nv1-1))]$color <- "green"
	V(gr)[c(nv1:(nv-1))]$color <- "red"
	mst <- minimum.spanning.tree(gr)
	
	sp <- apply(shortest.paths(mst), 1, max)
	radius <- min(sp)
	center <- which(sp == radius) 
	if (length(center)>1) center <- center[1]
	ranktree <- sort(shortest.paths(mst)[, center], decreasing = FALSE, index.return = TRUE)
	radial_ranking <- ranktree$ix
	permutations <- matrix(V(mst)$color, nv, 1)
	domain <- V(mst)$color
	
	#####################################################################################
	# if combinations > 3000, do an approximate permutation test (randomly picked       #
	# combinations withno replica). If nv1=nv/2, the condition is satisfied for nv>=14. #
	# if combinations < 3000, do full permutation test (all possible combinations with  #
	# no replica). Random search for 1000 distictive samples out of a total number of   #
	# samples which is little bit > 1000 will result in a prolonged loop execution	    #
	#####################################################################################
	if (combinations > 3000) {
		D <- array(0,c(1,number_perm))
		for (itr in 1:number_perm) 	{ 
			randperm <- sample(domain, replace = FALSE)
			while (sum(colSums(permutations == randperm) == nv) > 0)
			{
				randperm <- sample(domain, replace = FALSE)
			}
			permutations <- cbind(permutations,randperm)					
			D[itr] <- findTestStatRKS(radial_ranking,randperm,nv,nv1) 
		}
	} else {
		D <- array(0,c(1,combinations))
		co <- combn(c(1:nv), nv1)
		for (itr in 1:combinations) 	{ 
			mst2 <- mst
			V(mst2)$color <- "red"
			V(mst2)[co[,itr]-1]$color <- "green"
			randperm <- V(mst2)$color
			D[itr] <- findTestStatRKS(radial_ranking,randperm,nv,nv1) 
		}
	}
	
	alpha <- 0.05
	D_a <- findTestStatRKS(radial_ranking,domain,nv,nv1)
	p_value <- (sum(D > D_a) + 1) / (length(D) + 1)
	if (print_decision == TRUE) 	{
		if (p_value < alpha) print(paste("p_value = ", p_value, "      Reject Ho"), quote=FALSE) else print(paste("p_value = ", p_value, "      Fail to reject Ho"), quote=FALSE)
	}
	p_value
}



findTestStatKS <- function (KSranking,randperm,nv,nv1){
	ri <- 0
	si <- 0
	di <- array(0, c(1,nv))
	KSranking <- KSranking + 1
	for (i in 1:nv) {
		ri <- sum(randperm[KSranking[1:i]] == "green")
		si <- sum(randperm[KSranking[1:i]] == "red")
		di[i] <- (ri/nv1) - (si/(nv-nv1))
	}
	D <- sqrt((nv1 * (nv-nv1)) / (nv1 + (nv-nv1))) * max(abs(di))
	D
}


HDP.ranking <- function(mst,nv){
	rr <- farthest.nodes(mst, directed = FALSE, unconnected = TRUE)
	root <- floor(rr[1])
	terminal_nodes <- which(degree(mst) == 1)
	ltn <- length(terminal_nodes) - 1
	tn <- terminal_nodes - 1
	tn <- tn[tn != root]
	sp <- get.shortest.paths(mst, root, to = tn)
	
	path_len <- shortest.paths(mst)
	break_ties <- path_len[root+1, tn+1] / max(path_len)
	depth <- array(0, c(1,ltn))
	KSranks <- root
	for (k in 1:ltn)	{
		depth[k] <- length(sp[[k]])
	}
	md <- max(depth) 
	adjusted_depth <- depth + break_ties
	col_nodes <- array(0, c(1,ltn))
	alphabets <- rep("",ltn)
	for (col in seq(1,md,by=1)) 	{
		for (row in seq(1,ltn,by=1))  	{
			col_nodes[row] <- sp[[row]][col]
		}
		fcn <- factor(col_nodes)
		collevels <- levels(fcn)
		llev <- length(collevels)
		if (llev > 1) 	{
			mpg <- tapply(adjusted_depth,fcn,max)
			sortmpg <- sort(mpg, decreasing = FALSE, index.return = TRUE)
			smpg <- sortmpg$ix
			sorted_levels <- collevels[smpg]
			for (lind in seq(1,length(smpg),by=1)) 	{
				alphabets[which(col_nodes==sorted_levels[lind])]<- paste(alphabets[which(col_nodes==sorted_levels[lind])], letters[lind], sep="")
			}
		}
	}
	newranks <- sort(alphabets, decreasing = FALSE, index.return = TRUE)
	spm <- as.matrix(sp)
	sp_new <- spm[newranks$ix ,]
	sp_new <- as.matrix(sp_new)
	for (k in 1:ltn)	{
		len <- length(sp_new[[k]])
		for (u in 1:len)	{
			if (sum(KSranks == sp_new[[k]][u]) == 0){
				KSranks <- c(KSranks,sp_new[[k]][u]) }
		}
	}
	KSranks
}

#########################################################################
#									#
#	This function performs the multivariate Kolmogorov-Smirnov test	#
#									#
#########################################################################
# 'data' is a p-by-nv matrix that supply the gene expression profiles to the function
# 'nv' is the total number of samples in both groups
# 'nv1' is the number of samples of the first group. Obviously (nv-nv1) is the number of samples of the second group
# 'p' is the number of genes in the gene set

# The syntax for calling this function is:
# p_value <- MVKStest(data,nv1,print_decision=TRUE)	(print_decision can be either TRUE or FALSE)

MVKSHDPtest <- function (data,nv1,print_decision = FALSE) {
	library(igraph)
	library(combinat)
	dimensions <- dim(data)
	p <- dimensions[1]
	nv <- dimensions[2]
	number_perm <- 1000
	combinations <- factorial(nv) / (factorial(nv-nv1)*factorial(nv1))
	
	gt <- aperm(data, c(2,1))
	Wmat <- as.matrix(dist(gt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
	gr <- graph.adjacency(Wmat, weighted = TRUE, mode = "undirected")
	V(gr)[c(0:(nv1-1))]$color <- "green"
	V(gr)[c(nv1:(nv-1))]$color <- "red"
	mst <- minimum.spanning.tree(gr)
	KSranking <- HDP.ranking(mst,nv)
	
	permutations <- matrix(V(mst)$color, nv, 1)
	domain <- V(mst)$color
	
	#####################################################################################
	# if combinations > 3000, do an approximate permutation test (randomly picked       #
	# combinations withno replica). If nv1=nv/2, the condition is satisfied for nv>=14. #
	# if combinations < 3000, do full permutation test (all possible combinations with  #
	# no replica). Random search for 1000 distictive samples out of a total number of   #
	# samples which is little bit > 1000 will result in a prolonged loop execution	    #
	#####################################################################################
	if (combinations > 3000) {
		D <- array(0,c(1,number_perm))
		for (itr in 1:number_perm) 	{ 
			randperm <- sample(domain, replace = FALSE)
			while (sum(colSums(permutations == randperm) == nv) > 0)
			{
				randperm <- sample(domain, replace = FALSE)
			}
			permutations <- cbind(permutations,randperm)					
			D[itr] <- findTestStatKS(KSranking,randperm,nv,nv1) 
		}
	} else {
		D <- array(0,c(1,combinations))
		co <- combn(c(1:nv), nv1)
		for (itr in 1:combinations) 	{ 
			mst2 <- mst
			V(mst2)$color <- "red"
			V(mst2)[co[,itr]-1]$color <- "green"
			randperm <- V(mst2)$color
			D[itr] <- findTestStatKS(KSranking,randperm,nv,nv1) 
		}
	}
	
	alpha <- 0.05
	D_a <- findTestStatKS(KSranking,domain,nv,nv1)
	p_value <- (sum(D > D_a) + 1) / (length(D) + 1)
	if (print_decision == TRUE) 	{
		if (p_value < alpha) print(paste("p_value = ", p_value, "      Reject Ho"), quote=FALSE) else print(paste("p_value = ", p_value, "      Fail to reject Ho"), quote=FALSE)
	}
	p_value
}