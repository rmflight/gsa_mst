# simple test
require(igraph)
require(plyr)
nSamp <- 10 # n1 and n2, the number of samples, nv = nSamp * 2
nGene <- 40 # p, number of genes

grp1 <- cbind(rep(0, nGene), rep(0.1, nGene))

grp1Dat <- aaply(seq(1, nSamp), 1, function(x){rnorm(nGene, grp1[,1], grp1[,2])})

grp2 <- cbind(rep(1, nGene), rep(0.2, nGene))
grp2Dat <- aaply(seq(1, nSamp), 1, function(x){rnorm(nGene, grp2[,1], grp2[,2])})

geneProfile <- cbind(t(grp1Dat), t(grp2Dat))

rownames(geneProfile) <- paste("g", seq(1, nGene), sep=".")
colnames(geneProfile) <- c(paste("s1", seq(1, nSamp), sep="."), paste("s2", seq(1, nSamp), sep="."))

data <- geneProfile
nv1 <- nSamp

dimensions <- dim(data)
p <- dimensions[1]
nv <- dimensions[2]
number_perm <- 1000			# number of random permutations
combinations <- factorial(nv) / (factorial(nv-nv1)*factorial(nv1))
gt <- aperm(data, c(2,1)) # transpose the array?? -> yes, this is a simple transpose of the data
Wmat <- as.matrix(dist(gt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
gr <- graph.adjacency(Wmat, weighted = TRUE, mode = "undirected")

# this is the code supplied with the paper, and I assume used for analysis (see lines 32-33 of SupplementaryText.txt)
V(gr)[c(0:(nv1-1))]$color <- "green"
V(gr)[c(nv1:(nv-1))]$color <- "red"

mst <- minimum.spanning.tree(gr)



plot(gr)
par(ask=T)
plot(mst)

# this is how I think it was supposed to be done
gr2 <- gr
V(gr2)[c(1:(nv1))]$color <- "green"
V(gr2)[c(nv1+1:(nv))]$color <- "red"

mst2 <- minimum.spanning.tree(gr2)

plot(gr2)
plot(mst2)


# This script was originally sent to Y. Rahmatallah of University of Arkansas Medical School on Dec 6, 2012.