# Understanding GSA using MSTs


require(plyr)
# need genes, samples, and classes

# simplest case, all genes are diff expressed in the gene set

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

# from the function directly

dimensions <- dim(data)
p <- dimensions[1]
nv <- dimensions[2]
number_perm <- 1000			# number of random permutations
combinations <- factorial(nv) / (factorial(nv-nv1)*factorial(nv1))
gt <- aperm(data, c(2,1)) # transpose the array?? -> yes, this is a simple transpose of the data
Wmat <- as.matrix(dist(gt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
gr <- graph.adjacency(Wmat, weighted = TRUE, mode = "undirected")
V(gr)[c(1:(nv1))]$color <- "green"
V(gr)[c(nv1+1:(nv))]$color <- "red"

mst <- minimum.spanning.tree(gr)