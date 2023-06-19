# This script is intended to demonstrate the tools for causal discovery provided 
# by the 'pcalg' R package.

# installing pcalg
# install.packages('pcalg')


library(pcalg)
library(graph)
library(Rgraphviz)
library(RBGL)

set.seed(124)

# Exploring elementary objects
# Loading a data set with binary variables
data(gmB)

# The gmB was loaded in the workspace and contains two elements

# data (200 observations)
gmB$x[0:20,]

# model represeted by a DAG
gmB$g
# list of nodes and edges
gmB$g@nodes
gmB$g@edgeL

# plotting the DAG
plot(gmB$g)

# more data sets
# [B]inary, [D]iscrete, [G]aussian
par(mfrow=c(1,3))
data(gmD)
plot(gmD$g)

data(gmG)
plot(gmG$g)

# visualizing gmG distribution of some variables
hist(gmG$x[,2])

# creating a DAG
A <- rbind(c(0,1,0,0,1),
           c(0,0,0,1,1),
           c(1,0,0,1,0),
           c(1,0,0,0,1),
           c(0,0,0,0,0))

gA <- getGraph(A)
par(mfrow=c(1,2))
plot(gA)

## Working with simulated data
# generating a random DAG
myDAG1 <- randomDAG(n = 6, prob= 0.5, lB = 0.1, uB = 1)
myDAG2 <- randDAG(6,3)
par(mfrow=c(1,2))
plot(myDAG1,main="randomDAG")
plot(myDAG2,main="randDAG")
# listing edges
myDAG1@edgeL
myDAG2@edgeL

# myDAG2 is not topologically ordering, then we can use 'tsort()'
tsort(myDAG2)


# For now, we will use a DAG generate via randomDAG (topol. ordered)
par(mfrow=c(1,2))

p <- 8
pconn <- 0.25
myDAG1 <- randomDAG(p, prob = pconn, lB=0.1, uB=1)

plot(myDAG1, main = "randomDAG")

## Causal search without hidden and selection variables
# Estimate (Initial) Skeleton of a DAG using the PC / PC-Stable Algorithm
##################################################
## Using Gaussian Data
##################################################
data(gmG)
n <- nrow (gmG8$x)
V <- colnames(gmG8$x) # labels aka node names
# estimate Skeleton
skel.fit <- skeleton(suffStat = list(C = cor(gmG8$x), n = n),
                     indepTest = gaussCItest, ## (partial correlations)
                     alpha = 0.01, labels = V, verbose = TRUE)

# show estimated Skeleton
par(mfrow=c(1,2))
plot(skel.fit, main = "Estimated Skeleton")
plot(gmG8$g, main = "True DAG")
  
##################################################
## Using discrete data
##################################################
## Load data
data(gmD)
V <- colnames(gmD$x) # labels aka node names
## define sufficient statistics
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)
## estimate Skeleton
skel.fit <- skeleton(suffStat,
                     indepTest = disCItest, ## (G^2 statistics independence test)
                     alpha = 0.01, labels = V, verbose = FALSE)
## show estimated Skeleton
par(mfrow = c(1,2))
plot(skel.fit, main = "Estimated Skeleton")
plot(gmD$g, main = "True DAG")

##################################################
## Using binary data
##################################################
## Load binary data
data(gmB)
X <- gmB$x
## estimate Skeleton
skel.fm2 <- skeleton(suffStat = list(dm = X, adaptDF = FALSE),
                     indepTest = binCItest, alpha = 0.01,
                     labels = colnames(X), verbose = TRUE)

## show estimated Skeleton
par(mfrow = c(1,2))
plot(skel.fm2, main = "Binary Data 'gmB': Estimated Skeleton")
plot(gmB$g, main = "True DAG")

##################################################
## Using simulated data
##################################################
# Generating a DAG with node degrees >= 1
p <- 6
pconn <- 0.3
zerodegFlag <- TRUE

while (zerodegFlag ) {
  ## true DAG
  myDAG1 <- randomDAG(p, prob = pconn) 
  am <- (as(myDAG1,"matrix")!=0)*1L
  
  if (any(colSums(am + t(am)) == 0)) {
    zerodegFlag = TRUE
    next
  } else{
    zerodegFlag = FALSE
  }
}
plot(myDAG1, main = "randomDAG")


# Generating 1000 observations from myDAG1 using standard normal error distribution
n <- 1000
d.normMat <- rmvDAG(n, myDAG1, errDist="normal")

# Visualizing the distribution of the first four variables
par(mfrow=c(2,2))
for (i in 1:4){
  hist(d.normMat[,i])  
}

n <- nrow (d.normMat)
V <- colnames(d.normMat) # labels aka node names

skel.fit <- skeleton(suffStat = list(C = cor(d.normMat), n = n),
                     indepTest = gaussCItest, ## (partial correlations)
                     alpha = 0.01, labels = V, verbose = FALSE)
# show estimated Skeleton
par(mfrow=c(1,2))
plot(skel.fit, main = "Estimated Skeleton")
plot(myDAG1, main = "True DAG")

## PC
##################################################
## Using Gaussian Data
##################################################
## Load predefined data
n <- nrow (gmG8$ x)
V <- colnames(gmG8$ x) # labels aka node names
## estimate CPDAG
pc.fit <- pc(suffStat = list(C = cor(gmG8$x), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = FALSE)

## show estimated CPDAG
par(mfrow=c(1,2))
plot(pc.fit, main = "Estimated CPDAG")
plot(gmG8$g, main = "True DAG")

##################################################
## Using simulated data
##################################################
# Generating a DAG with node degrees >= 1
n <- nrow (d.normMat)
V <- colnames(d.normMat) # labels aka node names

pc.fit <- pc(suffStat = list(C = cor(d.normMat), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = FALSE)
# show estimated Skeleton
par(mfrow=c(1,2))
plot(pc.fit, main = "Estimated Skeleton")
plot(myDAG1, main = "True DAG")


##################################################
## FCI and RFCI
##################################################

## FCI with latent variables
data("gmL")

# visualizing the first 6 rows in the data set
# dimension of the data set (rows=observations, cols=Variables)
head(gmL$x)

# Number of nodes and edges in the true Model
gmL$g

suffStat1 <- list(C = cor(gmL$x), n = nrow(gmL$x))
pagfci.est <- fci(suffStat1, indepTest = gaussCItest, 
               p = ncol(gmL$x), alpha = 0.01, 
               labels = as.character(2:5))

par(mfrow = c(1, 2))
plot(pagfci.est)
plot(gmL$g, main = "True Model")

## RFCI
# same data set gML

suffStat1 <- list(C = cor(gmL$x), n = nrow(gmL$x))

pagrfci.est <- rfci(suffStat1, indepTest = gaussCItest, 
                p = ncol(gmL$x),alpha = 0.01, 
                labels = as.character(2:5))

par(mfrow = c(1, 3))
plot(pagfci.est)
mtext("FCI",side=3)
plot(pagrfci.est)
mtext('RFCI',side=3)
plot(gmL$g, main = "True Model")

### Comparing graphs via Structural Hamming Distance
Am <- rbind(c(0,1,0,0,1),
           c(0,0,0,1,1),
           c(1,0,0,1,0),
           c(1,0,0,0,1),
           c(0,0,0,0,0))

A <- as(getGraph(Am),"graphNEL")
A1 <- removeEdge(from="3",to="4",A)
A2 <- removeEdge(from="4",to="1",A1)

par(mfrow =c(1,3))
plot(A, main="A")
plot(A1, main="A1")
plot(A2, main="A2")

shdAA <- shd(A,A)
shdAA1 <- shd(A,A1)
shdAA2 <- shd(A,A2)

print(shdAA)
print(shdAA1)
print(shdAA2)

# normalized SHD value?

################## causal effects computation ida
set.seed(123)
p <- 10
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
## simulate Gaussian data from the true DAG
n <- 10000
dat <- rmvDAG(n, myDAG)
## estimate CPDAG and PDAG -- see help(pc)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p=p, alpha = 0.01)

par(mfrow=c(1,2))
plot(myDAG, main="True DAG")
plot(pc.fit, main="CPDAG")

# unfolding the CPDAG
amat <- as(pc.fit,"amat")
all.dags <- pdag2allDags(amat)

# for convenience a simple plotting function
# for the function output
plotAllDags <- function(res) {
  require(graph)
  p <- sqrt(ncol(res$dags))
  nDags <- ceiling(sqrt(nrow(res$dags)))
  par(mfrow = c(nDags, nDags))
  for (i in 1:nrow(res$dags)) {
    tmp <- matrix(res$dags[i,],p,p)
    colnames(tmp) <- rownames(tmp) <- res$nodeNms
    plot(as(tmp, "graphNEL"))
  }
}

plotAllDags(all.dags)

# computing the causal effects 
#ida(X, Y, mcov, CPDAG, method = c("local","global"))
ida(6,10, cov(dat), pc.fit@graph, method = "global", type = "cpdag")
ida(6,10, cov(dat), pc.fit@graph, method = "local", type = "cpdag")

#True causal effect
# causalEffect(DAG, Y , X)
causalEffect(myDAG, 10, 6)

sessionInfo()
