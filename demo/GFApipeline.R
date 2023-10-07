
require("GFA")
## Generating a simple data collection with 3 components:
##    2 dense ones and 1 sparse (bicluster).
## Aiming to show that GFA can infer this kind of structure accurately.
## The data sources, samples and features are named to mimick the multi source
## drug sensitivity study of Bunte et al.
set.seed(12345)
N <- 23; D <- c(6,25,25); K <- 3 #Data dimensions
X <- matrix(rnorm(N*K),N,K) #Latent variables
X[c(1:5,15:19),1] <- 0 #Sparsity for the first (bicluster) component
W <- matrix(rnorm(sum(D)*K),sum(D),K) #Loading matrix
W[c(3:4,10:25,32:sum(D)),1] <- 0 #Sparsity for k=1
W[1:6,3] <- 0 #k=3 inactive in the first view
W[2,1:2] <- W[6,1:2] + rnorm(2,0,0.2) #2nd and 6th feature similar
Y <- tcrossprod(X,W) + matrix(rnorm(N*sum(D),0,0.5),N,sum(D)) #Data generation
Y <- list(Y[,1:D[1]], Y[,1:D[2]+D[1]], Y[,1:D[3]+sum(D[1:2])]) #Splitting into groups
Y[[1]] <- sweep(Y[[1]],MARGIN=2,runif(D[1],2,4),"+") #Means for Y[[1]]
Y[[1]] <- sweep(Y[[1]],MARGIN=2,runif(D[1],0.5,2),"*") #Scale for Y[[1]]
#Data naming
names(Y) <- c("Drug sensitivity","Gene expression","Exome sequence")
for(m in 1:length(Y)) rownames(Y[[m]]) <- paste("Cell line",1:N)
colnames(Y[[1]]) <- paste("Drug",toupper(letters[1:D[1]]))
colnames(Y[[2]]) <- colnames(Y[[3]]) <- paste("Gene",1:D[2])
missing <- Y[[1]][7,2]
Y[[1]][7,2] <- NA #One missing value (to be predicted)

## DATA GENERATED, RUNNING THE MODEL PIPELINE
## Normalize the data - here we assume that every feature is equally important
norm <- normalizeData(Y,type="scaleFeatures")
## Get the model options - we want to detect bicluster structure
opts <- getDefaultOpts(bicluster=TRUE)
opts$convergenceCheck <- TRUE #We want to check for sampling chain convergence
## Giving a vague prior stating that roughly 1/4 of the variance is noise
opts <- informativeNoisePrior(Y,opts,noiseProportion=0.25,conf=0.1)
## Infer the model
res <- gfa(norm$train,opts=opts,K=6)
## Reconstruct the data based on inferred components (can be used to predict missing values)
rec <- reconstruction(res)
recOrig <- undoNormalizeData(rec, norm) #In the original data space
## Visualize the model
par(mfrow=c(2,4),mar=rep(4,4),oma=rep(0,4))
vis <- visualizeComponents(res,Y,norm,hclust=TRUE)

## FIGURE CAPTIONS 
# - Component activities: shows how strongly the 3 components are associated
#   with different data sources. White (or close to white) elements show that
#   a certain component is not associated with a certain data source.
#   In this case GFA is able to infer the correct component number, but in
#   general some components can be used to model (structured) noise as well.
# - Data reconstruction: the factorization XW^T, showing how the components
#   model the data collection.
# - Component X effect: The component active mostly in drug sensitivity and gene
#   expression, in a subset of the cell lines and drugs/genes, is a bicluster
#   component used in the data generation.

## PREDICTION OF A MISSING VALUE
cat(paste0("Missing drug sensitivity (cell line 7, drug B): ",round(missing,2),"\n",
           "Predicted drug sensitivity (cell line 7, drug B): ",round(recOrig[7,2],2),"\n"))

## LATENT SPACE DISTANCES TO DRUG B (F GENERATED AS MOST SIMILAR)
d <- as.matrix(dist(res$W))[2,1:6]; names(d) <- colnames(Y[[1]]); d <- d[-2]
print(sort(d))


