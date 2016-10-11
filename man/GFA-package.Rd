\name{GFA-package}
\alias{GFA-package}
\alias{GFA}
\docType{package}
\title{
Group factor analysis.
}
\description{
GFA does factor analysis for multiple data sets having matched observations, for exploratory or predictive data analysis.
}
\details{
The posterior distribution of GFA model parameters can be inferred with
function \code{\link{gfa}}, once the priors have been defined with
\code{\link{getDefaultOpts}}. The priors are widely customizable, with two
recommended setups: (i) dense group-sparse components (default; similar to
package CCAGFA that provides variational Bayesian inference for the same
model) and (ii) components interpretable as biclusters. It is recommended
to preprocess the data with function \code{\link{normalizeData}}.
Functions are provided for predicting missing data, choosing a prior for
the residual noise, identifying robust components and
visualizing the inferred model. A simple toy example of the pipeline
is provided as demo(GFApipeline), and a more elaborate one as
demo(GFAexample). Finally, the experiment presented in (Bunte, Leppaaho,
Saarinen and Kaski: Sparse group factor analysis for biclustering of
multiple data sources, Bioinformatics, 32(16):2457--2463, 2016) can
be replicated with demo(GFAdream).
}

\keyword{ package }
\examples{
#Data generation
X <- matrix(rnorm(20*3),20,3)                     #Latent variables
W <- matrix(rnorm(30*3),30,3)                     #Projection matrix
Y <- tcrossprod(X,W) + matrix(rnorm(20*30),20,30) #Observations
Y <- sweep(Y, MARGIN=2, runif(30), "+")           #Feature means
Y <- list(Y[,1:10], Y[,11:30])                    #Data grouping
#Model inference and visualization
norm <- normalizeData(Y, type="center")           #Centering
opts <- getDefaultOpts()                          #Model options
#Fast runs for the demo, default options recommended in general
opts[c("iter.burnin", "iter.max")] <- c(500, 1000)
res <- gfa(norm$train, K=5, opts=opts)            #Model inference
rec <- reconstruction(res)                        #Reconstruction
recOrig <- undoNormalizeData(rec, norm)           #... to original space
vis <- visualizeComponents(res, Y, norm)          #Visualization
}
