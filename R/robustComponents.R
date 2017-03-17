#' Robust GFA components
#'
#' \code{robustComponents} analyzes a collection of sampling chains and returns
#'   robust components.
#'
#' The function returns the effects (i.e. reconstructions) of robust components
#' to the data level. It is useful for a thorough model interpretation,
#' accumulating power over several sampling chains by
#' comparing them in the observation space (as opposed to the latent space).
#' The function is needed for this task, as the extreme multi-modality of factor
#' analysis prohibits efficient sampling techniques that would result in a
#' posterior estimate converging to the true posterior in practice.
#' The function uses a heuristic correlation-based procedure to analyze which
#' components occur frequently in GFA sampling chains.
#' 
#' @param models Either a vector containing the file names, where the models are
#'   saved as 'res', or a list containing the models.
#' @param corThr How close two components are required to be, in terms of
#'   correlation, in order to match them.
#' @param matchThr How big proportion of the chains need to contain the
#'   component to include it in the robust components.
#' @return A list with the following elements (when input data are paired in two
#'   modes, the returned list is of length 2, containing the following elements
#'   for each mode):
#' \item{Krobust}{The number of robust components found with the given
#'   thresholds.}
#' \item{effect}{The component effect in the data space; and array of size
#'   \eqn{N \times \sum_{m=1}^M D_m \times Krobust.}}
#' \item{indices}{The corresponding component indices; a \eqn{length(models)
#'   \times Krobust} matrix. Negative indices denote the closest component in
#'   the corresponding repetition, having no components above the threshold.}
#' \item{cor}{The correlations of the components matched to this robust
#'   component; a matrix of size \eqn{length(models) \times Krobust}. The
#'   correlations are reported relative to the first repetition with this
#'   component observed.}
#' @examples
#' X <- matrix(rnorm(10*2),10,2)
#' W <- matrix(rnorm(15*2),15,2)
#' Y <- tcrossprod(X,W) + matrix(rnorm(10*15),10,15)
#' opts <- getDefaultOpts() #Default options
#' #Fast runs for the demo, default options recommended in general
#' opts[c("iter.burnin", "iter.max")] <- c(500, 1000)
#' res <- list()
#' for(i in 1:4) res[[i]] <- gfa(list(Y[,1:6],Y[,7:15]),opts=opts,K=3)
#' rob <- robustComponents(res)
robustComponents <- function(models,corThr=0.90, matchThr=0.5) {
  if(length(models)<2) stop("Several sampling chains are needed in input parameter 'models'.")
  maxK <- c(0,0)
  comps <- list()
  reps <- length(models)
  dimnames <- list()
  for (rep in 1:reps) {
    if(typeof(models[[rep]])=="list")
      res <- models[[rep]]
    else
      load(models[rep])
    
    if(rep==1) {
      if(typeof(res$X)=="list") {
        V <- length(res$X); N <- sapply(res$X,nrow); D <- sapply(res$W,nrow)
        for(v in 1:V)
          dimnames[[v]] <- list(colnames(res$posterior$X[[v]]),colnames(res$posterior$W[[v]]))
      } else {
        V <- 1; N <- nrow(res$X); D <- nrow(res$W)
        dimnames[[1]] <- list(colnames(res$posterior$X),colnames(res$posterior$W))
      }
      comps <- list()
      for(v in 1:V) comps[[v]] <- list()
    }
    comps[[v]][[rep]] <- list()
    
    for(v in 1:V) {
      K <- res$K[v]
      maxK[v] <- max(K,maxK[v])
      for (k in 1:K) {
        if(V==1)
          comp <- crossprod(res$posterior$X[,,k], res$posterior$W[,,k]) #Sum over posterior samples
        else
          comp <- crossprod(res$posterior$X[[v]][,,k], res$posterior$W[[v]][,,k])
        comp <- comp / res$opts$iter.saved
        comps[[v]][[rep]][[k]] <- comp
      }
    }
  }
  
  # Matching and similarity
  rob <- list()
  for(v in 1:V) {
    rob[[v]] <- matchComponents(comps=comps[[v]], maxK=maxK[v], N=N[v], D=D[v],
                                dimnames=dimnames[[v]], corThr=corThr, matchThr=matchThr)
  }
  
  if(V==1) rob <- rob[[1]]
  return(rob)
}



## Internal function for matching components
matchComponents <- function(comps, maxK, N, D, dimnames, corThr, matchThr) {
  corList <- list()
  reps <- length(comps)
  rob <- list(Krobust=0,effect=c(),indices=matrix(NA,reps,0),cor=matrix(NA,reps,0))
  
  compStorage <- vector("list",length=reps) #Store the components that can still be used
  for(rep in 1:reps) compStorage[[rep]] <- which(sapply(comps[[rep]], function(x) {sum(abs(x)) > 0}))
  
  for (rep1 in 1:reps) {
    matching <- vector(length = maxK)
    sim <- array(NA, dim = c(maxK, reps))
    matchInd <- matrix(0, nrow = maxK, ncol = reps)
    
    for (k1 in compStorage[[rep1]]) {
      for (rep2 in which(sapply(compStorage, length)>0)) {
        #Correlation between the two components.
        #Note that we do not need to use the absolute values, as the correlation is computed in data space.
        d <- sapply(comps[[rep2]][compStorage[[rep2]]], function(x) cor(c(x), c(comps[[rep1]][[k1]])))
        sim[k1,rep2] <- max(d)
        matchInd[k1,rep2] <- compStorage[[rep2]][which.max(d)]
        if(sim[k1,rep2] < corThr) matchInd[k1,rep2] <- matchInd[k1,rep2]*-1 #Index saved for debugging purposes
      }
      
      if(sum(sim[k1,]>corThr, na.rm=T) >= matchThr*reps) { #Robust component found!
        # average over all similar components
        comp <- matrix(0, N, D)
        for(rep2 in which(matchInd[k1,] > 0)) {
          comp <- comp + comps[[rep2]][[matchInd[k1,rep2]]]
          compStorage[[rep2]] <- compStorage[[rep2]][!compStorage[[rep2]]==matchInd[k1,rep2]]
        }
        comp <- comp/sum(matchInd[k1,]>0)
        
        rob$Krobust <- rob$Krobust + 1
        rob$effect <- array(c(rob$effect, comp),dim=c(dim(comp),rob$K),
                            dimnames=list(dimnames[[1]],dimnames[[2]],
                                          paste0("K",1:rob$Krobust)))
        rob$indices <- cbind(rob$indices, matchInd[k1,])
        rob$cor <- cbind(rob$cor, sim[k1,])
      }
    }
  }
  rownames(rob$indices) <- rownames(rob$cor) <- paste0("rep",1:reps)
  colnames(rob$indices) <- colnames(rob$cor) <- paste0("K",1:rob$Krobust)
  return(rob)
}




