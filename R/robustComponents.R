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
#' @param corTrh How close two components are required to be, in terms of
#'   correlation, in order to match them.
#' @param matchTrh How big proportion of the chains need to contain the
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
robustComponents <- function(models,corTrh=0.90, matchTrh=0.5) {
  if(length(models)<2) stop("Several sampling chains are needed.")
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

    }
    comps[[rep]] <- list()
    for(v in 1:V) comps[[rep]][[v]] <- list()
    
    for(v in 1:V) {
      K <- res$K[v]
      maxK[v] <- max(K,maxK[v])
      for (k in 1:K) {
        comp <- matrix(0, nrow = N[v], ncol = D[v])
        for (post in 1:res$opts$iter.saved) {
          if(V==1)
            comp <- comp + tcrossprod(res$posterior$X[post,,k], res$posterior$W[post,,k])
          else
            comp <- comp + tcrossprod(res$posterior$X[[v]][post,,k], res$posterior$W[[v]][post,,k])
        }
        comp <- comp / res$opts$iter.saved
        comps[[rep]][[v]][[k]] <- comp
      }
    }
  }
  
  # matching and similarity
  rob <- list(); empty <- matrix(NA,reps,0)
  for(v in 1:V) rob[[v]] <- list(Krobust=0,effect=c(),indices=empty,cor=empty)
  
  for(v in 1:V) {
    compInc <- vector()
    corList <- list()
    for (rep1 in 1:reps) {
      matching <- vector(length = maxK[v])
      d <- array(NA, dim = c(maxK[v], reps))
      matchInd <- matrix(0, nrow = maxK[v], ncol = reps)
      for (k1 in 1:length(comps[[rep1]][[v]])) {
        comp1 <- comps[[rep1]][[v]][[k1]]
        #Component not used yet, and not empty
        if (!(paste(rep1,k1) %in% compInc) & sum(comps[[rep1]][[v]][[k1]]) != 0) {
          for (rep2 in 1:reps) {
            for (k2 in 1:length(comps[[rep2]][[v]])) { #Find component from rep2 best matching k1 of rep1
              if (!(paste(rep2,k2) %in% compInc)) { #Component not used yet
                comp2 <- comps[[rep2]][[v]][[k2]]
                distance <- cor(c(comp1), c(comp2))
                if (!is.na(distance) && (is.na(d[k1,rep2]) || (distance > d[k1,rep2]))) {
                  d[k1,rep2] <- distance
                  if(distance > corTrh)
                    matchInd[k1,rep2] <- k2
                  else
                    matchInd[k1,rep2] <- -k2 #Index saved for debugging purposes
                }
              }
            }
          }
        }
        
        ind <- which(d[k1,] > corTrh) #Which chains are closer than the threshold
        matching[k1] <- length(ind) #How many chains
        
        if(matching[k1] >= matchTrh*reps) {
          # average over all similar components
          robIndices <- robDist <- array(NA,dim=reps,dimnames=list(paste0("rep",1:reps)))
          
          Z_av <- vector(length = N[v])
          W_av <- vector(length = D[v])
          comp <- matrix(0, nrow = N[v], ncol = D[v])
          for(rep2 in 1:reps) {
            robDist[rep2] <- d[k1,rep2]
            robIndices[rep2] <- matchInd[k1,rep2]
            if(matchInd[k1,rep2] > 0) {
              comp <- comp + comps[[rep2]][[v]][[matchInd[k1,rep2]]]
              compInc <- c(compInc, paste(rep2,matchInd[k1,rep2]))
            }
          }
          comp <- comp/sum(matchInd[k1,]>0)
          rob[[v]]$Krobust <- rob[[v]]$Krobust + 1
          rob[[v]]$effect <- array(c(rob[[v]]$effect, comp),dim=c(dim(comp),rob[[v]]$K),
                                   dimnames=list(dimnames[[v]][[1]],dimnames[[v]][[2]],
                                                 paste0("K",1:rob[[v]]$Krobust)))
          rob[[v]]$indices <- cbind(rob[[v]]$indices, robIndices)
          rob[[v]]$cor <- cbind(rob[[v]]$cor, robDist)
        }
      }
    }
  }
  if(V==1) rob <- rob[[1]]
  return(rob)
}
