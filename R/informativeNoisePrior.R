#' Informative noise residual prior
#' 
#' \code{informativeNoisePrior} returns an informative noise prior for GFA, for
#' a given data collection. The function sets the noise residual parameters
#' such that the expected proportion of
#' variance explained is for all variables and groups (in contrast to being
#' proportional to their original scale). Recommended e.g. when the data is
#' 'small n, large p', and the standard prior from \code{\link{getDefaultOpts}}
#' seems to overfit the model by not shutting off any component with high
#' initial K.
#'
#' @param Y The data. For details, see function \code{\link{gfa}}.
#' @param opts Model options. See function \code{\link{getDefaultOpts}} for
#'   details. If option tauGrouped is TRUE (default), each data source is given
#'   equal importance (feature importance may vary within each
#'   source). If it is FALSE, each feature is given equal importance.
#' @param noiseProportion proportion of total variance to be
#'   explained by noise. Suggested to lie between 0.01 and 0.99.
#' @param conf Confidence in the prior, relative to confidence in the data.
#'   Suggested to lie between 0.01 and 100.
#' @return The input model options (opts) with an updated residual
#'   noise prior, corresponding to the elements prior.alpha_0t and
#'   prior.beta_0t. 
#' @examples
#' #Given data collection Y
#' opts <- getDefaultOpts()
#' \dontrun{opts <- informativeNoisePrior(Y,opts,0.2,1)}
#' \dontrun{res <- gfa(Y,opts=opts)}
informativeNoisePrior <- function(Y,opts,noiseProportion=0.5,conf=1) {
  if(typeof(Y[[1]])=="list") {
    V <- length(Y)
  } else {
    Y <- list(Y); V <- 1
  }
  
  prior.alpha_0t <- prior.beta_0t <- list()
  for(v in 1:V) {
    M <- length(Y[[v]])
    D <- unlist(lapply(Y[[v]],ncol))
    Ds <- c(1,cumsum(D)+1) ; De <- cumsum(D)
    gr <- vector("list")
    for(m in 1:M) {
      gr[[m]] <- Ds[m]:De[m]
    }
    N <- nrow(Y[[v]][[1]])
    prior.alpha_0t[[v]] <- rep(NA,sum(D))
    prior.beta_0t[[v]] <- rep(NA,sum(D))
    
    if(opts$tauGrouped) {
      for (m in 1:M) {
        prior.alpha_0t[[v]][gr[[m]]] <- conf*sum(!is.na(Y[[v]][[m]]))/2
        totalVariance <- sum(Y[[v]][[m]]^2,na.rm=T)/2
        sigma <- noiseProportion * totalVariance
        prior.beta_0t[[v]][gr[[m]]] <- conf*sigma
      }
    } else {
      if(V==2) stop("Tau needs to be grouped.")
      for (m in 1:M) {
        prior.alpha_0t[[v]][gr[[m]]] <- conf*colSums(!is.na(Y[[v]][[m]]))/2
        totalVariance <- colSums(Y[[v]][[m]]^2,na.rm=T)/2
        sigma <- noiseProportion * totalVariance
        prior.beta_0t[[v]][gr[[m]]] <- conf*sigma
      }
    }
  }
  opts$prior.alpha_0t <-  prior.alpha_0t
  opts$prior.beta_0t <-  prior.beta_0t
  return(opts)
}
