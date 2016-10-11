#' @title Full data reconstruction based on posterior samples
#'
#' @description
#' \code{reconstruction} returns the full data reconstruction based on given
#' posterior samples.
#' 
#' @param res The sampled model from function \code{\link{gfa}}
#' @param average If TRUE (default), averages the reconstruction over the
#'   posterior predictive samples. If set to
#'   FALSE, the output may require a large amount of memory. In case of large
#'   input data, we recommend acquiring the posterior predictive samples for
#'   subsets of data at a time, based on this implementation.
#' @return The data reconstruction, a numeric  \eqn{N \times \sum_{m=1}^M D_m} 
#'   matrix, if average is TRUE (default). Otherwise, the reconstruction is a
#'   \eqn{N \times \sum_{m=1}^M D_m \times Npost} array, with posterior samples
#'   in the third dimension. If the input data has been paired in two modes, the
#'   output will be a list of length 2, one element corresponding to each mode.

reconstruction <- function(res,average=TRUE) {
  if(res$opts$verbose>0 & is.null(res$posterior))
    print("Reconstructing data from the final posterior sample only.")
  if(typeof(res$X)=="list") { # Data paired in two modes
    Y <- list()
    for(v in 1:2) {
      if(is.null(res$posterior)) {
        Y[[v]] <- tcrossprod(res$X[[v]],res$W[[v]])
      } else {
        if(average) {
          Y[[v]] <- matrix(0,nrow(res$X[[v]]),nrow(res$W[[v]]))
          for(i in 1:res$opts$iter.saved)
            Y[[v]] <- Y[[v]] + tcrossprod(res$posterior$X[[v]][i,,], res$posterior$W[[v]][i,,])
          Y[[v]] <- Y[[v]]/res$opts$iter.saved
          
        } else {
          Y[[v]] <- array(NA,dim=c(nrow(res$X[[v]]),nrow(res$W[[v]]),res$opts$iter.saved))
          for(i in 1:res$opts$iter.saved)
            Y[[v]][,,i] <- tcrossprod(res$posterior$X[[v]][i,,], res$posterior$W[[v]][i,,])
        }
      }
    }
    #Effect from both the factorizations
    if(average) {
      Y[[1]][,res$groups[[1]][[1]]] <- Y[[1]][,res$groups[[1]][[1]]] + t(Y[[2]][,res$groups[[2]][[1]]])
      Y[[2]][,res$groups[[2]][[1]]] <- t(Y[[1]][,res$groups[[1]][[1]]])
    } else {
      for(i in 1:dim(Y[[1]])[3])
      Y[[1]][,res$groups[[1]][[1]],i] <- Y[[1]][,res$groups[[1]][[1]],i] + t(Y[[2]][,res$groups[[2]][[1]],i])
      Y[[2]][,res$groups[[2]][[1]],i] <- t(Y[[1]][,res$groups[[1]][[1]],i])
    }

    
  } else { # A single set of data sources, paired in one mode
    if(is.null(res$posterior)) {
      Y <- tcrossprod(res$X,res$W)
    } else {
      if(average) {
        Y <- matrix(0,nrow(res$X),nrow(res$W))
        for(i in 1:res$opts$iter.saved)
          Y <- Y + tcrossprod(res$posterior$X[i,,], res$posterior$W[i,,])
        Y <- Y/res$opts$iter.saved
        
      } else {
        Y <- array(NA,dim=c(nrow(res$X),nrow(res$W),res$opts$iter.saved))
        for(i in 1:res$opts$iter.saved)
          Y[,,i] <- tcrossprod(res$posterior$X[i,,], res$posterior$W[i,,])
      }

    }
  }
  return(Y)
}
