#' Sequential prediction of new samples from observed data views to unobserved
#'
#' \code{sequentialGfaPrediction} returns predictions for unobserved data
#' sources of a new set of data samples.
#' 
#' @param Y The new data samples, in a similar format as in function
#'   \code{\link{gfa}}
#' @param model The sampled model from function \code{\link{gfa}}, with
#'   model$opts$predict being a logical vector with the length matching the
#'   length of Y, describing which data views will be predicted (TRUE), and
#'   which have been observed (FALSE).
#' @return A list containing the predictions, with the observed
#'   data sources empty. Additionally, sampling MSE is given as element 'mse',
#'   and likelihood as 'cost'.
sequentialGfaPrediction <- function(Y, model) {
  opts <- model$opts
  if(typeof(Y[[1]])=="list") {
    V <- length(Y)
    if(V>1) stop("TODO: implementation for sequential learning or GFA paired in two modes.")
  } else {
    Y <- list(Y); V <- 1
    opts$prediction <- list(opts$prediction)
    model$groups <- list(model$groups)
    for(param in c("W","tau","beta","rz"))
      model$posterior[[param]] <- list(model$posterior[[param]])
  }
  opts$iter.max <- 20
  opts$iter.burnin <- 10
  opts$iter.saved <- 5
  if(length(opts$prediction)!=V) stop("Specify which data views should be predicted.")
  for(v in 1:V) if(length(opts$prediction[[v]])!=length(Y[[v]])) stop("Ill-specified 'model$opts$prediction'.")
  
  N <- rep(NA,V)
  X <- Y.prediction <- mse <- cost <- W <- beta <- tau <- rz <- list()
  for(v in 1:V) {
    N[v] <- nrow(Y[[v]][[which(!opts$prediction[[v]])[1]]]) # Number of samples in the prediction set
    
    X[[v]] <- vector("list",length=length(Y[[v]])) #Store the data to be predicted
    
    Y.prediction[[v]] <- vector(mode="list", length=length(model$groups))
    names(Y.prediction[[v]]) <- names(Y[[v]]) # names(De)
    for (mi in which(opts$prediction[[v]])) { # Go through all views that will be predicted.
      Y.prediction[[v]][[mi]] <- array(data=0, dim=c(N[v], length(model$groups[[v]][[mi]]))) # Initialize the prediction matrix for view 'mi'.
      rownames(Y.prediction[[v]][[mi]]) <- rownames(Y[[v]][[which(!opts$prediction[[v]])[1]]]) # Names of the prediction samples
      colnames(Y.prediction[[v]][[mi]]) <- names(model$groups[[v]][[mi]]) # Names of the variables in view 'mi', which will be predicted.
      X[[v]][[mi]] <- Y[[v]][[mi]]
      Y[[v]][[mi]] <- array(dim=dim(Y.prediction[[v]][[mi]]), dimnames=dimnames(Y.prediction[[v]][[mi]])) # Initialize view 'mi' as missing data. These will be predicted.
    }
    
    ##
    ## Prediction
    ##
    mse[[v]] <- matrix(NA,length(which(opts$prediction[[v]])),nrow(model$posterior$W[[v]]))
    cost[[v]] <- matrix(NA,opts$iter.max,nrow(model$posterior$W[[v]]))
  }
  
  for (ni in 1:nrow(model$posterior$W[[1]])) { # Go through all saved Gibbs samples.
    if(opts$verbose>1)
      print(paste0("Predicting, Gibbs sample ",ni))
    for(v in 1:V) {
      W[[v]] <- matrix(model$posterior$W[[v]][ni,,],ncol(model$posterior$W[[v]]),dim(model$posterior$W[[v]])[3])
      if(!opts$normalLatents)
        beta[[v]] <- matrix(model$posterior$beta[[v]][ni,,],ncol(model$posterior$beta[[v]]),dim(model$posterior$beta[[v]])[3])
      tau[[v]] <- model$posterior$tau[[v]][ni,]
      rz[[v]] <- model$posterior$rz[[v]][ni,]
    }
    prediction.ni <- gfa(Y=Y, K=NULL, opts=opts, projection=list(W=W, rz=rz, tau=tau, beta=beta))
    for(v in 1:V) {
      if(any(opts$prediction[[v]])) {
        for (mi in which(opts$prediction[[v]])) { # Go through the target views.
          Y.prediction[[v]][[mi]] = Y.prediction[[v]][[mi]] + prediction.ni[[v]][[mi]]/nrow(model$posterior$W[[v]])
          const <- nrow(model$posterior$W[[v]])/ni
          mse[[v]][which(mi==which(opts$prediction[[v]])),ni] <- mean((X[[v]][[mi]]-Y.prediction[[v]][[mi]]*const)^2)
          cost[[v]][,ni] <- prediction.ni[[v]]$cost
          if(opts$verbose>1)
            print(paste0("MSE at iteration ",ni,": ",mse[[v]][which(mi==which(opts$prediction[[v]])),ni]))
        }
      }
    }
  }
  for(v in 1:V) {
    Y.prediction[[v]]$mse <- mse
    Y.prediction[[v]]$cost <- cost
  }
  if(V==1) Y.prediction <- Y.prediction[[1]]
  
  return(Y.prediction)
  
}
