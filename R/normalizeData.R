#' Normalize data to be used by GFA
#'
#' \code{normalizeData} is used to transform a data collection into a normalized
#' form suitable for GFA.
#' This function does two things: 1. It centers each variable. GFA assumes
#' zero-mean data, as it models variances. 2. It normalizes the scales of
#' variables and/or variable groups. Features with higher variance will
#' affect the model structure more; if this is not desired, the normalization
#' should be done. In GFA it is additionally possible to normalize the
#' importance of variable groups (data sources), in addition or instead of
#' individual variables. Finally, the total variance of data is normalized for
#' numerical reasons. This is particularly important if no other normalization
#' is done. NOTE: the function assumes continuous-valued data.
#' If some features are e.g. binary with only a small portion of 1s, we do not
#' recommend centering them.
#' 
#' @param train a training data set. For a detailed description, see parameter Y
#'        in \code{\link{gfa}}.
#' @param test a test dataset. Should be provided if sequential prediction
#'        is used later.
#' @param type Specifies the type of normalization to do. Mean-centering of the
#'        features is performed in all the cases, and option "center"
#'        does not perform any scaling. Option "scaleOverall" (default) uses a
#'        single parameter to scale the variance of the
#'        whole data collection to 1, while "scaleSources" scales each data
#'        source to have variance 1. Finally, "scaleFeatures" performs
#'        z-normalization, i.e. assigns the variance of each feature to 1. 
#' @return A list containing the following elements:
#' \item{train}{Normalized training data.}
#' \item{test}{Normalized test data for sequential prediction (if provided as
#'         input).}
#' \item{trainMean}{Feature-wise means of the training data sources.}
#' \item{trainSd}{Feature-wise/overall standard deviations of the training data
#'         sources.}
#' 
normalizeData <- function(train, test=NULL,type="scaleOverAll") {
  if(typeof(train[[1]])=="list") { #Pairing in two modes
    V <- length(train)
    if(V==2)
      if(any(train[[1]][[1]]!=t(train[[2]][[1]]),na.rm=T))
        stop("The two sets of views must have the first view shared (transposed).")
  } else {
    train <- list(train); V <- 1
    if(!is.null(test)) test <- list(test)
  }
  
  trainMean <- trainSd <- list()
  for(v in 1:V) {
    trainMean[[v]] <- lapply(train[[v]],colMeans,na.rm=T)
    trainSd[[v]] <- lapply(train[[v]], function(x) rep(1,ncol(x)))
    for(m in V:length(train[[v]])) { #If V==2, skip m=1 for now
      #Center the data
      train[[v]][[m]] <- sweep(train[[v]][[m]], MARGIN=2, trainMean[[v]][[m]], "-")
      if(!is.null(test)) test[[v]][[m]] <- sweep(test[[v]][[m]], MARGIN=2, trainMean[[v]][[m]], "-")
      if(type=="center") {
        #All done
      } else if(type=="scaleFeatures") {
        trainSd[[v]][[m]] <- apply(train[[v]][[m]],2,sd,na.rm=T)
        trainSd[[v]][[m]][trainSd[[v]][[m]]==0] <- 1
      } else if(type=="scaleSources") {
        trainSd[[v]][[m]][] <- sd(c(train[[v]][[m]]),na.rm=T)
      } else if(type=="scaleOverall") { #First remove the means from each source
        
      } else {
        stop("Invalid normalization type specified.")
      }
      #Scale the data
      train[[v]][[m]] <- sweep(train[[v]][[m]], MARGIN=2, trainSd[[v]][[m]], "/")
      if(!is.null(test)) test[[v]][[m]] <- sweep(test[[v]][[m]], MARGIN=2, trainSd[[v]][[m]], "/")
    }
  }
  if(type=="scaleOverall") {
    s <- sd(unlist(train),na.rm=T)
    for(v in 1:V) {
      trainSd[[v]] <- lapply(train[[v]], function(x) rep(s,ncol(x)))
      for(m in 1:length(train[[v]]))
        train[[v]][[m]] <- train[[v]][[m]]/s
    }
  }
  if(V==2 & !type=="scaleOverall") { #The shared data source separately
    ov <- c(2,1) #The other mode
    #Centering
    for(v in 1:2) {
      trainMean[[v]][[1]] <- apply(train[[1]][[1]],ov[v],mean,na.rm=T)
      train[[1]][[1]] <- sweep(train[[1]][[1]], MARGIN=ov[v], trainMean[[v]][[1]], "-")
      if(!is.null(test)) test[[1]][[1]] <- sweep(test[[1]][[1]], MARGIN=ov[v], trainMean[[v]][[1]], "-")
    }
    #Scaling
    for(v in 1:2) {
      if(type=="scaleFeatures") {
        trainSd[[v]][[1]] <- apply(train[[1]][[1]],ov[v],sd,na.rm=T)
        trainSd[[v]][[1]][trainSd[[v]][[1]]==0] <- 1
      } else if(type=="scaleSources") {
        trainSd[[v]][[1]][] <- sd(c(train[[1]][[1]]),na.rm=T)
      }
      train[[1]][[1]] <- sweep(train[[1]][[1]], MARGIN=ov[v], trainSd[[v]][[1]], "/")
      if(!is.null(test)) test[[1]][[1]] <- sweep(test[[1]][[1]], MARGIN=ov[v], trainSd[[v]][[1]], "/")
    }
    train[[2]][[1]] <- t(train[[1]][[1]])
    if(!is.null(test)) test[[2]][[1]] <- t(test[[1]][[1]])
  }
  if(V==1) {
    params <- c("train","trainMean","trainSd")
    if(!is.null(test)) params <- c(params, "test")
    for(param in params)
      eval(parse(text=paste0(param," <- ",param,"[[1]]")))
  }

  return(list(train=train,test=test,trainMean=trainMean,trainSd=trainSd,
              type=type))
}


#' A function for returning predictions into the original data space
#'
#' \code{undoNormalizeData} returns the predictions on normalized data acquired
#' from \code{\link{normalizeData}} into the original data space.
#'
#' @param pred The predictions acquired from \code{\link{reconstruction}}.
#' @param normalization The output list obtained from
#'        \code{\link{normalizeData}}.
#' @return The predictions in the original data space.
undoNormalizeData <- function(pred,normalization) {
  mu <- normalization$trainMean
  sigma <- normalization$trainSd
  if(typeof("pred")!="list") { #One data collection paired in one mode
    if(length(dim(pred[[1]]))>2)
      stop("pred is not a matrix - average over the reconstruction.")
    pred <- sweep(pred, MARGIN=2, unlist(sigma), "*")
    pred <- sweep(pred, MARGIN=2, unlist(mu), "+")
    
  } else { #Two data collections paired in two modes
    V <- length(pred)
    if(length(dim(pred[[1]]))>2)
      stop("pred is not a matrix - average over the reconstruction.")
    
    ov <- c(2,1) #The other mode
    id <- list(1:length(mu[[2]][[1]]), 1:length(mu[[1]][[1]])) #Paired source
    for(v in V:1) { #Reverse order for the paired source
      pred[[1]][,id[[1]]] <- sweep(pred[[1]][,id[[1]]], MARGIN=ov[v], sigma[[v]][[1]], "*")
      pred[[v]][,-id[[v]]] <- sweep(pred[[v]][,-id[[v]]], MARGIN=2, unlist(sigma[[v]])[-id[[v]]], "*")
      pred[[v]][,-id[[v]]] <- sweep(pred[[v]][,-id[[v]]], MARGIN=2, unlist(mu[[v]])[-id[[v]]], "+")
    }
    for(v in V:1)
      pred[[1]][,id[[1]]] <- sweep(pred[[1]][,id[[1]]], MARGIN=ov[v], mu[[v]][[1]], "+")
  }
  return(pred)
}
