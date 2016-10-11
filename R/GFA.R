#' Gibbs sampling for group factor analysis
#'
#' \code{gfa} returns posterior samples of group factor analysis model.
#'
#' GFA allows factor analysis of multiple data sources (i.e. data sets).
#' The priors of the model can be set to infer bicluster structure
#' from the data sources; see \code{\link{getDefaultOpts}}.
#' Missing values (NAs) are inherently supported. They will not affect the model
#' parameters, but can be predicted with function \code{\link{reconstruction}},
#' based on the observed values of the corresponding sample and feature.
#' The association of a data source to each component is inferred based on
#' the data. Letting only a subset of the components to explain a data source
#' results in the posterior identifying relationships between any subset of the
#' data sources. In the extreme cases, a component can explain relationships
#' within a single data source only ("structured noise"), or across all the data
#' sources.
#'
#' @param Y Either \enumerate{
#'   \item{Data sources with co-occuring samples: a list of data
#'   matrices, where Y[[m]] is a numeric \eqn{N \times D_m} matrix, or}
#'   \item{Data sources paired in two modes (some data sources share the
#'   samples of the first data source, and some share its features): A list with
#'   two elements structured as 1. The data collections Y[[1]] and
#'   Y[[2]] should be connected by sharing their first data source, i.e.
#'   Y[[1]][[1]] should equal the transpose of Y[[2]][[1]].}
#'   }
#'   NOTE: The data features should have roughly zero mean and unit variance.
#'   If this is not the case, preprocessing with function
#'   \code{\link{normalizeData}} is recommended.
#' @param opts List of model options; see function \code{\link{getDefaultOpts}}.
#' @param K The number of components (i.e. latent variables). Recommended to be
#'   set somewhat higher than the expected component number, so that the sampler
#'   can determine the model complexity by shutting down excessive components.
#'   High values result in high CPU time. Default: half of the minimum of the
#'   sample size and total data dimensionality.
#' @param projection Fixed projections. Only intended for sequential prediction
#'   use via function \code{\link{sequentialGfaPrediction}}. Default: NULL.
#' @param filename A string. If provided, will save the sampling chain to this
#'   file every 100 iterations. Default "", inducing no saving.
#' @return A list containing the model parameters - in case of pairing in two
#'   modes, each element is a list of length 2; one element for each mode.
#'   For most parameters, the final posterior sample is provided to aid in
#'   initial checks; all the posterior samples should be used for model
#'   analysis. The list elements are:   
#' \item{W}{The loading matrix (final posterior sample); \eqn{D \times K}
#'          matrix.}
#' \item{X}{The latent variables (final sample); \eqn{N \times K} matrix.}
#' \item{Z}{The spike-and-slab parameters (final sample); \eqn{D \times K}
#'          matrix.}
#' \item{r}{The probability of slab in Z (final sample).}
#' \item{rz}{The probability of slab in the spike-and-slab prior of X
#'         (final sample).}
#' \item{tau}{The noise precisions (final sample); D-element vector.}
#' \item{alpha}{The precisions of the projection weights W (final sample);
#'         \eqn{D \times K} matrix.}
#' \item{beta}{The precisions of the latent variables X (final sample);
#'         \eqn{N \times K} matrix.}
#' \item{groups}{A list denoting which features belong to each data source.}
#' \item{D}{Data dimensionalities; M-element vector.}
#' \item{K}{The number of components inferred. May be less than the initial K.}
#' and the following elements:
#' \item{posterior}{the posterior samples of, by default, X, W and tau.}
#' \item{cost}{The likelihood of all the posterior samples.}
#' \item{aic}{The Akaike information criterion of all the posterior samples.}
#' \item{opts}{The options used for the GFA model.}
#' \item{conv}{An estimate of the convergence of the model's reconstruction
#'         based on Geweke diagnostic. Values significantly above 0.05 imply a
#'         non-converged model, and hence the need for a longer sampling chain.}
#' \item{time}{The CPU time (in seconds) used to sample the model.}
#'   
#' @examples
#' X <- matrix(rnorm(20*3),20,3) #Latent variables
#' W <- matrix(rnorm(30*3),30,3) #Loading matrix
#' Y <- tcrossprod(X,W) + matrix(rnorm(20*30),20,30) #Data generation
#' res <- gfa(list(Y[,1:10],Y[,11:30]), opts=getDefaultOpts(), K=5) #Inference
#' 
gfa <- function(Y, opts, K=NULL, projection=NULL, filename="") {
  if(file.exists(filename)) {
    print("Part of the sampling done already.")
    load(filename)
    print(paste("Continuing from iteration",iter))
    ptm <- proc.time()[1] - time #Starting time - time used so far
    
  } else {
    ## Initialization
    ptm <- proc.time()[1]
    if(typeof(Y[[1]])=="list") { #Pairing in two modes
      V <- length(Y)
      if(V==2)
        if(any(Y[[1]][[1]]!=t(Y[[2]][[1]]),na.rm=T))
          stop("The two sets of views must have the first view shared (transposed).")
      if(opts$verbose>0 & is.null(projection))
        print("Running GFA with data sources paired in two modes.")
    } else {
      Y <- list(Y); V <- 1
      if(opts$verbose>0) print(paste("Running GFA with",length(Y[[1]]),"data sources paired in one mode."))
    }
    
    # Store dimensionalities of the data sets
    M <- D <- N <- rep(NA,V)
    ov <- V:1 #The other view (in case of pairing in two modes)
    gr <- Yconst <- id <- list()
    for(v in 1:V) {
      M[v] <- length(Y[[v]])
      d <- sapply(Y[[v]],ncol)
      Ds <- c(1,cumsum(d)+1) ; De <- cumsum(d)
      gr[[v]] <- vector("list")
      for(m in 1:M[v]) {
        gr[[v]][[m]] <- Ds[m]:De[m]
        if (!is.null(colnames(Y[[v]][[m]]))) {
          names(gr[[v]][[m]]) <- colnames(Y[[v]][[m]])
        }
      }
      if(!is.null(names(Y[[v]]))) names(gr[[v]]) <- names(Y[[v]])
      Y[[v]] <- do.call(cbind,Y[[v]])
      D[v] <- ncol(Y[[v]])
      N[v] <- nrow(Y[[v]])
    }
    if(is.null(K) & is.null(projection)) {
      K <- round(apply(matrix(c(N,D),V),1,min)/2) #Half the components of a full initialization
      if(opts$verbose>0)
        print(paste("Initializing GFA with",paste0(K,collapse=","),"components."))
    }
    opts$K <- K
    
    if(any(K>200) & opts$verbose>0) {
      print("Using K>=200; computational cost may be high. This can be redeemed by 
             setting lower initial K, or if that does not allow automated model 
             complexity selection (does not shut down any components), with the help 
             of function informativeNoisePrior().")
    }

    
    alpha_0 <- opts$prior.alpha_0
    beta_0 <- opts$prior.beta_0
    alpha_0t <- opts$prior.alpha_0t
    beta_0t <- opts$prior.beta_0t
    #Make the alpha prior vector-valued, if it is not already
    if(typeof(alpha_0)!="list" && length(alpha_0)==1) {
      alpha_0 <- beta_0 <- list()
      for(v in 1:V) {
        if(opts$ARDW=="shared") Nrep <- 1
        if(opts$ARDW=="grouped") Nrep <- M[v]
        if(opts$ARDW=="element") Nrep <- D[v]
        alpha_0[[v]] <- rep(opts$prior.alpha_0, Nrep)
        beta_0[[v]] <- rep(opts$prior.beta_0, Nrep)
      }
    }
    if(length(opts$prior.alpha_0X)<V) {
      opts$prior.alpha_0X <- rep(opts$prior.alpha_0X, V)
      opts$prior.beta_0X <- rep(opts$prior.beta_0X, V)
    }
    #Make the tau prior vector-valued, if it is not already
    if(length(alpha_0t)==1 & typeof(alpha_0t)=="double") {
      alpha_0t <- list()
      for(v in 1:V) alpha_0t[[v]] <- rep(opts$prior.alpha_0t, D[v])
    }
    if(length(beta_0t)==1 & typeof(beta_0t)=="double") {
      beta_0t <- list()
      for(v in 1:V) beta_0t[[v]] <- rep(opts$prior.beta_0t, D[v])
    }
    
    const <- rep(NA,V)
    Yconst <- id <- b_tau <- Z <- alpha <- covZ <- r <- na.inds <- list()
    
    if (!is.null(projection)) { #Projection given as an argument.
      W <- projection$W; rz <- projection$rz
      tau <- projection$tau; beta <- projection$beta
      projection.fixed = TRUE
      for(v in 1:V)
        K[v] = ncol(W[[v]])
      
    } else {
      projection.fixed = FALSE
      W <- beta <- rz <- tau <- prediction <- list()
    }
    X <- XX <- prediction <- list()
    
    for(v in 1:V) {
      # Some constants for speeding up the computation
      #const <- - N[v]*Ds/2*log(2*pi) # Constant factors for the lower bound
      Yconst[[v]] <- colSums(Y[[v]]^2) #unlist(lapply(Y[[v]],function(x){sum(x^2)}))
      id[[v]] <- rep(1,K[v])              # Vector of ones for fast matrix calculations
      
      ##
      ## Initialize the model randomly
      ##
      # Other initializations could be done, but overdispersed random initialization is quite good.
      
      # Latent variables
      X[[v]] <- matrix(rnorm(N[v]*K[v],0,1),N[v],K[v])
      X[[v]] <- scale(X[[v]])
      covZ[[v]] <- diag(1,K[v]) # Covariance
      XX[[v]] <- crossprod(X[[v]]) # Second moments
      
      if(!opts$normalLatents & projection.fixed) {
        beta[[v]] <- matrix(colMeans(beta[[v]]),N,K,byrow=T) #Initialize as the mean over samples
      }else{
        beta[[v]] <- matrix(1,N[v],K[v])
      }
    }
    
    i.pred <- 0
    if (opts$iter.saved>0 & opts$iter.burnin<opts$iter.max) {
      mod.saved <- ceiling((opts$iter.max-opts$iter.burnin)/opts$iter.saved)
      #Save every 'mod.saved'th Gibbs sample.
      S <- floor((opts$iter.max-opts$iter.burnin)/mod.saved) # Number of saved Gibbs samples.
    }
    
    if (projection.fixed) { # Projection given as an argument to the sampler
      tmpW <- WtW <- WtWdiag <- list()
      for(v in 1:V) {
        obs = NULL # variable indices of the source views
        for (mi in 1:length(opts$prediction[[v]])) { # Go through all views.
          if (!opts$prediction[[v]][mi]) {
            obs = c(obs, gr[[v]][[mi]])
          }
        }
        covZ[[v]] <- diag(1,K[v]) + crossprod(W[[v]][obs,]*sqrt(tau[[v]][obs]))
        eS <- eigen(covZ[[v]])
        tmpW[[v]] <- eS$vectors*outer(id[[v]],1/sqrt(eS$values))
        covZ[[v]] <- tcrossprod(tmpW[[v]])
        if(!opts$normalLatents){
          WtW[[v]] <- crossprod(W[[v]][obs,]*sqrt(tau[[v]][obs]))
          WtWdiag[[v]] <- diag(WtW[[v]])
          diag(WtW[[v]]) <- 0
          if (is.null(dim(rz[[v]]))) {
            if (length(rz[[v]])==K[v]) {
              rz[[v]] = matrix(data=rz[[v]], nrow=N[v], ncol=K[v], byrow=TRUE)
            } else {
              stop(paste0("rz[[",v,"]] not of required length"))
            }
          }
        }
        if (any(opts$prediction[[v]])) {
          prediction[[v]] <- vector(mode="list", length=length(gr[[v]]))
          Y.true <- Y[[v]]
          for (m in which(opts$prediction[[v]])) {
            prediction[[v]][[m]] = matrix(0,N[v],length(gr[[v]][[m]]))
            Y[[v]][,gr[[v]][[m]]] <- 0
          }
        }
        #Added things
        const[v] <- 0
        for(m in which(opts$prediction[[v]]==F))
          const[v] <- const[v] - N[v]*length(gr[[v]][[m]])*log(2*pi)/2
        cost <- aic <- rep(sum(const),opts$iter.max)
        b_tau[[v]] <- beta_0t[[v]]
        #
      }
      
    } else { #Non-fixed projections
      for(v in 1:V) {
        WtW <- matrix(0,K[v],K[v])
        WtWdiag <- rep(0,K[v])
        
        const[v] <- - N[v]*D[v]*log(2*pi)/2
        tau[[v]] <- rep(opts$init.tau,D[v]) # The mean noise precisions
        b_tau[[v]] <- beta_0t[[v]]
        
        W[[v]] <- matrix(0,D[v],K[v])
        Z[[v]] <- matrix(1,D[v],K[v])
        covW <- diag(1,K[v])
        
        alpha[[v]] <- matrix(1,D[v],K[v])
        
        cost <- aic <- rep(sum(const), opts$iter.max)
        if(v==2) cost <- cost + N[v]*length(gr[[1]][[1]])*log(2*pi)/2
        r[[v]] <- matrix(0.5,D[v],K[v])
        rz[[v]] <- matrix(0.5,N[v],K[v])
      }
    }
    
    ## Missing values (NAs) in the data
    missingValues <- FALSE
    na.inds <- list()
    for(v in 1:V) {
      na.inds[[v]] <- which(is.na(Y[[v]]),arr.ind=T)
      if(V==2) na.inds[[2+v]] <- which(is.na(Y[[v]][,gr[[v]][[1]]]),arr.ind=T) #Shared data view
      if(length(na.inds[[v]])>0 & !projection.fixed) {
        missingValues <- TRUE
        if(opts$verbose>0)
          print("Missing Values detected, prediction using EM type approximation")
        #Update alpha_0t to take into account the observed data size
        if(V==1) a_tau <- alpha_0t[[1]] + colSums(!is.na(Y[[1]]))/2
        for(m in 1:M[v])
          alpha_0t[[v]][gr[[v]][[m]]] <- alpha_0t[[v]][gr[[v]][[m]]] + sum(!is.na(Y[[v]][,gr[[v]][[m]]]))/2
        
      } else {
        if(V==1) a_tau <- alpha_0t[[1]] + N/2 #Elementwise noise
        for(m in 1:M[v]) #View-wise noise
          alpha_0t[[v]][gr[[v]][[m]]] <- alpha_0t[[v]][gr[[v]][[m]]] + N[v]*length(gr[[v]][[m]])/2
      }
      for(m in 1:M[v])
        alpha_0t[[v]][m] <- alpha_0t[[v]][[gr[[v]][[m]][1]]]
      alpha_0t[[v]] <- alpha_0t[[v]][1:M[v]]
    }
    if(missingValues) {
      na.inds$Y <- list()
      for(v in 1:V)
        na.inds$Y[[v]] <- (!is.na(Y[[v]]))*1
      #Many matrix inversions are done - avoid duplicates by ordering
      na.inds$Nlist <- list()
      for(v in 1:V)
        na.inds$Nlist[[v]] <- hclust(dist(na.inds$Y[[v]]))$order
    }
    
    #Which X and W posterior samples to save in case of convergenceCheck?
    if(opts$convergenceCheck & opts$iter.saved>=8 & !projection.fixed) {
      ps <- floor(opts$iter.saved/4)
      start <- 1:ps
      end <- (-ps+1):0+opts$iter.saved
    } else {
      start <- end <- c()
    }
    
    ## W initialization
    if(!projection.fixed)
      for(v in 1:V) {
        covW <- diag(1,K[v]) + opts$init.tau*XX[[v]]
        eS <- eigen(covW)
        tmp <- eS$vectors*outer(id[[v]],1/sqrt(eS$values))
        covW <- tcrossprod(tmp)
        
        estimW = matrix(0,D[v],K[v])
        tmpY <- Y[[v]]
        if(v==2) { #Some variance explained by the first mode already
          tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
            t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        for(k in 1:K[v]){
          tmpY[na.inds[[v]]] <- 0
          estimW[,k] <- crossprod(tmpY,X[[v]][,k])
        }
        W[[v]] <- estimW%*%covW*opts$init.tau + matrix(rnorm(D[v]*K[v]),D[v],K[v])%*%tmp
      }
    
    posterior <- NULL
    iter <- 1
  }
  
  ##
  ## The main loop - Gibbs updates
  ##
  for(iter in iter:opts$iter.max) {
    if (!projection.fixed) {
      ##
      ## Sample the projections W
      ##
      for(v in 1:V) {
        XXdiag <- diag(XX[[v]])
        diag(XX[[v]]) <- 0
        if(missingValues) xx <- X[[v]]^2
        if(V==2) {
          tmpY <- Y[[v]]
          tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
            t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        for(k in 1:K[v]){
          if(missingValues) {
            if(opts$imputation=="Bayesian") { #||x||^2, based on observations only
              XXdiag <- apply(na.inds$Y[[v]], 2, function(y) sum(xx[,k]*y))
              lambda <- tau[[v]]*XXdiag + alpha[[v]][,k]
            } else if(opts$imputation=="conservative") {
              lambda <- tau[[v]]*XXdiag[k] + alpha[[v]][,k]
            } else {
              stop("Ill-defined option 'imputation'")
            }
            ss <- tcrossprod(X[[v]][,-k],W[[v]][,-k])
            if(V==1) tmp <- Y[[v]]-ss
            if(V==2) tmp <- tmpY-ss
            tmp[na.inds[[v]]] <- 0 #YX based on observations only
            mu <- tau[[v]]/lambda*crossprod(tmp,X[[v]][,k])
          } else {
            lambda <- tau[[v]]*XXdiag[k] + alpha[[v]][,k]
            if(V==1) mu <- tau[[v]]/lambda*c(crossprod(Y[[v]],X[[v]][,k]) - W[[v]]%*%XX[[v]][k,])
            if(V==2) mu <- tau[[v]]/lambda*c(crossprod(tmpY,X[[v]][,k]) - W[[v]]%*%XX[[v]][k,])
          }
          W[[v]][,k] <- mu + rnorm(D[v])*lambda^(-0.5)
          
          ##
          ## Updating the spike and slab prior, parameter Z (H in publications)
          ##
          if(iter > opts$sampleZ) {
            if(opts$spikeW!="group") {
              zone <- 0.5*( log(alpha[[v]][,k]) - log(lambda) + lambda*mu^2) +
                log(r[[v]][,k]) - log(1-r[[v]][,k])
              zone <- 1/(1+exp(-zone))
              zone <- as.double(runif(D[v]) < zone)
              Z[[v]][,k] <- zone
            } else {
              zone <- 0.5*(log(alpha[[v]][,k])-log(lambda) + lambda*mu^2)
              for(m in 1:M[v]) {
                logpr <- sum(zone[gr[[v]][[m]]]) + log(r[[v]][gr[[v]][[m]][1],k]) -
                  log( 1-r[[v]][gr[[v]][[m]][1],k] )
                logpr <- 1/(1+exp(-logpr))
                Z[[v]][gr[[v]][[m]],k] <- as.double(runif(1)<logpr)
              }
            }
          }
          W[[v]][,k] <- W[[v]][,k]*Z[[v]][,k]
        }
      }
    }
    
    ## Check if components can be removed.
    # If Gibbs samples are going to be saved, components are dropped only during the burn-in.
    # If the projection has been given as an argument (e.g. for prediction), components will not be dropped.
    if ((iter<=opts$iter.burnin | opts$iter.saved==0) & !projection.fixed) {
      for(v in 1:V) {
        keep <- which(colSums(Z[[v]])!=0 & colSums(abs(X[[v]]))>0)
        if(length(keep)==0) {
          print(paste0("All components shut down (mode ",v,"), returning a NULL model."))
          return(list())
        }
        if(length(keep)!=K[v]){
          K[v] <- length(keep)
          id[[v]] <- rep(1,K[v])
          alpha[[v]] <- alpha[[v]][,keep,drop=F]
          X[[v]] <- X[[v]][,keep,drop=F]
          W[[v]] <- W[[v]][,keep,drop=F]
          Z[[v]] <- Z[[v]][,keep,drop=F]
          r[[v]] <- r[[v]][,keep,drop=F]
          beta[[v]] <- beta[[v]][,keep,drop=F]
          rz[[v]] <- rz[[v]][,keep,drop=F]
        }
        mm <- ""; if(V==2) mm <- paste0(" (mode ",v,")")
        if(iter==opts$iter.burnin & K[v]<=opts$K[v]*.4 & opts$K[v]>10)
          warning(paste0("Over 60% of the initial components shut down",mm,". Excessive amount of
  initial components may result in a bad initialization for the model parameters.
  Rerun the model with K=",length(keep)+5,"."))
      }
    }
    
    ##
    ## sample r (the probability of slab)
    ##
    for(v in 1:V) {
      if(iter>opts$sampleZ & !projection.fixed){
        if(opts$spikeW=="shared"){
          zm <- colSums(Z[[v]])
          for(k in 1:K[v])
            r[[v]][,k] <- (zm[k]+opts$prior.beta[1])/(D[v]+sum(opts$prior.beta))
        }
        if(opts$spikeW=="element"){
          for(m in 1:M[v]){
            zm <- colSums(Z[[v]][gr[[v]][[m]],,drop=FALSE])
            for(k in 1:K[v])
              r[[v]][gr[[v]][[m]],k] <- (zm[k]+opts$prior.beta[1])/(length(gr[[v]][[m]])+sum(opts$prior.beta))
          }
        }
        if(opts$spikeW=="group"){
          zm <- rep(0,K[v])
          for(m in 1:M[v]){
            zm <- zm + Z[[v]][gr[[v]][[m]][1],]
          }
          for(m in 1:M[v]){
            for(k in 1:K[v]){
              r[[v]][gr[[v]][[m]],k] <- (zm[k]+opts$prior.beta[1])/(M[v]+sum(opts$prior.beta))
            }
          }
        }
      }
      
      ## 
      ## Sample the latent variables X
      ##
      if(iter>opts$sampleZ & !opts$normalLatents) {
        # Spike-and-slab prior for X (biclustering)
        if (!projection.fixed) {
          WtW <- crossprod(W[[v]]*sqrt(tau[[v]]))
          WtWdiag <- diag(WtW)
          diag(WtW) <- 0
          if(missingValues) ww <- W[[v]]^2*tau[[v]]
        }
        if(V==2) {
          tmpY <- Y[[v]]
          tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
            t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        for(k in 1:K[v]){
          if(missingValues) {
            if(opts$imputation=="Bayesian") { #||w||^2, based on observations only
              WtWdiag <- apply(na.inds$Y[[v]], 1, function(y) sum(ww[,k]*y))
              lambda <- WtWdiag + beta[[v]][,k]
            } else if(opts$imputation=="conservative") {
              lambda <- WtWdiag[k] + beta[[v]][,k]
            }
            ss <- tcrossprod(X[[v]][,-k],W[[v]][,-k])
            if(V==1) tmp <- Y[[v]]-ss
            if(V==2) tmp <- tmpY-ss
            tmp[na.inds[[v]]] <- 0 #YW based on observations only
            mu <- (tmp%*%(W[[v]][,k]*tau[[v]]))/lambda
            
          } else {
            if(projection.fixed) {
              lambda <- WtWdiag[[v]][k] + beta[[v]][,k]
              if(V==1) mu <- (Y[[v]]%*%(W[[v]][,k]*tau[[v]]) - X[[v]]%*%WtW[[v]][k,])/lambda
              if(V==2) mu <- (tmpY%*%(W[[v]][,k]*tau[[v]]) - X[[v]]%*%WtW[[v]][k,])/lambda
            } else {
              lambda <- WtWdiag[k] + beta[[v]][,k]
              if(V==1) mu <- (Y[[v]]%*%(W[[v]][,k]*tau[[v]]) - X[[v]]%*%WtW[k,])/lambda
              if(V==2) mu <- (tmpY%*%(W[[v]][,k]*tau[[v]]) - X[[v]]%*%WtW[k,])/lambda
            }
          }
          X[[v]][,k] <- mu + rnorm(N[v])*lambda^(-0.5)
          
          #Binary sparsity sampling
          zone <- 0.5*( log(beta[[v]][,k]) - log(lambda) + lambda*mu^2) +
            log(rz[[v]][,k]) - log(1-rz[[v]][,k])
          zone <- 1/(1+exp(-zone))
          zone <- as.double(runif(N[v]) < zone)
          X[[v]][,k] <- X[[v]][,k]*zone
          if (!projection.fixed)
            rz[[v]][,k] <- (sum(zone)+opts$prior.betaX[1])/(N[v]+sum(opts$prior.betaX))
        }
      } else {
        #Normal prior for X
        if (!projection.fixed) {
          covZ[[v]] <- diag(1,K[v]) + crossprod(W[[v]]*sqrt(tau[[v]]))
          eS <- eigen(covZ[[v]])
          tmp <- eS$vectors*outer(id[[v]],1/sqrt(eS$values))
          covZ[[v]] <- tcrossprod(tmp)
        }
        if(V==2) {
          tmpY <- Y[[v]]
          tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
            t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        if(missingValues) {
          tmpY[na.inds[[v]]] <- 0
          X[[v]] <- tmpY%*%(W[[v]]*tau[[v]])
          if(opts$imputation=="Bayesian") {
            #Missing values cause each sample to have a different covariance matrix
            #NOTE/TODO: Sampling components independently would be faster, but converge slower
            nn <- na.inds$Nlist[[v]] #Samples ordered to avoid duplicate covX computations
            compcovX <- TRUE
            for(i in 1:N[v]) {
              if(all(na.inds$Y[[v]][nn[i],]==1)) { #All features observed for this sample
                X[[v]][nn[i],] <- X[[v]][nn[i],,drop=F]%*%covZ[[v]] + matrix(rnorm(K[v]),1,K[v])%*%tmp
              } else {
                if(i>1) { #If different features observed, compute a new inverse (covX)
                  if(any(na.inds$Y[[v]][nn[i],]!=na.inds$Y[[v]][nn[i-1],])) compcovX <- TRUE
                  else compcovX <- FALSE
                }
                if(compcovX) {
                  covX <- diag(1,K[v]) + crossprod(W[[v]]*na.inds$Y[[v]][nn[i],]*sqrt(tau[[v]]))
                  eS <- eigen(covX) #Main source of computation!
                  tmp2 <- eS$vectors*outer(id[[v]],1/sqrt(eS$values))
                  covX <- tcrossprod(tmp2)
                }
                X[[v]][nn[i],] <- X[[v]][nn[i],,drop=F]%*%covX + matrix(rnorm(K[v]),1,K[v])%*%tmp2
              }
            }
          } else if(opts$imputation=="conservative") {
            X[[v]] <- X[[v]]%*%covZ[[v]] + matrix(rnorm(N[v]*K[v]),N[v],K[v])%*%tmp
          }
          
        } else {
          if(V==1) X[[v]] <- Y[[v]]%*%(W[[v]]*tau[[v]])
          if(V==2) X[[v]] <- tmpY%*%(W[[v]]*tau[[v]])
          X[[v]] <- X[[v]]%*%covZ[[v]]
          if(projection.fixed)
            X[[v]] <- X[[v]] + matrix(rnorm(N[v]*K[v]),N[v],K[v])%*%tmpW[[v]]
          else
            X[[v]] <- X[[v]] + matrix(rnorm(N[v]*K[v]),N[v],K[v])%*%tmp
        }
      }
      XX[[v]] <- crossprod(X[[v]])
      
      ##
      ## Update alpha and beta, the ARD parameters
      ##
      if(!opts$normalLatents) {
        # Spike-and-slab prior for X, including normal precision beta (ARD)
        if(opts$ARDLatent=="element") {
          beta[[v]] <- X[[v]]^2/2 + opts$prior.beta_0X[v]
          keep <- which(beta[[v]]!=0)
          for(n in keep) beta[[v]][n] <- rgamma(1,shape=opts$prior.alpha_0X[v]+0.5,rate=beta[[v]][n]) +1e-7
          beta[[v]][-keep] <- rgamma(length(beta[[v]][-keep]),
                                     shape=opts$prior.alpha_0X[v],rate=opts$prior.beta_0X) + 1e-7
        }
        if(opts$ARDLatent == "shared" & !projection.fixed) {
          tmpxx <- colSums(X[[v]]!=0)/2 + opts$prior.alpha_0X[v]
          xx <- colSums(X[[v]]^2)/2 + opts$prior.beta_0X[v]
          for(k in 1:K[v])
            beta[[v]][,k] <- rgamma(1,shape=tmpxx[k],rate=xx[k]) + 1e-7
        }
      }
      if (!projection.fixed) {
        if(opts$ARDW=="shared"){
          ww <- colSums(W[[v]]^2)/2 + beta_0[[v]]
          tmpz <- colSums(Z[[v]])/2 + alpha_0[[v]]
          for(k in 1:K[v])
            alpha[[v]][,k] <- rgamma(1,shape=tmpz[k],rate=ww[k]) + 1e-7
        }
        if(opts$ARDW=="grouped") {
          for(m in 1:M[v]){
            ww <- colSums(W[[v]][gr[[v]][[m]],,drop=FALSE]^2)/2 + beta_0[[v]][m]
            tmpz <- colSums(Z[[v]][gr[[v]][[m]],,drop=FALSE])/2 + alpha_0[[v]][m]
            for(k in 1:K[v])
              alpha[[v]][gr[[v]][[m]],k] <- rgamma(1,shape=tmpz[k],rate=ww[k]) + 1e-7
          }
        }
        if(opts$ARDW=="element") {
          alpha[[v]] <- W[[v]]^2/2 + beta_0[[v]]
          keep <- which(alpha[[v]]!=0)
          for(n in keep) alpha[[v]][n] <- rgamma(1,shape=alpha_0[[v]][n]+0.5,rate=alpha[[v]][n]) +1e-7
          alpha[[v]][-keep] <- rgamma( length(alpha[[v]][-keep]),shape=alpha_0[[v]][n],
                                       rate=beta_0[[v]][n]) + 1e-7
        }
        
        ##
        ## Update tau, the noise precisions
        ##
        tmpid <- 1:length(tau[[v]])
        if(V==2) tmpid <- tmpid[-gr[[v]][[1]]]
        
        XW <- tcrossprod(X[[v]],W[[v]])
        if(missingValues) {
          if(iter==opts$iter.max && !projection.fixed && missingValues) {
            id <- na.inds[[v]][,1] + (na.inds[[v]][,2]-1)*nrow(XW)
            if(length(id)>0 && (max(abs(XW)) > 2*max(abs(XW[-id])))) {
              vv <- if(V==1) "" else paste0(" (mode ",v,")")
              warning(paste0("Reconstruction of missing values resulted in larger scale than observed data",vv,".
  This may be a sign of overfitting to sparsely observed samples/columns. Consider
  using informative Gamma-priors for the parameters X and W to avoid this."))
            }
          }
          XW <- Y[[v]] - XW
          XW[na.inds[[v]]] <- 0
          b_tau[[v]][tmpid] <- colSums(XW[,tmpid]^2)/2
        } else{
          b_tau[[v]][tmpid] <- colSums((Y[[v]][,tmpid]-XW[,tmpid])^2)/2
        }
        if(opts$tauGrouped) {
          for(m in V:M[v]) { #For pairing in two modes (V==2), update m=1 later
            tau[[v]][gr[[v]][[m]]] <- rgamma(1, shape = alpha_0t[[v]][m],
                                             rate = beta_0t[[v]][gr[[v]][[m]][1]] +
                                               sum(b_tau[[v]][gr[[v]][[m]]]))
          }
        } else {
          if(V==2)
            stop("Tau needs to be grouped: otherwise the noise prior for the shared view is ill-defined")
          for(d in 1:D[[1]])
            tau[[1]][d] <- rgamma(1,shape=a_tau[d], rate= beta_0t[[1]][d]+b_tau[[1]][d])
        }
        # calculate likelihood.
        cost[iter] <- cost[iter] + N[v]*sum(log(tau[[v]][tmpid]))/2 -
          crossprod(b_tau[[v]][tmpid],tau[[v]][tmpid])
      }
      
    } #Loop over v (data collections) ends
    
    #Update tau for the paired view (pairing in two modes)
    if(!projection.fixed & V==2) {
      XW <- tcrossprod(X[[1]],W[[1]][gr[[1]][[1]],]) + t(tcrossprod(X[[2]],W[[2]][gr[[2]][[1]],]))
      if(missingValues && length(na.inds[[3]])) {
        XW <- Y[[1]][[1]][,gr[[1]][[1]]] - XW
        XW[na.inds[[3]]] <- 0
        b_tau[[1]][gr[[1]][[1]]] <- colSums(XW[,gr[[1]][[1]]]^2)/2
      } else {
        b_tau[[1]][gr[[v]][[1]]] <- colSums((Y[[1]][,gr[[1]][[1]]]-XW)^2)/2
      }
      
      if(opts$tauGrouped){
        tau[[1]][gr[[1]][[1]]] <- rgamma(1, shape = alpha_0t[[1]][1],
                                         rate = beta_0t[[1]][1] + sum(b_tau[[1]][gr[[1]][[1]]]))
        tau[[2]][gr[[2]][[1]]] <- tau[[1]][gr[[1]][[1]][1]]
      }else{
        stop("Tau needs to be grouped: otherwise the noise prior for the shared view is ill-defined")
      }
      cost[iter] <- cost[iter] + N[1]*sum(log(tau[[1]][gr[[1]][[1]]]))/2 -
        crossprod(b_tau[[1]][gr[[1]][[1]]],tau[[1]][gr[[1]][[1]]])
    }
    aic[iter] <- 2*cost[iter] #Akaike information criterion
    for(v in 1:V) aic[iter] <- aic[iter] - (D[v]*(K[v]+1)-K[v]*(K[v]-1)/2)*2
    
    
    if(filename != "" & iter%%100==0) { #Every 100 iterations, save the sampling
      time <- proc.time()[1]-ptm
      save(list=ls(),file=filename)
      if(opts$verbose>0)
        print(paste0("Iter ", iter, ", saved chain to '",filename,"'"))
    }
    
    # Calculate likelihood of observed views
    if(projection.fixed) {
      for(v in 1:V) {
        tautmp <- ( Yconst[[v]] + rowSums((W[[v]]%*%XX[[v]])*W[[v]]) -
                      2*rowSums(crossprod(Y[[v]],X[[v]])*W[[v]]) )/2
        for(m in which(opts$prediction[[v]]==F)) {
          cost[iter] <- cost[iter] + N[v]*sum(log(tautmp[gr[[v]][[m]]]))/2 -
            crossprod(b_tau[[v]][gr[[v]][[m]]],tautmp[gr[[v]][[m]]])
        }
      }
    }
    
    ##
    ## Prediction and collection of Gibbs samples
    ##
    for(v in 1:V) {
      if (any(opts$prediction[[v]])) { ## Prediction
        if (iter%%10==0 & opts$verbose>1) {
          print(paste0("Predicting: ",iter,"/",opts$iter.max))
        }
        if (iter>opts$iter.burnin & ((iter-opts$iter.burnin)%%mod.saved)==0) {
          i.pred <- i.pred+1/V
          for (m in which(opts$prediction[[v]])) { # Go through the views that will be predicted.
            prediction[[v]][[m]] <- prediction[[v]][[m]] + tcrossprod(X[[v]], W[[v]][gr[[v]][[m]],])
          }
        }
        
      } else if(!any(unlist(opts$prediction))) {
        if (((iter%%10==0 & opts$verbose>1) | (iter%%100==0 & opts$verbose>0)) & v==1) {
          print(paste0("Learning: ",iter,"/",opts$iter.max," - K=",paste0(K,collapse=",")," - ",Sys.time()))
        }
        if (opts$iter.saved>0) { ## Collection of Gibbs samples
          if (iter>opts$iter.burnin & ((iter-opts$iter.burnin)%%mod.saved)==0) { ## Save the Gibbs sample.
            if (iter-opts$iter.burnin==mod.saved) { # Initialize the list containing saved Gibbs samples.
              add <- ""; if(V==2) add <- paste0(" (mode ",v,")")
              if(K[v]==opts$K[v]) {
                if(opts$K[v]<min(N[v],D[v])) {
                  warning(paste0("Burn-in period finished, but no components have been shut down",add,",
  preventing automatic model complexity selection. Consider higher inital K, if
  computational resources permit (up to K=",min(N[v],D[v]),", default K=",
                                 round(min(N[v],D[v])/2),").
  Otherwise, refer to function informativeNoisePrior() for defining
  an informative prior for residual noise, which should enable
  automated model complexity selection."))
                } else {
                  warning(paste0("Burn-in period finished, but no components have been shut down",add,".
  The model seems to have a hard time determining suitable complexity 
  for this specific data set. Refer to function informativeNoisePrior() 
  for defining an informative prior for residual noise, which should enable
  automated model complexity selection."))
                }

              }
              if(v==1)
                posterior <- list(r=list(),rz=list(),beta=list()) # List containing saved Gibbs samples
              if ((opts$save.posterior$W | opts$convergenceCheck) & !any(opts$prediction[[v]])) {
                if(v==1) {posterior$W <- list()}
                posterior$W[[v]] <- array(dim=c(S, dim(W[[v]]))) # SxDxK[v]
                colnames(posterior$W[[v]]) <- colnames(Y[[v]])
              }
              posterior$rz[[v]] <- array(dim=c(S, ncol(rz[[v]])))
              # SxK[v] - Same 'rz' for all samples. Thus, save only a vector per iteration.
              posterior$r[[v]] <- array(dim=c(S, M[v],K[v]))
              # SxK[v]xM[v] - Same 'r' for all variables of a view.
              posterior$tau[[v]] <- array(dim=c(S, length(tau[[v]]))) # SxD
              colnames(posterior$tau[[v]]) = colnames(Y[[v]])
              if ((opts$save.posterior$X | opts$convergenceCheck) & !any(opts$prediction[[v]])) {
                if(v==1) {posterior$X <- list()}
                posterior$X[[v]] <- array(dim=c(S, dim(X[[v]]))) # SxN[v]xK[v]
                colnames(posterior$X[[v]]) <- rownames(Y[[v]])
              }
              posterior$beta[[v]] = array(dim=c(S,dim(beta[[v]])))
              if(v==1) {gr.start <- list()}
              gr.start[[v]] = vector(mode="integer", length=length(gr[[v]]))
              for (m in 1:length(gr[[v]])) { # Find the first index of each view.
                gr.start[[v]][m] = gr[[v]][[m]][1]
              }
            }
            
            s <- (iter-opts$iter.burnin)/mod.saved
            if (!projection.fixed) {
              if (opts$save.posterior$W | (opts$convergenceCheck & s%in%c(start,end)))
                posterior$W[[v]][s,,] <- W[[v]]
              posterior$tau[[v]][s,] <- tau[[v]]
              posterior$rz[[v]][s,] <- rz[[v]][1,] #'rz' is identical for all samples
              posterior$r[[v]][s,,] <- r[[v]][gr.start[[v]],]
              #'r' is identical for all variables within a view. Thus, save only these.
              posterior$beta[[v]][s,,] <- beta[[v]]
            }
            if (opts$save.posterior$X | (opts$convergenceCheck & s%in%c(start,end)))
              posterior$X[[v]][s,,] <- X[[v]]
          }
        }
      }
    }
  } ## Gibbs sampling done
  
  if(opts$convergenceCheck & opts$iter.saved>=8 & !projection.fixed) {
    #Estimate the convergence of the data reconstruction, based on the Geweke diagnostic
    if(opts$verbose>0)
      print("Starting convergence check")
    conv <- 0
    for(v in 1:V) {
      for(i in 1:N[v]) {
        if(opts$verbose>1) {
          if(i%%10==0)
            cat(".")
          if(i%%100==0)
            print(paste0(i,"/",N[v]))
        }
        y <- matrix(NA,0,ncol(Y[[v]]))
        for(ps in c(start,end)) {
          y <- rbind(y, tcrossprod(posterior$X[[v]][ps,i,],posterior$W[[v]][ps,,]))
        }
        foo <- rep(NA,ncol(Y[[v]]))
        for(j in 1:ncol(Y[[v]])) {
          if(sd(y[start,j])>1e-10 & sd(y[-start,j])>1e-10) {
            foo[j] <- t.test(y[start,j],y[-start,j])$p.value
          } else { #Essentially constant reconstruction
            if(abs(mean(y[start,j])-mean(y[-start,j]))<1e-10)
              foo[j] <- 1 #Same constant
            else
              foo[j] <- 0 #Different constant
          }
        }
        conv <- conv + sum(foo<0.05)/length(foo)/N[v]/V #check how many values are below 0.05
      }
      if(!opts$save.posterior$X)
        posterior$X[[v]] <- NULL
      if(!opts$save.posterior$W)
        posterior$W[[v]] <- NULL
      gc()
    }
    if(opts$verbose>0) {
      print(paste("Convergence diagnostic:",round(conv,4)))
      print("Values significantly greater than 0.05 imply a non-converged model.")
    }
    
  } else {
    conv <- NA
  }
  
  if(filename!="" & opts$iter.max>=10)
    file.remove(filename) #No need for temporary storage any more
  

  if(projection.fixed) {
    for(v in 1:V) {
      if(any(opts$prediction[[v]])) {
        for (m in which(opts$prediction[[v]])) # Go through the views that will be predicted.
          prediction[[v]][[m]] <- prediction[[v]][[m]]/i.pred
        prediction[[v]]$cost <- cost
      }
    }
    return(prediction)
    
  } else {
    for(v in 1:V) {
      d1 <- unlist(lapply(gr[[v]],function(x){x[1]}))
      if(opts$ARDW=="grouped")
        alpha[[v]] <- alpha[[v]][d1,]
      if(opts$spikeW=="group")
        Z[[v]] <- Z[[v]][d1,]
    }
    for(v in 1:V) {
      rownames(X[[v]]) <- rownames(Y[[v]])
      rownames(W[[v]]) <- colnames(Y[[v]])
    }
    if(V==1) { #List not needed for one set of data sources (paired in one mode)
      params <- c("W","X","Z","r","rz","tau","alpha","beta","gr")
      if(!is.null(posterior)) params <- c(params, paste0("posterior$",names(posterior)))
      for(param in params)
        eval(parse(text=paste0(param," <- ",param,"[[1]]")))
      D <- D[1]; K <- K[1]
    }
    
    time <- proc.time()[1] - ptm
    
    return(list(W=W, X=X, Z=Z, r=r, rz=rz, tau=tau, alpha=alpha, beta=beta, groups=gr, D=D,
                K=K, cost=cost, aic=aic, posterior=posterior, opts=opts, conv=conv, time=time))
  }
}


#' A function for generating the default priors of GFA model
#'
#' \code{getDefaultOpts} returns the priors of GFA
#' 
#' This function returns options defining the model's high-level structure
#' (sparsity priors) and the model's hyperparameters, defining uninformative
#' priors. We recommend keeping these as provided, with one exception: if the
#' uninformative prior of the noise residual (tau) seems to result in an overly
#' complex model (no components become shut down even if the initial K is set
#' high), risking overfitting, we recommend using
#' function \code{\link{informativeNoisePrior}} to adjust the priors.
#' 
#' @param bicluster Use binary sparsity priors in both the data modes? If FALSE
#'        (default), the components will be dense in the data sources, but
#'        group-sparse, i.e., each component is active in a (potentially
#'        different) subset of the data sources. If TRUE, binary sparsity is
#'        inferred for each data sample and feature, resulting in each component
#'        to be interpretable as a multi-source bicluster.
#' @return A list with the following model options:
#' \item{tauGrouped}{If TRUE (default), data views have separate noise
#'         precisions, otherwise each feature has.}
#' \item{normalLatents}{If TRUE, X will have a normal prior; if FALSE
#'  X, will have a spike-and-slab prior.}
#' \item{spikeW}{Sparsity prior of W. "group"=group sparsity, "element"=
#'  element-wise sparsity with shared hyperparameter across views, "shared"=
#'  element-wise sparsity with no grouping.}
#' \item{ARDW}{ARD prior type for W, determining the scale of the inferred
#'  components. "shared"=same scale for all the data sources, "grouped"
#'  (default)= separate scale for each data source, "element"=separate scale
#'  for each feature.}
#' \item{ARDLatent}{ARD prior type for X: "shared" (default)=shared scale for
#'   all the samples, "element"=separate scale for each sample.}
#' \item{imputation}{Missing value imputation type: "Bayesian" (default)=proper 
#'   Bayesian handling of missing values. "conservative"=missing values result
#'   in smaller parameter scale, which can be useful if tricky missing value
#'   structure causes exaggerated imputed values with the default setting
#'   (which can also be dealt with informative priors for alpha and beta).}
#' \item{iter.max}{The total number of Gibbs sampling steps (default 5000).}
#' \item{iter.saved}{The number of saved posterior samples (default 100).}
#' \item{iter.burnin}{The number of burn-in samples (default 2500).}
#' \item{init.tau}{The initial noise precision. High values imply initializing
#'   the model with an adequate number of components. Default 1000.}
#' \item{sampleZ}{When to start sampling spike and slab parameters (default:
#'  Gibbs sample 1).}
#' \item{prior.alpha_0t}{The shape parameter of tau's prior (default 10).}
#' \item{prior.beta_0t}{The rate parameter of tau's prior (default 10).}
#' \item{prior.alpha_0}{The shape parameter of alpha's prior (default 10).}
#' \item{prior.beta_0}{The rate parameter of alpha's prior (default 01).}
#' \item{prior.alpha_0X}{The shape parameter of beta's prior (default 10).}
#' \item{prior.beta_0X}{The rate parameter of beta's prior (default 1).}
#' \item{prior.beta}{Bernoulli prior for the spike-and-slab prior of W (counts
#'  for 1s and 0s; default c(1,1)).}
#' \item{prior.betaX}{Bernoulli prior for the possible spike-and-slab prior of
#'   X (default c(1,1)).}
#' \item{verbose}{The verbosity level. 0=no printing, 1=moderate printing,
#'  2=maximal printing (default 1).}
#' \item{convergenceCheck}{Check for the convergence of the data reconstruction,
#'  based on the Geweke diagnostic (default FALSE).}
#' \item{save.posterior}{A list determining which parameters' posterior samples are
#'   saved (default: X, W and tau).}
#'
#' @examples
#' #Given pre-specified data collection Y and component number K
#' opts <- getDefaultOpts(bicluster=FALSE)
#' opts$normalLatents <- FALSE #Binary sparsity for each sample and data source
#'  \dontrun{model <- gfa(Y,opts,K)}
getDefaultOpts <- function(bicluster=FALSE) {
  normalLatents <- TRUE
  spikeW <- "group"
  if(bicluster) {
    normalLatents <- FALSE
    spikeW <- "element"
  }
  
  return(list(tauGrouped=TRUE,normalLatents=normalLatents,spikeW=spikeW,ARDW="grouped",
              ARDLatent="shared",imputation="Bayesian",iter.max=5000,iter.saved=100,iter.burnin=2500,
              init.tau=10^3,sampleZ=1,prior.alpha_0t=10,prior.beta_0t=10,prior.alpha_0=10,prior.beta_0=1,
              prior.alpha_0X=10,prior.beta_0X=1,prior.beta=c(1,1),prior.betaX=c(1,1),
              verbose=1,convergenceCheck=FALSE,save.posterior=list(X=TRUE,W=TRUE,tau=TRUE)))
}
