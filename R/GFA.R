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
    ptm <- proc.time()[1] #timing
    
    ## Initialization
    init <- initializeParameters(Y=Y,opts=opts,K=K,projection=projection)
    for(param in init$unlist) { #Unlist the initialized model parameters (to be sampled)
      eval(parse(text=paste0(param," <- init$",param)))
      eval(parse(text=paste0("init$",param," <- NULL")))
    }
    gc()
    V <- length(Y) #Number of data modes
    ov <- V:1 #Indicator for the other data mode
    gr <- init$gr #Indicators pointing which features belong to which data views
    M <- sapply(gr, length) #Number of data views/sources in each mode
    N <- sapply(Y, nrow) #Number of samples per mode
    D <- sapply(Y, ncol) #Number of features per mode
    opts$K <- opts$initK <- init$K
    opts$missingValues <- init$missingValues
    opts$projection.fixed <- init$projection.fixed
    
    posterior <- NULL
    iter <- 1
  }
  
  ##
  ## The main loop - Gibbs updates
  ##
  for(iter in iter:opts$iter.max) {
    if (!opts$projection.fixed) {
      ##
      ## Sample the projections W
      ##
      for(v in 1:V) {
        tmpY <- Y[[v]]
        if(V==2) { #Substract the components of the other mode from the shared data view
          tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
            t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        groups <- if(opts$spikeW=="group") gr[[v]] else NULL #Group sparsity or not
        up <- updateSpikeAndSlab(Y=tmpY, A=W[[v]], Z=Z[[v]], mode=2, B=X[[v]], BB=XX[[v]], alpha=alpha[[v]],
                                 r=r[[v]], tau=tau[[v]], gr=groups, na.inds=list(V=init$na[[v]],Y=init$na$Y[[v]]),
                                 opts=opts, updateSpike=iter>opts$sampleZ)
        W[[v]] <- up$A
        Z[[v]] <- up$Z
      }
    }
    
    ## Check if components can be removed.
    # If Gibbs samples are going to be saved, components are dropped only during the burn-in.
    # If the projection has been given as an argument (e.g. for prediction), components will not be dropped.
    if ((iter<=opts$iter.burnin | opts$iter.saved==0) & !opts$projection.fixed) {
      for(v in 1:V) {
        keep <- which(colSums(Z[[v]])!=0 & colSums(abs(X[[v]]))>0)
        if(length(keep)==0) {
          print(paste0("All components shut down (mode ",v,"), returning a NULL model."))
          return(list())
        }
        if(length(keep)!=opts$K[v]) { #Prune out the shut down components
          opts$K[v] <- length(keep)
          for(param in c("X","W","Z","alpha","beta","r","rz")) {
            eval(parse(text=paste0(param,"[[v]] <- ",param,"[[v]][,keep,drop=FALSE]")))
          }
        }
        mm <- ""; if(V==2) mm <- paste0(" (mode ",v,")")
        if(iter==opts$iter.burnin & opts$K[v] <= opts$initK[v]*.4 & opts$initK[v]>10)
          warning(paste0("Over 60% of the initial components shut down",mm,". Excessive amount of
  initial components may result in a bad initialization for the model parameters.
  Rerun the model with K=",length(keep)+5,"."))
      }
    }
    

    for(v in 1:V) {
      if(iter>opts$sampleZ & !opts$projection.fixed){
        ## sample r (the probability of slab)
        r[[v]] <- updateBernoulli(r[[v]], Z[[v]], gr[[v]], opts) 
      }
      
      ## 
      ## Sample the latent variables X
      ##
      if(iter>opts$sampleZ & !opts$normalLatents) {
        # Spike-and-slab prior for X (biclustering)
        tmpY <- Y[[v]]
        if(V==2) { #Substract the components of the other mode from the shared data view
         tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
           t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        up <- updateSpikeAndSlab(Y=tmpY, A=X[[v]], Z=NULL, mode=1, B=W[[v]], BB=crossprod(W[[v]]*sqrt(tau[[v]])),
                                alpha=beta[[v]], r=rz[[v]], tau=tau[[v]], gr=NULL, na.inds=list(V=init$na[[v]], Y=init$na$Y[[v]]),
                                opts=opts, updateSpike=iter>opts$sampleZ)
        X[[v]] <- up$A
        if (!opts$projection.fixed) {
          ## sample rz (the probability of slab)
          for(k in 1:opts$K[v])
            rz[[v]][,k] <- (sum(up$Z[,k])+opts$prior.betaX[1])/(N[v]+sum(opts$prior.betaX))
        }
        
      } else {
        #Normal prior for X
        tmpY <- Y[[v]]
        if(V==2) {
          tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] - t(tcrossprod(X[[ov[v]]],W[[ov[v]]][gr[[ov[v]]][[1]],]))
        }
        X[[v]] <- updateNormal(Y=tmpY, W=W[[v]], tau=tau[[v]],
                               na.inds=list(V=init$na[[v]], Y=init$na$Y[[v]], Nlist=init$na$Nlist[[v]]), opts=opts)
          
      }
      XX[[v]] <- crossprod(X[[v]])
      
      ##
      ## Update alpha and beta, the ARD parameters
      ##
      if(!opts$normalLatents) {
        # Spike-and-slab prior for X, including normal precision beta (ARD)
        groups <- if(opts$ARDLatent=="shared") list(1:N[v]) else as.list(1:N[v])
        beta[[v]] <- updateGamma(alpha=beta[[v]], A=X[[v]],Z=(X[[v]]!=0)*1,alpha_0=init$prior.alpha_0X[v],
                                 beta_0=init$prior.beta_0X[[v]],gr=groups)

      }
      if (!opts$projection.fixed) {
        groups <- if(opts$ARDW=="shared") list(1:D[v]) else if(opts$ARDW=="grouped") gr[[v]] else as.list(1:D[v])
        alpha[[v]] <- updateGamma(alpha=alpha[[v]], A=W[[v]],Z=Z[[v]],alpha_0=init$alpha_0[[v]],beta_0=init$beta_0[[v]],gr=groups)

        
        ##
        ## Update tau, the noise precisions
        ##
        tmp <- updateTau(tau=tau[[v]], Y=Y[[v]], X=X[[v]], W=W[[v]], a_tau=init$a_tau, b_tau=b_tau[[v]],
                              alpha_0t=init$alpha_0t[[v]], beta_0t=init$beta_0t[[v]], gr=gr[[v]],
                              na.inds=init$na[[v]], opts, V=V, v=v, iter=iter)
        tau[[v]] <- tmp$tau
        b_tau[[v]] <- tmp$b_tau
        
        # Calculate the likelihood.
        cost[iter] <- cost[iter] + N[v]*sum(log(tau[[v]][tmp$id]))/2 -
          crossprod(b_tau[[v]][tmp$id],tau[[v]][tmp$id])
      }
      
    } #Loop over v (data collections) ends
    
    #Update tau for the paired view (pairing in two modes)
    if(!opts$projection.fixed & V==2) {
      XW <- tcrossprod(X[[1]],W[[1]][gr[[1]][[1]],]) + t(tcrossprod(X[[2]],W[[2]][gr[[2]][[1]],]))
      if(opts$missingValues && length(init$na[[3]])) {
        XW <- Y[[1]][,gr[[1]][[1]]] - XW
        XW[init$na[[3]]] <- 0
        b_tau[[1]][gr[[1]][[1]]] <- colSums(XW[,gr[[1]][[1]]]^2)/2
      } else {
        b_tau[[1]][gr[[v]][[1]]] <- colSums((Y[[1]][,gr[[1]][[1]]]-XW)^2)/2
      }
      
      if(opts$tauGrouped){
        tau[[1]][gr[[1]][[1]]] <- rgamma(1, shape=init$alpha_0t[[1]][1],
                                         rate=init$beta_0t[[1]][1] + sum(b_tau[[1]][gr[[1]][[1]]]))
        tau[[2]][gr[[2]][[1]]] <- tau[[1]][gr[[1]][[1]][1]]
      }else{
        stop("Tau needs to be grouped: otherwise the noise prior for the shared view is ill-defined")
      }
      cost[iter] <- cost[iter] + N[1]*sum(log(tau[[1]][gr[[1]][[1]]]))/2 -
        crossprod(b_tau[[1]][gr[[1]][[1]]],tau[[1]][gr[[1]][[1]]])
    }
    aic[iter] <- 2*cost[iter] #Akaike information criterion
    for(v in 1:V) aic[iter] <- aic[iter] - (D[v]*(opts$K[v]+1)-opts$K[v]*(opts$K[v]-1)/2)*2
    
    
    if(filename != "" & iter%%100==0) { #Every 100 iterations, save the sampling
      time <- proc.time()[1]-ptm
      save(list=ls(),file=filename)
      if(opts$verbose>0)
        print(paste0("Iter ", iter, ", saved chain to '",filename,"'"))
    }
    
    # Calculate likelihood of observed views
    if(opts$projection.fixed) {
      for(v in 1:V) {
        tautmp <- ( init$Yconst[[v]] + rowSums((W[[v]]%*%XX[[v]])*W[[v]]) -
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
        if (iter%%10==0 & opts$verbose>1)
          print(paste0("Predicting: ",iter,"/",opts$iter.max))
        if (iter>opts$iter.burnin & ((iter-opts$iter.burnin)%%init$mod.saved)==0) {
          init$i.pred <- init$i.pred+1/V
          for (m in which(opts$prediction[[v]])) # Go through the views that will be predicted.
            prediction[[v]][[m]] <- prediction[[v]][[m]] + tcrossprod(X[[v]], W[[v]][gr[[v]][[m]],])
        }
        
      } else if(!any(unlist(opts$prediction))) {
        if (((iter%%10==0 & opts$verbose>1) | (iter%%100==0 & opts$verbose>0)) & v==1) {
          print(paste0("Learning: ",iter,"/",opts$iter.max," - K=",paste0(opts$K,collapse=",")," - ",Sys.time()))
        }
        ## Collection of Gibbs samples
        if (opts$iter.saved>0 && iter>opts$iter.burnin && ((iter-opts$iter.burnin)%%init$mod.saved)==0) {
          if (iter-opts$iter.burnin==init$mod.saved) { # Initialize the list containing saved Gibbs samples.
            if(v==1) posterior <- NULL
            posterior <- initializePosterior(posterior=posterior,Y=Y[[v]],V=V,v=v,N=N[v],D=D[v],M=M[v],S=init$S,opts=opts)
            if(v==1) {gr.start <- list()}
            gr.start[[v]] = vector(mode="integer", length=length(gr[[v]]))
            for (m in 1:length(gr[[v]])) { # Find the first index of each view.
              gr.start[[v]][m] = gr[[v]][[m]][1]
            }
          }
          
          s <- (iter-opts$iter.burnin)/init$mod.saved
          if (!opts$projection.fixed) {
            if (opts$save.posterior$W | (opts$convergenceCheck && s%in%c(init$psStart,init$psEnd)))
              posterior$W[[v]][s,,] <- W[[v]]
            posterior$tau[[v]][s,] <- tau[[v]]
            posterior$rz[[v]][s,] <- rz[[v]][1,] #'rz' is identical for all samples
            posterior$r[[v]][s,,] <- r[[v]][gr.start[[v]],]
            #'r' is identical for all variables within a view. Thus, save only these.
            posterior$beta[[v]][s,,] <- beta[[v]]
          }
          if (opts$save.posterior$X | (opts$convergenceCheck && s%in%c(init$psStart,init$psEnd)))
            posterior$X[[v]][s,,] <- X[[v]]
        }
      }
    }
  } ## Gibbs sampling done
  
  if(opts$convergenceCheck && opts$iter.saved>=8 && !opts$projection.fixed) {
    #Estimate the convergence of the data reconstruction, based on the Geweke diagnostic
    conv <- checkConvergence(posterior=posterior, V=V, N=N, D=D, start=init$psStart,
                             end=init$psEnd, opts=opts)
    if(!opts$save.posterior$X)
      posterior$X[[v]] <- NULL
    if(!opts$save.posterior$W)
      posterior$W[[v]] <- NULL
    gc()
    
  } else {
    conv <- NA
  }
  
  if(filename!="" & opts$iter.max>=10)
    file.remove(filename) #No need for temporary storage any more
  

  if(opts$projection.fixed) {
    for(v in 1:V) {
      if(any(opts$prediction[[v]])) {
        for (m in which(opts$prediction[[v]])) # Go through the views that will be predicted.
          prediction[[v]][[m]] <- prediction[[v]][[m]]/init$i.pred
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
      D <- D[1]; opts$K <- opts$K[1]
    }
    
    time <- proc.time()[1] - ptm
    
    return(list(W=W, X=X, Z=Z, r=r, rz=rz, tau=tau, alpha=alpha, beta=beta, groups=gr, D=D,
                K=opts$K, cost=cost, aic=aic, posterior=posterior, opts=opts, conv=conv, time=time))
  }
}


## An internal function for updating the bernoulli spike and slab variable
# A is the factorization parameter being updated, and B the one fixed
updateSpikeAndSlab <- function(Y,A,Z,mode,B,BB,alpha,r,tau,gr,na.inds,opts,V,updateSpike) {
  BBdiag <- diag(BB)
  diag(BB) <- 0
  if(opts$missingValues) {
    xx <- B^2
    if(mode==1) xx <- xx*tau
  }
  
  tt <- if(mode==1) 1 else if(mode==2) tau
  for(k in 1:ncol(A)) {
    if(opts$missingValues) {
      if(opts$imputation=="Bayesian") { #||x||^2, based on observations only
        BBdiag <- apply(na.inds$Y, mode, function(y) sum(xx[,k]*y))
        lambda <- tt*BBdiag + alpha[,k]
      } else if(opts$imputation=="conservative") {
        lambda <- tt*BBdiag[k] + alpha[,k]
      } else {
        stop("Ill-defined option 'imputation'")
      }
      if(mode==1) ss <- tcrossprod(A[,-k],B[,-k])
      if(mode==2) ss <- tcrossprod(B[,-k],A[,-k])
      tmp <- Y-ss
      tmp[na.inds$V] <- 0 #YX based on observations only
      if(mode==1) mu <- tmp%*%(B[,k]*tau)/lambda
      if(mode==2) mu <- tau/lambda*c(crossprod(tmp,B[,k]))
    } else {
      lambda <- tt*BBdiag[k] + alpha[,k]
      if(mode==1) mu <- (Y%*%(B[,k]*tau) - A%*%BB[k,])/lambda
      if(mode==2) mu <- tau/lambda*c(crossprod(Y,B[,k]) - A%*%BB[k,])
    }
    A[,k] <- mu + rnorm(nrow(A))*lambda^(-0.5)
    
    ##
    ## Updating the spike and slab prior, parameter Z (H in publications)
    ##
    if(updateSpike) {
      if(is.null(Z)) Z <- matrix(NA, nrow(A), ncol(A))
      if(is.null(gr)) {
        zone <- 0.5*( log(alpha[,k]) - log(lambda) + lambda*mu^2) +
          log(r[,k]) - log(1-r[,k])
        zone <- 1/(1+exp(-zone))
        zone <- as.double(runif(nrow(Z)) < zone)
        Z[,k] <- zone
      } else {
        zone <- 0.5*(log(alpha[,k])-log(lambda) + lambda*mu^2)
        for(m in 1:length(gr)) {
          logpr <- sum(zone[gr[[m]]]) + log(r[gr[[m]][1],k]) -
            log( 1-r[gr[[m]][1],k] )
          logpr <- 1/(1+exp(-logpr))
          Z[gr[[m]],k] <- as.double(runif(1)<logpr)
        }
      }
    }
    A[,k] <- A[,k]*Z[,k]
  }
  
  return(list(A=A,Z=Z))
}


## An internal function for updating normally distributed X
updateNormal <- function(Y, W, tau, na.inds, opts) {
  N <- nrow(Y)
  K <- ncol(W)
  id <- rep(1,K)
  
  covZ <- diag(1,K) + crossprod(W*sqrt(tau))
  eS <- eigen(covZ)
  tmp <- eS$vectors*outer(id,1/sqrt(eS$values))
  covZ <- tcrossprod(tmp)
  
  if(opts$missingValues) {
    Y[na.inds$V] <- 0
    X <- Y%*%(W*tau)
    if(opts$imputation=="Bayesian") {
      #Missing values cause each sample to have a different covariance matrix
      #NOTE/TODO: Sampling components independently would be faster, but converge slower
      nn <- na.inds$Nlist #Samples ordered to avoid duplicate covX computations
      compcovX <- TRUE
      for(i in 1:N) {
        if(all(na.inds$Y[nn[i],]==1)) { #All features observed for this sample
          X[nn[i],] <- X[nn[i],,drop=F]%*%covZ + matrix(rnorm(K),1,K)%*%tmp
        } else {
          if(i>1) { #If different features observed, compute a new inverse (covX)
            if(any(na.inds$Y[nn[i],]!=na.inds$Y[nn[i-1],])) compcovX <- TRUE
            else compcovX <- FALSE
          }
          if(compcovX) {
            covX <- diag(1,K) + crossprod(W*na.inds$Y[nn[i],]*sqrt(tau))
            eS <- eigen(covX) #Main source of computation!
            tmp2 <- eS$vectors*outer(id,1/sqrt(eS$values))
            covX <- tcrossprod(tmp2)
          }
          X[nn[i],] <- X[nn[i],,drop=F]%*%covX + matrix(rnorm(K),1,K)%*%tmp2
        }
      }
    } else if(opts$imputation=="conservative") {
      X <- X%*%covZ + matrix(rnorm(N*K),N,K)%*%tmp
    }
    
  } else {
    X <- Y%*%(W*tau)%*%covZ
    X <- X + matrix(rnorm(N*K),N,K)%*%tmp
  }
  return(X)
}


## An internal function for updating a Gamma distribution
updateGamma <- function(alpha, A, Z, alpha_0, beta_0, gr) {
  for(m in 1:length(gr)) {
    alpha[gr[[m]],] <- rep(rgamma(ncol(alpha), shape=colSums(Z[gr[[m]],,drop=F])/2+alpha_0[m],
                                  rate=colSums(A[gr[[m]],,drop=F]^2)/2+beta_0[m]), each=length(gr[[m]]))
  }
  return(alpha+1e-7) #1e-7 added to precision to avoid ridiculously small values (huge variance)
}


## An internal function for updating a Bernoulli distribution
updateBernoulli <- function(r, Z, gr, opts) {
  if(opts$spikeW=="shared"){
    zm <- colSums(Z)
    for(k in 1:ncol(r))
      r[,k] <- (zm[k]+opts$prior.beta[1])/(nrow(r)+sum(opts$prior.beta))
  }
  if(opts$spikeW=="element"){
    for(m in 1:length(gr)) {
      zm <- colSums(Z[gr[[m]],,drop=FALSE])
      for(k in 1:ncol(r))
        r[gr[[m]],k] <- (zm[k]+opts$prior.beta[1])/(length(gr[[m]])+sum(opts$prior.beta))
    }
  }
  if(opts$spikeW=="group"){
    zm <- rep(0,ncol(r))
    for(m in 1:length(gr)){
      zm <- zm + Z[gr[[m]][1],]
    }
    for(m in 1:length(gr)){
      for(k in 1:ncol(r)){
        r[gr[[m]],k] <- (zm[k]+opts$prior.beta[1])/(length(gr)+sum(opts$prior.beta))
      }
    }
  }
  return(r)
}


## Internal function for updating noise residual tau
updateTau <- function(tau, Y, X, W, a_tau, b_tau, alpha_0t, beta_0t, gr, na.inds, opts, V, v, iter) {
  tmpid <- 1:length(tau)
  if(V==2) tmpid <- tmpid[-gr[[1]]]
  
  XW <- tcrossprod(X,W)
  if(opts$missingValues) {
    if(iter==opts$iter.max && !opts$projection.fixed && opts$missingValues) {
      id <- na.inds[,1] + (na.inds[,2]-1)*nrow(XW)
      if(length(id)>0 && (max(abs(XW)) > 2*max(abs(XW[-id])))) {
        vv <- if(V==1) "" else paste0(" (mode ",v,")")
        warning(paste0("Reconstruction of missing values resulted in larger scale than observed data",vv,".
    This may be a sign of overfitting to sparsely observed samples/columns. Consider
    using informative Gamma-priors for the parameters X and W to avoid this."))
      }
    }
    XW <- Y - XW
    XW[na.inds] <- 0
    b_tau[tmpid] <- colSums(XW[,tmpid]^2)/2
  } else{
    b_tau[tmpid] <- colSums((Y[,tmpid]-XW[,tmpid])^2)/2
  }
  if(opts$tauGrouped) {
    for(m in V:length(gr)) { #For pairing in two modes (V==2), update m=1 later
      tau[gr[[m]]] <- rgamma(1, shape=alpha_0t[m], rate=beta_0t[gr[[m]][1]]+sum(b_tau[gr[[m]]]))
    }
  } else {
    if(V==2)
      stop("Tau needs to be grouped: otherwise the noise prior for the shared view is ill-defined")
    for(d in 1:length(tau))
      tau[d] <- rgamma(1,shape=a_tau[d], rate= beta_0t[[1]][d]+b_tau[[1]][d])
  }
  return(list(tau=tau, b_tau=b_tau, id=tmpid))
}


## An internal function for intialization a list for storage of posterior samples
initializePosterior <- function(posterior,Y,V,v,N,D,M,S,opts) {
  add <- ""; if(V==2) add <- paste0(" (mode ",v,")")
  if(opts$K[v]==opts$initK[v]) {
    if(opts$initK[v]<min(N,D)) {
      warning(paste0("Burn-in period finished, but no components have been shut down",add,",
  preventing automatic model complexity selection. Consider higher inital K, if
  computational resources permit (up to K=",min(N,D),", default K=",
                     round(min(N,D)/2),").
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
    posterior$W[[v]] <- array(dim=c(S, D, opts$K[v]))
    colnames(posterior$W[[v]]) <- colnames(Y)
  }
  posterior$rz[[v]] <- array(dim=c(S, opts$K[v])) #Same 'rz' for all samples. Thus, save only a vector per iteration.
  posterior$r[[v]] <- array(dim=c(S, M, opts$K[v])) #Same 'r' for all variables of a view.
  posterior$tau[[v]] <- array(dim=c(S, D)) # SxD
  colnames(posterior$tau[[v]]) = colnames(Y)
  if ((opts$save.posterior$X | opts$convergenceCheck) & !any(opts$prediction[[v]])) {
    if(v==1) {posterior$X <- list()}
    posterior$X[[v]] <- array(dim=c(S, N, opts$K[v]))
    colnames(posterior$X[[v]]) <- rownames(Y)
  }
  posterior$beta[[v]] = array(dim=c(S, N, opts$K[v]))
  return(posterior)
}



## An internal function for checking the convergence of the sampling
checkConvergence <- function(posterior,V,N,D,start,end,opts) {
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
      y <- matrix(NA,0,D[v])
      for(ps in c(start,end)) {
        y <- rbind(y, tcrossprod(posterior$X[[v]][ps,i,],posterior$W[[v]][ps,,]))
      }
      foo <- rep(NA,D[v])
      for(j in 1:D[v]) {
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
  }
  if(opts$verbose>0) {
    print(paste("Convergence diagnostic:",round(conv,4)))
    print("Values significantly greater than 0.05 imply a non-converged model.")
  }
  return(conv)
}


## An internal function for initializing the GFA parameters
## note: no roxygen-comments, so that this function is undocumented in the manuals
initializeParameters <- function(Y,opts,K,projection) {
  init <- list()
  
  if(typeof(Y[[1]])=="list") { #Pairing in two modes
    V <- length(Y)
    if(V==2)
      if(!identical(Y[[1]][[1]], t(Y[[2]][[1]])))
        stop("The two sets of views must have the first view shared (transposed).")
    if(opts$verbose>0 & is.null(projection))
      print("Running GFA with data sources paired in two modes.")
  } else {
    Y <- list(Y); V <- 1
    if(opts$verbose>0) print(paste("Running GFA with",length(Y[[1]]),"data sources paired in one mode."))
  }
  
  # Store dimensionalities of the data sets
  M <- D <- N <- rep(NA,V)
  init$ov <- V:1 #The other view (in case of pairing in two modes)
  gr <- list()
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
  
  if(any(K>=200) & opts$verbose>0) {
    print("Using K>=200; computational cost may be high. This can be redeemed by 
             setting lower initial K, or if that does not allow automated model 
             complexity selection (does not shut down any components), with the help 
             of function informativeNoisePrior().")
  }
  
  for(p in c("Y","gr","K")) #Store the often used data dimensions in the init-list too
    eval(parse(text=paste0("init$",p," <- ",p)))
  
  init$alpha_0 <- opts$prior.alpha_0
  init$beta_0 <- opts$prior.beta_0
  init$alpha_0t <- opts$prior.alpha_0t
  init$beta_0t <- opts$prior.beta_0t
  #Make the alpha prior vector-valued, if it is not already
  if(typeof(init$alpha_0)!="list" && length(init$alpha_0)==1) {
    init$alpha_0 <- init$beta_0 <- list()
    for(v in 1:V) {
      if(opts$ARDW=="shared") Nrep <- 1
      if(opts$ARDW=="grouped") Nrep <- M[v]
      if(opts$ARDW=="element") Nrep <- D[v]
      init$alpha_0[[v]] <- rep(opts$prior.alpha_0, Nrep)
      init$beta_0[[v]] <- rep(opts$prior.beta_0, Nrep)
    }
  }
  init$prior.alpha_0X <- rep(opts$prior.alpha_0X, V)
  init$prior.beta_0X <- rep(opts$prior.beta_0X, V)
  #Make the tau prior vector-valued, if it is not already
  if(length(init$alpha_0t)==1 & typeof(init$alpha_0t)=="double") {
    init$alpha_0t <- list()
    for(v in 1:V) init$alpha_0t[[v]] <- rep(opts$prior.alpha_0t, D[v])
  }
  if(length(init$beta_0t)==1 & typeof(init$beta_0t)=="double") {
    init$beta_0t <- list()
    for(v in 1:V) init$beta_0t[[v]] <- rep(opts$prior.beta_0t, D[v])
  }
  
  const <- rep(NA,V)
  for(p in c("Yconst", "b_tau", "Z", "alpha", "r", "na",
             "W","beta","rz","tau","prediction","X","XX"))
    init[[p]] <- list()
  
  if (!is.null(projection)) { #Projection given as an argument.
    for(p in c("W","rz","tau","beta")) {
      init[[p]] <- projection[[p]]
    }
    init$projection.fixed <- TRUE
    for(v in 1:V) {
      K[v] <- init$K[v] <- ncol(init$W[[v]])
    }
    
  } else {
    init$projection.fixed <- FALSE
  }
  
  for(v in 1:V) {
    # Some constants for speeding up the computation
    init$Yconst[[v]] <- colSums(Y[[v]]^2) #unlist(lapply(Y[[v]],function(x){sum(x^2)}))
    
    ##
    ## Initialize the model randomly
    ##
    # Other initializations could be done, but overdispersed random initialization is quite good.
    
    # Latent variables
    init$X[[v]] <- matrix(rnorm(N[v]*K[v],0,1),N[v],K[v])
    init$X[[v]] <- scale(init$X[[v]])
    init$XX[[v]] <- crossprod(init$X[[v]]) # Second moments
    
    if(!opts$normalLatents & init$projection.fixed) {
      init$beta[[v]] <- matrix(colMeans(init$beta[[v]]),N,K,byrow=T) #Initialize as the mean over samples
    }else{
      init$beta[[v]] <- matrix(1,N[v],K[v])
    }
  }
  
  init$i.pred <- 0
  if (opts$iter.saved>0 & opts$iter.burnin<opts$iter.max) {
    init$mod.saved <- ceiling((opts$iter.max-opts$iter.burnin)/opts$iter.saved)
    #Save every 'mod.saved'th Gibbs sample.
    init$S <- floor((opts$iter.max-opts$iter.burnin)/init$mod.saved) # Number of saved Gibbs samples.
  }
  
  if(init$projection.fixed) { # Projection given as an argument to the sampler
    for(v in 1:V) {
      obs = NULL # variable indices of the source views
      for (mi in 1:length(opts$prediction[[v]])) { # Go through all views.
        if (!opts$prediction[[v]][mi]) {
          obs = c(obs, gr[[v]][[mi]])
        }
      }
      if(!opts$normalLatents && is.null(dim(init$rz[[v]]))) {
        init$rz[[v]] = matrix(data=init$rz[[v]], nrow=N[v], ncol=K[v], byrow=TRUE)
      }
      if (any(opts$prediction[[v]])) {
        init$prediction[[v]] <- vector(mode="list", length=length(gr[[v]]))
        for (m in which(opts$prediction[[v]])) {
          init$prediction[[v]][[m]] = matrix(0,N[v],length(gr[[v]][[m]]))
          init$Y[[v]][,gr[[v]][[m]]] <- 0
        }
      }
      #Added things
      const[v] <- 0
      for(m in which(opts$prediction[[v]]==F))
        const[v] <- const[v] - N[v]*length(gr[[v]][[m]])*log(2*pi)/2
      init$cost <- init$aic <- rep(sum(const),opts$iter.max)
      init$b_tau[[v]] <- init$beta_0t[[v]]
      #
    }
    
  } else { #Non-fixed projections
    for(v in 1:V) {
      const[v] <- - N[v]*D[v]*log(2*pi)/2
      init$tau[[v]] <- rep(opts$init.tau,D[v]) # The mean noise precisions
      init$b_tau[[v]] <- init$beta_0t[[v]]
      
      init$W[[v]] <- matrix(0,D[v],K[v])
      init$Z[[v]] <- matrix(1,D[v],K[v])
      
      init$alpha[[v]] <- matrix(1,D[v],K[v])
      
      init$cost <- init$aic <- rep(sum(const), opts$iter.max)
      if(v==2) init$cost <- init$cost + N[v]*length(gr[[1]][[1]])*log(2*pi)/2
      init$r[[v]] <- matrix(0.5,D[v],K[v])
      init$rz[[v]] <- matrix(0.5,N[v],K[v])
    }
  }
  
  ## Missing values (NAs) in the data
  init$missingValues <- FALSE
  for(v in 1:V) {
    init$na[[v]] <- which(is.na(init$Y[[v]]),arr.ind=T)
    if(V==2) init$na[[2+v]] <- which(is.na(Y[[v]][,gr[[v]][[1]]]),arr.ind=T) #Shared data view
    if(length(init$na[[v]])>0 & !init$projection.fixed) {
      init$missingValues <- TRUE
      if(opts$verbose>0)
        print("Missing Values detected, prediction using EM type approximation")
      #Update alpha_0t to take into account the observed data size
      if(V==1) init$a_tau <- init$alpha_0t[[1]] + colSums(!is.na(Y[[1]]))/2
      for(m in 1:M[v])
        init$alpha_0t[[v]][gr[[v]][[m]]] <- init$alpha_0t[[v]][gr[[v]][[m]]] + sum(!is.na(Y[[v]][,gr[[v]][[m]]]))/2
      
    } else {
      if(V==1) init$a_tau <- init$alpha_0t[[1]] + N/2 #Elementwise noise
      for(m in 1:M[v]) #View-wise noise
        init$alpha_0t[[v]][gr[[v]][[m]]] <- init$alpha_0t[[v]][gr[[v]][[m]]] + N[v]*length(gr[[v]][[m]])/2
    }
    for(m in 1:M[v])
      init$alpha_0t[[v]][m] <- init$alpha_0t[[v]][[gr[[v]][[m]][1]]]
    init$alpha_0t[[v]] <- init$alpha_0t[[v]][1:M[v]]
  }
  if(init$missingValues) {
    init$na$Y <- list()
    for(v in 1:V)
      init$na$Y[[v]] <- (!is.na(Y[[v]]))*1
    #Many matrix inversions are done - avoid duplicates by ordering
    init$na$Nlist <- list()
    for(v in 1:V)
      init$na$Nlist[[v]] <- hclust(dist(init$na$Y[[v]]))$order
  }
  
  #Which X and W posterior samples to save in case of convergenceCheck?
  if(opts$convergenceCheck & opts$iter.saved>=8 & !init$projection.fixed) {
    init$ps <- floor(opts$iter.saved/4)
    init$psStart <- 1:init$ps
    init$psEnd <- (-init$ps+1):0+opts$iter.saved
  } else {
    init$psStart <- init$psEnd <- c()
  }
  
  ## W initialization
  if(!init$projection.fixed) {
    for(v in 1:V) {
      covW <- diag(1,init$K[v]) + opts$init.tau*init$XX[[v]]
      eS <- eigen(covW)
      tmp <- eS$vectors*outer(rep(1,init$K[v]),1/sqrt(eS$values))
      covW <- tcrossprod(tmp)
      
      estimW = matrix(0,D[v],K[v])
      tmpY <- init$Y[[v]]
      if(v==2) { #Some variance explained by the first mode already
        tmpY[,gr[[v]][[1]]] <- tmpY[,gr[[v]][[1]]] -
          t(tcrossprod(init$X[[init$ov[v]]], init$W[[init$ov[v]]][gr[[init$ov[v]]][[1]],]))
      }
      for(k in 1:init$K[v]){
        tmpY[init$na[[v]]] <- 0
        estimW[,k] <- crossprod(tmpY,init$X[[v]][,k])
      }
      init$W[[v]] <- estimW%*%covW*opts$init.tau + matrix(rnorm(D[v]*K[v]),D[v],K[v])%*%tmp
    }
  }
  
  #Which parameters to unlist from init?
  init$unlist <- c("Y","X","XX","W","Z","alpha","beta","r","rz","tau","b_tau","prediction","aic","cost")
  
  return(init)
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
