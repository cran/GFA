#' Visualize GFA components
#'
#' \code{visualizeComponents} illustrates the factorization inferred by GFA,
#' averaging over the posteriors of the parameters, if they have been stored.
#'
#' @param model The learned GFA model.
#' @param Y The used input data to be plotted, if supplied. Default NULL.
#' @param norm The normalization acquired from \code{\link{normalizeData}}, if
#'        applied. If provided, the reconstruction is shown in the original data
#'        space. Default NULL.
#' @param mode Determines which mode to visualize in case of pairing in two
#'        modes (default: 1).
#' @param showAll Show the full predictions and factorizations? May be
#'        cumbersome for large data. Default TRUE.
#' @param hclust Order features and samples based on hierarchical clustering?
#'        Default FALSE.
#' @param topK Number of strongest components visualized in the data space.
#'        Default 3.
#' @param topFeatures How many most relevant features to show for the data space
#'        visualizations? Default NA, showing all the features.
#' @param topSamples How many most relevant samples to show for the data space
#'        visualizations? Default NA, showing all the samples.
#' @return A list containing the matrices that have been visualized.
visualizeComponents <- function(model,Y=NULL,norm=NULL,mode=1,showAll=TRUE,hclust=FALSE,
                                topK=3,topFeatures=NA,topSamples=NA) {
  if(showAll) {
    rec <- reconstruction(model, average=TRUE)
    if(!is.null(norm)) rec <- undoNormalizeData(rec, norm)
  }
  
  V <- 1 #Number of modes the data is paired in
  if(typeof("pred")=="list") {
    V <- 2
    pred <- pred[[mode]]
    model$groups <- model$groups[[mode]]
    model$W <- model$W[[mode]]
    model$X <- model$X[[mode]]
    model$K <- model$K[mode]
    if(!is.null(Y)) Y <- Y[[mode]]
  }
  K <- model$K
  gr <- model$groups; M <- length(gr)
  if(is.null(names(gr))) names(gr) <- paste("Source",1:M)
  if(!is.null(Y)) Y <- do.call(cbind,Y)
  
  Yk <- vector("list",length=K) #Component effects sizes
  
  if(!is.null(model$posterior$X) & !is.null(model$posterior$W)) { #Posterior available
    if(V==2) {
      model$posterior$X <- model$posterior$X[[mode]]
      model$posterior$W <- model$posterior$W[[mode]]
    }
    N <- nrow(model$X); D <- nrow(model$W)
    ps <- nrow(model$posterior$X)
    X <- array(0,dim=c(N,K),dimnames=list(rownames(model$X),c()))
    W <- array(0,dim=c(D,K),dimnames=list(rownames(model$W),c()))
    for(p in 1:ps) {
      X <- X + model$posterior$X[p,,]/ps
      W <- W + model$posterior$W[p,,]/ps
    }
    
  } else {
    ps <- NA
    X <- model$X; W <- model$W
  }
  kEff <- colMeans(abs(X))*colMeans(abs(W))
  if(topK==0) kShow <- c()
  else kShow <- order(kEff,decreasing=T)[1:min(topK,K)]

  if(is.na(ps)) {
    for(k in kShow)
      Yk[[k]] <- outer(X[,k],W[,k])
    
  } else {
    for(k in kShow) {
      Yk[[k]] <- matrix(0,N,D)
      for(p in 1:ps)
        Yk[[k]] <- Yk[[k]] + outer(model$posterior$X[p,,k], model$posterior$W[p,,k])/ps
    }
  }
  for(k in kShow) {
    rownames(Yk[[k]]) <- rownames(X)
    colnames(Yk[[k]]) <- rownames(W)
  }

  
  vis <- function(x,title,shown,gr=NULL,hclust=c(FALSE,FALSE)) {
    N <- nrow(x); D <- ncol(x)
    
    if(hclust[1]) {
      if(is.null(rownames(x))) rownames(x) <- paste0("n",1:N)
      x <- x[hclust(dist(x))$order, ]
    }
    if(hclust[2]) {
      if(is.null(gr)) stop("Clustering of features TRUE, but feature groups not provided.")
      id <- c(); coln <- c()
      for(m in 1:(length(gr)-1)) {
        dd <- (gr[m]+1):gr[m+1]
        id <- c(id, hclust(dist(t(x[,dd])))$order+gr[m])
        coln <- c(coln, paste0(names(gr)[m]," - ","d",1:length(dd)))
      }
      if(is.null(colnames(x))) colnames(x) <- coln
      x <- x[,id]
    }
    
    mar <- c(4,4,4,6)
    for(i in 1:2) if(!is.null(dimnames(x)[[c(2,1)[i]]])) mar[i] <- mar[i]+3
    #if(!is.null(dimnames(x)[[1]])) mar[4] <- mar[4]+2
    par(mar=mar)
    cols <- colorRampPalette(c("orange","red","white","blue","cyan"))(101)
    if(any(is.na(x))) cols <- colorRampPalette(c("orange","red","#DDDDDD","blue","cyan"))(101)
    M <- max(abs(x),na.rm=T)
    breaks <- seq(-M,M,length=102)
    s <- names(shown)
    shown[[length(shown)+1]] <- x
    names(shown) <- c(s,title[1])
    
    image(1:D,1:N,t(x[N:1,,drop=F]),col=cols,breaks=breaks,axes=F,main=title[1],
          xlab="",ylab="")
    title(xlab=title[3],line=mar[1]-1)
    title(ylab=title[2],line=mar[2]-1)
    box()
    for(i in 1:2) { #Data modes
      at <- 1:c(N,D)[i]
      if(is.null(dimnames(x)[[i]])) {
        par(las=1)
        while(length(at)>25) {
          if(length(at)<70) axis(i,at,rep(NA,length(at)),tck=-0.01)
          at <- at[seq(10,length(at),by=10)]
        }
        if(i==1) axis(c(2,1)[i],rev(at),at)
        else axis(c(2,1)[i],at,at)
      } else { #Named axes
        par(las=2)
        if(i==1) axis(c(2,1)[i],rev(at),dimnames(x)[[i]], cex.axis=N^(-1/5))
        else axis(c(2,1)[i],at,dimnames(x)[[i]],cex.axis=D^(-1/5))
      }
    }
    #Grouping
    par(xpd=T)
    if(!is.null(gr)) {
      mu <- gr[-1]/2+gr[-length(gr)]/2
      for(i in 1:length(mu)) {
        if(i!=length(mu)) lines(rep(gr[i+1]+1/2,2), c(.5, N*1.03+.5), lwd=2)
        text(mu[i],N*1.03+.5,names(gr)[i])
      }
    }
    #Colorbar
    n <- length(cols)
    cba <- D + 1/2 + D/60; cbw <- D/40
    for(i in 1:n){
      polygon(c(0,cbw,cbw,0)+cba, c(0,0,N/n,N/n)+N*(i-1)/n+1/2,
              col=cols[i], border=NA)
    }
    #Colorbar: axis
    lines(rep(cba+cbw,2),c(0,N)+1/2)
    m <- 10^floor(log10(M)); m <- floor(M/m)*m
    for(l in c(-m,0,m)) {
      ly <- N*(l/M/2+.5)+1/2
      lines(cba+cbw-c(cbw,-cbw)/5, rep(ly,2))
      text(cba+cbw*2.5,ly,l)
    }
    par(xpd=F)
    return(shown)
  }
  
  # Plot view-component activities
  viewComp <- matrix(NA,M,K); rownames(viewComp) <- names(gr)
  for(m in 1:length(gr))
    for(k in 1:K)
      viewComp[m,k] <- mean(abs(W[gr[[m]],k]))
  
  shown <- list()
  shown <- vis(viewComp, c("Component activities","","Components"), shown)
  
  #Plot the full data, and the full factorization
  gr1 <- c(0,cumsum(sapply(gr,length))); names(gr1) <- c(names(gr),"NA")
  if(showAll) {
    #par(mfrow=c(1,1))
    if(!is.null(Y)) vis(Y, c("Observed data","Samples","Features"), shown, gr1)
    shown <- vis(rec, c("Data reconstruction","Samples","Features"), shown, gr1)
    #par(mfrow=c(1,2))
    shown <- vis(X, c("Latent variable X","Samples","Components"), shown,
                 hclust=c(hclust,F))
    shown <- vis(t(W), c("Projection matrix W^T","Components","Features"),
                 shown, gr1, hclust=c(F,hclust))
  }
  
  #Highlight the components explaning most variance - possibly for subsets of
  #samples and features
  for(k in kShow) {
    grk <- gr1
    #par(mfrow=c(1,1))
    lab1 <- "Samples"
    lab2 <- "Features"
    if(!is.na(topSamples)) {
      if(is.null(rownames(Yk[[k]])))
        rownames(Yk[[k]]) <- paste0("n",1:nrow(Yk[[k]]))
      keep <- order(rowSums(abs(Yk[[k]])),decreasing=T)
      Yk[[k]] <- Yk[[k]][keep[1:min(length(keep),topSamples)],]
      lab1 <- paste0(lab1," (top",topSamples,")")
    }
    if(!is.na(topFeatures)) {
      coln <- c(); grTop <- list()
      for(m in 1:length(gr)) {
        coln <- c(coln, paste0(names(gr)[m],", d",1:length(gr[[m]])))
        keep <- gr[[m]][order(colSums(abs(Yk[[k]][,gr[[m]]])),decreasing=T)]
        grTop[[m]] <- keep[1:min(length(keep),topFeatures)]
      }
      if(is.null(colnames(Yk[[k]]))) colnames(Yk[[k]]) <- coln
      grk <- c(0,cumsum(sapply(grTop,length))); names(grk) <- names(gr1)
      Yk[[k]] <- Yk[[k]][,unlist(grTop)]
      lab2 <- paste0(lab2," (top",topFeatures,")")
      
    }
    shown <- vis(Yk[[k]], c(paste("Component",k,"effect"),lab1,lab2),
                 shown, grk, hclust=c(hclust,hclust))
  }
  return(shown)
}


