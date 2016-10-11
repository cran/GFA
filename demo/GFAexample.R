
require("GFA")

#Function for generating toy data (input: pairing in V modes)
genData <- function(V=1) {
  if(V==1) print("Generating three data views paired in one mode with 6 group sparse components.")
  if(V==2) {
    print(paste("Generating data with pairing in two modes, with a total of 4",
                "data views. 4 bicluster components."))
  }
  N <- c(200,100) #Number of samples
  D <- list(c(100,50,60),c(200,70)) #Dimensions
  K <- c(2,2); if(V==1) K <- 6 #Component numbers
  Nzero <- list(list(1:100,101:200),list(1:50,51:100))
  if(V==1) Dzero <- list(list(c(),1:100,101:150,151:210,1:150,101:210))
  else Dzero <- list(list(1:100,c(1:30,151:210)),list(1:200,c(1:140,250:270)))
  
  Y <- y <- X <- W <- list()
  for(v in 1:V) {
    X[[v]] <- matrix(rnorm(N[v]*K[v]),N[v],K[v]) #Generate the latents
    W[[v]] <- matrix(rnorm(sum(D[[v]]*K[v])),sum(D[[v]]),K[v]) #Generate the projections
    if(V==1) { #Varying scale to components
      W[[v]] <- scale(W[[v]]); W[[v]] <- sweep(W[[v]],MARGIN=2,c(.8, 2.2, 1.4, 0.7, 1.7, 0.5),"*") 
    }

    for(k in 1:K[v]) { #Sparsity
      if(V==2) X[[v]][Nzero[[v]][[k]],k] <- 0 #Sample sparsity only with biclustering
      if(length(Dzero[[v]][[k]])>0) W[[v]][Dzero[[v]][[k]],k] <- 0
    }
    
    tau <- 1 #tau <- 0.2; if(V==1) tau <- 2
    Y[[v]] <- tcrossprod(X[[v]],W[[v]]) + matrix(rnorm(N[v]*sum(D[[v]]),0,tau),N[v],sum(D[[v]]))
    y[[v]] <- list()
    d <- 0
    for(m in 1:length(D[[v]])) {
      y[[v]][[m]] <- Y[[v]][,1:D[[v]][m]+d]
      if(v==2) { #Effect of both component sets to the shared view
        kk <- 1:K[v]
        y[[2]][[1]] <- t(y[[1]][[1]]) + tcrossprod(X[[v]][,kk],W[[v]][1:D[[2]][1],kk])
        y[[1]][[1]] <- t(y[[2]][[1]])
      }
      d <- d + D[[v]][m]
    }
  }
  groups <- list()
  for(v in 1:V) {
    groups[[v]] <- list()
    d <- c(0, cumsum(D[[v]]))
    for(m in 1:(length(d)-1)) groups[[v]][[m]] <- (d[m]+1):d[m+1]
  }

  yNew <- YNew <- NA
  ntest <- 1:20 # First 20 samples for test use
  if(V==1) {
    YNew <- Y[[1]][ntest,]; Y[[1]] <- Y[[1]][-ntest,]
    yNew <- list(list())
    for(m in 1:length(y[[1]])) {
      yNew[[1]][[m]] <- y[[1]][[m]][ntest,]; y[[1]][[m]] <- y[[1]][[m]][-ntest,]
    }
  }
  
  return(list(y=y,Y=Y,X=X,W=W,N=N,D=D,K=K,groups=groups,yNew=yNew,YNew=YNew))
}

#Function for plotting a data collection paired in two modes
showComp <- function(toy,res=NULL) { #Y,main,D1,D2) {
  D1 <- cumsum(toy$D[[1]])
  D2 <- cumsum(toy$D[[2]])
  Y <- list()
  
  if(is.null(res)) {
    for(v in 1:2) {
      for(k in 1:toy$K[v]) {
        Y[[length(Y)+1]] <- tcrossprod(toy$X[[v]][,k],toy$W[[v]][,k])
      }
    }
    main <- "True bicluster components"
    cat("True bicluster components:
  similar to a toy data experiment in Bunte et al., gray rectangle denote data
  sources. Various colors denote different components, each active in a subset
  of the samples and the features (bicluster).")
    
  } else {
    for(v in 1:2) {
      for(k in 1:res$K[v]) {
        X <- matrix(0,nrow(res$X[[v]]),nrow(res$W[[v]]))
        for(i in 1:opts$iter.saved)
          X <- X + tcrossprod(res$posterior$X[[v]][i,,k], res$posterior$W[[v]][i,,k])
        Y[[length(Y)+1]] <- X/opts$iter.saved
      }
    }
    main <- "Inferred components"
    cat("Inferred components:
  the decomposition GFA has inferred, matching closely to the true structure.
  The components have not been ordered, so the colors may vary.")
  }
  
  K <- length(Y)
  X <- matrix(0,max(D2),max(D1))
  for(k in 1:K) {
    N <- nrow(Y[[k]])
    D <- ncol(Y[[k]])
    x <- abs(Y[[k]])+k*10
    if(N==D2[1]) {
      id <- which(x%%10 >= X[1:N,1:D]%%10,arr.ind=T)
      X[id] <- x[id] #Showing only the strongest component
    } else {
      id <- which(t(x)%%10 >= X[1:D,1:N]%%10)
      X[id] <- t(x)[id]
    }
  }
  X[is.na(X)] <- 10 #Base level (gray)
  X[-1:-(D2[1]),-1:-(D1[1])] <- NA #Top right white (no data)
  
  #Separate the data sets
  s <- 20 #Whitespace
  Y <- X[1:(D2[1]*1),1:(D1[1]*1)]
  for(i in 2:length(D1))
    Y <- cbind(Y,Y <- matrix(NA,D2[1]*1,s),X[1:(D2[1]*1),(D1[i-1]*1+1):(D1[i]*1)])
  Y <- rbind(Y,matrix(NA,(D2[2]-D2[1])*1+s,ncol(Y)))
  Y[(D2[1]*1+1):(D2[2]*1)+s,1:(D1[1]*1)] <- X[(D2[1]*1+1):(D2[2]*1),1:(D1[1]*1)]
  
  cc <- c("red","blue","green","orange","cyan","brown","yellow","black")
  cols <- "grey"
  for(k in 1:K) cols <- c(cols,rep(cc[k],2),"grey")
  cols <- colorRampPalette(cols)(50*K)
  at <- c()
  for(k in 1:K)
    at <- c(at,seq(0,6,length=50)+k*10)
  at <- c(at,(K+1)*10)
  
  image(t(Y),col=cols,breaks=at,axes=F,main=main)
  return(X)
}

#Plotting function for group sparse data collection, paired in one mode
plotModel <- function(toy,res,predY) {
  ptve <- c()
  for(m in 1:2) {
    ptve[m] <- (1-mean((toy$yNew[[1]][[m]]-predY[[m]])^2)/mean(toy$yNew[[1]][[m]]^2))*100
    print(paste0("PTVE in predicting view ",m," is ",round(ptve[m]),"%"))
  }
  
  par(mfrow=c(2,3), mar=rep(4,4),oma=c(0,4,0,0))
  cols <- colorRampPalette(c("white","black"))(100)
  
  s <- list()
  s[[1]] <- lapply(toy$groups[[1]],function(x) {apply(toy$W[[1]][x,],2,var)}) #True structure
  w <- apply(res$posterior$W,2:3,mean)
  s[[2]] <- lapply(toy$groups[[1]],function(x) {apply(w[x,],2,var)}) #Inferred structure
  mains <- paste("Component strength plot",c("(true)","(inferred)"))
  ord <- list()
  for(i in 1:2) {
    s[[i]] <- matrix(unlist(s[[i]]),ncol=length(s[[i]]))
    ord[[i]] <- order(rowSums(s[[i]]),decreasing=TRUE)
    s[[i]] <- s[[i]][ord[[i]],]
    
    image(s[[i]][,3:1],col=cols,xlab="Component",ylab="Data source",axes=F,main=mains[i])
    axis(1,(1:nrow(s[[i]])-1)/(nrow(s[[i]])-1),1:nrow(s[[i]]))
    axis(2, (1:ncol(s[[i]])-1)/(ncol(s[[i]])-1), ncol(s[[i]]):1)
    text(0,1,'*',col=c("red","blue")[i])
    if(i==2) { #Add PTVEs to the plot
      for(m in 1:2) {
        lines(c(1.1,1.1)-m/10,c(0,(3-m)/2)); lines(c(1.08,1.1,1.12)-m/10,c(-.05,0,-.05)+(3-m)/2)
        txt <- paste0("PTVE=",round(ptve[m]),"%")
        if(m==1) txt <- paste("Predicting source 1 | 3,",txt)
        text(1.13-m/10,0,txt,srt=270,adj=1)
      }
      cat("Component strength plot (inferred):
  clearly GFA can infer the correct component affiliations accurately. Components 7 and 8
  have fitted to residual noise, and are significantly weaker than the 'true components'.
  Additionaly PTVEs are shown for doing sequential prediction from new observations of
  data source 3 to sources 1 and 2.\n")
    } else {
      cat("Component strength plot (true):
  dark blocks denote strong activation in component-source pair, whereas white denotes
  no activation. Hence components active in multiple data sources imply dependencies
  between the data sources; these can be used to do prediction from a data source to
  another. E.g. component 2 can be used to predict data source 1 from 3 and vice versa.\n")
    }
  }
  plot(toy$W[[1]][,ord[[1]][1]],xlab="Feature",ylab="Weight",
       main="Strongest component highlighted",col="red",type="l")
  lines(w[,ord[[2]][1]],col="blue",lty=2)
  cc <- cor(toy$W[[1]][,ord[[1]][1]], w[,ord[[2]][1]])
  legend("bottomleft",c("True",paste0("Inferred, cor=",round(cc,2))),col=c("red","blue"),lty=c(1,2),bty="n")
  cat("Strongest component highlighted:
  the loadings of component 1 are shown for each feature. The values for data source 1
  are correctly identified as zeros.\n")
  
  mtext(side=2,"Group sparse factorization",outer=T,line=1, at=0.75, las=3)
  mtext(side=2,"Biclustering",outer=T,line=1, at=0.25, las=3)
}

## DEMO 1: showing that GFA can infer a group sparse factorization correctly,
##         and how it can be used for sequential prediction
set.seed(12345)
toy <- genData(V=1) #Data generation
opts <- getDefaultOpts(bicluster=FALSE) #Model options
res <- gfa(toy$y[[1]],opts=opts,K=10) #Model inferece

## Sequential prediction for the two latter data views,
## based on new observations in the third data view
res$opts$prediction <- c(TRUE,TRUE,FALSE)
predY <- sequentialGfaPrediction(Y=toy$yNew[[1]], model=res)
plotModel(toy,res,predY) #Visualization



## DEMO 2: Inferring sprase decomposition (biclustering structure),
##         where the data collection contains pairings in two data modes
##         (Bottom figures.)
toy <- genData(V=2) #Generate the toy data
X <- showComp(toy) #Plot true bicluster-components

# Set missing values to a part of the data: prediction with reconstruction
omit <- unique(cbind(sample(151:160,30,TRUE),sample(21:30,30,TRUE)))
trueY <- toy$y[[1]][[2]][omit]; toy$y[[1]][[2]][omit] <- NA

opts <- getDefaultOpts(bicluster=TRUE) #Model options
res <- gfa(toy$y,opts=opts,K=c(6,6)) #Model inference
X <- showComp(toy,res) #Plot the inferred components

#Predict the missing values
predY <- reconstruction(res)[[1]] #Reconstruction based on the component structure
predY <- predY[cbind(omit[,1], toy$D[[1]][1]+omit[,2])] #Data sources concatenated

plot(c(trueY),c(predY),xlab="True missing values",ylab="Predicted values",
     main=paste("Prediction: correlation",round(cor(c(trueY),c(predY)),2),"with true values"))
