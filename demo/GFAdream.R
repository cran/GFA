
require("GFA")
#' Preprocess DREAM7 challenge data sources
#'
#' \code{dreamExp} returns a list containing the DREAM7 challenge data
#' used in Bunte et al. (2016).
#' 
#' @param path A string containing the path where the data sources have been
#' stored. The data collection can be downloaded from:
#' https://www.synapse.org/#!Synapse:syn2785778/wiki/70254
#' @return Y A list containing the following DREAM7 challenge data sources:
#'   drug sensitivity, gene expression, methylation, RPPA and exome sequence.
#'   
dreamExp <- function(path="data/") {
  D <- 500 #Number of genes used
  
  #Separate processing script for exome sequence
  processExomeSeqDataset <- function(ex) {
    ex <- ex[ex$Avg.Mismatch.alt. <= 0.5 & ex$MismatchQualitySum.alt. <= 7.5, ]
    cellLines <- unique(ex[,1])
    count <- 0
    exomeStatList <- list()
    geneName <- (ex[which(cellLines[1] == ex[,1]),7])
    for(i in 2:length(cellLines))
      geneName <- union(geneName, (ex[which(cellLines[i] == ex[,1]),7]))
    
    genes <- unique(geneName)
    ex.processed <- matrix(0, nrow=length(cellLines),ncol=length(genes))
    rownames(ex.processed) <- cellLines
    colnames(ex.processed) <- genes
    class(ex.processed) <- "numeric"
    
    for(k in 1:length(cellLines)) {
      dataSubSet <- (ex[which(cellLines[k] == ex[,1]),])
      length(dataSubSet[,7] %in% genes)
      alt.ref.count <- vector('numeric',length=length(dataSubSet[,7]))
      names(alt.ref.count) <- dataSubSet[,7]
      
      for(i in 1:dim(dataSubSet)[1])
        alt.ref.count[i] <- (dataSubSet[i,20]/(dataSubSet[i,20]+dataSubSet[i,19]))
      
      duplicateGenes <- unique(names(alt.ref.count)[duplicated(names(alt.ref.count))])
      for(d in 1:length(duplicateGenes)) {
        index <- which(duplicateGenes[d] == names(alt.ref.count))
        maxValue <- alt.ref.count[index[1]]
        for(a in 2:length(index)) {
          if(alt.ref.count[index[a]] >= maxValue) {
            maxValue <- alt.ref.count[index[a]]
            alt.ref.count[index[a-1]] <- 0
          } else {
            alt.ref.count[index[a]] <- 0
          }
        }
      }
      maxValue2 <- max(alt.ref.count)
      if(maxValue2 < 0.5)
        alt.ref.count[which(maxValue2==alt.ref.count)] <- 0
      ex.processed[which(as.character(cellLines[k]) == rownames(ex.processed)),names(alt.ref.count)] <- alt.ref.count
    }
    return(t(ex.processed))
  }
  
  path <- paste0(path,"DREAM7_DrugSensitivity1_")
  if(!file.exists(paste0(path,"Predictions.csv"))) {
    stop("NCI-DREAM challenge data needs to be downloaded to 'data/' from
  https://www.synapse.org/#!Synapse:syn2785778/wiki/70254 including:
   - DREAM7_DrugSensitivity1_Predictions.csv
   - DREAM7_DrugSensitivity1_Drug_Response_Training.txt
   - DREAM7_DrugSensitivity1_GeneExpression.txt
   - DREAM7_DrugSensitivity1_Methylation.txt
   - DREAM7_DrugSensitivity1_RPPA.txt
   - DREAM7_DrugSensitivity1_Exomeseq.txt.
   - DREAM7_DrugSensitivity1_RNAseq_expressed_calls.txt
   - DREAM7_DrugSensitivity1_SNP6_gene_level.txt ")
    
  } else {
    #Data sources to read
    files <- c("Drug_Response_Training","GeneExpression","Methylation","RPPA",
               "Exomeseq","RNAseq_quantification","SNP6_gene_level")
    names(files) <- c("dss","ge","met","rppa","exome","rna","cnv")
    #All cell lines with observations
    pred <- read.table(paste0(path,"Predictions.csv"),header=T, sep=',',row.names=1)
    cellLines <- rownames(pred)
    Y <- list()
    for(f in 1:length(files)) {
      if(!file.exists(paste0(path,files[f],".txt"))) {
        stop(paste("File",paste0(path,files[f],".txt"),"needs to be downloaded to 'path' from
  https://www.synapse.org/#!Synapse:syn2785778/wiki/70254"))
      }
      header <- T; rown <- 1
      set <- names(files)[f]
      if(set%in%c("cnv")) {header <- F; rown <- NULL}
      if(set%in%c("rppa","rna")) {header <- T; rown <- NULL}
      # Data reading
      if(set=="exome") {
        tmp <- read.delim(paste0(path,files[f],".txt"))
        tmp <- processExomeSeqDataset(tmp)
      } else {
        tmp <- read.table(paste0(path,files[f],".txt"),header=header, row.names=rown)
      }
      if(f==1) #Training cell lines (with observed drug responses) first
        cellLines <- c(cellLines[cellLines%in%rownames(tmp)], cellLines[!cellLines%in%rownames(tmp)])
      else tmp <- t(tmp)
      #Row and column names to correct places
      if(set=="cnv") {cl <- tmp[-1:-2,1]; cn <- tmp[2,-1]; tmp <- tmp[-1:-2,-1]
                      colnames(tmp) <- cn; rownames(tmp) <- cl
                      tmp[c("21NT","184B5","600MPE","21MT1"),] <- NA}
      if(set=="met") {colnames(tmp) <- tmp[1,]; tmp <- tmp[-1:-3,]}
      if(set=="rna") {colnames(tmp) <- tmp[1,]; tmp <- tmp[-1:-2,]}
      if(set=="rppa") {tmp <- tmp[,tmp["FullyValidated",]=='Yes']; colnames(tmp) <- tmp[1,];
                       tmp <- tmp[-1:-2,]}
      if(f>1) {class(tmp) <- "numeric"} else {tmp <- as.matrix(tmp)}
      #Full matrices, non-observed cell lines as NAs
      Y[[set]] <- array(NA, dim=c(length(cellLines),ncol(tmp)), dimnames=list(cellLines, colnames(tmp)))
      Y[[set]][intersect(cellLines,rownames(tmp)),] <- tmp[intersect(cellLines,rownames(tmp)),]
    }
  }
  N <- nrow(Y$dss)
  
  for(set in c("met","rna","cnv")) {
    Y[[set]] <- Y[[set]][,colnames(unique(as.matrix(Y[[set]]), MARGIN=2))]
    Y[[set]] <- unique(Y[[set]], MARGIN=2)
  }
  
  # Dimensionality reduction: common genes, chosen by biggest average
  # (over the views) standard deviations
  for(i in 2:length(Y)) Y[[i]] <- Y[[i]][,sort(colnames(Y[[i]]))]
  fullGenes <- c("ge","exome","met","rna","cnv")
  genes <- Reduce(intersect, lapply(Y[fullGenes],colnames))
  sdev <- array(NA, dim=c(length(fullGenes), length(genes)), dimnames=list(fullGenes,genes))
  for(set in fullGenes) {
    Y[[set]] <- Y[[set]][,genes]
    sdev[set,] <- apply(Y[[set]], 2, sd, na.rm=T)
    sdev[set,] <- sdev[set,] / max(sdev[set,])
  }
  sdev <- apply(sdev, 2, mean, na.rm=T)
  genes <- order(-sdev)[1:D]
  # Include all the genes present in RPPA
  genes <- union(which(colnames(Y$ge) %in% colnames(Y$rppa)), genes)[1:D]
  
  for(set in fullGenes) Y[[set]] <- Y[[set]][,genes]
  common <- which(colnames(Y$rppa) %in% colnames(Y$ge))
  uncommon <- setdiff(1:ncol(Y$rppa), common)
  Y$rppa <- Y$rppa[,c(common,uncommon)]
  
  # Scale DSS, GE, MET and RPPA
  for(set in c("dss","ge","met","rppa")) Y[[set]] <- scale(Y[[set]])
  Y$rna <- NULL; Y$cnv <- NULL #RNA and CNV omitted from the data collection
  
  return(Y)
}


dreamScore <- function(models) {
  reps <- length(models)
  load("data/Dream_DrugSensitivity_allCelllines.RData")
  
  dssPred <- matrix(0, nrow = nrow(dream7.drugResponse.mat), ncol = ncol(dream7.drugResponse.mat))
  for (rep in 1:reps) {
    rec <- reconstruction(models[[rep]])[,1:31]
    dssPred <- dssPred + rec
  }
  ranks <- apply(-dssPred, 2, rank)
  
  if(!file.exists("dreamscore/")) {
    stop("The scoring scripts need to be downloaded to directory 'dreamscore' from
  https://www.synapse.org/#!Synapse:syn2785778/wiki/70258")
  }
  
  txt <- "dream7_drugsensitivity1_predictions_"
  file <- paste0("dreamscore/",txt,"reps",reps,".csv")
  write.csv(ranks,file=file,row.names = rownames(dream7.drugResponse.mat),quote=F)
  #print(file)
  ds <- readLines(file)
  ds[1] <- paste0("DrugAnonID,",paste0(paste0("Drug",1:31),collapse=","))
  writeLines(ds,con=file)
  
  setwd("dreamscore/")
  perlscript <- "weighted_average_concordance_index.pl"
  fname <- paste0(txt,"reps",reps,".csv")
  tmp <- system(paste0("perl ", perlscript," ",fname),intern=T,ignore.stderr=T)
  x <- as.numeric(substr(tmp[38],80,86))
  
  setwd("../")
  print(paste("DREAM score:",x))
  return(x)
}

#Load and preprocess the data
Y <- dreamExp(path="data/")

#Model options (including informative noise prior)
opts <- getDefaultOpts(bicluster=TRUE)
opts$iter.max <- 12000
opts$iter.burnin <- 10000
opts$iter.saved <-  100
opts$verbose <- 1
opts <- informativeNoisePrior(Y,noiseProportion=0.5,conf=1,opts)
#Conservative missing value imputation was used in Bunte et al., as the
#missing value structure in combination with the Bayesian imputation resulted in
#some exaggerated predictions from the sparse exome sequence to gene expression,
#methylation, and drug sensitivity. Rectifying this while keeping the Bayesian
#imputation would have required carefully defined informative priors for the
#model parameters alpha and beta, discouraging overfitting to the exome sequence.
opts$imputation <- "conservative"

cat("Demonstrating with one model run only. Adjust 'reps' to average
prediction over sampling chains.\n")
models <- list()
scores <- c()
reps <- 1
for(rep in 1:reps) {
  print(paste("Running GFA with 60 components; seed",rep))
  set.seed(rep)
  models[[rep]] <- gfa(Y, opts=opts, K=60)
  
  cat(paste0("Single model (rep ",rep,")"))
  scores[rep] <- dreamScore(models[rep])
}
#NOTE: the replication is not exact, as GFA code was improved after publication of Bunte et al.
#to correctly account for the missing values when setting the informative noise prior.
cat(paste0("Averaged model (over ",reps," repetitions)"))
x <- dreamScore(models)

#Illustrate the model (rep 1); show only top 30 features in the component-effect plots,
#as the 500 genes per data source would overcrowd the plot
par(mfrow=c(2,4))
vis <- visualizeComponents(models[[1]],Y,hclust=T,topK=3,topFeatures=30)





