
summariseLoci <- function(cD, perReplicate = TRUE)
  {
    if(perReplicate) {
      colSums(exp(cD@locLikelihoods))
    } else {
      sum(1 - exp(rowSums(log(1-exp(cD@locLikelihoods)))))
    }
  }

.controlFDR <- function(likes, FDR) {
  selnum <- max(which(cumsum((1 - sort(likes, decreasing = TRUE)) / 1:length(likes)) < FDR))
  if(selnum > 0) sellikes <- sort(order(likes, decreasing = TRUE)[1:selnum]) else sellikes <- integer()
  sellikes
}

.controlFWER <- function(likes, FWER) {
  llsum <- likes
  selnum <- max(which(1 - cumprod(sort(llsum, decreasing = TRUE)) < FWER))
  if(selnum > 0) sellikes <- sort(order(llsum, decreasing = TRUE)[1:selnum]) else sellikes <- integer()
  sellikes
}

selectLoci <- function(cD, likelihood, FDR, FWER, perReplicate = TRUE#) {  
                                        , returnBool = FALSE) {

#  returnBool <- FALSE
  if(!missing(likelihood)) {
      if(perReplicate) {
          selLoc <- cD@locLikelihoods > log(likelihood)               
          if(returnBool) return(selLoc) else return(which(rowSums(selLoc) > 0))
      } else {
          selLoc <- 1 - exp(rowSums(log(1 - exp(cD@locLikelihoods)))) > likelihood
          if(returnBool) return(selLoc) else return(rowSums(selLoc) > 0)
      }
  } else {
    if(!missing(FDR)) {
      controlFunction <- .controlFDR
      controlCrit <- FDR
    } else if(!missing(FWER)) {
      controlFunction <- .controlFWER
      controlCrit <- FWER
    } else stop ("No criterion for locus selection given.")    
    if(perReplicate) {
      selRep <- lapply(1:ncol(cD@locLikelihoods), function(jj) controlFunction(exp(cD@locLikelihoods[,jj]), controlCrit))
      if(returnBool) {
        bool <- do.call("cbind", lapply(1:length(selRep), function(ii) {          
          selBool <- rep(FALSE, nrow(cD))
          if(length(selRep[[ii]]) > 0) selBool[selRep[[ii]]] <- TRUE
          selBool
        }))
        colnames(bool) <- colnames(cD@locLikelihoods)
        return(bool)
      }
      selLoc <- sort(unique(unlist(selRep)))
    } else {
      selLoc <- controlFunction(1 - exp(rowSums(log(1 - exp(cD@locLikelihoods)))), controlCrit)
      if(returnBool) {
        bool <- rep(FALSE, nrow(cD))
        bool[selLoc] <- TRUE
        return(bool)
      }
    }
  }
  if(length(selLoc) == 0) stop("No loci found for the given selection criterion.")
  cD[selLoc,]
}
