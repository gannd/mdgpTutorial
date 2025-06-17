# =================================================================================================
# Partial Credit Weighted Overall Accuracy
#	- estimate weighted overall, user's and producer's accuracies based on error severity determined by label agreement
#	- estiamte standard errors and confidence intervals for weighted confusion matrix
#			
# 	Equations derived from: --> add reference
#			 
# packages: {terra},{sf}
# =================================================================================================

# calculate weights for partial credit based on label agreement
prtlCrdtWghts <- function(r_nme,c_nme){
  
  # extract numbers from strings
  require(stringr)
  numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
  }
  
  if (grepl('_x_',r_nme) == TRUE){
    r_clsLst <- unlist(strsplit(as.character(r_nme),'_x_'))
  } else {r_clsLst <- as.character(r_nme)}
  
  if (grepl('_x_',c_nme) == TRUE){
    c_clsLst <- unlist(strsplit(as.character(c_nme),'_x_'))
  } else {c_clsLst <- as.character(c_nme)}
  
  clsNmePrtLst <- as.data.frame(unique(c(gsub('[[:digit:]]','',r_clsLst),gsub('[[:digit:]]','',c_clsLst))))
  names(clsNmePrtLst) <- 'clsCmps'
  
  cls1 <- as.data.frame(as.numeric(numextract(r_clsLst)))
  rownames(cls1) <- gsub('[[:digit:]]','',r_clsLst)
  names(cls1) <- 'prc1'
  cls2 <- as.data.frame(as.numeric(numextract(c_clsLst)))
  rownames(cls2) <- gsub('[[:digit:]]','',c_clsLst)
  names(cls2) <- 'prc2'
  
  # merge class components
  clsNmePrtLst <- merge(clsNmePrtLst,cls1, by.x='clsCmps', by.y='row.names', all=TRUE)
  clsNmePrtLst <- merge(clsNmePrtLst,cls2, by.x='clsCmps', by.y='row.names', all=TRUE)
  
  clsNmePrtLst[is.na(clsNmePrtLst)] <- 0  
  clsNmePrtLst <- transform(clsNmePrtLst, min = pmin(prc1,prc2))
  
  # calculate partial credit weight
  w <- sum(clsNmePrtLst$min)/100
  
  return(w)
}

# calculate weight adjusted accuracy
wghtAdjAcc <- function(cm, wm) {
  
  # convert confusion and weight tables to matrices
  cMtx <-as.matrix(cm); wMtx <- as.matrix(wm)
  
  if (any(dim(cMtx) != dim(wMtx)) == TRUE){ print("Confusion Matrix and Weight Matrix are not of Equal Size") }
  
  #summarize cmx
  rowSum <-apply(cMtx,1,sum); colSum <- apply(cMtx,2,sum)
  
  # confusion matrix and margin proportions
  prp <- cMtx/sum(cMtx)
  colPrp <- colSum/sum(cMtx)
  rowPrp <- rowSum/sum(cMtx)
  
  # check if sums equal 1
  if ( round((sum(rowPrp) + sum(colPrp))/2, 2) != 1){ print("Marginal Proportions don't add to 1") }
  if ( round(sum(prp), 2) != 1){ print("Table Proportions don't add to 1") }
  
  # expected proportions
  prpExp <- rowPrp %o% colPrp
  
  # weighted weights
  rowWghts <- wMtx%*%colPrp
  colWghts <- t(t(wMtx)%*%rowPrp)
  
  # user accuracy - rows
  ua <- apply(wMtx*prp,1,sum)/rowPrp
  uaSD <-sqrt(ua*(1-ua)/rowSum)
  
  # producer accuracy - columns
  pa <- apply(wMtx*prp,2,sum)/colPrp
  paSD <- sqrt(pa*(1-pa)/colSum)
  
  naiveWght <- sum(sum(prp * wMtx))
  naiveWghtVar <- ((naiveWght*(1-naiveWght))/sum(cMtx))
  
  # add references for equations
  naiveWghtExp <- sum(sum(prpExp * wMtx))
  kpWght <- (naiveWght-naiveWghtExp)/(1-naiveWghtExp)
  naiveWghtCmp <- 1 - naiveWght
  naiveWghtExpCmp <- 1 - naiveWghtExp
  thWght <- 0; for (i in 1:nrow(cMtx)) for (j in 1:ncol(cMtx))
    thWght <- thWght + (prp[i,j]*((wMtx[i,j]*naiveWghtExpCmp - (rowWghts[i]+colWghts[j]) * naiveWghtCmp)^2 ))
  kpWghtVar <- (thWght - (naiveWght*naiveWghtExp - 2*naiveWghtExp + naiveWght)^2) / (sum(cMtx) * naiveWghtExpCmp^4)
  
  return(list(	sum.n=sum(cMtx),
               sum.kappa=kpWght, sum.kvar=kpWghtVar, theta=c(naiveWght,naiveWghtExp,thWght),
               sum.naive=naiveWght, sum.var=naiveWghtVar,
               user.wa=ua, prod.wa=pa,
               user.wsd=uaSD, prod.wsd=paSD,
               weights.row=rowWghts, weights.col=colWghts, expected=prpExp))
}

# summarize weight adjusted accuracy estimates
summary.wghtAdjAcc <- function(wghtAdjAcc, alpha=0.05) {
  
  ciw <-function(var, n) {
    qnorm(1-(alpha/2))*sqrt(var) + (1/(2*n))
  }
  print(paste("Number of observations:", wghtAdjAcc$sum.n), quote=F)
  print(paste("Sum of weighted sum of row, column weights:",
              round(sum(wghtAdjAcc$weights.row), 2), ",",
              round(sum(wghtAdjAcc$weights.col), 2) ), quote=F)
  print("Summary of weighted naive statistics", quote=F)
  print(paste(
    "Overall accuracy, stdev, CV%:",
    round(wghtAdjAcc$sum.naive, 4), ",", round(sqrt(wghtAdjAcc$sum.var), 4), ",",
    round((sqrt(wghtAdjAcc$sum.var)/wghtAdjAcc$sum.naive)*1000,0)/10),
    quote=F)
  w <- ciw(wghtAdjAcc$sum.var, wghtAdjAcc$sum.n)
  print(paste(
    round((1-alpha)*100,0),"% confidence limits for accuracy:",
    round((wghtAdjAcc$sum.naive-w),4), "...",
    round((wghtAdjAcc$sum.naive+w),4), sep=""), quote=F)
  print("User weighted accuracy", quote=F)
  print(round(wghtAdjAcc$user.wa,4));
  print("Producer weighted reliability:", quote=F)
  print(round(wghtAdjAcc$prod.wa,4));
  print("Summary of weighted kappa statistics", quote=F)
  print(paste("Overall weighted kappa, stdev, & CV%:",
              round(wghtAdjAcc$sum.kappa,4), ",",
              round(sqrt(wghtAdjAcc$sum.kvar),4), ",",
              round((sqrt(wghtAdjAcc$sum.kvar)/wghtAdjAcc$sum.kappa)*1000,0)/10), quote=F)
  w<-ciw(wghtAdjAcc$sum.kvar, wghtAdjAcc$sum.n)
  print(paste(
    round((1-alpha)*100,0),"% confidence limits for weighted kappa:",
    round((wghtAdjAcc$sum.kappa-w),4), "...",
    round((wghtAdjAcc$sum.kappa+w),4), sep=""), quote=F)
}

# confidence interval function
ci <-function(var,n,alpha) {
  qnorm(1-(alpha/2))*sqrt(var) + (1/(2*n))
}