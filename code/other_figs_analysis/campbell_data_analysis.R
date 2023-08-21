
#campbell_data_analysis.R: code for fly data analysis: Shyam Srinivasan (C), shyam@snl.salk.edu
#this file contains functions to analyze the data provided by Rob Campbell, which provides 
#a summary of data from his paper

#to do: a function which gets all the significant responses and the significant cells for each odor
#peak responses instead of mean responses.
#the no of singicant cells can be a df, while the list of significant cells could be a list of vectors
#to calculate number of SDs above mean that give you alpha use qnorm, and for calculating alpha for given number of SDs use pnorm

#this function reads in a descriptor file 
#
#dirname: current directory name
#pattern: that will differentiate the header from other files
readCampDesc <- function(foldname='.',pattern="*"){
  #just look in this directory and read the files which fit the description of the 
  #for now match those files that are headers: contain the word head and are csv files
  pattern = '(head)+[a-z.]*csv'
  headnames <- traversedir(foldnam = foldname,searchstr = pattern)
  if(length(headnames) == 0) return(headnames)
  #now read these files
  headfiles <- lapply(headnames,function(x) read.csv(file = x,header = T))
  headfiles
}

#this function gets all the header and data files in this directory and constructs the response matrix
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for (header,data) files
readCampMat <- function(foldname='.',fpatterns=c('tseries.*head.*csv','tseries.*data[0-9]+.*csv'),op=1){
  #just look in this directory and read the files which fit the description of the 
  #for now match those files that are headers: contain the word head and are csv files
  #cat('readCampMat')
  headerfile <- list.files(path = foldname,pattern = fpatterns[1])
  datafiles <- list.files(path = foldname,pattern = fpatterns[2])
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile)
  res.ls <- lapply(1:length(datafiles), function(x) {
    data <- read.csv(datafiles[x],header = F)
    tmp <- getCampDataTrialStats(data,header = headerfile,trialno = x)
  })
  res.df <- transposeDF(convertNestedListsDF(res.ls))
  res.df <- insertColDf(res.df,headerfile[,1])
  rownames(res.df) <- 1:length(datafiles)
  names(res.df) <- c('odor',1:(length(res.df[1,])-1))
  #cat(rownames(res))
  res.df
}

#this function checks if all the entries of the header match each other
#returns T if they do, F otherwise
#header: the header structure
isCampHeadSame <- function(header,op=1){
  #except the name all the other columns should be similar
  headermean <- apply(header[,-1], 2, mean) #get the mean of all columns except the 1st
  #compare it to the first row, if its not equal return F
  if ( mean(header[1,-1] == headermean) < 1 ) return(F) 
  T
}

#thsi function gets stats on the header information
#header: the header structure
#stats: 1 - mean
getCampHeadStats <- function(header,stats=1,op=1){
  #except the name all the other columns should be similar
  headermean <- apply(header[,-1], 2, mean) #get the mean of all columns except the 1st
  headermean
}

#this function extracts the odors names or any other entity from the header file
#header: the header structure
#col: by default col 1 which is the odor column
getCampHeadOdors <- function(header,col=1,op=1){
  as.character(header[,col]) #de-factor the odor columsn
}

#this function collapses the time series responses into one quantity like the mean, peak etc
#data: the ddata file
#trialno: vector containing extratime,stimlat,stimdur,offset for this particular trial
#stats: the kind of stats needed, 1 - mean, 2 - peak
#op: 1- return signal, 2 - return base_mean, base_sd, thresh.base 3 - signal, base.mean, base.sd, thresh.base
getCampDataTrialStats <- function(data,header,trialno=1,stats=1,alpha=.01,op=1){
  #figure out the stats for the mean the same way Campbell did
  assignVarVals(names(header[trialno,]),header[trialno,])
  #now you can figure out the period over which you should average
  #based on campbell's code in responseperiodframes and ROI-responsematrix
  st <- (stimlat+offset)/fp
  fin <- st + (stimdur+extratime+offset)/fp
  st <- ceiling(st)
  fin <- ceiling(fin)
  endb <- st - 1 #endof baseline
  base.mean <- round(apply(data[,c(1:endb)],1,mean),4)
  base.sd <- round(apply(data[,c(1:endb)],1,sd),4)
  nosds <- qnorm(1-alpha) 
  thresh.base <- base.mean + nosds*base.sd #the threshold for baseline
  #now calculate the mean/max response
  sig.fn <- switch(stats,mean,max)
  res <- round(apply(data[,c(st:fin)],1,sig.fn),4) #the data signal
  #return the following
  #print(cbind.data.frame(base.mean,base.sd,thresh.base))
  switch(op,res,list(base.mean,base.sd,thresh.base),list(res,base.mean,base.sd,thresh.base))
}


#given the data, header, and trial no, will calculate the base mean, base sd, and signal mean/max, and return them
#the calling function can then make a decision about how this information will be used.



#this function calculates the significance threshold for each of the cells in the odor matrix
#data: the ddata file
#trialno: vector containing extratime,stimlat,stimdur,offset for this particular trial
#alpha is used for calculating the SDs above mean where the significance should be fixed
#meanmax: mean or max signal to be used for calculating significance thresholds 1 - max, 2 - mean
#sigval: the signal value to be returned is based on mean =1  or max=2
#rank: basis for ranking significant responses. 1 - max sig -threshold, 2 - z score, 3 - max signal signal, 4 - mean signal
#op: 1 - dont smooth, 2 - smooth the signals
calcCampDataSig <- function(data,header,trialno=1,alpha=0.01,rank=2,maxmean=1,sigval=1,op=1){
  #figure out the significance the same way Campbell did
  assignVarVals(names(header[trialno,]),header[trialno,])
  #now you can figure out the period over which you should average
  #based on campbell's code in responseperiodframes and ROI-responsematrix
  #cat(stimlat,stimdur,extratime,offset,fp)
  #cat((stimlat+offset),(stimdur+extratime+offset),'\n')
  st <- ceiling((stimlat+offset)/fp)
  stb <- st-1 #for baesline
  fin <- ceiling(st + (stimdur+extratime+offset)/fp)
  #calculate sd, mean, and thresh, and max signal
  resm <- apply(data[,c(1:stb)],1,mean)
  ressd <- apply(data[,c(1:stb)],1,sd)
  nosds <- qnorm(1-alpha) 
  thresh <- resm + nosds*ressd
  sigmax <- apply(data[,c(st:fin)],1,max)
  sigmean <- apply(data[,c(st:fin)],1,mean) #this is for sig.norm

  #now, calculate the cells that cross the threshold
  sign.smooth <- sapply(1:getLength(data), function(x){
    vec <- unlist(data[x,])
    #cat('\n',x,':',vec)
    res.smooth <- lowess(vec,y=NULL,f=5/length(vec))
    res.smooth <- res.smooth$y
    #cat('\nsmoothed: ',res.smooth)
    #cat('\npval:',res.smooth[c(st:fin)],'\t,',1-pnorm(res.smooth[c(st:fin)],mean = resm[x],sd = ressd[x]))
    max(res.smooth[c(st:fin)]) #the data is already smoothed, so we might not need to do this
  })
  res.smooth <- ifelse(sign.smooth >= thresh,1,0) #the data is already smoothed, so we might not need to do this
  #figure out the significant cells
  sig <- switch(maxmean,sigmax,sigmean) #signal for determining significance
  ret.sig <- switch(sigval,sigmean,sigmax) #signal that needs to be returned
  res <- ifelse(sig >= thresh,1,0) #if the max between st:fin exceeds thresh
  #rank these cells by their activity
  sig.norm <- switch(rank,sigmax-thresh,(sigmax-resm)/ressd,sigmax,sigmean)
  rank.sig.norm <- getMatRank(sig.norm,decrease = T)
  #print(cbind.data.frame(resm,ressd,thresh,sigmax,sigmean,res,data[,c(st:fin)]))
  #print(cbind.data.frame(resm,ressd,thresh,max.sig,res,res*rank.sig,max.sig.norm,res*rank.sig.norm))
  #return the ranks and the values used for making those rankings
  list(rank.sig.norm*res,ret.sig*res)
  #list(which(res.smooth>0),length(which(res.smooth>0)))
}

#this function goes through a bunch of data matrices and gives the stats for the different
#ways for ranking signigicantly responding cells.
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for (header,data) files
#alpha is used for calculating the SDs above mean where the significance should be fixed
#rank: is the kind of ranking used for calculating the significant cell ranking
#meanmax: mean or max signal to be used for calculating significance thresholds 1 - max, 2 - mean
#sigval: the signal value to be returned is based on mean =1  or max=2
#rank: basis for ranking significant responses. 1 - max sig -threshold, 2 - z score, 3 - max signal signal
#op: types of return functions
calcCampDataSigList <- function(foldname='.',fpatterns=c('tseries.*head.*csv','tseries.*data[0-9]+.*csv'),alpha=0.01,
                                rank=2,maxmean=1,sigval=1,op=1){
  #just look in this directory and read the files which fit the description of the 
  #for now match those files that are headers: contain the word head and are csv files
  #cat('calcCampMatSigRanks')
  headerfile <- list.files(path = foldname,pattern = fpatterns[1])
  datafiles <- list.files(path = foldname,pattern = fpatterns[2])
  #cat('here',datafiles,list.files())
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile)
  res.ls <- lapply(1:length(datafiles), function(x) {
    data <- read.csv(datafiles[x],header = F)
    tmp <- calcCampDataSig(data,header = headerfile,trialno = x,alpha = alpha,rank = rank,maxmean = maxmean,sigval = 1)
  })
  names(res.ls) <- headerfile[,1] #set the names of various lists to their odor names
  #res.df <- transposeDF(convertNestedListsDF(res.ls))
  #res.df <- insertColDf(res.df,headerfile[,1])
  #rownames(res.df) <- 1:length(datafiles)
  #names(res.df) <- c('odor',1:(length(res.df[1,])-1))
  #cat(rownames(res))
  res.ls
}

#prepare thge data set for linear or other classifiers. The result is to produce the odors
#where each odor is along the row, and the last column is the name of the other odor, And, order
#the rows by bunching all the same odor trials together.
#train: the percentage of data that are to be training samples
#trainop: the option for how training data are chosen, 1 - train % of every odor, 2 - trials are chosen
#randomly
#celltype: types of cells in the DF, 1, significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - specified cells,  
#6 - cells with a certain reliability specified by ctypepars
#ctypepars: the parameter that specifies the reliabilities that we want 
#filterop: which columns to filter: 0 dont filter, 1 - filter out columns that have the same value throughout
#op:1 - just the data in the form above with the rows sorted by the odors
#2: 2 : break them into training and test samples, based on train
#3: break them into training and test data based on loocv
#loocv : specifies the number of the fold that should be chosen for testing, so all other folds for training
campPrepDataSig <- function(data.sig,train=80,trainop=1,filterop=1,celltype=4,ctypepar=c(),loocv=1,op=1){
  #get the data as lists first, then convert to data frame
  resp.lst <- campGenRUCellMat(alldata = data.sig,op=celltype,cellop = 2,params = ctypepar)
  #resp.lst is a list of data frames, where each df represents all the trials for an odor
  tresp.lst <- lapply(1:length(resp.lst),function(i) {
    res.df <- cbind.data.frame(transposeDF(resp.lst[[i]]),rep(names(resp.lst)[i],ncol(resp.lst[[i]])) )
  })
  resp.df <- joinListDFs(tresp.lst)
  #assign the last column as odor, and change all the cell labels to have the letter 'c' so its not just a number
  names(resp.df) <- paste('c',names(resp.df),sep = '')
  names(resp.df)[ncol(resp.df)] <- 'odor'
  #also, the labels have to be factors.
  resp.df <- ConvertDfCols(resp.df,cols = ncol(resp.df),op=3) #set the odor col to be factors
  #Data is ready, now prepare it for the algorithm
  res <- prepDataAlgo(dat.df = resp.df,train = train,filterop = filterop,op=op,loocv = loocv)
  res
}



#given data.lst, the return of campDATasiglist, will add extra trials to it, adding
#op: 1 - sig. cells, 2 - reliable cells, 3 - unreliable cells, 4 all cells, 6 - cells with a specific reliability specified by ctypepar
#ctypepar: specifies the reliabilities that we want to add some noise to.
#notrials: no of new trials to add for every odor
#noispar: the noise that is added to the new cells in the trials
#op: 1 - return added trials, 2 - orig. list and new trials
campAddTrials <- function(data.lst,notrials=1,ctypepar=c(),noisepars=c(0,0.05),op=1){
  #returns the cell positions that are in the specified category
  numtrials <- campGetTrialDetails(alldata = data.lst,op=2)[[1]] #get number of trials
  cells.relresp <- campGetReliabResp(data.lst = data.lst,respop = 1) #get reliability and average value
  cells <- list(as.matrix(cells.relresp[[1]]),as.matrix(cells.relresp[[2]]))
  lenvec <- length(cells[[1]][,1]) #length of the odor vector
  #cat('\ncells',str(cells))
  #makes a data frame with the amount of noise that should be added to each cell, based on noisepars
  cells.new <- lapply(colnames(cells[[1]]),function(x) {
    #get the posns
    posns <- unlist(sapply(ctypepar,function(i) which(cells[[1]][,x] == i)))
    vec.ones <- genVecOnes(lenvec,posns = posns)
    #now determine the trials to add
    reliab.vec <- cells[[1]][,x] * vec.ones #the posns that are going to be generated for every trial
    others.vec <- cells[[2]][,x] * as.numeric(!vec.ones) #the other types of cells that are not going to change
    #cat('\nreliab',cells[[1]][,x],'\n',reliab.vec/numtrials)
    vecs <- lapply(1:notrials,function(i) {#
      reliab.pvec <- sapply(1:length(reliab.vec),function(i) rbinom(1,1,reliab.vec[i]/numtrials[x]))
      reliab.vec <- (cells[[2]][,x] * reliab.pvec)
      reliab.vec <- reliab.vec * (1 + rnorm(length(reliab.vec),mean = noisepars[1],sd = noisepars[2]))
      reliab.vec <- reliab.vec + others.vec
      list(setAboveThreshOne(reliab.vec),reliab.vec)
    })
    #names(vecs) <- rep(x,notrials)
    vecs
  })
  names(cells.new) <- colnames(cells[[1]])
  #ok now add the trials. Each new trial consist of the centroid of the other trials with some noise added to the specified cell
  #also making sure that the noise is not negative 
  res <- flattenLists(cells.new,level = 1)
  res <- switch(op,res,c(data.lst,res))
}

#similar to add trials, but here you modify trials so that the the unreliable cells are closer to reliable cells using the 
#min-max normalization and scaling procedure.
#params: c(min,max) for the mimnimum and mximum values possible 
campScaleTrials <-function(dat.lst,params=c(0,1),op=1){
  #go through each trial, and just scale the values. Only consider those  cells that have a signigicant response.
  res <- lapply(dat.lst,function(x){
    vec <- filterWithVec(x[[2]],x[[1]]) #filter the response vector with the significance vector
    newvec <- normalizeVec(vec,params = params,op=4)
    newvec[newvec==params[1]] <- 0
    list(x[[1]],newvec)
  })
  res
}



# a wrapper function around the classifyDAta function
# dat.ls in the form of a return from campDataSigList
# celltype: 
#train: the percentage of data that are to be training samples, if more than one element, we have to do training on all of them
#trainop: the option for how training data are chosen, 1 - train % of every odor, 2 - trials are chosen
#randomly
#celltype: types of cells in the DF, 1, significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - specified cells
#6 - cells with a certain reliability specified by ctypepars
#ctypepars: the parameter that specifies the reliabilities that we want 
#filterop: which columns to filter: 0 dont filter, 1 - filter out columns that have the same value throughout
#op: the type of learning algorithm to be used. 1-3 - linear: lda, qda, mda, 10 - knn, 20 : binary logistic regression, 
#21 -multinomial logistic regression, 30 - SVM
#resop: the kind of results to be returned. 1: the % of correct, 2 - list(pred_labels,actual values), 3 - list(probability matrix,actual values)
#loocv: 0 - do not do leave one out cross validation, 1 - do it by choosing test stimulation sequentially like 1,2; 3,4; ...
# 2 - choose all possible combinations of test stimuli  
# #normop: 0 - do no normalization, 1 - standard z-score normalization, 2- do min-max normalization
campClassifyData <- function(dat.ls,celltype=4,ctypepar=c(),train=80,trainop=1,filterop=1,resop=1,loocv=0,normop=0,op=10){
  #get the data in DF format  
  data.df <- campPrepDataSig(dat.ls,celltype = celltype,ctypepar = ctypepar,op=1) 
  #if just one training/test ratio
  if(length(train)==1) res <- classifyData(dat.df = data.df,train = train,trainop = trainop,filterop = filterop,op=op,
                                           resop = resop,loocv = loocv,normop = normop)
  else {#iterate through all the training parameter values
    res <- lapply(train,function(x) classifyData(dat.df = data.df,train = x,trainop = trainop,filterop = filterop,op=op,
                                                 resop = resop,loocv = loocv,normop = normop) )
    names(res) <- train
  }
  res
}


#finds out if those odors that are correlated are also hard to tell apart, and vice-versa
#dat.lst: the return of campDAtaSiglist
#classifier: the type of learning algorithm to be used. 1-3 - linear: lda, qda, mda, 10 - knn, 20 : binary logistic regression, 
#21 -multinomial logistic regression, 30 - SVM
#celltype: types of cells in the DF, 1, significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - specified cells
#6 - cells with a certain reliability specified by ctypepars
#ctypepar: the parameter that specifies the reliabilities that we want 
#train: % of data that needs to be used for training
#trainop: the option for how training data are chosen, 1 - train % of every odor, 2 - trials are chosen
#randomly
#filterop: which columns to filter: 0 dont filter, 1 - filter out columns that have the same value throughout
#resop: the kind of results to be returned. 1: the % of correct, 2 - list(pred,true), 3 - list(results,true)
#loocv: 0 - do not do leave one out cross validation, 1 - do it by choosing test stimulation sequentially like 1,2; 3,4; ...
# 2 - choose all possible combinations of test stimuli  
#normop: 0 - do no normalization, 1 - standard z-score normalization, 2- do min-max normalization
#op:1 
campDetCorVsClass <- function(dat.lst,classifier=10,train=80,trainop = 1,filterop=1,celltype=4,ctypepar=c(),loocv=0,normop=0,op=1){
  res.class <- campClassifyData(dat.ls = dat.lst,op=classifier,resop = 3,celltype = celltype,ctypepar = ctypepar,
                                train = train,trainop = trainop,loocv = loocv,normop = normop)
  cat('\ncampDetCorVsC',str(res.class))
  print(res.class)
  res.proc <- processClassResults(res.class)
  res.probmat <- meanDF(res.class[[1]],col = 0) #the prob matrices
  res.corr <- computeListMatCorr(dat.lst,matop = 1)
  #now get the correlation results: fixed so that the results have names
  #correlation vs accuracy
  coracc <- sapply(1:nrow(res.corr),function(i) cor(unlist(res.probmat[i,]),res.corr[i,]) ) 
  names(coracc) <- rownames(res.corr)
  coracc <- cleanNA(coracc)
  #correaltion vs AUC
  corauc <- sapply(1:nrow(res.corr),function(i) cor(unlist(res.proc[[6]][i,]),res.corr[i,]) ) 
  names(corauc) <- rownames(res.corr)
  corauc <- cleanNA(corauc)
  res <- list(coracc,corauc,res.proc,res.probmat,res.corr)
  names(res) <- c('coracc','corauc','res.proc','res.probmat','res.corr')
  res
}
#res8 <- campCorrClassPairs(tst,reliability = list(4:6,1:3,1:6),classifier = 21,simcor = 0.45,dissimcor = 0.2,loocv = 1,normop = 0)


#for a data set gets the correlation vs the classification results for odor-pairs 
#dat.lst: the return of campDAtaSiglist
#classifier: the type of learning algorithm to be used. 1-3 - linear: lda, qda, mda, 10 - knn, 20 : binary logistic regression, 
#21 -multinomial logistic regression, 30 - SVM
#celltype: types of cells in the DF, 1, significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - specified cells
#6 - cells with a certain reliability specified by ctypepars
#ctypepar: the parameter that specifies the reliabilities that we want 
#train: % of data that needs to be used for training
#trainop: the option for how training data are chosen, 1 - train % of every odor, 2 - trials are chosen
#randomly
#filterop: which columns to filter: 0 dont filter, 1 - filter out columns that have the same value throughout
#simcor: the filter for similarity
#dissimcor: the filter for dissimilarity
#loocv: 0 - do not do leave one out cross validation, 1 - do it by choosing test stimulation sequentially like 1,2; 3,4; ...
# 2 - choose all possible combinations of test stimuli  
#normop: 0 - do no normalization, 1 - standard z-score normalization, 2- do min-max normalization
campCorrClassPairs <-function(dat.lst,reliability=c(),classifier=10,train=80,trainop = 1,filterop=1,simcor=0.55,
                              dissimcor=0.15,loocv=0,normop=0,op=1){
  res <- lapply(reliability,function(x){
    corclass <- campDetCorVsClass(dat.lst=dat.lst,celltype = 6,ctypepar = x,classifier = classifier,train = train,
                                  trainop = trainop,filterop = filterop,loocv = loocv,normop=normop)
    cors <- corclass[[5]][lower.tri(corclass[[5]],diag = T)]
    classes <- corclass[[3]][[6]][lower.tri(corclass[[3]][[6]],diag = T)]
    dissimpairs <- which(cors <= dissimcor)
    simpairs <- which(cors > simcor)
    res.dissim <- list(cors[dissimpairs],classes[dissimpairs])
    res.sim <- list(cors[simpairs],classes[simpairs])
    list(res.sim,res.dissim)
  })
  names(res) <- reliability
  res
}



#functions for classification in the figures
classificationcCommands <- function(){
  model <- lda(Species ~ ., data = res2[[1]])
  predictions <- predict(model,res2[[2]])
  
  
  temp1 <- campPrepDataSig(tmp);temp3 <- prepDataAlgo(temp1,op=2)
  model <- mda(odor ~ .,data = temp3[[1]]);predictions <- predict(model,temp3[[2]])
  
  temp4 <- scale(temp1[,-ncol(temp1)],center = T,scale = F)
  temp5 <- cbind.data.frame(temp4,temp1[,ncol(temp1)])
  names(temp5)[ncol(temp5)] <- 'odor'
  temp6 <- prepDataAlgo(temp5,op=2)
  
  #svm
  #mouse:
  res5 <- campPrepDataSig(tst,celltype = 3)
  model_svm <- svm(odor ~ ., data = res5, kernel = "linear", cost = 10, scale = FALSE)
  pred_svm <- predict(model_svm,res5[,-ncol(res5)])
  confusionMatrix(pred_svm,res5[,ncol(res5)])  
  # Accuracy : 0.2             
  # 95% CI : (0.0573, 0.4366)
  
  res6 <- campPrepDataSig(tst,celltype = 4,op=2)
  model_svm <- svm(odor ~ ., data = res6[[1]], kernel = "linear", cost = 10, scale = FALSE)
  pred_svm <- predict(model_svm,res6[[2]][,-ncol(res6[[2]])])
  confusionMatrix(pred_svm,res6[[2]][,ncol(res6[[2]])]) 
  # Accuracy : 0.4             
  # 95% CI : (0.1912, 0.6395)
  
  res6 <- campPrepDataSig(tst,celltype = 2,op=2)
  model_svm <- svm(odor ~ ., data = res6[[1]], kernel = "linear", cost = 10, scale = FALSE)
  pred_svm <- predict(model_svm,res6[[2]][,-ncol(res6[[2]])])
  confusionMatrix(pred_svm,res6[[2]][,ncol(res6[[2]])]) 
  # Accuracy : 0.55            
  # 95% CI : (0.3153, 0.7694)
  
  #fly
  res6 <- campPrepDataSig(tmp,celltype = 4,op=2)
  model_svm <- svm(odor ~ ., data = res6[[1]], kernel = "linear", cost = 10, scale = FALSE)
  pred_svm <- predict(model_svm,res6[[2]][,-ncol(res6[[2]])])
  confusionMatrix(pred_svm,res6[[2]][,ncol(res6[[2]])]) 
  # Accuracy : 0.4286          
  # 95% CI : (0.1766, 0.7114)
  
  res6 <- campPrepDataSig(tmp,celltype = 2,op=2)
  model_svm <- svm(odor ~ ., data = res6[[1]], kernel = "linear", cost = 10, scale = FALSE)
  pred_svm <- predict(model_svm,res6[[2]][,-ncol(res6[[2]])])
  confusionMatrix(pred_svm,res6[[2]][,ncol(res6[[2]])]) 
  # Accuracy : 0.6429          
  # 95% CI : (0.3514, 0.8724)
  
  res6 <- campPrepDataSig(tmp,celltype = 3,op=2)
  model_svm <- svm(odor ~ ., data = res6[[1]], kernel = "linear", cost = 10, scale = FALSE)
  pred_svm <- predict(model_svm,res6[[2]][,-ncol(res6[[2]])])
  confusionMatrix(pred_svm,res6[[2]][,ncol(res6[[2]])]) 
  # Accuracy : 0.1429          
  # 95% CI : (0.0178, 0.4281)
  
  #mouse and flies: nearest neighbour, first is the version from the package. Second is your code.
  res <- campClassifyData(tst,celltype = 3,train = seq(30,90,5),trainop = 1)
  unlist(res)
  # 30         35         40         45         50         55         60         65         70         75         80         85         90 
  # 0.08333333 0.10000000 0.16000000 0.16000000 0.10000000 0.12500000 0.15000000 0.13333333 0.10000000 0.10000000 0.10000000 0.10000000 0.10000000 
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1000000 0.1000000 0.1200000 0.1200000 0.1500000 0.1500000 0.1500000 0.1333333 0.1333333 0.1000000 0.1000000 0.1000000 0.1000000 
  res <- campClassifyData(tst,celltype = 4,train = seq(30,90,5),trainop = 1)
  unlist(res)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1833333 0.1166667 0.2200000 0.2400000 0.2250000 0.2250000 0.2500000 0.2666667 0.3333333 0.3000000 0.3000000 0.3000000 0.4000000 
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  #0.2500000 0.2500000 0.2200000 0.2200000 0.2500000 0.2500000 0.2500000 0.2666667 0.2666667 0.2500000 0.2500000 0.2500000 0.5000000 
  res <- campClassifyData(tst,celltype = 2,train = seq(30,90,5),trainop = 1)
  unlist(res)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.3666667 0.2500000 0.3800000 0.4000000 0.3500000 0.3500000 0.4500000 0.4333333 0.4333333 0.3500000 0.4000000 0.2500000 0.4000000 
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.2833333 0.2833333 0.4000000 0.4000000 0.3750000 0.3750000 0.3750000 0.4333333 0.4333333 0.4500000 0.4500000 0.4500000 0.4000000 
  res <- campClassifyData(tmp,celltype = 2,train = seq(30,90,5),trainop = 1)
  unlist(res)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1714286 0.2142857 0.3214286 0.2857143 0.3333333 0.4285714 0.3809524 0.3809524 0.5714286 0.5714286 0.5714286 0.2857143 0.4285714 
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1428571 0.3571429 0.3571429 0.3571429 0.5238095 0.5238095 0.5238095 0.5238095 0.6428571 0.6428571 0.6428571 0.4285714 0.4285714 
  res <- campClassifyData(tmp,celltype = 3,train = seq(30,90,5),trainop = 1)
  unlist(res)
  # 30         35         40         45         50         55         60         65         70         75         80         85         90 
  # 0.11428571 0.25000000 0.14285714 0.14285714 0.19047619 0.19047619 0.09523810 0.14285714 0.07142857 0.14285714 0.14285714 0.14285714 0.00000000 
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1428571 0.1428571 0.1428571 0.1428571 0.2380952 0.2380952 0.2380952 0.2380952 0.2142857 0.2142857 0.2142857 0.0000000 0.0000000 
  res <- campClassifyData(tmp,celltype = 4,train = seq(30,90,5),trainop = 1)
  unlist(res)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1142857 0.2857143 0.2142857 0.3571429 0.2857143 0.2857143 0.2857143 0.2380952 0.3571429 0.3571429 0.3571429 0.2857143 0.2857143 
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1428571 0.3928571 0.3928571 0.3928571 0.4285714 0.4285714 0.4285714 0.4285714 0.3571429 0.3571429 0.3571429 0.2857143 0.2857143 
  
  #SVM
  res1 <- campClassifyData(tmp,celltype = 4,train = seq(30,90,5),trainop = 1,op=30)
  unlist(res1)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.4857143 0.6428571 0.6428571 0.6428571 0.5238095 0.5238095 0.5238095 0.5238095 0.4285714 0.4285714 0.4285714 0.4285714 0.4285714 
  res2 <- campClassifyData(tmp,celltype = 4,op=30,resop = 3,train = 45)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.75
  #
  res1 <- campClassifyData(tmp,celltype = 2,train = seq(30,90,5),trainop = 1,op=30)
  unlist(res1)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.6857143 0.6785714 0.6785714 0.6785714 0.6666667 0.6666667 0.6666667 0.6666667 0.6428571 0.6428571 0.6428571 0.8571429 0.8571429 
  res2 <- campClassifyData(tmp,celltype = 2,op=30,resop = 3,train = 85)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.8095238
  #
  res1 <- campClassifyData(tmp,celltype = 3,train = seq(30,90,5),trainop = 1,op=30)
  unlist(res1)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.2571429 0.2500000 0.2500000 0.2500000 0.2857143 0.2857143 0.2857143 0.2857143 0.1428571 0.1428571 0.1428571 0.1428571 0.1428571 
  res2 <- campClassifyData(tmp,celltype = 3,op=30,resop = 3,train = 65)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.61
  #
  #mouse
  res1 <- campClassifyData(tst,celltype = 4,train = seq(30,90,5),trainop = 1,op=30)
  unlist(res1)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.3333333 0.3333333 0.3600000 0.3600000 0.3250000 0.3250000 0.3250000 0.5333333 0.5333333 0.4000000 0.4000000 0.4000000 0.4000000 
  res2 <- campClassifyData(tst,celltype = 4,op=30,resop = 3,train = 65)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.8395062
  svm.results <- cbind(c('fly','svm','all cells',0.64,0.75,45),c('fly','svm','reliable',0.85,0.8,85),
                       c('fly','svm','unreliable',0.29,0.61,65),
                       c('mouse','svm','all cells',0.53,0.83,65),c('mouse','svm','reliable',0.63,0.94,65),
                       c('mouse','svm','unreliable',0.23,0.51,65))
  
  res1 <- campClassifyData(tst,celltype = 2,train = seq(30,90,5),trainop = 1,op=30)
  unlist(res1)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.5166667 0.5166667 0.6200000 0.6200000 0.5750000 0.5750000 0.5750000 0.6333333 0.6333333 0.5500000 0.5500000 0.5500000 0.6000000 
  res2 <- campClassifyData(tst,celltype = 2,op=30,resop = 3,train = 65)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.945679
  #
  res1 <- campClassifyData(tst,celltype = 3,train = seq(30,90,5),trainop = 1,op=30)
  unlist(res1)
  # 30        35        40        45        50        55        60        65        70        75        80        85        90 
  # 0.1500000 0.1500000 0.1800000 0.1800000 0.1750000 0.1750000 0.1750000 0.2333333 0.2333333 0.2000000 0.2000000 0.2000000 0.1000000 
  res2 <- campClassifyData(tst,celltype = 3,op=30,resop = 3,train = 65)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.5111111
  
  #logit
  res <- campClassifyData(tmp,celltype = 4,train = seq(30,90,10),op = 21)
  unlist(res)
  # 30        40        50        60        70        80        90 
  # 0.4571429 0.6428571 0.4285714 0.4285714 0.5000000 0.5000000 0.4285714 
  res2 <- campClassifyData(tmp,celltype = 4,op=21,resop = 3,train = 40)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.90625
  
  res <- campClassifyData(tmp,celltype = 2,train = seq(30,90,10),op = 21)
  unlist(res)
  # 30        40        50        60        70        80        90 
  # 0.6571429 0.6785714 0.5714286 0.5714286 0.7857143 0.7857143 0.7142857 
  res2 <- campClassifyData(tmp,celltype = 2,op=21,resop = 3,train = 80)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.9464286
  
  res <- campClassifyData(tmp,celltype = 3,train = seq(30,90,10),op = 21)
  unlist(res)
  # 30        40        50        60        70        80        90 
  # 0.2285714 0.2500000 0.2857143 0.2857143 0.2142857 0.2142857 0.2857143 
  res2 <- campClassifyData(tmp,celltype = 3,op=21,resop = 3,train = 60)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.7619048
  logit.results <- cbind(c('fly','svm','all cells',0.64,0.9,40),c('fly','svm','reliable',0.78,0.95,80),
                         c('fly','svm','unreliable',0.29,0.76,60),
                         c('mouse','svm','all cells',NA,NA,NA),c('mouse','svm','reliable',0.57,0.7,70),
                         c('mouse','svm','unreliable',NA,NA,NA))
  
  # mouse
  res <- campClassifyData(tst,celltype = 4,train = seq(30,90,10),op = 21)
  unlist(res)
  #no result
  res <- campClassifyData(tst,celltype = 2,train = seq(30,90,10),op = 21)
  unlist(res)
  # 30        40        50        60        70        80        90 
  # 0.1666667 0.2400000 0.4000000 0.4000000 0.5666667 0.3500000 0.4000000 
  res2 <- campClassifyData(tst,celltype = 2,op=21,resop = 3,train = 70)
  res3 <- processClassResults(res2)
  res3[[5]]
  #[1] 0.708642
  
  res <- campClassifyData(tst,celltype = 3,train = seq(30,90,10),op = 21)
  unlist(res)
  #no result
  
  lda.results <- cbind(c('fly','svm','all cells',0.64,0.9,40),c('fly','svm','reliable',0.78,0.95,80),
                         c('fly','svm','unreliable',0.29,0.76,60),
                         c('mouse','svm','all cells',NA,NA,NA),c('mouse','svm','reliable',0.57,0.7,70),
                         c('mouse','svm','unreliable',NA,NA,NA))
  
  #LDA
  res <- campClassifyData(tst,celltype = 2,train = seq(30,90,10),op = 1,resop=1)
  unlist(res)
  # 30   40   50   60   70   80   90 
  # 0.55 0.58 0.25 0.25 0.30 0.60 0.40 
  res <- campClassifyData(tst,celltype = 4,train = seq(30,90,10),op = 1,resop=1)
  unlist(res)
  # 30        40        50        60        70        80        90 
  # 0.1833333 0.1400000 0.1250000 0.1250000 0.2000000 0.4000000 0.2000000 
  res <- campClassifyData(tst,celltype = 3,train = seq(30,90,10),op = 1,resop=1)
  unlist(res)
  # 30        40        50        60        70        80        90 
  # 0.1333333 0.0600000 0.1250000 0.1250000 0.1000000 0.2000000 0.0000000 
  res <- campClassifyData(tmp,celltype = 2,train = seq(40,90,10),op = 1,resop=1)
  unlist(res)
  # 40        50        60        70        80        90 
  # 0.6428571 0.6666667 0.6666667 0.5714286 0.5714286 0.5714286 
  res <- campClassifyData(tmp,celltype = 3,train = seq(40,90,10),op = 1,resop=1)
  unlist(res)
  # 40        50        60        70        80        90 
  # 0.2500000 0.2857143 0.2857143 0.2857143 0.2857143 0.2857143 
  res <- campClassifyData(tmp,celltype = 4,train = seq(40,90,10),op = 1,resop=1)
  unlist(res)
  # 40        50        60        70        80        90 
  # 0.3928571 0.5238095 0.5238095 0.5000000 0.5000000 0.4285714 
  
  #doing AUC and other results
  res2 <- campClassifyData(tmp,op=21,resop = 3)
  res3 <- processClassResults(res2)
  
  #svm
  res2 <- campClassifyData(tst,op=30,resop = 3,celltype = 4,train = 80)
  res3 <- processClassResults(res2)
  res3[[5]]
  
  #plotting the the accuracy for each set of classifiers as a barplot
  #f stands for flies, m is mouse
  temp1 <- rbind.data.frame(c(0.6,0.2,0.4),c(0.57,0.28,0.5),c(0.4,0.1,0.5),c(0.64,0.21,0.36),c(0.55,0.2,0.4),c(0.64,0.14,0.42))
  #f: flies, m:mouse
  rownames(temp1) <- c('lda.m','lda.f','knn.m','knn.f','svm.m','svm.f')
  colnames(temp1) <- c('Rel.','Unrel.','all')
  stackedHorzBarPlot(temp1,horz = T,spaces = 0.2,ticknox = 5,fontsize = 0.9,sepwidth = 1)
  #lda_knn_svm_horizontal
  #lda_knn_svm_class_legend

  #plotting the AUC
  temp2 <- rbind.data.frame(c(0.8,0.59,0.64),c(0.72,0.72,0.83),c(0.74,0.39,0.59),c(0.71,0.36,0.65),c(0.88,0.7,0.84),c(0.93,0.66,0.82))
  rownames(temp2) <- c('lda.m','lda.f','knn.m','knn.f','svm.m','svm.f')
  colnames(temp2) <- c('Rel.','Unrel.','all')
  stackedHorzBarPlot(temp2,horz = T,spaces = 0.2,ticknox = 5,fontsize = 0.9,sepwidth = 1,fixx = c(0,2.5)) 
  #stacked_barplot_auc_svm_lda_knn
  #both results together
  res.class <- cbind(temp1,temp2)
  xtable::xtable(res.class) #matrix in latex
  
  #join the campClassifyData, wherein similar odors, their correspnding probabilities are averaged, and then compare this 
  #to the correlation matrix returned by computeListMatCorr
  res2 <- campClassifyData(tst,op=30,resop = 3,celltype = 4,train = 80)
  res3 <- processClassResults(res2)
  res4 <- meanDF(res2[[1]],col = 0)
  res5 <- computeListMatCorr(tst,matop = 1)
  sapply(1:nrow(res4),function(i) cosine(unlist(res4[i,]),res5[i,]) )
  
  #computing the classifier performance of accuracy and AUC versus correlation: Figure S3D: figure classifier AUC validation. Figure XX
  res6 <- campDetCorVsClass(tst,classifier = 10)
  cor(res6[[5]][lower.tri(res6[[5]],diag = T)],res6[[3]][[6]][lower.tri(res6[[3]][[6]],diag = T)])
  
  #comparing how odors that are similar (correlation) have low AUC and vice versa to show AUC's ability to distinguish odors.
  tst1 <- res6[[5]][lower.tri(res6[[5]],diag = T)] # corr
  tst2 <- res6[[3]][[6]][lower.tri(res6[[3]][[6]],diag = T)] #auc
  tst3 <- which(tst2 < 0.25) #get all the AUCs less thatn 0.25, i.e., similar ones
  getStatsLst(list(tst1[tst3],tst2[tst3]))
  #$mean
  #[1] 0.49167220 0.02678571
  fstripchartvecs(list(tst1[tst3],tst2[tst3]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #auc_corr_similar_comparison
  
  tst3 <- which(tst2 > 0.75)
  getStatsLst(list(tst1[tst3],tst2[tst3]))
  # $mean
  # [1] 0.2413992 0.9545455
  fstripchartvecs(list(tst1[tst3],tst2[tst3]),tickno = 2,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #auc_corr_discriminating_comparison
 
  #now compare reliable and unreliable for similar and dissimilar odors
  res6r <- campDetCorVsClass(tst,celltype = 2,classifier = 10)
  tst4 <- res6r[[5]][lower.tri(res6r[[5]],diag = T)]
  tst5 <- res6r[[3]][[6]][lower.tri(res6r[[3]][[6]],diag = T)]
  #dissimilar odors
  tst6 <- which(tst4 < 0.25)
  fstripchartvecs(list(tst4[tst6],tst5[tst6]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #0.07745502 0.85416667
  #knn_dissimilar_reliable_auc_cor
  
  #similar odors
  tst6 <- which(tst4 > 0.55)
  fstripchartvecs(list(tst4[tst6],tst5[tst6]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #0.6329533 0.3416667
  #
  
  #now compare reliable and unreliable for similar and dissimilar odors
  res6ur <- campDetCorVsClass(tst,celltype = 3,classifier = 10)
  tst4 <- res6r[[5]][lower.tri(res6r[[5]],diag = T)]
  tst5 <- res6r[[3]][[6]][lower.tri(res6r[[3]][[6]],diag = T)]
  #dissimilar odors
  tst6 <- which(tst4 < 0.25)
  fstripchartvecs(list(tst4[tst6],tst5[tst6]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #0.07745502 0.85416667
  
  #similar odors
  tst6 <- which(tst4 > 0.55)
  fstripchartvecs(list(tst4[tst6],tst5[tst6]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #0.6329533 0.3416667
  
  #unreliable
  res6ur <- campDetCorVsClass(tst,celltype = 3,classifier = 10)
  tst7 <- res6ur[[5]][lower.tri(res6ur[[5]],diag = T)]
  tst8 <- res6ur[[3]][[6]][lower.tri(res6ur[[3]][[6]],diag = T)]
  tst9 <- which(tst7 < 0.25)
  fstripchartvecs(list(tst7[tst9],tst8[tst9]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #knn_dissimilar_unreliable_auc_cor
  
  tst9 <- which(tst7 > 0.55)
  fstripchartvecs(list(tst7[tst9],tst8[tst9]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  #knn_similar_unreliable_auc_cor
  #0.6329533 0.2666667
  
  #cells with reliability 3,4
  res6.34 <- campDetCorVsClass(tst,celltype = 6,ctypepar = c(3,4),classifier = 10)
  tst7.34 <- res6.34[[5]][lower.tri(res6.34[[5]],diag = T)]
  tst8.34 <- res6.34[[3]][[6]][lower.tri(res6.34[[3]][[6]],diag = T)]
  tst9.34 <- which(tst7.34 < 0.25)
  fstripchartvecs(list(tst7.34[tst9.34],tst8.34[tst9.34]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  
  tst9.34 <- which(tst7.34 > 0.55)
  fstripchartvecs(list(tst7[tst9],tst8[tst9]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  
  #cells with reliability 5,6
  res6.56 <- campDetCorVsClass(tst,celltype = 6,ctypepar = c(5,6),classifier = 10)
  tst7.56 <- res6.56[[5]][lower.tri(res6.56[[5]],diag = T)]
  tst8.56 <- res6.56[[3]][[6]][lower.tri(res6.56[[3]][[6]],diag = T)]
  tst9.56 <- which(tst7.56 < 0.25)
  fstripchartvecs(list(tst7.56[tst9.56],tst8.56[tst9.56]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  
  tst9.56 <- which(tst7.56 > 0.55)
  fstripchartvecs(list(tst7.56[tst9.56],tst8.56[tst9.56]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  
  #cells with reliability 7,8
  res6.78 <- campDetCorVsClass(tst,celltype = 6,ctypepar = c(7,8),classifier = 10)
  tst7.78 <- res6.78[[5]][lower.tri(res6.78[[5]],diag = T)]
  tst8.78 <- res6.78[[3]][[6]][lower.tri(res6.78[[3]][[6]],diag = T)]
  tst9.78 <- which(tst7.78 < 0.25)
  fstripchartvecs(list(tst7.78[tst9.78],tst8.78[tst9.78]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  
  tst9.78 <- which(tst7.78 > 0.55)
  fstripchartvecs(list(tst7.78[tst9.78],tst8.78[tst9.78]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  
  
  
  
  #corr_vs_auc_dissimilar_odors_reliabilty_series_svm
  #corr_vs_auc_dissimilar_odors_reliabilty_series_knn
  #corr_vs_auc_dissimilar_odors_reliabilty_series_svm_fly
  #corr_vs_auc_dissimilar_odors_reliabilty_series_knn_fly
  
  res7 <- campCorrClassPairs(tst,reliability = as.list(8:3),classifier = 30,simcor = 0.45,dissimcor = 0.2)
  fstripchartvecs(c(res7[[1]][[1]],res7[[2]][[1]][2],res7[[3]][[1]][2],res7[[4]][[1]][2],res7[[5]][[1]][2],res7[[6]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  #corr_vs_auc_similar_odors_reliabilty_series_svm
  
  res6 <- campCorrClassPairs(tmp,reliability = as.list(6:1),classifier = 30,simcor = 0.35,dissimcor = 0.2)
  fstripchartvecs(c(res6[[1]][[1]],res6[[2]][[1]][2],res6[[3]][[1]][2],res6[[4]][[1]][2],res5[[5]][[1]][2],res5[[6]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  #corr_vs_auc_similar_odors_reliabilty_series_svm_fly
  
  #with loocv
  res2 <- campCorrClassPairs(tst,reliability = as.list(8:3),classifier = 30,simcor = 0.45,dissimcor = 0.2,loocv = 1)
  fstripchartvecs(c(res2[[1]][[1]],res2[[2]][[1]][2],res2[[3]][[1]][2],res2[[4]][[1]][2],res2[[5]][[1]][2],res2[[6]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  
  #grouping cells and added trials
  #auc_fly_added_trials_12_vs_56
  #auc_fly_added_trials_12_vs_56_change_allcells
  #auc_fly_added_trials_12_vs_56_change_12
  res5 <- campAddTrials(tmp,ctypepar = c(1:3),notrials = 4,op=2) #adding trials
  res8 <- campCorrClassPairs(res5,reliability = list(5:6,2:1),classifier = 30,simcor = 0.45,dissimcor = 0.2,loocv = 1)
  fstripchartvecs(c(res8[[1]][[1]],res8[[2]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  
  #auc_similar_odors_rel_unrel_logit
  res5 <- campAddTrials(tst,ctypepar = c(3:4),notrials = 1,op=2)
  res8 <- campCorrClassPairs(res5,reliability = list(1:4,5:8),classifier = 10,simcor = 0.45,dissimcor = 0.2,loocv = 1)
  fstripchartvecs(c(res8[[1]][[1]],res8[[2]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  
  res5 <- campAddTrials(tst,ctypepar = c(1:8),notrials = 4,op=2)
  res8 <- campCorrClassPairs(tst,reliability = list(5:8,1:4),classifier = 10,simcor = 0.45,dissimcor = 0.2,loocv = 1)
  fstripchartvecs(c(res8[[1]][[1]],res8[[2]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  
  #auc_similar_odors_rel_unrel_knn_mouse
  #auc_similar_odors_rel_unrel_knn_mouse_extra4trials
  #auc_similar_odors_rel_unrel_logit_mouse_extra4trials
  
  res8 <- campCorrClassPairs(tst,reliability = list(5:8,1:4,1:8),classifier = 21,simcor = 0.45,dissimcor = 0.2,loocv = 1)
  fstripchartvecs(c(res8[[1]][[1]],res8[[2]][[1]][2],res8[[3]][[1]][2]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 2,methodstr = 'overplot',semthick = 1.2)
  #auc_similar_odors_rel_unrel_all_logit_mouse;auc_similar_odors_rel_unrel_all_logit_fly
  res8 <- campCorrClassPairs(tmp,reliability = list(4:6,1:3,1:6),classifier = 21,simcor = 0.45,dissimcor = 0.2,loocv = 1)
  
  #accuracy
  HorzBarPlot(temp1[c(1,3,5),c(2,1,3)],ticknox = 2,fixx = c(0,0.6),rndfact = 1,sepwidth = 1,horz = T,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_horz_group_mouse_accuracy
  HorzBarPlot(temp1[c(1,3,5),c(2,1,3)],ticknox = 2,fixx = c(0,0.6),rndfact = 1,sepwidth = 1,horz = F,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_vert_group_mouse_accuracy
  HorzBarPlot(temp1[c(2,4,6),c(2,1,3)],ticknox = 2,fixx = c(0,0.6),rndfact = 1,sepwidth = 1,horz = T,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_horz_group_fly_accuracy
  HorzBarPlot(temp1[c(2,4,6),c(2,1,3)],ticknox = 2,fixx = c(0,0.6),rndfact = 1,sepwidth = 1,horz = F,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_vert_group_fly_accuracy

  #AUC
  HorzBarPlot(temp2[c(2,4,6),c(2,1,3)],ticknox = 2,fixx = c(0,0.8),rndfact = 1,sepwidth = 1,horz = T,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_horz_group_fly_auc
  HorzBarPlot(temp2[c(2,4,6),c(2,1,3)],ticknox = 2,fixx = c(0,0.8),rndfact = 1,sepwidth = 1,horz = F,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_vert_group_fly_auc
  HorzBarPlot(temp2[c(1,3,5),c(2,1,3)],ticknox = 2,fixx = c(0,0.8),rndfact = 1,sepwidth = 1,horz = T,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_horz_group_mouse_auc
  HorzBarPlot(temp2[c(1,3,5),c(2,1,3)],ticknox = 2,fixx = c(0,0.8),rndfact = 1,sepwidth = 1,horz = F,barwidth = 40,op=3,spaces = .1,font=c(2,2),grplabel = T)
  #barplot_vert_group_mouse_auc
  
  
}

#just want to get to the end of the previous function which is too long
classificationCommandsEnd <- function(){

}  

#remove those data cols whose values are constant across rows
#op= 0: nothing to do, 1: return the filtered data frame
#2: return the cols that were removed
filterDataCols <-function(dat.df,op=1){
  if(nrow(dat.df)==1 || op==0) return(dat.df) #just one row, nothing to do
  #get all the cols with constant values
  cols <- sapply(1:ncol(dat.df),function(i){
    vals <- unique(dat.df[,i])
    length(vals)>1
  })
  res.df <- dat.df[,cols]
  #cat('\nsiim cols',which(cols==F),':',ncol(res.df))
  switch(op,res.df,which(cols==F) )
}

#does nearest neighbor classification on a bunch of training vectors, and then can be 
#tested on the test vectors
#dat.ls:lkist(train.data,test.data)
#traininig and testing data is in the form of a dF where rows are the vectors and columns are vector features. last col is the label
#for the vector
nearestNeighborClass <-function(dat.ls,op=1){
  train.dat <- dat.ls[[1]] #the training data
}



#********************************Start here ****************************** Liam

#the subsequent functions are the ones that should be used to understand the campbell data, preferebaly, as they use 
#the significance values from campbell's own mat files. 

#sigpval: gives the value of the significant responses for each cell across all trials
#tseries_dataX.csv: the responses for trial X
#

#this function goes through the tseries_sigpvals file and gets the significant cells for each of the trials
#returns a a list where each element contains the following vectors, the sigcells, the response values, and basemean,base sd, and z-scores
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for the sigpvals file
#alpha is used for calculating the SDs above mean where the significance should be fixed
#rank: is the kind of ranking used for calculating the significant cell ranking
#meanmax and sigval and ranking too are not needed here
#op: types of return functions
getCampDataSigList.old <- function(foldname='.',fpatterns=c('tseries.*head.*csv','tseries.*sigpvals.*csv','tseries.*data[0-9]+.*csv'),alpha=0.01,
                                op=1){
  #just look in this directory and read the files which fit the description of the 
  #for now match those files that are headers: contain the word head and are csv files
  #cat('calcCampMatSigRanks'),
  headerfile <- list.files(path = foldname,pattern = fpatterns[1])
  sigpfile <- list.files(path = foldname,pattern = fpatterns[2])
  datafiles <- list.files(path = foldname,pattern = fpatterns[3])
  #cat('here',datafiles,list.files())
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile)
  if(length(sigpfile) == 1) sigpfile <- read.csv(file = sigpfile,header = T)
  #cat(str(sigpfile),str(headerfile))
  #cat(datafiles,'\n',sortStringNos(datafiles))
  datafiles <- sortStringNos(datafiles)
  datalen <- checkCampMinLength(sigpfile,datafiles)
  #and get the datarows that have data for all timeframes
  res.ls <- lapply(1:getLength(headerfile), function(x){
  #res.ls <- lapply(3:4, function(x){
    #ok, get all the elements below alpha, and then rank them
    vec <- unlist(sigpfile[x,])
    res <- ifelse(vec <= alpha,1,0) #gets the sig. cells
    rank.sig <- getMatRank(1-vec,decrease = T) #have to reverse tank as no. closer to 1 is higher sign.
    names(rank.sig) <- 1:length(rank.sig)
    names(res) <- 1:length(rank.sig)
    #cat('\n',x,'datafile',datafiles[x],'\n')
    data.vec <- getCampDataTrialStats(data=read.csv(datafiles[x],header = F),headerfile,trialno=x,stats=1,op=3)
    #adjust -ve values of significant cells with the thresh. in data.vec
    data.adj <- adjustCampNegValues(data.vec,res,datalen) #adjusts negatvive values for those sig. cells whose response is -ve
    #cat('data.vec',str(data.vec),'data.adj',str(data.adj))
    #list(res,data.vec[[1]],data.vec[[2]],data.vec[[3]],data.vec[[4]],(data.vec[[1]]-data.vec[[2]])/data.vec[[3]]) #the sig. cells, response values, basemean,base sd, sig. as z-score
    list(res[1:datalen],data.adj[[1]],data.adj[[2]],data.adj[[3]],(data.adj[[1]]-data.adj[[2]])/data.adj[[3]]) #the sig. cells, response values, basemean,base sd, sig. as z-score
  })
  cat(dim(res.ls))
  names(res.ls) <- headerfile[data.entries,1] #set the names of various lists to their odor names
  res.ls
}

#this function goes through the tseries_sigpvals file and gets the significant cells for each of the trials
#returns a a list where each element contains the following vectors, the sigcells, the response values, and basemean,base sd, and z-scores
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for the sigpvals file
#alpha is used for calculating the SDs above mean where the significance should be fixed
#rank: is the kind of ranking used for calculating the significant cell ranking
#meanmax and sigval and ranking too are not needed here
#op: types of return functions
#dirop: whether the foldname is op=1, '.', 2 : relative, or 3: absolute e.g., /home/shyam....
getCampDataSigList <- function(foldname='.',fpatterns=c('tseries.*head.*csv','tseries.*sigpvals.*csv','tseries.*data[0-9]+.*csv'),alpha=0.01,
                               op=1,dirop=1){
  #just look in this directory and read the files which fit the description of the 
  #for now match those files that are headers: contain the word head and are csv files
  headerfile <- list.files(path = foldname,pattern = fpatterns[1])
  sigpfile <- list.files(path = foldname,pattern = fpatterns[2])
  datafiles <- list.files(path = foldname,pattern = fpatterns[3])
  if(foldname != '.' || dirop==2) {#we are not in the directory, so save the current one and set the path
    currentdir <- getwd()
    newdir <- paste(currentdir,'/',foldname,'/',sep = '') #go to the actual directory
    setwd(newdir)
  }
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile)
  if(length(sigpfile) == 1) sigpfile <- read.csv(file = sigpfile,header = T)
  datafiles <- sortStringNos(datafiles)
  #in some files, not all tseries_data have data for all time frames. So, pick the ones that have the 
  #highest number of timeframes and only analyze those
  datalen <- checkCampMinLength(sigpfile,datafiles) #datalen is the size of the timeframe
  #and get the datarows that have data for all timeframes
  #old deprecated: datarows <- sapply(datafiles, function(x) nrow(data.table::fread(x, select = 1L)))
  datarows <- sapply(datafiles, function(x) nrow(read.csv(x,header = F) ))
  if( checkVectorSimilar(datarows,op=3) < getLength(headerfile) ) data.entries <- checkVectorSimilar(datarows,op=2)
  else data.entries <- 1:getLength(headerfile)
  res.ls <- lapply(data.entries, function(x){
    #ok, get all the elements below alpha, and then rank them
    vec <- unlist(sigpfile[x,])
    res <- ifelse(vec <= alpha,1,0) #gets the sig. cells
    rank.sig <- getMatRank(1-vec,decrease = T) #have to reverse rank as no. closer to 1 is higher sign.
    names(rank.sig) <- 1:length(rank.sig)
    names(res) <- 1:length(rank.sig)
    #cat('\n',x,'datafile',datafiles[x],'\n')
    data.vec <- getCampDataTrialStats(data=read.csv(datafiles[x],header = F),headerfile,trialno=x,stats=1,op=3)
    #adjust -ve values of significant cells with the thresh. in data.vec
    data.adj <- adjustCampNegValues(data.vec,res,datalen) #adjusts negatvive values for those sig. cells whose response is -ve
    #cat('data.vec',str(data.vec),'data.adj',str(data.adj))
    #list(res,data.vec[[1]],data.vec[[2]],data.vec[[3]],data.vec[[4]],(data.vec[[1]]-data.vec[[2]])/data.vec[[3]]) #the sig. cells, response values, basemean,base sd, sig. as z-score
    list(res[1:datalen],data.adj[[1]],data.adj[[2]],data.adj[[3]],(data.adj[[1]]-data.adj[[2]])/data.adj[[3]]) #the sig. cells, response values, basemean,base sd, sig. as z-score
  })
  names(res.ls) <- headerfile[data.entries,1] #set the names of various lists to their odor names
  if(foldname != '.' || dirop==2)  setwd(currentdir) #set the directory back
  res.ls
}

#this function checks between sigpfile and datafiles to see which one has fewer entures. The 
#for a couple of the files, sigpfile has fewer entries than the datafile and vice versa
checkCampMinLength <- function(sigdata,datafiles,op=1){
  #cat('checkCampMinLength',str(sigdata))
  if(length(sigdata) == 0) return(0)
  #gets the size of each data row
  #cat('\ndatafiles',str(datafiles))
  #old deprecatd: datarows <- sapply(datafiles, function(x) nrow(data.table::fread(x, select = 1L)))
  datarows <- sapply(datafiles, function(x) nrow(read.csv(x,header = F) ))
  sigsize <- length(sigdata[1,])
  mindata <- checkVectorSimilar(datarows) #min(datarows)
  #cat('checkcampminlength ',mindata,dim(sigdata),'\n',datarows,'\n')
  ifelse(mindata <= sigsize,mindata,sigsize)
}

#this function takes in the data vector and adjusts it so that if the cell is significant, the mean response is at least at the threshold level
#also afjust it so that all NAs are set to 0s, or better yet, just truncate the matrices accordingly
#op =1 truncate the matrices to the one with the shortest length, 2 - set NAs to 0
#datalen: the minimum data length
adjustCampNegValues <- function(datavec,sig.cells,datalen,op=1){
  #cat(length(datavec[[1]]),length(datavec[[4]]),length(sig.cells))
  #go through all sig. cells or data vectors, whichever are shorter
  res <- sapply(1:length(sig.cells), function(x){
    #if a sig. cells response is below threshold, set it to thresh
    #cat('\t',datavec[[1]][x],'\t',datavec[[4]][x],'\t',sig.cells[x],'\n')
    if ((sig.cells[x]==1) & (datavec[[1]][x] < datavec[[4]][x]) ) tmp <- datavec[[4]][x]
    else tmp <- datavec[[1]][x]
    tmp
  })
  #cat('\n',res)
  posns <- 1:datalen
  list(res[posns],datavec[[2]][posns],datavec[[3]][posns],datavec[[4]][posns])
}

#this function given a data set, the return of getCAmpDataSigList, will give you two lists
#list(all significant cell response means, all significant cells base mean)
campGetDataSigBase <- function(dat,op=1){
  #get the reliable and unreliable cells for the data set
  dat.sigcalls <- campFreqCellsAnalog(data = dat)
  #list 2 gives the frequency
  #names(dat gives the names of the odors)
  #so get all cells >= 1 and do their resp and base mean: go nuts
  trial.names <- names(dat) #get all the trial names
  #each trial, go through, get the signal mean and base mean as two separate lists
  res.ls <- lapply(trial.names, function(x){
    #at this point just get the list of all significant cells
    #if you want reliable vs unreliable, do that later
    sig.resp <- dat[[x]][[1]]
    sig.signal <- dat[[x]][[2]]
    sig.base <- dat[[x]][[3]]
    #cat('\n',x,str(dat[[x]]))
    signals <- sig.signal[which(sig.resp==1)]
    basemeans <- sig.base[which(sig.resp==1)]
    list(signals,basemeans)
  })
  res.ls
}


#this function takes the return of getCampDataSigList, and gets #signal cells per trial and #reliable cells and #non realiable cells
#data: the return of getCampDataSigList
#op=2, data frame, 1 -  means
getCampDataRel <- function(data,op=1){
  res <- getGroupedListsDf(data,sel = 1,op=1) #sel =1 is fgrouped based on names
  res.ls <- lapply(res, function(x){
    notrials <- length(x)/2
    #sig cells per trial
    sig.cells <- apply(x, 2, sum)
    #reliable cells, majority of trials here is > notrials/2 i.e., only 4-6 for 6 trials
    rel.cells <- apply(x, 1, sum)
    nonrel.cells <- which(rel.cells <= notrials & rel.cells >0 )
    rel.cells <- which(rel.cells > notrials)
    c(mean(sig.cells),length(rel.cells),length(nonrel.cells))
  })
  #cat('getRel')
  #print(res.ls)
  res.df <- convertNestedListsDF(res.ls)
  res.df <- transposeDF(res.df)
  row.names(res.df) <- names(res.ls)
  names(res.df) <- c('sig.','rel.','unrel.')
  switch(op,apply(res.df,2,mean),res.df)
}

#function that gets the frequency list, reliable cells, and unreliable cells for each odor
#cell is reliable if it has more than trials/2 no of significant responses
#data: the return of getCampDataSigList
#op, 1 - look at the directory
#op, 2 - get the headder stuff from the data
campClassifyCells <- function(data,foldname='.',op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  #incorrect, no of trials cann differ across odors
  if(op==1) notrials <- campGetTrialDetails(foldname = foldname)[[1]]#length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  if(op==2) notrials <- campGetTrialDetails(alldata = data,op=2)[[1]]
  #cat('\ncampClass',notrials)
  #have to calculate the freq. list for all the odors, which is a single list
  res.ls <- lapply(names(res.sig), function(x){
    #first the freq cell list
    freq.cells <- apply(res.sig[[x]], 1, sum)
    #reliable and non-reliable cells. This is the place to code what majority of trials
    #majority of trials is > notrials/2
    nonrel.cells <- which(freq.cells <= notrials[x]/2 & freq.cells >0 )
    rel.cells <- which(freq.cells > notrials[x]/2)
    list(freq.cells,rel.cells,nonrel.cells)
  })
  freq.lst <- lapply(res.ls,'[[',1) #each list item contains a frequency vector, and vector names are the cell no
  rel.lst <- lapply(res.ls,'[[',2) #each list item contains the reliable cells for that odor
  nonrel.lst <- lapply(res.ls,'[[',3)
  res <- list(freq.lst,rel.lst,nonrel.lst)
  res
}

#function that gets the same results as campClassifyCells in addition to a list of DFs of reliable cells, and unreliable cells for each odor trial
#cell is reliable if it has more than trials/2 no of significant responses
#data: the return of getCampDataSigList
#trials: the number of trials for every odor. If a single number all the odors have the same number of trials.
#trialop: tells you where you should get the trial data from 1- from the file, 2 - from the getCampDataSiglist data structure
#op, 1 - 
campClassifyTrialCells <- function(data,foldname='.',trials=c(),trialop=2,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  trialdata <- campGetTrialDetails(foldname = foldname,alldata = data,op=trialop)
  #incorrect, no of trials cann differ across odors
  if(length(trials)>0) {#the trial vector listing nu,ber of trials for every odor
    if(length(trials)==1) notrials <- rep(trials,length(names(res.sig)))
    else notrials <- trials
    #cat('\ntrials',notrials)
  } 
  else notrials <- trialdata[[1]]#no of trials, assumes that all odors have same # of Trials
  names(notrials) <- names(res.sig)
  nocells <- trialdata[[3]]
  #have to calculate the freq. list for all the odors, which is a single list
  res.ls <- lapply(names(res.sig), function(x){
    #first the freq cell list
    freq.cells <- apply(res.sig[[x]], 1, sum)
    #reliable and non-reliable cells. This is the place to code what majority of trials
    #majority of trials is > notrials/2
    nonrel.cells <- which(freq.cells <= notrials[x]/2 & freq.cells >0 )
    rel.cells <- which(freq.cells > notrials[x]/2)
    #cat('\nrel',freq.cells,'\t',notrials,':',rel.cells,'rel.test',notrials[x]/2,x,':',which(freq.cells>3))
    #get the number of rel. and unrel. cells in each trial, percentage wise
    trialres <- sapply(1:ncol(res.sig[[x]]),function(y){
      #cat(str(trial.rel.cells),str(nocells))
      trial.rel.cells <- sum(res.sig[[x]][rel.cells,y])
      trial.nonrel.cells <- sum(res.sig[[x]][nonrel.cells,y])
      trial.freq.cells <- sum(res.sig[[x]][freq.cells,y])
      #cat('\nEach cell',c(trial.rel.cells,trial.nonrel.cells),';;',res.sig[[x]][nonrel.cells,y],';;',res.sig[[x]][rel.cells,y])
      c(trial.rel.cells,trial.nonrel.cells)/(nocells/100) #make it a percentage
    })
    #cat('\ntrialres',as.matrix(trialres))
    trialres <- t(trialres)
    colnames(trialres) <- c('rel.','unrel.')
    list(freq.cells,rel.cells,nonrel.cells,trialres)
  })
  freq.lst <- lapply(res.ls,'[[',1) #each list item contains a frequency vector, and vector names are the cell no
  rel.lst <- lapply(res.ls,'[[',2) #each list item contains the reliable cells for that odor
  nonrel.lst <- lapply(res.ls,'[[',3) #unreliable cells for that odor
  trialres.lst <- lapply(res.ls,'[[',4)#number of rel. and unrel. cells in each trial, percentage wise
  #cat('\n',names(res.ls),names(res.sig))
  names(freq.lst) <- names(res.sig)
  names(rel.lst) <- names(res.sig)
  names(nonrel.lst) <- names(res.sig)
  names(trialres.lst) <- names(res.sig)
  res <- list(freq.lst,rel.lst,nonrel.lst,trialres.lst)
  names(res) <- c('freq.list','list rel.','list unrel.','percent rel-unrel/trial')
  res
  #returns a frequency list, the reliable cells for each odor, the unreliable cells, and the percent of rel. and unrel. cells per trial
}

#function that gets 3 lists ripe for fitting with a Gamma dist.: response magnitudes, reliability, and overlap
#dat.ls: the getCampDataSiglList like strucrure
#op:1
campGetFitVecsData <- function(dat.ls,op=1){
  notrials <- max(campGetTrialDetails(alldata = dat.ls,op=2)[[1]])
  
  #responses
  res.resp <- getAboveThresh(unlist(lapply(dat.ls,'[[',2)))
  
  #reliability
  cells.info <- campClassifyTrialCells(data = dat.ls)[[1]]  
  res.rel <- getAboveThresh(unlist(cells.info)/notrials)
  
  #overlap
  overlap <- campComputeOverlap(dat.ls,op=2)
  res.overlap <- joinListDFs(overlap)[,2]
  
  list(res.resp,res.rel,res.overlap)
}

#counts the average number of cells in each reliability class given the getCampDataSigList structure
#returns op=1: list(the average number of cells in each reliability class as a vector, no trials)
#2: list(avg. no od cells,sem of #cells,notrials)
countNoRelClass <- function(dat.ls,op=1){
  #cells.info contains a list of odor vectors, where each vector is the reliability level of that cell
  cells.info <- campClassifyTrialCells(data = dat.ls)[[1]]  
  trialinfo <- campGetTrialDetails(alldata = dat.ls,op=2)
  notrials <- max(trialinfo[[1]])
  noodors <- length(trialinfo[[1]])
  nocells <- trialinfo[[3]]
  percentfactor <- 100/nocells
  # cat('\ncell.infor',str(cells.info))
  # cat('\nfirst cell info',cells.info[[1]])
  res <- sapply(1:notrials, function(i) mean(sapply(cells.info, function(x) length(which(x==i)) ) ) )
  ressem <- sapply(1:notrials, function(i) sd(sapply(cells.info, function(x) length(which(x==i)) ) )/sqrt(noodors) )
  switch(op,list(res*percentfactor,notrials),list(res*percentfactor,ressem*percentfactor,notrials))
}



#comapres the overlap between cells of a certain reliability class.
campCompRelOverlap <- function(dat.ls,op=1){
  cells.freq <- campClassifyTrialCells(data = dat.ls)[[1]]  
  trialdata <- campGetTrialDetails(alldata = dat.ls,op=2)
  notrials <- max(trialdata[[1]])
  nocells <- trialdata[[3]]
  
  res.ls <- campCalcOverlap(cells.freq = cells.freq,notrials = notrials,nocells = nocells,op=1)
  pred.ls <- campCalcOverlap(cells.freq = cells.freq,notrials = notrials,nocells = nocells,op=2)
  
  res.comp <- sapply(1:notrials, function(i) sum(upper.tri(res.ls[[1]])*res.ls[[i]]))/sapply(1:notrials, function(i) sum(upper.tri(res.ls[[1]])*pred.ls[[i]]))
  
  list(res.ls,pred.ls,res.comp)
  
}

#given the cells freq (a vector of the freq of response of that cell), notirals, and no cells will 
#give you either the actual overlap or theoretical overlap based on a random shuffle
#op: 1, actual overlap, 2 - theoretical overlap
campCalcOverlap <- function(cells.freq,notrials,nocells,op=1){

  fnop <- switch(op, function(a,b) length(intersect(a,b)),
                 function(a,b) length(a)/nocells * length(b)/nocells * nocells )
  
  res.ls <- lapply(1:notrials, function(i){
    res <- sapply(1:length(cells.freq), function(j)
      over <- sapply(1:length(cells.freq), function(k){
        rel1 <- which(cells.freq[[j]]==i)
        rel2 <- which(cells.freq[[k]]==i)
        over.no <- fnop(rel1,rel2) #length(intersect(rel1,rel2))
        #theoretical prediction: p(cell/reli) * p(cell/relia) * nocells
        #pred.no <- length(rel1)/nocells * length(rel2)/nocells * nocells 
      })
    )
  })
  res.ls
}


#function that returns a list of odors, where each odor instead of a frequency list contains
#a list whose length is the number of trials, and for each of these entries specifies 
#the cells that had that many significant responses, e.g., 1 would be all cells with sig. responses for one trial
#cell is reliable if it has more than trials/2 no of significant responses
#data: the frequnecy list for all odors. for each odor the cell and the number of sig. responses
#the return of campClassifyCells
#data: results of getCAmpDataSigList
#op, 1 - 
campFreqCellsAnalog <- function(data,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  #not correct, have to use the new fn campgettrialsdetails
  notrials <- campGetTrialDetails(alldata = data,op=2)[[1]]#length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  #have to calculate the freq. list for all the odors, which is a single list
  freq.lst <- lapply(res.sig, function(x){
    #first the freq cell list
    freq.cells <- apply(x, 1, sum)
  })
  #cat('campFreq',str(freq.lst[names(notrials)[1]]) )
  #ok, lets go through the freqnecy list and get freq of cells in each reliability class
  res.ls <- lapply(names(freq.lst), function(x){
    #cat(str(freq.lst),'\n',notrials[x],x,str(notrials))
    tmp.ls <- lapply(1:notrials[x], function(y) which(freq.lst[x][[1]]==y))
  })
  names(res.ls) <- names(freq.lst)
  list(res.ls,freq.lst)
}

#result: from an observation of the frequency list, where there are multiple cells that are active for all trials, it is obvious that either
#this is artifact or the firing of kenyon cells is not random. On way to check between the two would be to see if these cells are acctive for
#multuple odors. If they are, it is probably artifact, and those cells should be excluded. If not,maybe things are not random. 

#function that gets the overlap for odor pairs for each of the reliability scores
#data: is the return of getcampDataSiglist
#skip: the odor name to skip
#op, 1 - 
campOverlapAnalog <- function(data,skip='paraffin',op=1){
  tmp <- campFreqCellsAnalog(data = data)
  #get the analog counts and freq of responses
  analogcells <- tmp[[1]] #for each significance levels, gives the cell numbers that are at that level for each odor
  freq <- tmp[[2]]
  odors <- names(analogcells) #names of the odors
  odorno <- length(odors) #no of odors
  notrials <- campGetTrialDetails(data)[[1]] #a vector of no of trials for each odor
  odorpairs <- lapply(1:odorno^2, function(x){
    c(floor((x-1)/odorno)+1,(x %% odorno+1))
  })
  #now, go through the list while skipping paraffin oil
  res.ls <- lapply(odorpairs, function(x){
    campOverlapPairAnalog(odors[x[1]],odors[x[2]],analogcells,freq,notrials)
  })
  #just get the odor pair names
  res.names <- sapply(odorpairs, function(x){
    paste(odors[x[1]],odors[x[2]],sep = ',')
  })
  names(res.ls) <- res.names
  #if you want to take out paraffin, this is the point to do it.
  res.ls
}

#given a pair of odors, will get the overlap for cells within each class of 
#significance
#odor1,odor2: the index of the odor pairs to consider in the 
#analogcells and frequency list lists
#notrials:the total no of trials
campOverlapPairAnalog <- function(odor1,odor2,analogcells,freq,notrials,op=1){
  #go through each level of significance and calculate the overlap
  #cat('here')
  #cat(odor1,odor2,analogcells[[odor1]][[1]])
  #cat(str(analogcells),str(freq))
  #this a symmetric function even if the odors have different numbers of trials
  #e.g., if odor1 has 3 trials and odor2 has 6 trials, we have to calculate the significance
  #for 3 levels of significance, or if it's the other way around. Doesnt make sense to calculate
  #for greater than 3 as the other cell is not that active
  res.ls <- sapply(1:min(notrials[c(odor1,odor2)]), function(x){
    #first make sure that neither of the odors have no cell of this significance
    if(length(analogcells[[odor1]][[x]])>0 && length(analogcells[[odor2]][[x]])>0){
      overlapcells <- intersect(analogcells[[odor1]][[x]],analogcells[[odor2]][[x]])
      #cat('\novverlap cells',overlapcells)
      overlap <- sapply(overlapcells, function(y){
        prob <- freq[[odor1]][y]*freq[[odor2]][y]/(notrials[odor1]*notrials[odor2])#odorfreq
      })
      overlap <- sum(unlist(overlap))
    }
    else overlap <- 0
    overlap
  })
  res.ls
}

#function goes through all the cells for any pair of odors and gets their reliability
#scores
#data: the returns of the getCampDataSigList function
#nolevels: the number of reliability levels = notrials*nolevels
#relscore: 1, average of both odors, 2 - reliability of first odor
#trialop: the op for campGetTrialDetails, 1 - go through the directory, 2 - use data to get trials. default = 2
#skip: odors to skip
#odorpairs: from the other side, only do these odor pairs
#op=1, overlap; 2, difference, 3 - probability difference
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#6 - given a cell, tells you if it is active for both odors of the odor pair
#7 - gives you the overlap for all cells that have a positive overlap score
#8 - gives you the probabilities difference, i.e., weighted by prob(odor1)
#9 - probability overlap, all cells,
#10 - 1 * prob of resp for odor2
#11 - resp. is freq_i*resp_i, and overlap is r_i*r_j/(r_i + r_j)^2
#12 - probability overlap, all cells, report overlap for 0 reliability or silent cells too
#13 - dot product of the two vectors
#14 - cosine distance
#15 - fractional hamming distance distance
#16 - normed hamming distance
#17 - euclidean distance
#18 - pearson's correlation
#cells: selectively do overlap only for the specified cells. 1 - default do all cells, if cells = list(c(cell1,cell2,...),,..). Specifies the cell numbers
#or whatever index is used to differentiatte cells
campOverlapAllAnalog <- function(data,nolevels=2,skip=c('paraffin'),odorpairs=c(),relscore=2,trialop=2,cells=1,op=1){
  notrials <- campGetTrialDetails(alldata = data,op=trialop)[[1]] #no of trials, assumes that all odors have same # of Trials
  #first, get the reliability and response rates as two different lists
  resprel <- campGetReliabResp(data.lst = data,op=2)
  assignVarVals(c('freq.lst','resp.lst'),resprel)
  odors <- getSkipStrIndex(names(freq.lst),skip = skip,op = 2) #get odors without skip
  if(length(odorpairs)>0) {#if odorpairs are already chosen, then pick the corresponding numbered pairing
      pair.t <- lapply(convertMatVec(length(odors)),function(i) paste(odors[i[1]],odors[i[2]],sep = ',') %in% odorpairs  )
      pairs <- convertMatVec(length(odors))[which(unlist(pair.t) == T)]
  }
  else  pairs <- convertMatVec(length(odors)) #get the pairing of odors 
  res.ls <- lapply(pairs, function(x){#calculate overlap/difference for all odorpairs
    campOverlapAllPairAnalog(odor1=odors[x[1]],odor2=odors[x[2]],freq=freq.lst,resp = resp.lst,cells=cells,
                             notrials=notrials,nolevels=nolevels,relscore = relscore,op=op)
  })
  #just get the odor pair names
  res.names <- sapply(pairs, function(x){
    paste(odors[x[1]],odors[x[2]],sep = ',')
  })
  names(res.ls) <- res.names
  #cat('\nAllanalog',names(freq.lst),'\n',names(res.ls))
  res.ls
}

#given a pair of odors, will get the overlap for cells within each class of 
#significance
#odor1,odor2: the odor pairs to consider in the 
#frequency list lists
#response rate lists
#notrials:the total no of trials
#nolevels: the number of reliability levels = notrials*nolevels
#relscore: 1, average of both odors, 2 - reliability of first odor
#op=1, overlap; 2, difference, 3 - probability difference, 
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#6 - given a cell, tells you if it is active for both odors of the odor pair
#7 - gives you the overlap for all cells that have a positive overlap score
#8 - gives you the probabilities difference, i.e., weighted by prob(odor1)
#9 - probability overlap, all cells,
#10 - 1 * prob of resp for odor2
#11 - resp. is freq_i*resp_i, and overlap is r_i*r_j/(r_i + r_j)^2
#12 - probability overlap, all cells, report overlap for 0 reliability or silent cells too
#13 - dot product of the two vectors
#14 - cosine distance
#15 - fractional hamming distance distance
#16 - normed hamming distance
#17 - euclidean distance
#18 - pearson's correlation
#cells: selectively do overlap only for the specified cells. 1 - default do all cells, if cells = list(odor1,odor2) with odor=c(cell1,cell2,...). 
#Specifies the cell numbers or whatever index is used to differentiatte cells
campOverlapAllPairAnalog <- function(odor1,odor2,freq,resp,notrials,nolevels,cells=1,
                                     relscore=1,op=1){
  #prepare the vectors
  tmp <- campPrepVecs(odor1 = odor1,odor2 = odor2,freq = freq,resp = resp,cells = cells,op = op)
  assignVarVals(c('vec1','vec2','params'),tmp)
  #calculate the reliability score
  if(relscore==1) rel <- sapply(1:length(vec1), function(x){
    floor(((vec1[x]+vec2[x])/2)*10*nolevels)/(10*nolevels)
  })
  else rel <- sapply(1:length(vec1), function(x) vec1[x])
  #cat('\nvectors',vec1,'2:',vec2,'\n',odor1,odor2,str(freq),'op',op,'\n')
  #make sure you pass only those notrials for the two pairs of odors
  over.res <- campOverlapFunction(vec1,vec2,notrials[c(odor1,odor2)],cells=cells,params=params,op=op)
  #cat('\ncampOverlapfn ',str(over.res))
  #process the results of campOverlap function: do any population level processing you need to do
  res <- campProcessResults(over.res = over.res,rel = rel,odor1 = odor1,odor2 = odor2,freq = freq,
                            resp = resp,vec1=vec1,vec2=vec2,cells = cells,op = op)
  #cat('\nres',str(res))
  res
}

#prepares the vectors for the campOverlapFunction function. For argument definitions, refer to campOverlapAllPairAnalog
#returns a list of the vector pair
campPrepVecs <-function(odor1,odor2,freq,resp,cells=1,op=1){
  params <- 0 #extra params that might have to be passed.
  vec1 <- freq[[odor1]]
  vec2 <- freq[[odor2]]
  #further processing of vecs needed
  if(op==11 || op>=13){#the vector is avg. response rate
    # vec1 <- vec1*resp[[odor1]]
    # vec2 <- vec2*resp[[odor2]]
    vec1 <- resp[[odor1]]
    vec2 <- resp[[odor2]]
    params <- list(mean(resp[[odor1]][which(resp[[odor1]]>0)]),mean(resp[[odor2]][which(resp[[odor2]]>0)])) #the mean of responses to both odors, for op=15 e.g,
    if(op==18) params <- c(params,list(sd(resp[[odor1]][which(resp[[odor1]]>0)]),sd(resp[[odor2]][which(resp[[odor2]]>0)]) ))
    #cat('\n',unlist(params),resp[[odor1]] )
  }
  #filtering cells
  if(isDataType(cells)==4){#if the cells list exists, get only the selected cells
    cellpair <- list(cells[[odor1]],cells[[odor2]])
    if(op>=12) allcells <- union(cellpair[[1]],cellpair[[2]]) #we want cells from both odors
    else allcells <- cells[[odor1]] #we only want cells from the first odor
    vec1 <- vec1[allcells]
    vec2 <- vec2[allcells]
    #update params too, mean and sd should only be for selected cells
    # cat('\ncOAPA',odor1,odor2,length(allcells),':',unlist(params),c(mean(vec1),mean(vec2),sd(vec1),sd(vec2)),':',cor(vec1,vec2),
    #     '\n1:',sort(vec1,decreasing = T),'\n2:',sort(vec2,decreasing = T),'\n',sortStringNos(allcells))
    params <- list(mean(vec1),mean(vec2),sd(vec1),sd(vec2))
  }
  #cat('\ncprep',vec1,allcells,resp[[odor1]][allcells])
  list(vec1,vec2,params)
}

#processes the results of the campOverlapFunction function. For argument definitions, refer to campOverlapAllPairAnalog
#returns the processed result
campProcessResults <- function(over.res,rel,odor1,odor2,freq,resp,vec1,vec2,cells=1,op=1){
  over <- over.res
  #cat('\n',str(over))
  if(op==12 || op==13 || op==15 || op==16){#this a population level activity, so sum it
    over <- sum(over.res)
    #cat('\n',odor1,odor2,over)
  }
  if(op==14){#population cosine
    #cat('\n',sum(over),sqrt(sum(vec1*vec1)),sum(vec2*vec2))
    over <- sum(over.res)/(sqrt(sum(vec1*vec1))*sqrt(sum(vec2*vec2)) )
    #cat(':over',over)
  }
  if(op==17){#population euclidean
    over <- sqrt(sum(over.res))
  }
  if(op==18){#pearsons correlation
    over <- sum(over.res)/(length(over.res)-1) #ivide by number of entries - 1
    #cat('\n18',str(over),str(over.res) )
  }
  #cat('\nsum',over)
  if(op >= 12) {
    res.df <- cbind(mean(rel),over) #population level computation, only one number
    #cat('\nhere',str(res.df),'\nover',str(over))
  }
  else res.df <- cbind(rel,over) #make the results into a matrix
  #cat('\ncampprocess',str(getThreshDF(res.df,op=1)))
  #cat('\nCampProcess op',op,':',over,'\nres.df',str(res.df))
  #and sort the matrix by relliability, and sum all the overlaps with the same reliability
  switch(op,sumDF(getThreshDF(res.df,op=1)),sumDF(getThreshDF(res.df,op=1)),
         sumDF(getThreshDF(res.df,op=1)),getThreshDF(res.df,op=1),getThreshDF(res.df,op=1),getThreshDF(res.df,op=1),
         res.df[which(res.df[,2]>0),],res.df[which(res.df[,2]!=0),],getThreshDF(res.df,op=1),getThreshDF(res.df,op=1), #7,8,9,10
         getThreshDF(res.df,op=1),res.df,res.df,res.df,res.df,res.df,res.df,#11,12,13,14,15,16,17 #old 12: getThreshDF(res.df,thresh = -1,op=1)
         res.df) #18
}


#this applies the overlap function operation that needs to be performed
#returns a vector of results
#vec1 and vec2 are the cells frequency response to the two odors
#notrials: the no of trials for the two pairs of odors
#ec50: the saturation coefficient , default value .02-.04
#op=1, overlap; 2, difference, 3 - probability difference, 
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#6 - given a cell, tells you if it is active for both odors of the odor pair
#7 - gives you the overlap for all cells that have a positive overlap score
#8 - gives you the probabilities difference, i.e., weighted by prob(odor1)
#9 - overlap all cells
#10 - 1 * prob of resp for odor2
#11 - resp. is freq_i*resp_i, and overlap is r_i*r_j/(r_i + r_j)^2
#12 - probability overlap, all cells, report overlap for 0 reliability or silent cells too
#13 - dot product of the two vectors
#14 - cosine distance
#15 - fractional hamming distance distance
#16 - normed hamming distance
#17 - euclidean distance
#18 - pearson's correlation
#cells: selectively do overlap only for the specified cells. 1 - default do all cells, if cells = list(odor1,odor2) with odor=c(cell1,cell2,...). 
#Specifies the cell numbers or whatever index is used to differentiatte cells
#params: the extra parameters that might be needed for each of the option that might require population level moments like the mean for option 15 
campOverlapFunction <- function(vec1,vec2,notrials,ec50=0.02,cells=1,params=params,op=1){
  #calculate the overlap or difference
  #bug fix: fixed the notrials per odor being different situation
  #cat('\noverlap',vec1[1],vec2[1],(notrials[1]*notrials[2]))
  #cat('\nparams',unlist(params))
  over <- switch(op,sapply(1:length(vec1), function(x){#1. overlap
    vec1[x]*vec2[x]/(notrials[1]*notrials[2])
  }),
  sapply(1:length(vec1), function(x){#2. difference
    abs((vec1[x]/notrials[1])-(vec2[x]/notrials[2])) #asymmetric notrials
  }),
  sapply(1:length(vec1), function(x){#3. prob. difference is not symmetric, more A, summed
    (vec1[x]/notrials[1])*((vec1[x]/notrials[1])-(vec2[x]/notrials[2])) #asymmetric notrials
  }),
  sapply(1:length(vec1), function(x){#4. prob. difference, not symmpetric, weights are weighted
    (vec1[x]/notrials[1])*((vec1[x]/notrials[1])-(vec2[x]/notrials[2])) #asymmetric notrials
  }),
  sapply(1:length(vec1), function(x){#5. prob. difference, not symmpetric, weights are saturated/weighted
    (vec1[x]/notrials[1])*(satFn(vec1[x]/notrials[1],ec50 = ec50)-satFn(vec2[x]/notrials[2],ec50 = ec50) )
  }),
  sapply(1:length(vec1), function(x){#6. whether the cell is active for both odors
    as.numeric(vec1[x]*vec2[x] > 0)
  }),
  sapply(1:length(vec1), function(x){#7. overlap for all cells
    vec1[x]*vec2[x]/(notrials[1]*notrials[2])
  }),
  sapply(1:length(vec1), function(x){#8. prob. difference is not symmetric, more A, summed
    (vec1[x]/notrials[1])*((vec1[x]/notrials[1])-(vec2[x]/notrials[2])) #asymmetric notrials
  }),
  sapply(1:length(vec1), function(x){#9. prob. overlap
    vec1[x]*vec2[x]/(notrials[1]*notrials[2])
  }),
  sapply(1:length(vec1), function(x){#10. prob. of odor 2
    vec2[x]/(notrials[2])
  }),
  sapply(1:length(vec1), function(x){#11. overlap=r_i*r_j/(r_i + r_j)^2
    vec1[x]*vec2[x]/((vec1[x]+vec2[x])/2)^2
    #vec1[x]*vec2[x]/(notrials[1]*notrials[2])
  }),  
  sapply(1:length(vec1), function(x){#12. prob. overlap, and reporting for silent cells too, op=12
    vec1[x]*vec2[x]/(notrials[1]*notrials[2])
  }),
  sapply(1:length(vec1), function(x){#13. dot product of the two population vectors
    vec1[x]*vec2[x]
  }),
  sapply(1:length(vec1), function(x){#14. cosine distance of the two population vectors
    vec1[x]*vec2[x]
  }),
  sapply(1:length(vec1), function(x){#15. fractional hamming distance
    no <- countAboveThresh(c(vec1[x],vec2[x])) #no of entries > 0
    #cat('\t',(vec1[x]-params[[1]]),(vec2[x]-params[[2]]),(no * 2 -3) * (((vec1[x]-params[[1]])+(vec2[x]-params[[2]]))/no),'::')
    res <- ((no * 2) -3) * abs(((vec1[x]-params[[1]])+(vec2[x]-params[[2]]))/no) #non-overlapping subtract, and overlapping add to score
    cleanNAVec(res)
  }),
  sapply(1:length(vec1), function(x){#16. normed hamming distance
    abs(vec1[x] - vec2[x]) #hamming but without the binary bit
  }),
  sapply(1:length(vec1), function(x){#17. euclidean distance
    (vec1[x] - vec2[x])^2 #euclidean distance
  }),  
  sapply(1:length(vec1), function(x){#18. pearson's correlation
    #cat('\n',((vec1[x]-params[[1]])/params[[3]]),((vec2[x]-params[[2]])/params[[4]]) )
    ((vec1[x]-params[[1]])/params[[3]]) * ((vec2[x]-params[[2]])/params[[4]]) 
  }))
  over
}



#function that prepares the output of campOverlapAllAnalog in a form ready for 
#campMatchHallemRob. The idea is that each entry will have nolevels*notrials number of 
#ntries with blank ones set to 0
#might not be needed, why plot the 0s, just do it
campPrepareAnalog <- function(dat.ls,nolevels=2,op=1){
  notrials <- max(campGetTrialDetails()[[1]])#kludgy fix for now, have to figure this out
  #cat(notrials,nolevels)
  levels <- c(1:(nolevels*notrials))/nolevels #determine ideal no of levels
  #cat(levels,str(dat.ls))
  res.ls <- lapply(dat.ls, function(x){
    odorlevels <- x[,1] 
    newlevels <- setdiff(levels,odorlevels) #get missing levels
    #now insert these levels with overlap of 0
    new.df <- cbind.data.frame(newlevels,rep(0,length(newlevels)))
    #cat(str(x),str(new.df),str(joinDF(x,new.df)))
    res <- getMatSortPosns(joinDF(x,new.df),col = 1,op=3)
    #cat(newlevels,'\t',odorlevels)
    unlist(res[,2]) #just return the overlap in proper sorted order now
  })
  res.ls
}

#gets the number of trials and the frequency of each odor in this data set, as well as the number of cells
#data: get details of the trails and the odornames
#op:1= notrials and trial names
#op:2= notrials and trial names without using foldname
campGetTrialDetails <- function(foldname='.',alldata=c(),op=1){
  #cat('\ncGET',op)
  if(op==2){#without using foldname
    #cat('\ngetrial',str(alldata))
    odornames <- unique(names(alldata))
    odorfreq.tab <- table(names(alldata))
    odorfreq <- as.vector(odorfreq.tab)
    names(odorfreq) <- names(odorfreq.tab)
    nocells <- length(alldata[[1]][[1]])
    #cat('\nhere',odorfreq,odornames)
  }
  #print(datahead)
  if(op==1){#with using foldname
    datahead <- campGetDataDetails(foldname=foldname) #get header data and data frames with valid data
    if( length(searchString(names(datahead[[1]])[1],'mouse')) > 0 ){#mnouse, so different type of processing
      details <- simReadTrialDetailsFile(foldname = foldname)
      #odors <- simReadTrialOdorMap(foldname = foldname) #get all the odors
      odorfreq <- details[,2]
      odornames <- details[,1]
      names(odorfreq) <- odornames
      nocells <- details[1,3]
    }
    else {#the fly  
      #read in a datafile and count the number of rows to get no of cells
      data.df <- read.csv(paste(foldname,names(datahead[[1]])[1],sep = '/'),header = F)
      nocells <- nrow(data.df)
      tmp <- datahead[[2]][datahead[[1]],] #only choose those rows.trials that have valid data
      odors <- split(tmp[,1],tmp[,1]) #get the odors and list them based on odor names
      odornames <- names(odors)
      odorfreq <- sapply(odors, function(x) length(x)) #no of trials for each odor
      trial.freq <- tabulate(odorfreq)
      notrials <- which(trial.freq == max(trial.freq)) 
    }
  }
  # the number of trials is frequency of trials that are the highest, e.g. most odors have 4 trials
  #cat('trial freq',trial.freq,odorfreq,max(trial.freq),'\n')
  #cat(notrials)
  list(odorfreq,odornames,nocells) #freq of each odor, names of odors
}

#this function gets details of the header and which data frames actually havve valid data in them
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for the sigpvals file
#alpha is used for calculating the SDs above mean where the significance should be fixed
#meanmax and sigval and ranking too are not needed here
#op: types of return functions
#change: changed tseries to just 'series' to align it with simon's files
campGetDataDetails <- function(foldname='.',fpatterns=c('*series.*head.*csv','*series.*sigpvals.*csv','*series.*data[0-9]+.*csv'),alpha=0.01,
                               op=1){
  #just look in this directory and read the files which fit the description of the 
  #for now match those files that are headers: contain the word head and are csv files
  headerfile <- list.files(path = foldname,pattern = fpatterns[1])
  sigpfile <- list.files(path = foldname,pattern = fpatterns[2])
  datafiles <- list.files(path = foldname,pattern = fpatterns[3])
  #cat('here',datafiles,list.files())
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile)
  if(length(sigpfile) == 1) sigpfile <- read.csv(file = sigpfile,header = T)
  #cat(str(sigpfile),str(headerfile))
  #cat('\ndata',datafiles)
  #cat('\n',datafiles,'\n',sortStringNos(datafiles))
  datafiles <- sortStringNos(datafiles)
  #in some files, not all tseries_data have data for all time frames. So, pick the ones that have the 
  #highest number of timeframes and only analyze those
  #cat(str(datafiles),'done\n')
  datalen <- checkCampMinLength(sigpfile,datafiles) #datalen is the size of the timeframe
  #and get the datarows that have data for all timeframes
  datarows <- sapply(datafiles, function(x) nrow(data.table::fread(x, select = 1L)))
  if( checkVectorSimilar(datarows,op=3) < getLength(headerfile) ) data.entries <- checkVectorSimilar(datarows,op=2)
  else data.entries <- 1:getLength(headerfile)
  list(data.entries,headerfile) #gives the timeframes or trials that are valid
}


#given the hallem dataset will get the the vectors of interest
#hallem: the hallem data set, the rows contain the vectors of interest
#odors: the names of the odors
#op=1, on the rob odors, 2 - all the odors
campGetHallemOdors <- function(hallem,odors,op=1){
  #go through find the odors of interest and then get the corresponding vectors
  res <- lapply(odors, function(x){
    index <- which(hallem[,1]==x)
    vec <- cleanNA(as.numeric(hallem[index,2:ncol(hallem)]))
    vec
  })
  #the odors that are not present hacve to be taken out, get their indices
  skipindex <- sapply(res, function(x){
    if(length(x) == 0) F
    else T
  })
  res <- res[skipindex] #and take the indicated odors out
  names(res) <- odors[skipindex]
  #cat('campGetHallemOdors',odors,'::',names(res),'\n')
  res
}

#gets the correlations between the odors and return a matrix
#op=1, on the rob odors, 2 - all the odors
#topn=the topn percentile, by default=100
campGetHallemCorr <- function(hallem,odors,topn=100,op=1){
  odorvecs <- campGetHallemOdors(hallem = hallem,odors = odors,op=op)
  #cat('campgethallemcorr',odors,str(odorvecs),length(odorvecs),1:0)
  if(length(odorvecs)==0) return(NULL) #no mtch so return nothing
  res.ls <- lapply(1:length(odorvecs),function(x){
    res <- sapply(1:length(odorvecs),function(y){
      #cor(odorvecs[[x]],odorvecs[[y]]) #old
      #cat('\ncampgetHallemCorr',x,y,'\n',odorvecs[[x]],';',odorvecs[[y]],';',topn)
      compareTopNPercentiles(odorvecs[[x]],odorvecs[[y]],topn = topn,op = 1)
    })
  })
  res.df <- convertNestedListsDF(res.ls)
  rownames(res.df) <- names(odorvecs)
  colnames(res.df) <- names(odorvecs)
  res.df
}

#function that comparse and matches hallem and rob data sets
#thresh: only consider those odors with a correlation higher than this
#match: only pick those odors that have a correlation within match, default gets all odors
#op=1, do hallem, 2 - rob, 3 - do both using match, 4 - get it for odorpairs specified in match
campCompareHallRobCor <- function(hallem,data.lst,skip=c(),topn=100,thresh=0,match=1,op=1){
  robcor <- computeListMatCorr(data.lst,matop = 2,params = c(topn,topn),op=1,celltype = 5)
  hallcor <- campGetHallemCorr(hallem = hallem,odors = unique(names(data.lst)),topn = topn)
  #screen for match and skip odors in skip
  robcor.vec <- unlist(flattenDFtoList(robcor))
  hallcor.vec <- unlist(flattenDFtoList(hallcor))
  common <- intersect(names(robcor.vec),names(hallcor.vec))
  res <- cbind.data.frame(hallcor.vec[common],robcor.vec[common])
  rownames(res) <- names(robcor.vec[common])
  #now go through the options op
  if(op==1){#hallem set
    #cat('\n',which(res[,1]>thresh),res[,1],'\n',res[which(res[,1]>thresh),1],'\n')
    res <- res[which(res[,1]>thresh),]
  }
  if(op==2){#rob set
    res <- res[which(res[,2]>thresh),]
  }
  if(op==3){#hallem and rob set
    res <- res[which(res[,1]>thresh & res[,2]>thresh),]
  }
  if(op==4){#get odorpairs specified by match
    res <- res[match,]  
  }
  match.vec <- abs(res[,2]-res[,1])
  res <- res[which(match.vec < match),]
  res <- res[getIndexSameSubStr(rownames(res),op=2),] #filter out the same odor-pair
}

#same as above except this does it for a range of topn values
campCompareHallRobCorTopn <- function(hallem,data.lst,skip=c(),topn=c(50,100),thresh=0,match=1,op=1){
  #first get the list for all topn and and then apply threshold to see which ones you want
  orig <- campCompareHallRobCor(hallem = hallem,data.lst = data.lst,skip = skip,topn = 100,thresh = thresh,match = 2,op=op)
  odors <- rownames(orig)
  #now for these odor pairs get the percentile scores specified for the topn in topn
  res.ls <- lapply(topn,function(x){
    hallrob <- campCompareHallRobCor(hallem = hallem,data.lst = data.lst,skip = skip,topn = x,
                                     thresh = -1,match = odors,op=4)
    
  })
  names(res.ls) <- topn
  res.ls
}


#function that compares hallem and robs data for a particular reliability level
#goto fix this function, what are the inputs? and forget about nolevels
#op=1,analog, 2 - digital- splot into reliable and
campMatchHallemRob <- function(roboverlap,hallemcorr,notrials=6,nolevels=2,op=1){
  #make the hallem corr into a flat list
  hall.ls <- flattenDFtoList(hallemcorr)
  totlevels <- nolevels*notrials
  #go in and get the rob overlap for the different 
  res.ls <- lapply(names(hall.ls), function(x){
    #cat(str(roboverlap[x][[1]]))
    robover <- sapply(1:totlevels, function(y) roboverlap[x][[1]][y])
    if(op==2) robover <- c(sum(robover[1:(totlevels/2)]),
                           sum(robover[(1+(totlevels/2)):totlevels]))
    c(hall.ls[[x]],robover)
  })
  #cat(str(res.ls))
  res.df <- transposeDF(convertNestedListsDF(res.ls))
  rownames(res.df) <- names(hall.ls)
  res.df
}

#matches the Hallem and Rob correlations and gets the useful difference for those odor pairs that are within 0.1 correlations
#according to both Hallem and Rob
#hallem: the hallem data set
#rob: correlations between Rob's odors
#odors: the odor names for doing correlations
#optype: 1 - overlap; 2, difference, 3 - probability difference, 
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#9 - overlap, all cells
#hallrobmatch: how to balance the two sets: 1 - only odors present in both, no similar odor pairs 
#2 - only odors present in both, with similar odor pairs,
#3 -rob odors, similar odor pairs ok, 4 - rob odors, similar odor pairs ok
#also in all cases i
#celltype: type of cells for which we should get correlations, 1 - all significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells
#op=1, the sum of reliability cells of the same type, 2 - mean of all reliability cells of the same type, 
#todo: stuff to weed out empty and mismatched and same odor pairs
campMatchHalRobCor <- function(hallem,skip=c('paraff','empty'),matchlevel = 0.1,halrobmatch=1,celltype=4,optype=4,op=1){
  #getall the data, and then rob's odor correlations
  alldata <- getCampDataSigList() #gets the original data
  notrials <- campGetTrialDetails()[[1]] #get no of trials
  #cat('\ncMHRC',length(notrials))
  if (length(notrials) < 2) return(NULL) #less than 2 odors, nothing to compare
  #cat('\ncMHRC',length(notrials[[1]]))
  #get the Rob correlations
  rob <- computeListMatCorr(alldata,celltype =celltype)
  #get the hallem correlations for the odors from rob
  hallemdf <- campGetHallemCorr(hallem,rownames(rob))
  #calcaulte tjje overlap/difference for each reliability level
  overlap.cols <- campOverlapAnalogRel(alldata = alldata,notrials=notrials,skip = skip,relscore = 2,optype = optype,op=op)
  if(halrobmatch<=2){#odors have to be present in both
    odors <- intersect(rownames(rob),rownames(hallemdf))
    res <- pickMatchedVals(hallemdf[odors,odors],rob[odors,odors],prec = matchlevel)
    if (length(res)==0) return(NULL) #no match, return NULL
    res <- t(res) #transpose so that the odors are in rows
    colnames(res) <- c('match','hallem','rob')
  }
  else{#only doing the rob thing
    res <- cbind.data.frame(unlist(flattenDFtoList(rob))) #make it into a DF
    rownames(res) <- names(flattenDFtoList(rob))
    colnames(res) <- c('rob')
  }
  odors <- getSkipStrIndex(rownames(res),skip=skip,op=2) #odors without skip
  #now, take out similar odor pairs
  if(halrobmatch==1 || halrobmatch==3){#take out similar odor pairs
    odors <- odors[getIndexSameSubStr(odors,op=2)] #op=2 gets strings with dissimilar substrings
  }  
  #now, we have to match Rob and Hallem according to the conditions
  #cat('\ncMHRC ',setdiff(odors,rownames(res)),'\n',odors)
  res <- cbind.data.frame(res[odors,],overlap.cols[odors,])
  names(res)[1] <- 'rob' #just have to reset the robname for cases of halrobmatch 3 and 4
  res
}

#given the campsigdatalist structure will get alll those odor-pairs with correlation over a particular threshold
# and get the overlap vs reliability for all the cells of these odor-pairs
#thresh.win - lower and upper bounds of the correlations for specifying similar or dissimilar odors
#overop: the kind of overlap option you want: 9 - overlap of cells, 10 - 1 * freq of odor2, 11- measure
#op: 1 - list of overlaps, 2 - df of all cells with reliabiltiy in 1st col and overlap in second, 3: correlation between the two, 
# and equation ofthe regression
campGetOverlapCells <- function(data.lst,thresh.win=c(0,1),overop=9,op=1){
  #get the overlap for all cells of each odor pair  
  #overlap.allcells <- campComputeOverlap(datasig = data.lst,op=2)
  #cat('\ncampGetOverlapCells')
  overlap.allcells <- campOverlapAllAnalog(data = data.lst,skip = c(),trialop = 2,op = overop)
  #cat(str(overlap.allcells))
  #now, get the odor pairs whose correlations are within the thresh.win, and get the overlaps of the cells only for these pairs
  odor.corr <- computeListMatCorr(alldata = data.lst,celltype=4) #correlations of all odor pairs
  odor.corr <- flattenDFtoList(odor.corr)
  #cat('\nCGOC')
  odorpair.overlap <- joinListDFs(overlap.allcells[names(odor.corr)[which(odor.corr>=thresh.win[1] & odor.corr<=thresh.win[2])]])
  odor.over.lst <- lapply(unique(odorpair.overlap[,1]),function(i) odorpair.overlap[which(odorpair.overlap[,1]==i),2] )
  #assign the reliability labels to names, and make sure they are sorted
  names(odor.over.lst) <- unique(odorpair.overlap[,1])
  odor.over.lst <- odor.over.lst[getSortPosns(names(odor.over.lst),op=2)]
  #cat('\nCGOC',unique(odorpair.overlap[,1]))
  #calculate the correlations and line fit
  overrel.corr <- cor(odorpair.overlap[,1],odorpair.overlap[,2])
  overrel.fit <- fitlinearlaw(odorpair.overlap[,1],odorpair.overlap[,2])
  switch(op,odor.over.lst,odorpair.overlap,list(overrel.corr,overrel.fit))
}

#given the campsigdatalist structure will get the overlap for all cells that are in the top N percentile
# and get the overlap vs reliability for all the cells of these odor-pairs
#thresh.win - lower and upper bounds of the correlations for specifying similar or dissimilar odors
#op: 1 - list of overlaps, 2 - df of all cells with reliabiltiy in 1st col and overlap in second, 3 - vector of overlaps  
#topn: the top Nth percentile. 0 means all cells and 5 would mean top 95 percentile of cells
#overop: the kind of overlap option you want: 9 - overlap of cells, 10 - 1 * freq of odor2, 11- measure, 12 - prob. overlap but report 
#for silent cells, too.
#symop: 1 - asymmetric for the first of the odor pair, 2 - symmetric as in both of the odor pairs are considered for topn cells
#sameop: same odors, 1 - yes, 2 - no
#deciles - c(), not using deciles just topN percentile, number. e.g., 10 signifies the decile from topN to topN+decile
campGetOverlapTopN <- function(data.lst,thresh.win=c(0,1),topn=0,deciles=c(),symop=1,overop=9,sameop=1,op=1){
  sig.topn <- campGetTopNCells(data.lst = data.lst,topn = topn,deciles = deciles,op=2) #get the topN neurons for all odors
  #get the overlap for all cells of each odor pair
  #overlap.allcells <- campComputeOverlap(datasig = data.lst,overop = overop,op=2)
  overlap.allcells <- campOverlapAllAnalog(data = data.lst,skip = c(),op = overop,trialop = 2,cells = sig.topn)
  #cat('\n',str(overlap.allcells),'\n',str(sig.topn))
  #Algo:what we have to do is for each odor pair, designate the first of the pair as the primary odor, then
  #get the topn neurons for that odor, then go back to the overlap list entry and only get those cells and their overlap.
  #now, get the odor pairs whose correlations are within the thresh.win, and get the overlaps of the cells only for these pairs
  odor.corr <- computeListMatCorr(alldata = data.lst,celltype = 4) #correlations of all odor pairs
  odor.corr <- flattenDFtoList(odor.corr)
  odorpairs <- names(odor.corr)[which(odor.corr>=thresh.win[1] & odor.corr<=thresh.win[2])] #names of odor pairs that satisfy the criterion
  overlap.target <- overlap.allcells[odorpairs] #the overlap numbers for the targetted population
  # res <- lapply(names(overlap.target),function(x){
  #   #choose the appropriate top N cells from one or both odors, based on symop
  #   odors <- strsplit(x=x,split = ',')[[1]]
  #   # if(symop==1) odornames <- names(sig.topn[odors[1]][[1]])
  #   # else odornames <- union(names(sig.topn[odors[1]][[1]]),names(sig.topn[odors[2]][[1]]))
  #   overlap.target[[x]][odornames,] #get the particular cells specific in sig.topn by using their names.
  # })
  res <- overlap.target
  names(res) <- names(overlap.target)
  if(sameop==2){#take out similar odor pairs
    res <- res[getIndexSameSubStr(names(res),op=2)] #op=2 gets strings with dissimilar substrings
  } 
  #cat('\nCGOY',str(res))
  res.df <- joinListDFs(res)
  #list(names(sig.topn),names(overlap.allcells),overlap.allcells[odorpairs],res,res.df )
  res <- switch(op,res,res.df,res.df[,2])
}

#thresh.win - lower and upper bounds of the correlations for specifying similar or dissimilar odors
#op: 1 - list of overlaps, 2 - df of all cells with reliabiltiy in 1st col and overlap in second, 3 - vector of overlaps 
#topn: the top Nth percentile. 0 means all cells and 5 would mean top 95 percentile of cells
#symop: 1 - asymmetric for the first of the odor pair, 2 - symmetric as in both of the odor pairs are considered for topn cells
#overop: the kind of overlap option you want: 9 - overlap of cells, 10 - 1 * freq of odor2, 11- measure, 12 - prob. overlap but report 
#for silent cells, too.
#sameop: same odors, 1 - yes, 2 - no
campGetOverlapTopNSeq <- function(data.lst,thresh.win=c(0,1),topnseq=seq(10,50,10),deciles=c(),overop=9,sameop=1,symop=1,op=1){
  overlap.res <- lapply(topnseq,function(i) campGetOverlapTopN(data.lst = data.lst,thresh.win = thresh.win,topn = i,
                                                               overop=overop,symop=symop,sameop = sameop,deciles = deciles,op=op) )
  #res <- lapply(overlap.res,function(x) x[,2])
  names(overlap.res) <- topnseq
  res <- rev(overlap.res)
  res
}



#given the campsigdatalist structure will get the overlap for all cells that are in the top N percentile
# and get the overlap vs reliability for all the cells of these odor-pairs at the population level
#thresh.win - lower and upper bounds of the correlations for specifying similar or dissimilar odors
#op: 1 - list of overlaps, 2 - df of all cells with reliabiltiy in 1st col and overlap in second, 3 - vector of overlaps  
#topn: the top Nth percentile. 0 means all cells and 5 would mean top 95 percentile of cells
#overop: the kind of overlap option you want: 9 - overlap of cells, 10 - 1 * freq of odor2, 11- measure, 12 - prob. overlap but report 
#for silent cells, too, 13 - dot product too, 14 - cosine distance
#symop: 1 - asymmetric for the first of the odor pair, 2 - symmetric as in both of the odor pairs are considered for topn cells
#sameop: same odors, 1 - yes, 2 - no
#deciles - c(), not using deciles just topN percentile, number. e.g., 10 signifies the decile from topN to topN+decile
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campGetOverlapTopNPop <- function(data.lst,thresh.win=c(0,1),topn=0,deciles=c(),symop=1,overop=9,sameop=1,respop=1,op=1){
  sig.topn <- campGetTopNCells(data.lst = data.lst,topn = topn,deciles = deciles,respop = respop,op=2) #get the topN neurons for all odors
  #get the overlap for all cells of each odor pair
  #overlap.allcells <- campOverlapAllAnalog(data = data.lst,skip = c(),op = overop,trialop = 2,cells = sig.topn)
  #Algo:what we have to do is for each odor pair, designate the first of the pair as the primary odor, then
  #get the topn neurons for that odor, then go back to the overlap list entry and only get those cells and their overlap.
  #now, get the odor pairs whose correlations are within the thresh.win, and get the overlaps of the cells only for these pairs
  odor.corr <- computeListMatCorr(alldata = data.lst,celltype=4,matop=2) #correlations of all odor pairs
  odor.corr <- flattenDFtoList(odor.corr)
  odorpairs <- names(odor.corr)[which(odor.corr>=thresh.win[1] & odor.corr<=thresh.win[2])] #names of odor pairs that satisfy the criterion
  #cat('\ncor',str(odor.corr),'\n',odorpairs)
  overlap.cells <- campOverlapAllAnalog(data = data.lst,skip = c(),op = overop,trialop = 2,cells = sig.topn,odorpairs = odorpairs)
  #overlap.target <- overlap.allcells[odorpairs] #the overlap numbers for the targetted population
  #res <- overlap.target
  res <- overlap.cells
  if(sameop==2){#take out similar odor pairs
    res <- res[getIndexSameSubStr(names(res),op=2)] #op=2 gets strings with dissimilar substrings
  }
  res.df <- joinListDFs(res)
  vec <- res.df[,2]
  names(vec) <- rownames(res.df)
  #cat('\ntop',str(res.df))
  res <- switch(op,res,res.df,vec)
}

#do it for a seq of topn values, but here the data points are sums for the whole pair
#thresh.win - lower and upper bounds of the correlations for specifying similar or dissimilar odors
#op: 1 - list of overlaps, 2 - df of all cells with reliabiltiy in 1st col and overlap in second, 3 - vector of overlaps 
#topn: the top Nth percentile. 0 means all cells and 5 would mean top 95 percentile of cells
#symop: 1 - asymmetric for the first of the odor pair, 2 - symmetric as in both of the odor pairs are considered for topn cells
#overop: the kind of overlap option you want: 9 - overlap of cells, 10 - 1 * freq of odor2, 11- measure, 12 - prob. overlap but report 
#for silent cells, too.
#sameop: same odors, 1 - yes, 2 - no
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campGetOverlapTopNPopSeq <- function(data.lst,thresh.win=c(0,1),topnseq=seq(10,50,10),deciles=c(),overop=9,sameop=1,symop=1,respop=1,op=1){
  overlap.res <- lapply(topnseq,function(i) campGetOverlapTopNPop(data.lst = data.lst,thresh.win = thresh.win,topn = i,
                                                                  overop=overop,symop=symop,sameop = sameop,deciles = deciles,respop = respop,op=op) )
  #res <- lapply(overlap.res,function(x) x[,2])
  names(overlap.res) <- topnseq
  res <- rev(overlap.res)
  res
}


campCorCosTest <- function(no=1,sds=1,op=1){
  res1 <- c(rep(10,24),11)
  res2 <- c(rep(9,24),12)
  res3 <- lapply(1:10,function(i) c(res1,rnorm(length(res1)*no,10,sds)))
  res4 <- lapply(1:10,function(i) c(res2,rnorm(length(res2)*no,10,sds)))
  res.cor <- sapply(1:10,function(i) cor(res3[[i]],res4[[i]]))
  res.cos <- sapply(1:10,function(i) cosine(res3[[i]],res4[[i]]))
  c(mean(res.cor),mean(res.cos))
 
  res3 <- lapply(1:10,function(i) c(res1,rep(7.5,12),rnorm(13,7.5,4),
                                    rep(5,6),rnorm(19,5,3),
                                    rep(2.5,3),rnorm(22,2.5,1))  )
  res4 <- lapply(1:10,function(i) c(res1,rep(7.5,12),rnorm(13,7.5,4),
                                    rep(5,6),rnorm(19,5,3),
                                    rep(2.5,3),rnorm(22,2.5,2))  )
  
  res.cor <- sapply(1:10,function(i) cor(res3[[i]],res4[[i]]))
  res.cos <- sapply(1:10,function(i) cosine(res3[[i]],res4[[i]]))
  
  c(mean(res.cor),mean(res.cos))
  
  
  #c(cor(res3,res4),cosine(res3,res4))
}

# result: the email explanation, plus it is possible that when you have the vecror as a whole, if you have smaller numbers in addition to bigger ones
# that tends to drive correlation


#prepares the data needed for the similarity vs overlap with different lines for each topk plot
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campOverSimRankData <- function(alldata,decile=25,overop=14,respop=1,op=1){
  #get the correlations first
  res.cor <- computeListMatCorr(alldata,celltype = 1,matop = 2)
  res.cor <- sort(unlist(flattenDFtoList(res.cor)))
  res.cor <- res.cor[getIndexSameSubStr(names(res.cor),op=2)] #op=2 gets strings with dissimilar substrings
  #
  overlap <- campGetOverlapTopNPopSeq(data.lst = alldata,thresh.win = c(-1,1),topnseq = seq(decile,100,decile),
                                      sameop = 2,overop=overop,op=3,deciles = decile,respop = respop)
  
  #now go through the top k's and for each one, make a pairinh
  res.ls <- lapply(overlap, function(x){
    overvals <- sort(x) #sort the overlap vals 
    #cat('\n',str(x))
    res <- cbind.data.frame(res.cor,x[names(res.cor)])
  })
  res.ls
}


#gets the odor trial responses as a list, with only the signifant responses havving values, others are 0
#data.sig: the data
#op
campGetOdorTrials <- function(data.sig,op=1){
  len <- length(data.sig)
  res.lst <- lapply(data.sig,function(x){
    resp <- x[[1]]*x[[2]]
  })
  res.lst
}

#gets the topN percentile of cells that respond. Topn is on the basis of response rate
#topop: 1 - do it based on average activity, 2 - highest activity in any of the trials
#op=1, get their values and names are posns , 2 - get their positions
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campGetTopNCells <- function(data.lst,topn=5,deciles=c(),respop=1,topop=1,op=1){
  res.sig <- getGroupedListsDf(data.lst,op=1) #gets the significant response Dfs
  res.val <- getGroupedListsDf(data.lst,op=2) #gets the significant response Dfs
  notrials <- campGetTrialDetails(alldata = data.lst,op=2) #get the number for trials for each odor
  #now get the topn for every odor as a list
  res.avgres <- lapply(1:length(res.sig),function(i){
    #the average response of only the cells that respond significantly
    res.filter <- res.val[[i]] * res.sig[[i]]
    #calculate top based on average or max activity in significant trials
    res.avg <- switch(topop,switch(respop,apply(res.filter,1,sum)/apply(res.sig[[i]],1,sum),
                                   apply(res.filter,1,sum)/notrials[[1]][i] ),
                      apply(res.filter,1,max))
    #cat('\nGettopn',str(res.avg),str(res.filter))
    res.avg <- cleanNA(res.avg) #get rid of NAs
    if(length(deciles)>0) getTopNDecile(res.avg,topn = topn,decile = deciles,zero = 2,op=op+1)
    else getTopNPercentile(res.avg,topn = topn,op=op+1,zero = 2)
  })
  names(res.avgres) <- names(res.sig)
  res.avgres
}

#gets the percentile of cells that respond within percentile range perange
#op=1, get their values and names are posns , 2 - get their positions
# campGetPercNCells <- function(data.lst,perange=c(10,5),deciles=c(),op=1){
#   top1 <- campGetTopNCells(data.lst = data.lst,topn = perange[1],deciles = topn[2]-topn[1],op=op)
# }



#gets the reliability and response as DFs. where the reliability of a cell and mean response are in the columns
#data.lst: return of getCampDataSigList()
#op:1 - as dfs 2 - as lists
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campGetReliabResp <- function(data.lst,respop=1,op=1){
  res.sig <- getGroupedListsDf(data.lst,op=1) #gets the significant response Dfs
  res.val <- getGroupedListsDf(data.lst,op=2) #gets the response Dfs
  notrials <- campGetTrialDetails(alldata = data.lst,op=2) #get the number for trials for each odor
  #now get the topn for every odor as a list
  res.avgres <- lapply(1:length(res.sig),function(i){
    #the average response of only the cells that respond significantly
    res.filter <- res.val[[i]] * res.sig[[i]]
    res.avg <- switch(respop,cleanNAVec( apply(res.filter,1,sum)/apply(res.sig[[i]],1,sum)),
                      cleanNAVec( apply(res.filter,1,sum)/notrials[[1]][i]) )
  })
  res.rel <- lapply(res.sig,function(x) apply(x,1,sum))
  res.rel.df <- convertNestedListsDF(res.rel)
  res.avgres.df <- convertNestedListsDF(res.avgres)
  names(res.avgres) <- names(res.rel)
  names(res.avgres.df) <- names(res.sig)
  names(res.rel.df) <- names(res.sig)
  rownames(res.avgres.df) <- 1:nrow(res.avgres.df)
  rownames(res.rel.df) <- 1:nrow(res.rel.df)
  switch(op,list(res.rel.df,res.avgres.df),list(res.rel,res.avgres) )  
}


#function compares the top n % of cells and gets the number of cells that are in each reliability 
#class for each of the odors
#topn: the topn decile/percentile to consider
#trialop: get trial detsils from the data.lst data structure
#deciles: if we want deciles, specify the decile size and it goes from top to topn+decile
#op:1 - abolute number, 2 - percentage
#topop: the top cells are decided on the basis of average response (1) or max response (2) across trials
#retop: 1- a list with the number of cells of that reliability for each odor, 2 - a list with 
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campCompTopRelCells <- function(data.lst,topn=10,trialop=2,deciles=1,retop=1,topop=1,respop=1,op=1){
  #get the top cells
  if (length(deciles) > 0) topcells <- campGetTopNCells(data.lst = data.lst,topn=topn,deciles = deciles,topop = topop,respop = respop)
  else topcells <- campGetTopNCells(data.lst = data.lst,topn = topn,topop = topop,respop = respop) 
  #get cells according to their reliability
  relcells <- campGetReliabResp(data.lst = data.lst)
  odors <- names(topcells)
  notrials <-campGetTrialDetails(alldata = data.lst,op = trialop)[[1]]
  #for each odor, the reliability of the top responding cells
  odor.rel <- lapply(1:length(odors),function(i) {
    res <- relcells[[1]][odors[i]][names(topcells[[odors[i]]]),]
  })
  # for each reliability level, get their frequencies per odor
  reltop.lst <- lapply(1:max(notrials),function(i) {
    freq <- sapply(1:length(odor.rel),function(j) switch(op,countElem(vec=odor.rel[[j]],elem=i), #absolute freq
                                                         countElem(vec=odor.rel[[j]],elem=i)/length(odor.rel[[j]]) ) ) #relative freq
    names(freq) <- 1:length(odor.rel)
    freq
  })
  names(reltop.lst) <- 1:max(notrials)
  reltop.lst <- lapply(reltop.lst,function(x) getAboveThresh(x))
  reltop.lst
}


#comparese the reliability of population deciles across the whole range
#decsize: size of the deciles. Has to be a perfect divisor of 100 as it has to split all 100 %
#topop: the top cells are decided on the basis of average response (1) or max response (2) across trials
#op:1: reliability make up, just the return of campCompTopRelCells, 2 = absolute number of cells of each reliability class,
#3 - normalized number of each reliability class
campRelDeciles <- function(data.lst,decsize,topop=1,op=1){
  topn.size <- floor(decsize)
  topn.no <- floor(100/topn.size)
  topn.range <- sapply(topn.no:1,function(i) (i*topn.size) )
  if((100-topn.range[1]) < (topn.size/10)) topn.range <- c(100,topn.range[-1],0)
  else topn.range <- c(100,topn.range,0)
  #cat('\n',topn.range)
  reldeciles <- lapply(1:(length(topn.range)-1),function(i){
    #cat('\ntopn',topn.range[i],topn.range[i]-topn.range[i+1])
    res <- campCompTopRelCells(data.lst = data.lst,topn = topn.range[i],deciles = topn.range[i]-topn.range[i+1],topop = topop)
    #return reliability make up, absolute number of cells of each reliability class, normalized number of each reliability class
    switch(op,res,sapply(res, mean),100*sapply(res, mean)/sum(sapply(res, mean)))
  })
  names(reldeciles) <- topn.range[-length(topn.range)]
  
  reldeciles
}

#thigns to fix: make outputs as percentages, then fix the paper too, Fits can be power law or exponential, might be better than linear fits
#gets the statistics of reliable and unreliable cells: Mainly the ones in Fig. 1 siuch as 
#R and U cells per odor, per trial, and correlation and eqn fit between reliability and response rate
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
campGetReliabStats <- function(data.ls,skip=c(),respop=1,op=1){
  #get details of this set: noodor, notrials
  trial.det <- campGetTrialDetails(alldata = data.ls,op=2)
  noodors <- length(trial.det[[2]])
  notrials <- max(trial.det[[1]])
  nocells <- trial.det[[3]]
  odor.dets <- list(noodors,notrials,nocells)
  names(odor.dets) <- c('# odors','# trials','#cells')
  percent.factor <- nocells/100 #the factor for correcting numbers so that they are in percent terms
  #get the reliability and response rates
  res.resp <- campGetReliabResp(data.ls,respop = respop) 
  #cor.relresp <- sapply(1:ncol(res.resp[[1]]),function(i) cor(res.resp[[1]][,i],res.resp[[2]][,i]) )
  #now get the average response rate for each reliability level, here average is E(x) or average based on respop
  rel.vec <- as.vector(as.matrix(res.resp[[1]]))
  resp.vec <- as.vector(as.matrix(res.resp[[2]]))
  #cat('\n',str(rel.vec),str(resp.vec))
  resp.trial <- sapply(1:notrials,function(i) mean(resp.vec[which(rel.vec==i)])) 
  resp.trial.sem <- sapply(1:notrials,function(i) sd(resp.vec[which(rel.vec==i)])/sqrt(length(which(rel.vec==i))) ) 
  #get the average number of rel.,unrel., and sig. cells across trials
  res.trials <- campClassifyTrialCells(data = data.ls)
  res.trials.df <- joinListDFs(res.trials[[4]]) #gets all the odor trials in a single DF
  res.trials.df <- cbind(res.trials.df,apply(res.trials.df, 1, sum)) #add a col giving the number of sig. cells
  names(res.trials.df)[3] <- 'sig'
  #odors, sig. cells then rela and unrel cells
  res.odors.no <- list(sapply(res.trials[[2]],length)/percent.factor,sapply(res.trials[[3]],length)/percent.factor ) # no of rel and unrel cells, & sig. cells per odor
  #now, the stats for these things, linear fit first, then exponential fit
  #cat('\n',str(resp.trial)) 
  resprel.stats <- list(cor(1:notrials,resp.trial),fitlinearlaw(1:notrials,resp.trial),
                        cor(1:notrials,log10(resp.trial)),fitexplaw(1:notrials,resp.trial))
  names(resprel.stats) <- c('cor.linear','linear fit','cor.exp','exp fit')

  trial.stats <- list(apply(res.trials.df,2,mean),apply(res.trials.df,2,sd),apply(res.trials.df,2,sd)/sqrt(nrow(res.trials.df)))
  names(trial.stats) <- c('trial means','trial sd','trial sem')
  odor.stats <- list(mean(res.odors.no[[1]]),mean(res.odors.no[[2]]),sd(res.odors.no[[1]]),sd(res.odors.no[[2]]),
                     sd(res.odors.no[[1]])/sqrt(length(res.odors.no[[1]])),sd(res.odors.no[[2]])/sqrt(length(res.odors.no[[1]])) )
  names(odor.stats) <- c('tot. rel. no/odor mean','tot. unrel. no/odor mean','tot. rel. no/odor sd','tot. unrel. no/odor sd',
                         'tot. rel. no/odor sem','tot. unrel. no/odor sem')
  resdat.ls <- list(cbind(1:notrials,resp.trial,resp.trial.sem),res.trials.df,res.odors.no[[1]],res.odors.no[[2]])   
  names(resdat.ls) <- c('rel. resp. data','trial data','total rel. odor','total unrel. odor')
  res <- list(resprel.stats,trial.stats,odor.stats,resdat.ls,odor.dets)  
  names(res) <- c('resp. rel. stats','trial stats','odor stats','all the data','odor details')
  res
}

#gets the statistics of reliable and unreliable cells for a list of data sets: Mainly the ones in Fig. 1 siuch as 
#R and U cells per odor, per trial, and correlation and eqn fit between reliability and response rate
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
#op:1: get a list of results, 2: get , 3- only the odor details
campGetReliabStatsLists <- function(data.ls,skip=c(),respop=1,op=1){

  res.ls <- lapply(data.ls, function(x){
    dat <- campGetReliabStats(data.ls = x,skip = skip,respop = respop)
    res1 <- c(dat$`resp. rel. stats`$cor.linear,dat$`resp. rel. stats`$`linear fit`,
              dat$`resp. rel. stats`$cor.exp,dat$`resp. rel. stats`$`exp fit`) #the resp rel eqn fits
    names(res1) <- c('lin. fit','a','b','exp. fit','a','b')
    res2 <- c(dat$`trial stats`$`trial means`,dat$`odor stats`$`tot. rel. no/odor mean`,dat$`odor stats`$`tot. unrel. no/odor mean`)
    names(res2) <- c('#rel.trials','#unrel.trials','#sig.trials','#rel.odors','#unrel.odors')
    res3 <- unlist(dat$`odor details`)
    #cat('\nres',dat$`resp. rel. stats`$cor.linear)
    #res2 <- dat[[2]][[1]][1:2]
    #res3 <- dat[[3]][1:2]
    switch(op,dat,c(res1,res2,res3))
  })
  res.ls
}  

#gets the general correlation between all odor-pairs for reliable and unreliblae cells and the whole cells popultion
#thresh.rel - correlactions above this amount are considered similar
#thresh.unrel - correlations below this amount aee considered dissimilar
campCompCorCells <- function(dat.ls,thresh.rel=0.5,thresh.unrel=0.15,op=1){
  res.all <- computeListMatCorr(alldata = dat.ls,matop = 2,celltype = 4,op=op)
  res.rel <- computeListMatCorr(alldata = dat.ls,matop = 2,celltype = 2,op=op)
  res.unrel <- computeListMatCorr(alldata = dat.ls,matop = 2,celltype = 3,op=op)
  
  #odor.corr <- computeListMatCorr(alldata = data.lst,celltype=4,matop=2) #correlations of all odor pairs
  res.all.flat <- flattenDFtoList(res.all)
  res.rel.flat <- flattenDFtoList(res.rel)
  res.unrel.flat <- flattenDFtoList(res.unrel)
  odorpairs <- names(res.all.flat)[which(res.all.flat>=thresh.rel)] #names of odor pairs that satisfy the criterion
  #cat('\n',odorpairs)
  odorpairs <- odorpairs[getIndexSameSubStr(odorpairs,op=2)]
  
  res <- c(cor(res.all[upper.tri(res.all)],res.rel[upper.tri(res.rel)]),
    cor(res.all[upper.tri(res.all)],res.unrel[upper.tri(res.unrel)]),
    mean(res.all[upper.tri(res.all)]),mean(res.rel[upper.tri(res.rel)]),mean(res.unrel[upper.tri(res.unrel)]) )
  names(res) <- c('cor.all_rel','cor.all_unrel','cor.all','cor.rel','cor.unrel')
  res
  c(cor(unlist(res.rel.flat[odorpairs]),unlist(res.all.flat[odorpairs])),
    cor(unlist(res.unrel.flat[odorpairs]),unlist(res.all.flat[odorpairs])) )
  list(unlist(res.rel.flat[odorpairs]),unlist(res.unrel.flat[odorpairs]),unlist(res.all.flat[odorpairs]) )
  #todo: do cosines and correlations, also maybve just get average correlations. Flip side do dissimilar odors, too.
}

#Result: ok. for each percentile we calculate the cells that are in common, and their product. That's what contributes to the numerator.
#then see with increasing percentiles, how this number changes. This would also require a quick calculation of how the distance metric
#of A and B changes with percentile.

#function that calculate how often a reliable cell is reliable for another odor, and how often it is unreliable for another odor
#the number of times it is silent is noodors - (sum of the first two)
campCalcRelFreq <- function(alldata,op=1){
  #trial details
  trialdets <- campGetTrialDetails(alldata = alldata,op = 2)
  noodors <- length(trialdets[[2]])
  rel <- campGenRUCellMat(alldata,op=2,retype = 5) #rel. cells
  unrel <- campGenRUCellMat(alldata,op=3,retype = 5) # unrel cells
  rel.set <- unique(unlist(rel))
  unrel.set <- unique(unlist(unrel))
  res.rel <- sapply(rel.set,function(x){
    res <- sum(sapply(1:noodors,function(i) x %in% rel[[i]] ))
  })
  names(res.rel) <- rel.set
  res.rel
}

#an alis function to be backward compatible.
campCompPercentRanges <- function(...){
  campCompPercentOdors(...)
}


#function that gets the overlap for different deciles
#odorpairs: do the overlap for these odor pairs across different deciles
#op: 1 - do overlap, 2 - do nonoverlap, 3 - percentage overlap, 4 - percent nonoverlap, 5 - ratio of continuous overlap/non-overlap, 6 - overlap/nonoverlap
#simop: the option for measuring similarity, 1 = corre, 2 - cosine
#trialop: 1 get the result by computing it across all trial combinations, and averaging, leave out same-trial pairs
#trialop: 2 get the result by computing result for averaged trial response
campComputeDecileOver <- function(dat.ls,odorpairs=c(),thresh=0.5,deciles=seq(25,100,25),decsize=25,simop=1,trialop=1,op=1){
  #get the correlation for all odors
  res.all <- computeListMatCorr(alldata = dat.ls,matop = 2,celltype = 4,op=simop)
  
  #get the flattened vectors
  res.all.vec <- unlist(flattenDFtoList(res.all))
  res.upper <- unlist(flattenDFtoList(upper.tri(res.all)))
  
  #now, get the overlap for each decile or percentile
  overop <- switch(op,1,5,6,7,8,9) # the kind of overlap or non-overlap you want
  res.over.ls <- lapply(deciles,function(x){
    #res.over <- campCompPercentRanges(dat.ls,params = c(x,decsize),op=overop)
    res.over <- campCompPercentOdors(dat.ls,params = c(x,x),trialop = trialop,op=overop)
    unlist(flattenDFtoList(res.over[[1]]))
  })
  res.over.df <- (convertNestedListsDF(res.over.ls) )
  colnames(res.over.df) <- deciles
  #cat('\n',str(res.over.df),str(res.all.vec))
  res.df <- cbind.data.frame(res.upper,names(res.all.vec),res.all.vec,res.over.df)
  posns <- which(res.df[,1]==T) #only get the upper triangular matrix
  res.df <- res.df[posns,-1]
  res.df
}


#compares (the statistics of) the overlap between different percentile ranges
#poutput: the number of overlapping cells that don't include similar elements
#matop:
#params:
#retype:
#op:1 - the no of overlapping cells, 2 - dot product of overlapping cells, 3 - relative dot product, i.e., dot product/norm(a)*norm(b)
# 4 - norm product for the percentile, 5 - the non-overlap between the two odors, 6 - percentage overlap, 7 - percent nonoverlap
# 8 - continuous overlap/non-overlap, #9 - overlap/non-overlap, 10 - total number of responding cells, 11 - continuous overlap
#12 - continuous non-overlap, 13 - asymmetric continuous non-overlap,
#trialop is 2 for op=2,3,4.
#trialop: 1 get the result by computing it across all trial combinations, and averaging, leave out same-trial pairs
#trialop: 2 get the result by computing result for averaged trial response
campCompPercentOdors <- function(alldata,cells=c(),matop=2,params=c(10,10),trialop=1,op=1){
  dat.ls <- campGenRUCellMat(alldata,op=5,matop = 1,params = params,retype = 3) #select the kind of cells we want to get correlations for
  #cat('\n',str(dat.ls))
  if(length(dat.ls)<2) return(NULL) #if there is only one or no odors, there is nothing to compare it with
  notrials <- campGetTrialDetails(alldata = alldata,op=2)[[1]] #get no of trials for each odor
  stim.norms <- sapply(names(dat.ls[[2]][[2]]),function(x) lpnorm(unlist(dat.ls[[2]][[2]][x])) )
  #iterate through odor pairs, get the number of cells in common, and also get their dot. product.
  res.ls <- sapply(names(dat.ls[[1]]),function(i){
    res.in <- sapply(names(dat.ls[[1]]),function(j){
      #get number of cells in common
      #cat('\n',i,j,':')
      over <- campComparePercentile(dat1.df = dat.ls[[1]][[i]], dat2.df = dat.ls[[1]][[j]],
                            dat1.avg = dat.ls[[2]][[2]][i],dat2.avg = dat.ls[[2]][[2]][j],trialop = trialop,op=op)
      over
    })
  })
  #print(res.ls)
  res.no <- sum(res.ls[upper.tri(res.ls,diag = F)])
  list(res.ls,res.no)
}


#function that looks at two stimulus responses and then gives a percentile breakup of their overlap and other 
#chatracteristics. The two stimuli are in the form of data frames or matrices, with columns giving the response values across trials, 
# and rownames giving the position or cell no that is providing the responses. Zero values are left out.
#datx.df: is the response data with responding cells along rows and trials along columns
#datx.avg: is the average trial response across all cells.
#op:1 - the no of overlapping cells, 2 - dot product of overlapping cells, 3 - relative dot product, i.e., dot product/norm(a)*norm(b)
# 4 - norm product for the percentile, 5 - the non-overlap between the two odors, 6 - percentage overlap, 7 - percent nonoverlap
# 8 - continuous overlap/non-overlap, #9 - overlap/non-overlap, 100 - test condition, 10 - total number of responding cells, 11 - continuous overlap
#12 - continuous non-overlap,13 - asymmetric continuous non-overlap,
#trialop is 2 for op=2,3,4.
#trialop: 1 get the result by computing it across all trial combinations, and averaging, leave out same-trial pairs
#trialop: 2 get the result by computing result for averaged trial response
campComparePercentile <- function(dat1.df,dat2.df,dat1.avg,dat2.avg,trialop=1,op=1){
  #get number of cells in common
  stim1 <- rownames(dat1.df) #cells for stim 1
  stim2 <- rownames(dat2.df) #cells for stim2
  total <- length(union(stim1,stim2)) #total number of cells from both stimuli
  overlap <- intersect(stim1,stim2) #cells in common
  #cat('\toverlap',length(overlap),'total:',total,'nonoverlap',length(union(setdiff(stim1,overlap),setdiff(stim2,overlap))) )
  nonoverlap <- setdiff(union(stim1,stim2),intersect(stim1,stim2)) #cells that are not in common
  if(op <= 9 ) choice <- switch(op,1,2,3,4,2,3,4,5,6)
  else choice <- op
  #all these have to be done over all trials, and then averaged or done on the average of all trials
  if(trialop==1 && !(op %in% c(2:4)) ){#second condn: for op=2:4, we look at trial averaged response only
    res.mn <- campComparePercentileTrial(dat1.df = dat1.df,dat2.df = dat2.df,dat1.avg = dat1.avg,dat2.avg = dat2.avg,choice=choice,op=op)
  }
  if(trialop==2){#we are looking at average response comparisons
    res.mn <- switch(choice,
                     length(overlap),
                     length(c(setdiff(stim1,overlap),setdiff(stim2,overlap))),
                     length(overlap)/total,
                     length(c(setdiff(stim1,overlap),setdiff(stim2,overlap)))/total
    )
    if(choice==5){
      common <- intersect(stim1,stim2)
      #cat('\ncommon',common,':',stim1,':',stim2,str(dat1.avg) )
      common.over <- sum(unlist(sapply(common,function(s) min(dat1.avg[s,1],dat2.avg[s,1]) ) ) )
      common.nonover <- sum(abs(dat1.avg[common,1]-dat2.avg[common,1]) )
      noncommon <- sum(dat1.avg[setdiff(stim1,stim2),1]) + sum(dat2.avg[setdiff(stim2,stim1),1])
      if(identical(dat1.df,dat2.df) ) res.mn <- 1
      else res.mn <- common.over/(common.nonover+noncommon) 
      #cat('\ncommon',common,':',common.over,':',common.nonover,':',noncommon )
    }
    if(choice==6) 
      #cat('\n',length(c(setdiff(stim1,overlap),setdiff(stim2,overlap))) )
      res.mn <- length(overlap)/length(c(setdiff(stim1,overlap),setdiff(stim2,overlap)))
    if(choice==100) {#the test condition
      tmp <- total
    }
  }
  #cat('\tres:',res.mn)
  res.mn
}

#this is the option where we do comparisons between every trial for two odors and then average the results
campComparePercentileTrial <- function(dat1.df,dat2.df,dat1.avg,dat2.avg,choice,op=1){
  #go through all the trials for both odors
  res.ls <- sapply(1:ncol(dat1.df),function(i){
    res.en <- sapply(1:ncol(dat2.df),function(j){
      #cat('\n',str(dat1.df[setdiff(stim1,overlap),i]),':',stim1,':',overlap,':',setdiff(stim1,overlap),'\n',choice )
      if(choice==3 || choice ==4){
        total <- length(union(names(dat1.df)[which(dat1.df[,i]>0)],names(dat2.df)[which(dat2.df[,j]>0)]))
      }
      tmp <- switch(choice,
                    countNonZeroes(dat1.df[overlap,i]*dat2.df[overlap,j]), #gets the overlap of cell between 2 trials       
                    #the cells of stim1 not overlapping and stim 2 not overlapping
                    countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i]),
                    countNonZeroes(dat1.df[overlap,i]*dat2.df[overlap,j])/total, #gets the overlap of cell between 2 trials       
                    (countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i]))/total
      )
      if(choice==5){#continuous overlap
        stim1 <- rownames(dat1.df)[which(dat1.df[,i]>0)]
        stim2 <- rownames(dat2.df)[which(dat2.df[,j]>0)]
        common <- intersect(stim1,stim2)
        common.over <- sum(unlist(sapply(common,function(s) min(dat1.df[s,i],dat2.df[s,j]) ) ) )
        common.nonover <- sum(abs(dat1.df[common,i]-dat2.df[common,j]) )
        noncommon <- sum(dat1.df[setdiff(stim1,stim2),i]) + sum(dat2.df[setdiff(stim2,stim1),j])
        if(identical(dat1.df,dat2.df) ) tmp <- 1
        else tmp <- common.over/(common.nonover+noncommon)  
      }
      if(choice==6) {
        tmp <- countNonZeroes(dat1.df[overlap,i]*dat2.df[overlap,j])/(countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i]))
        #cat('\t',countNonZeroes(dat1.df[setdiff(stim1,overlap),i]),countNonZeroes(dat2.df[setdiff(stim2,overlap),i]) )
      }
      if(choice==100) {#the test condition
        tmp <- total
      }
      if(choice==10) {#total
        tmp <-  length(union(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)])) 
      }
      if(choice==11){#continuous overlap
        common <-  intersect(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)])
        tmp <- sum(unlist(sapply(common,function(s) min(dat1.df[s,i],dat2.df[s,j]) ) ) )
      }
      if(choice==12){#continuous non-overlap
        noncommon1 <-  setdiff(rownames(dat1.df)[which(dat1.df[,i]>0)],
                              intersect(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)]) )
        noncommon2 <-  setdiff(rownames(dat2.df)[which(dat2.df[,j]>0)],
                               intersect(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)]) )
        common <-  intersect(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)])
        tmp <- sum(unlist(sapply(common,function(s) abs(dat1.df[s,i]-dat2.df[s,j]) ) ) )
        tmp <- tmp + sum(dat1.df[noncommon1,i]) + sum(dat2.df[noncommon2,j])
      }
      if(choice==13){#continuous non-overlap asymmetric
        noncommon1 <-  setdiff(rownames(dat1.df)[which(dat1.df[,i]>0)],
                            intersect(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)]) )
        common <-  intersect(rownames(dat1.df)[which(dat1.df[,i]>0)],rownames(dat2.df)[which(dat2.df[,j]>0)])
        tmp <- sum(unlist(sapply(common,function(s) abs(dat1.df[s,i]-dat2.df[s,j]) ) ) )
        tmp <- tmp + sum(dat1.df[noncommon1,i])
      }
      
      tmp
    })
    #cat('\n',res.en)
    #cat('\n',(countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i])) )
    res.en #overlap of dat1 col i with dat2 all cols
  })
  #cat('\n',res.ls)
  res.mn <- mean(cleanNAVec(res.ls))
  
}

#this is when we are going to look at or compare between averages of responses
campComparePercentileAvg <- function(dat1.df,dat2.df,dat1.avg,dat2.avg,choice,op=1){
  res.mn <- switch(choice,
                   length(overlap),
                   length(c(setdiff(stim1,overlap),setdiff(stim2,overlap))),
                   length(overlap)/total,
                   length(c(setdiff(stim1,overlap),setdiff(stim2,overlap)))/total
  )
  if(choice==5){
    common <- intersect(stim1,stim2)
    #cat('\ncommon',common,':',stim1,':',stim2,str(dat1.avg) )
    common.over <- sum(unlist(sapply(common,function(s) min(dat1.avg[s,1],dat2.avg[s,1]) ) ) )
    common.nonover <- sum(abs(dat1.avg[common,1]-dat2.avg[common,1]) )
    noncommon <- sum(dat1.avg[setdiff(stim1,stim2),1]) + sum(dat2.avg[setdiff(stim2,stim1),1])
    if(identical(dat1.df,dat2.df) ) res.mn <- 1
    else res.mn <- common.over/(common.nonover+noncommon) 
    #cat('\ncommon',common,':',common.over,':',common.nonover,':',noncommon )
  }
  if(choice==6) 
    #cat('\n',length(c(setdiff(stim1,overlap),setdiff(stim2,overlap))) )
    res.mn <- length(overlap)/length(c(setdiff(stim1,overlap),setdiff(stim2,overlap)))
  if(choice==100) {#the test condition
    tmp <- total
  }
  
  
}



#gets the statistics on the reliability cells of each degree
#alldata: the return of getCampDataSigList
#relscore: 1, average of both odors, 2 - reliability of first odor
#optype=1, overlap; 2, difference, 3 - probability difference
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#9 - overlap, all cells
#op=1, the sum of reliability cells of the same type, 2 - mean of all reliability cells of the same type, 
campOverlapAnalogRel <- function(alldata,notrials,skip=c(),relscore=2,optype=4,op=1){
  #compute the usuefull difference for all cells for odors pairs
  alloverlap <- campOverlapAllAnalog(alldata,relscore = relscore,skip = skip,op = optype)
  #now, compute the sums of differences for each cell
  overlap.cols <- switch(op,lapply(alloverlap, function(x) sumDF(x)),lapply(alloverlap, function(x) meanDF(x)) )
  #cat('\ncOAR',str(overlap.cols),str(alloverlap))
  #make a data frame with overlap/difference showing for each reliability cell level
  overlap.cols <- transposeDF(combineListOfUnevenDfs(overlap.cols,factorcol = 1,cols = 2,padval = 0))
  names(overlap.cols) <- overlap.cols[1,] #the first row, which holds the reliability numbers
  overlap.cols <- overlap.cols[-1,] #we dont need the first row, which holds the reliability numbers
  #gets the number of odor trials for the first odor of each pair
  odortrials <- sapply(rownames(overlap.cols), function(x) notrials[strsplit(x,split = ',')[[1]][1]])
  names(odortrials) <- rownames(overlap.cols) #gotta reassign the names
  #cat('\ncOARodor trials',rownames(overlap.cols),'\not',odortrials,'\n',notrials,str(notrials),'\n',rownames(res))
  trialno <- odortrials[rownames(overlap.cols)] #no of trials for the first odor of the pair
  res.df <- cbind.data.frame(overlap.cols,trialno)
  res.df
}

#do a sexcond data frame that contains the fit, contribution of reliable and unreliable cells.
#idea: have to introduce a correlation calculation with the overlap function for each pair of odors
#could be done in the pairanalog function above. just do a correlation between reliability and usefulness or difference, i.e., col 1 and 2
#also calcualte it at the end with whatever is calling the campOverlapAnalogAll function. Something like campAnalogStats

#function that compares the contributes of reliable and unreliable cells with and without saturation
#matchlevel: amout of correlation match between HAllem and rob that is required
#halrobmatch: how to balance the two sets: 1 - only odors present in both, no similar odor pairs 
#2 - only odors present in both, with similar odor pairs,
#3 -rob odors, similar odor pairs ok, 4 - rob odors, similar odor pairs ok
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#op=1, summed data, 2 - mean data
campCompareSatNorm <- function(hallem,matchlevel=0.1,halrobmatch=1,skip = c('empty','paraff'),op=1){
  alldata <- getCampDataSigList() #get the data
  #get the correlations for hallem and rob
  rescor <- campMatchHalRobCor(hallem,matchlevel = matchlevel,halrobmatch = halrobmatch,skip = skip)  # will get all correaltions, we are interested in the first 2 columms
  if (length(rescor)==0) return(list(NULL,NULL,NULL,NULL,NULL)) #nothing more to do, all null
  rescor.sat <- campMatchHalRobCor(hallem,matchlevel = matchlevel,halrobmatch = halrobmatch,skip=skip,optype = 5)  # now, saturation
  odors <- rownames(rescor) #we want to examine odors common to rescor and res33: rescor subset of res33
  #get non-saturated data
  res33 <- campOverlapAllAnalog(alldata,skip = skip,op = 4)#op=4, normal
  res33 <- campAnalogStats(res33,dat.ls = alldata,op=5) #get all the stats
  res33 <- res33[odors,]
  res35 <- campOverlapAllAnalog(alldata,skip = skip,op = 6)#op=6, which cells are active for both odors
  res35 <- campAnalogStats(res35,dat.ls = alldata,op=6) #get the cell overlap stats
  res35 <- res35[odors,]
  #saturated data
  res34 <- campOverlapAllAnalog(alldata,skip = skip,op = 5)#op=4, saturated 
  res34 <- campAnalogStats(res34,dat.ls = alldata,op=5) #get all the stats, slope, difference ...
  res34 <- res34[odors,]
  #put it all togetehr, odor similarity(hallem and rob), slope, unreliable and reliable cell contributtions: normal first then saturated
  #cat('\nnoyo1',str(res33),str(rescor),str(res34),str(res35))
  res <- cbind(rescor[odors,c(2,3)],res33[odors,2:4],res34[odors,2:4],res33[odors,c(7,8)],rescor[odors,ncol(rescor)],res35[odors,c(1,2)])
  names(res)[3:13] <- c('slope.nor','ursum.nor','rsum.nor','slope.sat','ursum.sat','rsum.sat','ur.cells','r.cells','trials','urel.act','rel.act') 
  list(res,res33,res34,rescor,rescor.sat)
}

#perform analysis on the return of the campNormSatNorm function
#similar: the amount of similarity threshold above which an odor pair is considered similar
#dissimilar: the amount of similarity threshold below which an odor pair is considered dissimilar
#op: 1 - get stats , 2: get the plots for unreliable and reliable, for normal and saturated in that order
campGetSatNormStats <- function(dat.ls=c(),hallem=c(),similar=0.3,dissimilar=0.15,op=1){
  if (length(dat.ls) == 0) data.ls <- campCompareSatNorm(hallem = hallem)
  else data.ls <- dat.ls
  dat.anal <- data.ls[[1]]
  #now, get similar and dissimilar, reliable vs unreliable,
  res.dissimilar <- dat.anal[which(dat.anal[,2] < dissimilar),] #filter all cells that are below the dissiilar threshold
  res.similar <- dat.anal[which(dat.anal[,2] > similar),] #filter all cells that are above the similar threshold
  stats.dis <- sapply(c(4,7,5,8), function(x) getStats(res.dissimilar[,x]))
  colnames(stats.dis) <- c('unrel.norm','unrel.sat','rel.norm','rel.sat')
  row.names(stats.dis) <- c('mean','sd','sem') 
  stats.sim <- sapply(c(4,7,5,8), function(x) getStats(res.similar[,x]))
  colnames(stats.sim) <- c('unrel.norm','unrel.sat','rel.norm','rel.sat')
  row.names(stats.sim) <- c('mean','sd','sem') 
  switch(op,list(stats.dis,stats.sim,res.dissimilar,res.similar),list(list(res.dissimilar[,4],res.dissimilar[,7]),list(res.dissimilar[,5],res.dissimilar[,8]),
                                                                      list(res.similar[,4],res.similar[,7]),list(res.similar[,5],res.similar[,8])) )
}

#function that gets the stats that are required on the output of campOverlapAllAnalog.
#data.ls: all the odor pairs
#dat.ls: the return of campGEtDataSigList to get the number of trials
#notrials: the total number of trials
#samee: F - take out all odor pairs that are similar, T - keep similar odor pairs
#op=1, does a reliability vs overlap measure correlation and slope
#op=2, contribution of unreliable cells, contribution of reliable cells (summed)
#op=3, contribution of unreliable cells, contribution of reliable cells (mean)
#op=4, # of unreliable cells, # of reliable cells
#op=5, all the above results
#op=6, contribution of reliable and unreliable cells (summed) and their number
campAnalogStats <- function(data.ls,dat.ls=c(),notrials=c(),same=T,op=1){
  if(length(data.ls) == 0) return(NULL) #nothing to do
  if(length(notrials)>0) notrials <- notrials
  else if(length(dat.ls)>0) notrials <- campGetTrialDetails(alldata = dat.ls,op=2)[[1]]
       else return(NULL)
  #cat('\nnotrials: ',notrials)
  samestr <- getIndexSameSubStr(names(data.ls))
  if(same==F) dat.ls <- data.ls[-samestr]
  else dat.ls <- data.ls
  #cat('\ncAS',names(dat.ls),'\nop = ',op)
  #skip elements or odor pairs that are the same
  res <- switch(op,
  sapply(names(dat.ls), function(x){#1.go through the data list and get slope and fit slope
    corr <- cor(dat.ls[[x]][,1],dat.ls[[x]][,2])
    fitslope <- fitlinearlaw(dat.ls[[x]][,1],dat.ls[[x]][,2])
    c(corr,fitslope[2])
  }),
  sapply(names(dat.ls), function(x){#2.get contribution of unreliable and reliable cells
    odor <- getSubString(x,split=',',sub = 1)
    unrel <- sum(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel <- sum(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    c(unrel,rel)
  }),
  sapply(names(dat.ls), function(x){#3.get contribution of unreliable and reliable cells
    odor <- getSubString(x,split=',',sub = 1)
    unrel <- mean(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel <- mean(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    c(unrel,rel)
  }),
  sapply(names(dat.ls), function(x){#4.get contribution of unreliable and reliable cell numbers
    #cat('\nop',4)
    odor <- getSubString(x,split=',',sub = 1)
    unrel <- length(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel <- length(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    c(unrel,rel)
  }),
  sapply(names(dat.ls), function(x){#5.all together, 
    odor <- getSubString(x,split=',',sub = 1)
    #todo: there is no overlap between these two odors. So return 0? Check the function which returns
    #these results: campOverlapAllAnalog.
    #cat('\ncAS',x,identical(dat.ls[[x]][,1],dat.ls[[x]][,2]),dat.ls[[x]][,1],'\t',dat.ls[[x]][,2])
    #cat('\n',str(dat.ls[[x]]))
    if(elemsIdentical(dat.ls[[x]][,1]) || elemsIdentical(dat.ls[[x]][,2]) || nrow(dat.ls[[x]])==0) corr <- 0
    else corr <- cor(dat.ls[[x]][,1],dat.ls[[x]][,2])
    fitslope <- fitlinearlaw(dat.ls[[x]][,1],dat.ls[[x]][,2])
    unrel.sum <- sum(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel.sum <- sum(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    unrel.mean <- mean(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel.mean <- mean(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    unrel.no <- length(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel.no <- length(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    tmp <- c(corr,fitslope[2],unrel.sum,rel.sum,unrel.mean,rel.mean,unrel.no,rel.no)
    names(tmp) <- c('rel.res','slope','urel.sum','rel.sum','urel.mn','rel.mn','urel.no','rel.no')
    #cat('\ncas:',tmp)
    tmp
  }),
  sapply(names(dat.ls), function(x){#all together, 
    odor <- getSubString(x,split=',',sub = 1)
    unrel.sum <- sum(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel.sum <- sum(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    unrel.no <- length(dat.ls[[x]][which(dat.ls[[x]][,1] <= floor(notrials[odor]/2)),2])
    rel.no <- length(dat.ls[[x]][which(dat.ls[[x]][,1] > floor(notrials[odor]/2)),2])
    tmp <- c(unrel.sum,rel.sum,unrel.no,rel.no)
    names(tmp) <- c('urel.sum','rel.sum','urel.no','rel.no')
    tmp
  }))
  res <- cleanNAMat(res) #turn NAs to 0, and we're off!
  t(res)
}


#gets various statistics given: Actually not sure what it does. 
campPrepareStatsDf <- function(dat.df){
  #get rid of the elements where the overlap is 0
  res.cor <- sapply(3:ncol(dat.df), function(x){
    tmp <- dat.df[,c(1,x)]
    res <- getThreshDF(tmp,col = 2,op=3)
    cor(res[,1],res[,2])
  })
  res.fit <- sapply(3:ncol(dat.df), function(x){
    tmp <- dat.df[,c(1,x)]
    res <- getThreshDF(tmp,col = 2,op=3)
    fitlinearlaw(res[,1],res[,2])
  })
  res.mn <- sapply(3:ncol(dat.df), function(x){
    tmp <- dat.df[,c(1,x)]
    res <- getThreshDF(tmp,col = 2,op=3)
    mean(res[,2])
  })
  
  list(res.fit,res.cor,res.mn)
}

#gets stats on the match of Rob;s and Hallem's results
#dat.df: The return of the campMatchHallemRob function
#relno: the number of relibability leve;ls
campHallRobGetStats <- function(dat.df,relno=6*2,op=1){
  #get the fit
  res1 <- sapply(2:(1+relno), function(x) fitlinearlaw(dat.df[,1],dat.df[,x]))
  #get the correlation
  res2 <- sapply(2:(1+relno), function(x) cor(dat.df[,1],dat.df[,x]))
  #get fit rsquare
  #tmp <- lapply(2:7, function(x) fitlinearlaw(res3[,1],res3[,x],op=3))
  #res3 <- sapply(tmp, function(x) x$r.squared)
  list(res1,res2)
}

#given a data set will generate a lst of dataframes, where each data frame is just the values of the 
#reliable or unreliable cells
#op=1, significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - specified cells, 6 - cells with a specifed reliability
#where reliability is in params, Could be one number or multiple numbers
#headerop: how we get header information, 1 = look through the folder, 2 - look through alldata
#retype: 1- datframe, 2 - list, 3, 4 - same as 1,2 except it includes list(results,campgetreiableResp), 5 - posns of the cells
#topparams: for op=5, params(topnstart,size) topnstart gives the start of the percentile range, and size gives the size, e.g., (50,5) is top 50 to 45
#matop: if you need to return the cell responses as 1 - a list of dfs or 2 - a list of vectors: averages of all of trials 
#respop: for the average response: 1 - average of all the trials it responds in, 
#2 - average based on all trial responses, even those where it does not respond
#cellop: whether you return responses of all cells or just responding cells. 1 - only responding cells, 2 - all cells
campGenRUCellMat <- function(alldata,headerop=1,retype=1,matop=1,params=c(),respop=1,cellop=1,op=1){
  #get trial details and nocells
  assignVarVals(c('trials','nocells'),campGetTrialDetails(alldata = alldata,op=2)[c(1,3)])
  #cat('\ntrials dets',str(trials))
  allnames <- names(trials) #get the names of the odors in the order of the trials
  if(matop==1){
    dat.ls <- getGroupedListsDf(alldata,op=2) #get response values
    rel.ls <- getGroupedListsDf(alldata,op=1) #get reliability values
  } 
  res <- campGetReliabResp(data.lst = alldata,respop = respop) #get reliability and response rates for the averages
  if(op==5) cells <- campGetTopNCells(data.lst = alldata,topn = params[1],deciles = params[2],op=2)
  if(op==6) {
    if(length(params)==0) pars <- 1:trials
    else pars <- params
  }
  #now, get the posns for the options specified
  posns.lst <- switch(op,
         lapply(allnames,function(i){#=1, significant cells
           posns <- which(res[[1]][i]>0)
         }),
         lapply(allnames,function(i){#2 - reliable cells
           posns <- which(res[[1]][i]>(trials[i]/2))
         }),
         lapply(allnames,function(i){#3 - unreliable cells
           posns <- which(res[[1]][i]<=(trials[i]/2) & res[[1]][i]>0)
         }),
         lapply(allnames,function(i){#4 - all cells
           posns <- 1:length(alldata[[1]][[1]]) #no of cells in the first trial
         }),
         lapply(allnames,function(i){#5 - specified cells
           posns <- cells[[i]]
         }),
         lapply(allnames,function(i){#6 - cells with specified reliability
           #cat('\n',str(res[[1]][i]),';',pars)
           posns <- which(res[[1]][[i]] %in% pars) #for each elem in params, gets their corresponding positions in res[[1]][i]
         }))
  names(posns.lst) <- allnames #name results appropriately
  #cat('\nposns',str(trials),'\n')
  if(matop==2) cells <- lapply(allnames,function(i) {#now pick all the cells for each  odor
    #cat('\n',i,str(res[[2]]),unlist(res[[2]][,i]),'\n',posns.lst[[i]] )
    trialavg <- unlist(res[[2]][i])[posns.lst[[i]] ]  #pick the appropritate col. and posn cells
    names(trialavg) <- posns.lst[[i]]
    #return all cells or just responding cells, and if there no cells of this type, return 0 vector
    if(cellop==2 || length(trialavg)==0) trialavg <- fillVector(vec=trialavg,pos=posns.lst[[i]],len=nocells)
    switch(retype,as.data.frame(trialavg),trialavg,as.data.frame(trialavg),trialavg)
  })
  else cells <- lapply(allnames,function(i){#now pick all the cells for each trial of every odor
    #cat('\ntrials',i,str(dat.ls[[i]]),trials[i] )
    trial.cells <- lapply(1:trials[i],function(j) {
      tmp <- dat.ls[[i]][,j]*rel.ls[[i]][,j]
      resp <- tmp[posns.lst[[i]] ]
      names(resp) <- posns.lst[[i]] #the names of each vector carry that position
      #return all cells or just responding cells, and if there no cells of this type, return 0 vector
      if(cellop==2 || length(resp)==0) resp <- fillVector(vec=resp,pos=posns.lst[[i]],len=nocells)
      resp
    })
    #cat('\ntrial',str(trial.cells))
    trial.cells <- switch(retype,convertNestedListsDF(trial.cells),trial.cells,convertNestedListsDF(trial.cells),trial.cells)
  })
  names(cells) <- allnames
  switch(retype,cells,cells,list(cells,res),list(cells,res),posns.lst)
}


#function to compute the effect of alpha, the significance, on the number of reliable and unreliable cells
#foldname: the folder that holds the dataset
#op = 1 get the number of reliable and unreliable cells.
# 2: return the classified data, the getCampDAtaSigList structures
# 3: list of reliable and unreliable cells
# alpharange: the range or values of alphas for which the function is evaluated
campGetAlphaEffectRUCells <- function(foldname=".",alpharange=c(0.05,0.01,0.005,0.001,0.0005),op=1){
  #cat('\ncampGetAlpha')
  res.ls <- lapply(alpharange,function(x){
    #cat('\nalpha: ',x)
    dat.sig <- getCampDataSigList(foldname = foldname,alpha = x )
    trial.class <- campClassifyTrialCells(dat.sig)
    cells <- joinListDFs(trial.class[[4]])/length(dat.sig[[1]][[1]]) * 100
    res <- switch(op,apply(cells,2,mean), #reliable and unreliable cells
                  dat.sig,                #the data set
                  list(trial.class[[2]],trial.class[[3]]))    
  })
  names(res.ls) <- alpharange #set the names to alpha value (significance) that is being served
  res.ls
}

#function that takes the contents of the campGetAlphaEffectRUCells op =3, i.e., the reliable and unreliable cells list
#then tells you the cells or number of cells that move from being reliable to unreliable, and the number that drop out of the unreliable
#list. Assumption is that the list is ordered from the least significance to most significance
#dat.ls : the list of RU cells for each alpha
#compst: the starting index to use as a comparison point
#op: 1 : return the cells
#2: return the number of cells
campCompareAlphaRUcells <- function(dat.ls,compst=1,compend=length(dat.ls),op=1){
  noalphs <- length(dat.ls)
  odors <- names(dat.ls[[1]][[1]])
  #cat('odors ',odors)
  #first find the reliable cells that drop off into the unreliable lot, comparison is between compst and compend
  rel2ur.drop <- lapply(1:length(odors),function(i) {
    noint <- intersect(dat.ls[[compst]][[1]][[i]],dat.ls[[compend]][[2]][[i]])
  })
  names(rel2ur.drop) <- odors
  #next find the loss from reliable cells overall
  rel.drop <- lapply(1:length(odors),function(i) {
    noint <- setdiff(dat.ls[[compst]][[1]][[i]],dat.ls[[compend]][[1]][[i]])
  })
  names(rel.drop) <- odors
  #finally, number lost from the unreliable cells group
  unrel.drop <- lapply(1:length(odors),function(i) {
    noint <- setdiff(dat.ls[[compst]][[2]][[i]],dat.ls[[compend]][[2]][[i]])
  })
  names(unrel.drop) <- odors
  res <- switch(op,list(rel2ur.drop,rel.drop,unrel.drop),
                list( sapply(rel2ur.drop, function(x) length(x)),sapply(rel.drop, function(x) length(x)),
                      sapply(unrel.drop, function(x) length(x)) ) )
  names(res) <- c('rel2ur.drop','rel.drop','unrel.drop')
  res
}

#compares the loss in the number of cells from the reliable to unreliable cell group, and from the unreliable cell group
#compares losses between adjacent levels of significance
#dat.ls: is the list of data frames of the responses for different levels of significance
#op =1 : compare for group alpha and alpha + 1
#op =2 : give the number for each alpha
campCompareAlphaGroups <- function(dat.ls,op=1){
  alpharange <- names(dat.ls) 
  res.ls <- lapply(1:(length(dat.ls)-1), function(i){
    res <- campCompareAlphaRUcells(dat.ls = dat.ls,compst = i,compend = i+1,op=2)
    sapply(res,mean)
  })
  names(res.ls) <- alpharange[1:(length(alpharange)-1)]
  res.ls
}

#computes the useful discrimination given the a data set, i.e., the output of campgetdatasiglist
#This does it for a single data set
#meas.sim: measures the amount of similarity and dissimilarity c(dissim,similarity)
#dat.ls: the he output of campgetdatasiglist
#if op=2, list contains 2 vectors : list(dissimilar odors, similar odors)
#op:1 - the odors are determined here, 2 - the odors similar and dissimilar are predetermined.
#retop: 1 - normal return UD results, 2 - the dissim and sim odors, 
#3 - the dissim and sim odors as odor names, 4 - 1 and 2 as lists
campComputeUD <- function(dat.ls,meas.sim=c(0.15,0.5),odors=list(),retop=1,op=1){
  #get the correaltion to determine sim and dissim odors, then compare the UD for normal and extended training
  res.corr <- unlist(flattenDFtoList(computeListMatCorr(dat.ls,matop = 1)))
  #calculate the reliable and unreliable cell differences at the start and at saturation for all odors
  res.start <- campAnalogStats(campOverlapAllAnalog(dat.ls,op = 4),dat.ls = dat.ls,op=5)
  start.arrange <- insertColDf(res.start,newcol = res.corr[rownames(res.start)],colname = 'cor.')
  res.sat <- campAnalogStats(campOverlapAllAnalog(dat.ls,op = 5),dat.ls = dat.ls,op=5)
  sat.arrange <- insertColDf(res.sat,newcol = res.corr[rownames(res.sat)],colname = 'cor.')
  
  #choose dissimilar and similar odors
  if (op!=2){
    res.notsame <- getIndexSameSubStr(rownames(res.start),op=2)
    res2 <- intersect(which(start.arrange[,1]<meas.sim[1]),res.notsame)  #dissimilar
    res3 <- intersect(which(start.arrange[,1]>meas.sim[2]),res.notsame)   #similar
  } else {#op==2 for now, need it to make sure that two consecuitve calls in computeUDRel check the same odors
    res2 <- odors[[1]];res3 <- odors[[2]]
  } 

  #now, get the stats on similar and dissimilar, not saturated and saturated
  #cleanNAVec gives zeroes instead of NA: happens when there are no similar odors
  res.dis <- cleanNAVec(apply(res.start[res2,],2,mean)) #dissimilar
  res.dissat <- cleanNAVec(apply(res.sat[res2,],2,mean)) #dissimilar saturated
  res.sim <- cleanNAVec(apply(res.start[res3,],2,mean)) #similar
  res.simsat <- cleanNAVec(apply(res.sat[res3,],2,mean)) #similar saturated
  #assemble the results
  res <- list(res.dis,res.dissat,res.sim,res.simsat)
  names(res) <- c('dissim.st','dissim.sat','sim.st','sim.sat')
  res.odors <- switch(retop,list(res2,res3),list(res2,res3),
                      list(rownames(start.arrange)[res2],rownames(start.arrange)[res3]),list(res2,res3))
  names(res.odors) <- c('dissim','sim')
  switch(retop,res,res.odors,res.odors,list(res,res.odors) )
}

#function that calculates the UD for similar and dissimilar odors for a range of alphas, i.e., significance/
# alpharange: the range or values of alphas for which the function is evaluated
# foldname: the folder that holds the dataset
# meas.sim: measures the amount of similarity and dissimilarity c(dissim,similarity)
# cols: the cols that you are interested in. empty is get everything
campComputeUDAlphas <- function(foldname='.',alpharange=c(0.05,0.01,0.005,0.001,0.0005),cols=c(),meas.sim=c(0.15,0.5),op=1){
  #get the different campDatasiglist with different alphas
  datsets.ls <- campGetAlphaEffectRUCells(foldname = foldname,alpharange = alpharange,op=2)
  res.ls <- lapply(datsets.ls, function(x){
    res <- campComputeUD(x,meas.sim = meas.sim)
    if (length(cols)==0) res.df <- convertNestedListsDF(res)
    else res.df <- convertNestedListsDF(lapply(res,'[',cols))
    res.df
  })
  #cat('\n',str(res.ls))
  #res.ls
  res.df <- joinListDFs(res.ls)
  res.df  
}

#this function calculates the effect of any particular group on discrimination with 
#normal and extended training. Works like this. Do normal and extended training with normal group of 
#cells. Then, repeat with the reliability group of cells taken out
#dat.ls: the data that we are considering return of getCampDataSigList
#rel: the reliability level that we want to get the data for
#meas.sim: the measures for judging odors to be dissimilar or similar
#op: 1- two lists, 2 - a data frame first 4 cointrol, next 4 reliabiltyu out results, 3 - comparison
#of all the results
#odorop: determines whether you used the same odors for the comparisons (2) or determine the odors
#within the correlations for each data set (1) 
campComputeUDReliability <- function(dat.ls,meas.sim=c(0.15,0.45),rel=1,odorop=1,op=1){
  #determine the odors, and thats what you should use for both comparisons
  odors <- campComputeUD(dat.ls = dat.ls,meas.sim = meas.sim,retop = 2)
  contUD <- campComputeUD(dat.ls = dat.ls,meas.sim = meas.sim,op = odorop,odors = odors)
  #now, generate the data set with all rel cells removed and do the computation again
  dat.ls.rel <- campDataRemoveRel(dat.ls = dat.ls,rel = rel)
  contUD.rel <- campComputeUD(dat.ls = dat.ls.rel,meas.sim = meas.sim,op = odorop,odors = odors)
  res <- list(contUD,contUD.rel)
  #we can also just get the UD components
  #dissim
  res.df <- as.data.frame(res)
  dissim.st.ud <- sum(res.df[c('urel.sum','rel.sum'),'dissim.st'])-sum(res.df[c('urel.sum','rel.sum'),'dissim.st.1'])
  dissim.sat.ud <- sum(res.df[c('urel.sum','rel.sum'),'dissim.sat'])-sum(res.df[c('urel.sum','rel.sum'),'dissim.sat.1'])
  sim.st.ud <- sum(res.df[c('urel.sum','rel.sum'),'sim.st'])-sum(res.df[c('urel.sum','rel.sum'),'sim.st.1'])
  sim.sat.ud <- sum(res.df[c('urel.sum','rel.sum'),'sim.sat'])-sum(res.df[c('urel.sum','rel.sum'),'sim.sat.1'])
  res.cont.ud <- c(dissim.st.ud,dissim.sat.ud,sim.st.ud,sim.sat.ud)
  switch(op,res,res.df[c('urel.sum','rel.sum'),],res.cont.ud)
}


#calculates the reliability contributions for all levels of reliability. It fixes the similar and 
#dissimilar odors right at the start based on the original getCampDataSigList without any cells removed
#dat.ls: the data that we are considering return of getCampDataSigList
#rel: the reliability level that we want to get the data for
#meas.sim: the measures for judging odors to be dissimilar or similar
#op: 1-lists, 2 - a data frame 
campCompUDSingleRel <- function(dat.ls,meas.sim=c(0.15,0.45),rel=1,op=1){
  #determine the odors, and thats what you should use for both comparisons
  odors <- campComputeUD(dat.ls = dat.ls,meas.sim = meas.sim,retop = 2)
  #cat('\nodors',str(odors))
  #now, generate the data set with all rel cells removed and do the computation again
  dat.ls.rel <- campDataRemoveRel(dat.ls = dat.ls,rel = rel)
  contUD.rel <- campComputeUD(dat.ls = dat.ls.rel,meas.sim = meas.sim,op = 2,odors = odors)
  switch(op,contUD.rel,as.data.frame(contUD.rel))
}

#do all the reliabilities
#op: list of DFs, or 2 - df
campCompUDAllSingleRel <- function(dat.ls,meas.sim=c(0.15,0.45),odorop=1,op=1){
  notrials <- round(mean(campGetTrialDetails(alldata = dat.ls,op=2)[[1]]))
  res.ls <- lapply(1:notrials, function(i){
    res <- campCompUDSingleRel(dat.ls = dat.ls,meas.sim = meas.sim,rel = setdiff(1:notrials,i),op=2)
    res <- res[c('urel.sum','rel.sum'),]
  })
  switch(op,res.ls,convertListToDF(res.ls))
}


#this function calculates campComputeUDReliability for all reliabilty levels, and gives you
#an answer in terms of the trend across reliability levels
#odorop: determines whether you used the same odors for the comparisons (2) or determine the odors
#within the correlations for each data set (1) 
#op: 1- two lists, 2 - a data frame first 4 cointrol, next 4 reliabiltyu out results, 3 - comparison
#of all the results
campCompUDAllRel <- function(dat.ls,meas.sim=c(0.15,0.45),odorop=1,op=1){
  trials <- campGetTrialDetails(alldata = dat.ls,op=2)[[1]]
  notrials <- min(trials) #if all get equal no of trials, this sovers everything, else lowest #
  res.ls <- lapply(1:notrials, function(i){
    res <- campComputeUDReliability(dat.ls = dat.ls,meas.sim = meas.sim,rel = i,op=op,odorop = odorop)
  })
  names(res.ls) <- 1:notrials
  #cat('\nnames',names(res.ls))
  res.ls
}

#function takes in the return of getCampSigDataList and gives you the same data without 
#cells at that reliability level
#dat.ls: the data that we are considering, return of getCampDataSigList
#rel: the reliability level that we want to get the data for, can be a single rel level or multiple
campDataRemoveRel <- function(dat.ls,rel=1,op=1){
  #now take out all cells with the particular reliability
  rel.ls <- campFreqCellsAnalog(data = dat.ls)
  #get a list with the reliable rel cells for each odor
  cell.ls <- lapply(rel.ls[[2]], function(x){
    cells <- lapply(1:length(rel),function(i) {#get cells at every reliability level
      which(x==rel[i])
    })
    res <- unlist(cells) #make it into a vector
  })
  #now take these out of of dat.ls
  trial.odors <- names(dat.ls)
  res.ls <- lapply(1:length(trial.odors), function(i){
    cells <- cell.ls[[trial.odors[i]]]
    #go through each component vector of the list, and set the rel cells to 0
    odor.vals <- lapply(dat.ls[[i]], function(x){
      x[cells] <- 0
      x
    })
    odor.vals
  })
  names(res.ls) <- names(dat.ls)
  res.ls
}

#function that gets the correlation between a list of matrices
#dat.ls: the list of matrices or data frames
#sel=1, correlation between every vector column across the two matrices and 
#the number is their average, sel=2, correlation between selected vectors in each matrix,
#selected col specified in by 'col'
#celltype: 1 - all sig. cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - get correlation for the prescribed cells, 6 - cells with a specifed reliability
#where reliability is in params, Could be one number or multiple numbers
#cells: for option 5, they are in the form of a list which gives the cells for each odor
#matop: correlation based on 1- averages of all correlations, 2 - based on average response rate comparison
#params: for celltype=5, params(topnstart,size) topnstart gives the start of the percentile range, and size gives the size, e.g., (50,5) is top 50 to 45
#op: 1 - pearson, 2 - cosine
#respop: the average response: 1 - average of all the trials it responds in, 2 - average based on all trial responses, even those where it does not respond
computeListMatCorr <- function(alldata,sel=1,col=1,cells=c(),celltype=4,matop=1,params=c(10,10),respop=1,op=1){
  dat.ls <- campGenRUCellMat(alldata,op=celltype,matop = matop,params = params,respop = respop) #select the kind of cells we want to get correlations for
  #cat('\n',str(dat.ls))
  if(length(dat.ls)<2) return(NULL) #if there is only one or no odors, there is nothing to compare it with
  notrials <- campGetTrialDetails(alldata = alldata,op=2)[[1]] #get no of trials for each odor
  #iterate through odor pairs, get the matrix Unions for the celltypes, and compute
  res.ls <- sapply(names(dat.ls),function(i){
    res.in <- sapply(names(dat.ls),function(j){
      if(nrow(dat.ls[[i]])==0 || nrow(dat.ls[[j]])==0 ) cor.mat <- 0
      else{
        #cat('\n',i,j,':',unlist(dat.ls[[i]][,1]),'\n',dat.ls[[j]][,1],'\n' )
        matunion <- getMatUnion(dat.ls[[i]],dat.ls[[j]])
        cor.mat <- computeCor(matunion[[1]],matunion[[2]],matop = matop,op=op)
      }
      cor.mat
    })
  })
  res.ls
}

#computes the correlation between two matrices mat1 and mat2
#op: 1 - cor, 2 - linear fit, 3- the two vectors for plotting as alist
computeMatCor <- function(mat1,mat2,op=1){
  #do it row or column wise, actually do both and check if there is a difference. If there is, 
  #switch to using the average.
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  res <- switch(op,cor(vec1,vec2),fitlinearlaw(vec1,vec2),list(vec1,vec2))
  res
}

#function to compute the correlation between vectors that make up a list
computeListCorr <- function(dat.ls,op=1){
  novecs <- length(dat.ls)
  #algo: do the correlation between all (i,j) pairs and then name them at the end
  res <- sapply(1:novecs,function(i){
    res.in <- sapply(1:novecs,function(j){
      cor(dat.ls[[i]],dat.ls[[j]])
    })
  })
  row.names(res) <- names(dat.ls)
  colnames(res) <- names(dat.ls)
  res
}

#computes the correlation 
#matop: specifies the correlation type if it is a matrix, 
#op: 1 - Pearson, 2 - cosine
computeCor <- function(mat1,mat2,matop=1,op=1){
  if(length(dim(mat1)) == 2){#matrices, so do it here
    if(matop==1) {#average of all correlations
      # matvecs1 <- elimZeroCols(mat1[,1:ncol(mat1)])
      # matvecs2 <- elimZeroCols(mat2[,1:ncol(mat2)])
      matvecs1 <- elimZeroCols(mat1)
      matvecs2 <- elimZeroCols(mat2)
      cor.mat <- mean( switch(op,cor(matvecs1,matvecs2),cosine(matvecs1,matvecs2)) )
    }
    if(matop==2){#average all the matrix columsn and then correlate, mat[,1] already contains the average
      mat1vec <- mat1[,1]#cleanNAVec(sapply(1:nrow(mat1),function(i) sum(mat1[i,])/countNonZeroes(mat1[i,]) ) )
      mat2vec <- mat2[,1]#cleanNAVec(sapply(1:nrow(mat2),function(i) sum(mat2[i,])/countNonZeroes(mat2[i,]) ) )
      #cat('\n1:',sort(mat1vec,decreasing = T),'\n2:',sort(mat2vec,decreasing = T),':length:',length(mat1vec))
      cor.mat <- switch(op,cor(mat1vec,mat2vec),cosine(mat1vec,mat2vec) )
    }
    return(cor.mat)
  }
  switch(op,cor(mat1,mat2),cosine(mat1,mat2))
}

#checks if any of the matrix columns are all zeroes, if they are just skips them
#mat: matrix to be checked
#val: default 0 that has to be checked
elimZeroCols <- function(mat,val=0,op=1){
  mcols <- sapply(1:ncol(mat),function(i) length(which(mat[,i]==val)) == length(mat[,i]) )
  mcols <- !mcols
  mat[,mcols]
}

#takes a union of the two matrix rows, so that they have the same number of row
#returns a list of two matrices or DFs, where the first one has the missing rows of the second, though 
#all entries are 0, and vice-versa
#at1 and mat2: thw two matrices
getMatUnion <- function(mat1,mat2,op=1){
  #do the matrices in turn
  mat1.new <- fillMatMissingRows(mat1,mat2)
  mat2.new <- fillMatMissingRows(mat2,mat1)
  list(mat1.new,mat2.new)
}


#fills the first matrix, with the rows that are missing from the second matrix, i.e., the rows of matrix
#2 that are not in matrix 1.
#mat1: the matrix to change, mat2: the matrix to scour for the missing rows
fillMatMissingRows <- function(mat1,mat2,op=1){
  #do the matrices in turn
  row1 <- rownames(mat1)
  row2 <- rownames(mat2)
  #cat('\nfill',row1,':',row2,'\n')
  r21 <- setdiff(row2,row1) #the rows of 2 that are not in 1
  vec1 <- as.numeric(row1)
  mat1.new <- mat1
  emptyrow <- rep(0,ncol(mat1))
  for(i in as.numeric(r21)){
    #see where each row should be inserted
    pos <- findPosVal(i,vec1 )
    #use vec1 to update the rownames
    vec1 <- insertValPos(i,vec1,pos)
    #cat('\n',i,str(mat1.new),pos,emptyrow)
    #print(mat1.new)
    mat1.new <- insertRowDf(mat1.new,newrow = emptyrow,posn = pos)
    rownames(mat1.new) <- vec1
  }
  mat1.new #return this new matrix
}

#gets the statistics in euclidean distance terms the way glenn had requested
#dat.ls: return of getCampDataSigList
#op: 1- filter with sig cells, 2 - do not filter, 
#3 - compare trial1-2, then 2-3 so on with filtering
#rettype: 1 - list, 2 - DF
campGetTrialDistances <- function(dat.ls,rettype=1,op=1){
  sig.dat <- getGroupedListsDf(dat.ls,op=1) #sig. cells
  resp.dat <- getGroupedListsDf(dat.ls,op=2) #responses
  res.ls <- lapply(names(sig.dat), function(x){#filter one vector by another
    sig <- sig.dat[[x]]
    resp <- resp.dat[[x]]
    res <- switch(2-(op %% 2),filterWithMat(resp,sig),resp) #choose sig cell or resp
    #cat(str(sig),str(resp),str(res))
    #now get the euclidean distance of the odor vector for every trial
    res.vec <- switch(floor(op/2)+1,
                      sapply(1:ncol(res), function(y) lpnorm(res[,y],p=2)), #distance from origin
                      sapply(2:ncol(res), function(y) lpnorm(res[,y]-res[,y-1],p=2))  #distance between t1-t2, ...
                     )
  })
  if(rettype==1){#list return
    names(res.ls) <- names(sig.dat)
    res <- res.ls
  }
  if(rettype==2){#DF return
    res <- transposeDF(convertNestedListsDF(res.ls));
    colnames(res) <- 1:ncol(res);
    row.names(res) <- names(sig.dat);
  }
  res
}

#given the campDataSiglist call, will return data for comparin no of trial responses and average response
#dat.ls: return of getCampDataSigList
campTrialNoVsResp <-function(dat.ls,op=1){
  hab.dat <- campHabituation(dat.ls,op=4) #get all the responses for each cell across odors
  hab.df <- joinUnevenDFs(hab.dat) #join them into a single DF or matrix
  #cat(str(hab.df))
  res.freq <- getRowsDFNonzerosAll(hab.df,op=4) #get frequency responses for each possibility   
  #get the avg. response for each specific number of trials, and report that with trial number
  #cat(str(res.freq))
  res.ls <- lapply(1:length(res.freq), function(x){
    if (length(res.freq[[x]]) != 0){#get all cells and their with exactly x responses
      no <- nrow(res.freq[[x]])
      #cat(x)
      vals <- sapply(1:no, function(y) mean(res.freq[[x]][y,]))
      tmp <- cbind.data.frame(rep(x,no),vals)
      rownames(tmp) <- rownames(res.freq[[x]]) 
    }
    else tmp <- cbind.data.frame(x,0)
    tmp
  })
  res <- switch(op,res.ls,joinUnevenDFs(res.ls))
  res
}

#result: The explanation to Saket
#if the sign response of an unreliable cell is in an eaerlier trial, it is more likely to respond a second time
#if the onoe signigifcant response is in the first 6 vs 3 trials, there are likely to be fewer cells that respond a second time

#results: 
#The number of unreliable cells goes up with the number of 
#trials. The correlation is around 0.52. The unerliable and reliable sum 
#however is not very correlated. It is about -0.02. I think this is because 
#when there are more unreliable cells, the overlap between cells increases, 
#so the per cell usefulness goes down, but since there are more cells, you 
#still get similar amounts of useful discrimination. Got to make sure of 
#this by measuring amount of overlap.


#function that takes the output of the campAllDirDiscScores and getsthe output in one of various forms
#specified by op
#dfno: the dfno to process and the processing to be done is specified by op
#norm:0 - do not normalize, 1 or higher - normalize before joining the dfs; 1 - normalize by mean
#op=1, get the output of one of the five data frames as an aggregated data frame, 
#2: get the usefdul D list for each class of odor-pairs, dissimilar, midsimilar and similar
#3: get the slope change for each class of odor-pairs, dissimilar, midsimilar and similar
#4, get the averages of each experimental set or all together
campProcessDirSets <- function(dat.ls,norm=0,dfno=1,dissimilar=0.15,similar=0.5,op=1){
  data.ls <- lapply(dat.ls, '[[',dfno)
  data.ls <- campNormalizeOdorSets(data.ls,op=norm) #normalize
  #just join everything so we can plot it
  #old agg.df <- convertNestedListsDF(data.ls)
  agg.df <- joinListDFs(data.ls)
  #have to get the names right too.
  names.df <- lapply(names(data.ls), function(x) paste(x,rownames(data.ls[[x]]),sep = '/'))
  rownames(agg.df) <- unlist(names.df) #set the row and col names
  names(agg.df) <- names(data.ls[[1]])
  if(op==1) res <- cleanNAMat(agg.df) #get rid of NA
  if(op==2){#get the useful D for each odor-pair class
    tmp6 <- agg.df[which(agg.df[,1]< dissimilar),]
    dis <- list(tmp6[,4],tmp6[,7],tmp6[,5],tmp6[,8],tmp6[,'ur.cells'],tmp6[,'r.cells'],tmp6[,'trials'])
    tmp6 <- agg.df[which(agg.df[,1]> similar),]
    sim <- list(tmp6[,4],tmp6[,7],tmp6[,5],tmp6[,8],tmp6[,'ur.cells'],tmp6[,'r.cells'],tmp6[,'trials'])
    tmp6 <- agg.df[which(agg.df[,1] >= dissimilar & agg.df[,1] <= similar),]
    mid <- list(tmp6[,4],tmp6[,7],tmp6[,5],tmp6[,8],tmp6[,'ur.cells'],tmp6[,'r.cells'],tmp6[,'trials'])
    res <- list(dis,mid,sim)
  }
  if(op==3){#get the change in slope for each odor-pair class
    tmp61 <- agg.df[which(agg.df[,1]< dissimilar),]
    dis <- abs(tmp61[,'slope.nor']-tmp61[,'slope.sat'])
    tmp62 <- agg.df[which(agg.df[,1]> similar),]
    sim <- abs(tmp62[,'slope.nor']-tmp62[,'slope.sat'])
    tmp63 <- agg.df[which(agg.df[,1] >= dissimilar & agg.df[,1] <= similar),]
    mid <- abs(tmp63[,'slope.nor']-tmp63[,'slope.sat'])
    res <- list(list(dis,mid,sim),list(tmp61[,'slope.nor'],tmp61[,'slope.sat'],
                                       tmp62[,'slope.nor'],tmp62[,'slope.sat'],tmp63[,'slope.nor'],tmp63[,'slope.sat']))
  }
  res
  
  
}

#normalize each of the dfs in the list 
#op: specifies the normalization, 0 - do not, 1 - normalize by mean
campNormalizeOdorSets <- function(dat.ls,op=1){
  norm <- op
  if(norm>0){
    if(norm == 1){
      #cols to normalize 
      datnorm.ls <- lapply(dat.ls, function(x){
        colavg <- x[,c('rob','ursum.nor','rsum.nor','ursum.sat','rsum.sat')]
        colvag.mn <- apply(colavg[,-1], 2, mean)
      }) 
    }
    #now generate the normalized DFs
    normdfs <- lapply(names(dat.ls), function(x){
      #noramlize the first column and then choose the norm for the other columns to be 
      #the real ratio between them
      tmp <- dat.ls[[x]] #get the df
      tmp[,4] <- tmp[,4]/datnorm.ls[[x]][1]
      tmp[,-c(1:4,6)] <- tmp[,-c(1:4,6)]/datnorm.ls[[x]][-1] *(datnorm.ls[[x]][-1]/datnorm.ls[[x]][1])
      tmp
    })
  }
  else normdfs <- dat.ls
  normdfs
}


#uses a cutoff as a paarameter to calculate normal verssus fine grained discrimination. How many odor pairs 
#satsfy either requirement
#dat.df: the return of campProcessDirSets with option 1. Basically the return of of campSatnorm and all the experimental
#sets joined together.
#cutoff: the cutoff
#op:1 - get the number of odor pairs that cross the cutoff for normal,saturated, and for both w/o ur
#op=2: for a given cutoff get the fraction in each of the three categories of odor similarity classes
campCalcCutoff <- function(dat.df,cutoff,op=1){
  if(op==1){
    #get the cells that are above the cutoff for normal
    posns.nor <- which((dat.df[,'ursum.nor']+dat.df[,'rsum.nor'])>cutoff)
    posnsr.nor <- which((dat.df[,'rsum.nor'])>cutoff) #purely r cells alone
    #get the cells that are above the cutoff for saturation
    posns.sat <- which((dat.df[,'ursum.sat']+dat.df[,'rsum.sat'])>cutoff)
    posnsr.sat <- which((dat.df[,'rsum.sat'])>cutoff) #purely r cells alone
    res <- c(length(posns.nor),length(posns.sat),length(posnsr.nor),length(posnsr.sat))
    #cat('cCC:',res,'\t:')
  }
  if(op==2){
    dis <- which(dat.df[,'rob'] < 0.15)
    mid <- which(dat.df[,'rob'] > 0.15 & dat.df[,'rob'] < 0.5)
    sim <- which(dat.df[,'rob'] > 0.5)
    #cat('\ncCC',(dis))
    dis.no <- length(which((dat.df[dis,'ursum.nor']+dat.df[dis,'rsum.nor'])>cutoff) )/length(dis)
    mid.no <- length(which((dat.df[mid,'ursum.nor']+dat.df[mid,'rsum.nor'])>cutoff) )/length(mid)
    sim.no <- length(which((dat.df[sim,'ursum.nor']+dat.df[sim,'rsum.nor'])>cutoff) )/length(sim)
    res <- c(dis.no,mid.no,sim.no)
    dis.no <- length(which((dat.df[dis,'ursum.sat']+dat.df[dis,'rsum.sat'])>cutoff) )/length(dis)
    mid.no <- length(which((dat.df[mid,'ursum.sat']+dat.df[mid,'rsum.sat'])>cutoff) )/length(mid)
    sim.no <- length(which((dat.df[sim,'ursum.sat']+dat.df[sim,'rsum.sat'])>cutoff) )/length(sim)
    res <- c(res,dis.no,mid.no,sim.no)
  }
  res
}

#calculates the cutoff for a range of cutoff parameters
#dat.df: the return of campProcessDirSets with option 1. Basically the return of of campSatnorm and all the experimental
#sets joined together.
#cutop:1 - get the number of odor pairs that cross the cutoff for normal,saturated, and for both w/o ur
#cutop=2: for a given cutoff get the fraction in each of the three categories of odor similarity classes
#op=1: get the cutoof for ptrrange equally spread between highest and lowest useful D parameters
#op=2, get the cutoff from the cutrange specificatoins
#op=3, get the cutoff from cut off vector range
campCalcCutoffRange <- function(dat.df,ptsrange=10,cutrange=seq(0,1,0.5),cutoffvec=c(),cutop=1,op=1){
  d.vec <- dat.df[,'ursum.nor'] + dat.df[,'rsum.nor']
  #determine the number of cutoffs based on op
  switch(op,pts <-cutrange,pts <- seq(min(d.vec),max(d.vec),(max(d.vec)-min(d.vec))/ptsrange),
         pts <- cutoffvec)
  res <- sapply(pts, function(x){
    tmp <- campCalcCutoff(dat.df = dat.df,cutoff = x,op = cutop)
    #cat(x,'. ',tmp,'\n')
    tmp <- c(x,tmp)
  })
  res.df <- as.data.frame(t(res))
  res.df
}

#function that calculates the odor habituation responses over a number of trials
#op=1, for every odor, for each cell across all trials, get the trial where it hits its peak
#op=2, now, across odor see if its always the first trial that is the peak response
#op=3, across odors, gets the significant responses for each cell, if not all cells, padded out
#op-4, just get all the responses for each odor across trials as df
campHabituation <- function(dat.ls,op=1){
  res.rel <- getGroupedListsDf(dat.ls,op = 1) #reliable cells
  res.sig <- getGroupedListsDf(dat.ls,op = 2) #the responses of the cells
  #use the first list to filter the second for all the odors
  res.lsts <- lapply(1:length(res.rel), function(x){
    res <- filterListWithList(vec1.lst = res.sig[[x]],vec2.lst = res.rel[[x]],thresh = 0)
  }) 
  #now, lets get the stats that we want
  #first trials vs no of cells that have the peak for this odor
  res.ls <- switch(op,
  lapply(res.lsts, function(x.df){
    res.df <- convertNestedListsDF(x.df)
    ranks <- getDfColPeak(dat.df = res.df,op=2)
  }),
  #go through each cell and figure out if the 
  lapply(res.lsts, function(x.df){
    res.df <- convertNestedListsDF(x.df)
    ranks <- getDfColDecPeak(dat.df = res.df,op=2)
  }),
  #go through each cell see how the first and next significant response are related 
  lapply(res.lsts, function(x.df){
    res.df <- convertNestedListsDF(x.df)
    vals.df <- transposeDF(getDfColNVals(dat.df = res.df,pad=0,op=3)) #gets the first two sig responses
    #vals.df <- convertNestedListsDF(vals)
  }),
  lapply(res.lsts, function(x.df){#just get the responses for each odor:op = 4
    res.df <- convertNestedListsDF(x.df)
  }) )
  res.ls
}



#process the different data sets of odors and their cell differences: the result of campOverlapAnalogAll 
#op=1, the earlier one, which analyses a list of DFs 
campProcessOdorSets <- function(lstdfs,op=1){
  odorsets <- joinUnevenDFs(lstdfs)
  #average of rob and hallem
  overlap <- apply(odorsets[,c(2,3)],1,mean) #get the overlap for these odor sets:mean of Rob and Hallem
  #gets the scores for each reliable population set
  res <- lapply(1:length(overlap), function(x){
    tmp <- odorsets[x,4:(3+odorsets[x,ncol(odorsets)])]
    #cat(tmp,rep(overlap[x],length(tmp)),'\n')
    cbind.data.frame(rep(overlap[x],length(tmp)),tmp)
  })  
  res
}


#this function gets the correlation between mean response levels and # of significant responses
#data: the return of getCampDataSigList
# op = 1, return corr, 2 - the data as a data frame for plotting, 3 - SD corr, 4 - SD data frame for plotting, 5 - SD data frame, that is scaled by the mean
getCampResponseSigCor <- function(data,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  res.data <- getGroupedListsDf(data,op=2) #gets the response DF
  #cat(str(res.sig),str(res.data))
  #cat(str(res.data),str(res.sig),' fin\n')
  res.ls <- lapply(1:length(res.sig), function(x){ #go throughh all the odors
    notrials <- length(res.sig[[x]]) # no of cols in the DF
    #reliable cells
    nosig.cells <- apply(res.sig[[x]], 1, sum)
    #mean response levels
    #got to get the mean of only those cells that are significant
    sigresp.cells <- res.data[[x]] * res.sig[[x]]
    #y is the index or cell no
    resp <- switch(op,sapply(1:getLength(nosig.cells),function(y) sum(sigresp.cells[y,])/nosig.cells[y]),
                   sapply(1:getLength(nosig.cells),function(y) sum(sigresp.cells[y,])/nosig.cells[y]),
                   sapply(1:getLength(nosig.cells),function(y) sd(sigresp.cells[y,][sigresp.cells[y,]>0])),
                   sapply(1:getLength(nosig.cells),function(y) sd(sigresp.cells[y,][sigresp.cells[y,]>0])),
                   sapply(1:getLength(nosig.cells),function(y) sd(sigresp.cells[y,][sigresp.cells[y,]>0])/mean(sigresp.cells[y,][sigresp.cells[y,]>0])),
                   sapply(1:getLength(nosig.cells),function(y) mean(sigresp.cells[y,][sigresp.cells[y,]>0])) )
    resp <- cleanNAVec(resp)
    #cat('\n',op,';',resp,'\nsig:',nosig.cells,'\n' )
    #print(res.data)
    #cat(cor(nosig.cells,resp))
    posns <- which(nosig.cells > 0)
    switch(op,cor(nosig.cells[posns],resp[posns]),cbind.data.frame(nosig.cells,resp),
           cor(nosig.cells[posns],resp[posns]),cbind.data.frame(nosig.cells,resp),
           cbind.data.frame(nosig.cells,resp),cbind.data.frame(nosig.cells,resp) )
  })
  names(res.ls) <- names(res.sig)
  switch(op,unlist(res.ls),res.ls,unlist(res.ls),res.ls,res.ls,res.ls)
}


#signal noise comparisons
#result: scalar number, 2 = vector of values
#op=1, cor sig vs sd, 2 cv, 3 - ff
campSigNoiseComp <-function(data,result=1,op=1){
  res.sig <- getCampResponseSigCor(data,op=2)
  res.sd <- getCampResponseSigCor(data,op=4)
  res <- lapply(1:length(res.sig),function(i){
    posns <- which(res.sig[[i]][,1]>1)
    res.cor <- switch(result,
                      switch(op,cor(res.sig[[i]][posns,2],res.sd[[i]][posns,2]),
                             mean(res.sd[[i]][posns,2]/res.sig[[i]][posns,2]),
                             mean(res.sd[[i]][posns,2]^2/res.sig[[i]][posns,2]) ),
                      switch(op,cbind(res.sig[[i]][posns,2],res.sd[[i]][posns,2]),
                             c(res.sd[[i]][posns,2]/res.sig[[i]][posns,2]),
                             c(res.sd[[i]][posns,2]^2/res.sig[[i]][posns,2]) ) )
  })  
  if(result==1) res <- unlist(res)
  #cat('\n',length(res),names(res.sig))
  names(res) <- names(res.sig)
  res
}

#signal noise comparisons
#result: scalar number, 2 = vector of values
#op=1, the df of all cell responses, 2 = separated by freq./reliability, 3 - fF
campSigNoiseFF <-function(data,result=1,op=1){
  #res.sig <- getCampResponseSigCor(data,op=2)
  #res.sd <- getCampResponseSigCor(data,op=4)
  res.scalesd <- getCampResponseSigCor(data,op=5)
  scalesd.df <- joinListDFs(res.scalesd)
  scalesd.df <- scalesd.df[which(scalesd.df[,1]>1),] #you can't have an SD for trials <= 1
  #get a vector of the different frequencies
  freq.levels <- unique(scalesd.df[,1])
  #cat('\nfreq',sort(freq.levels,decreasing = F))
  scalesd.lst <- lapply(sort(freq.levels,decreasing = F),function(i) scalesd.df[which(scalesd.df[,1]==i),2])
  scaleff.lst <- lapply(sort(freq.levels,decreasing = F),function(i) scalesd.df[which(scalesd.df[,1]==i),2]^2)
  names(scalesd.lst) <- sort(freq.levels,decreasing = F)
  names(scaleff.lst) <- sort(freq.levels,decreasing = F)
  
  switch(op,scalesd.df,scalesd.lst,scaleff.lst)
}



#this function gets the probability overlap between odors for reliable and non-reliable cells
#data: the return of getCampDataSigList
#compare = 1: computer the overlap, 2 - computer the difference in probability b/w the odors
computeCampOverlap <- function(data,compare=1,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  notrials <- length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  #have to calculate the freq. list for all the odors, which is a single list
  res.ls <- lapply(res.sig, function(x){
    #first the freq cell list
    freq.cells <- apply(x, 1, sum)
    #reliable and non-reliable cells. This is the place to code what majority of trials
    #majority of trials is > notrials/2
    nonrel.cells <- which(freq.cells <= notrials/2 & freq.cells >0 )
    rel.cells <- which(freq.cells > notrials/2)
    list(freq.cells,rel.cells,nonrel.cells)
  })
  notrials <- campGetTrialDetails()[[1]]#length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  freq.lst <- lapply(res.ls,'[[',1) #each list item contains a frequency vector, and vector names are the cell no
  rel.lst <- lapply(res.ls,'[[',2) #each list item contains the reliable cells for that odor
  nonrel.lst <- lapply(res.ls,'[[',3)
  #cat(str(freq.lst),str(rel.lst)) #print out the freq.lst descriptions
  #now call the overlap calculating function
  #cat(str(freq.lst),str(rel.lst))
  overlap.rel <- getCampOverlapScores(freq.lst,rel.lst,notrials = notrials)
  overlap.nonrel <- getCampOverlapScores(freq.lst,nonrel.lst,notrials = notrials)
  overlap.relnonrel <- getCampOverlapScores(freq.lst,rel.lst,nonrel.lst,notrials = notrials,op=2)
  res <- list(overlap.rel,overlap.nonrel,overlap.relnonrel)
  res
}

#this function gets the probability difference between odors for reliable and non-reliable cells
#data: the return of getCampDataSigList
#probno: whether the result shoudl be a probability or numbers
#op, 1 - computer the difference in probability b/w the odors for the whole population
computeCampDiff <- function(data,probno=1,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  notrials <- length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  #have to calculate the freq. list for all the odors, which is a single list
  res.ls <- lapply(res.sig, function(x){
    #first the freq cell list
    freq.cells <- apply(x, 1, sum)
    #reliable and non-reliable cells. This is the place to code what majority of trials
    #majority of trials is > notrials/2
    nonrel.cells <- which(freq.cells <= notrials/2 & freq.cells >0 )
    rel.cells <- which(freq.cells > notrials/2)
    list(freq.cells,rel.cells,nonrel.cells)
  })
  notrials <- length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  freq.lst <- lapply(res.ls,'[[',1) #each list item contains a frequency vector, and vector names are the cell no
  rel.lst <- lapply(res.ls,'[[',2) #each list item contains the reliable cells for that odor
  nonrel.lst <- lapply(res.ls,'[[',3)
  #cat(str(freq.lst),str(rel.lst)) #print out the freq.lst descriptions
  #now call the overlap calculating function
  #cat(str(freq.lst),str(rel.lst))
  #four types of differences
  dif.rel <- getCampDiffScores(freq = freq.lst,cell.lst = rel.lst,cell2.lst = rel.lst,notrials=notrials,probno=probno,op=1)
  dif.nonrel <- getCampDiffScores(freq = freq.lst,cell.lst = nonrel.lst,cell2.lst = nonrel.lst,notrials=notrials,probno=probno,op=1)
  dif.relnonrel <- getCampDiffScores(freq = freq.lst,cell.lst = rel.lst,cell2.lst = nonrel.lst,notrials=notrials,probno=probno,op=2)
  dif.sigsilent <- getCampDiffScores(freq = freq.lst,cell.lst = rel.lst,cell2.lst = nonrel.lst,notrials=notrials,probno=probno,op=3)
  res <- list(dif.rel,dif.nonrel,dif.relnonrel,dif.sigsilent)
  res
}

#will get the cells ranked by increasing reliability versus the differences in their probabilities for both odors
#data: the return of getCampDataSigList
#odorpair: the odor pair for which we are going to calculate this, e.g., odor 1 and 2
#op, 1 - computer the difference in probability b/w the odors
campGetReliableDiff <- function(data,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  notrials <- length(res.sig[[1]]) #no of trials, assumes that all odors have same # of Trials
  #have to calculate the freq. list for all the odors, which is a single list
  res.ls <- lapply(res.sig, function(x){
    
  })
  #now go through both odors and get the difference for those two cells positions.
}

#function that gets the stats for the functions campGetReliableDiff and computeCampOverlap
#
#dat.ls: the return of getCampDataSigList
#anal: the type of analysis required. 1 - compare overlap vs differences for the different populations
campOverlapDiffStats <- function(dat.ls,anal=1,op=1){
  #get overlap and difference in terms of probabilities
  overlap <- computeCampOverlap(data = dat.ls)
  dif <- computeCampDiff(data = dat.ls,probno = 1) 
  dif.no <- computeCampDiff(data = dat.ls,probno = 2) #get the numbers for different populations 
  res.compare <- sapply(1:length(overlap), function(x) {
    #print(diag(dif[[x]]))
    tmp <- overlap[[x]]
    tmp1 <- dif[[x]]
    cor(as.vector(tmp),as.vector(tmp1))
  })
  dif.avg <- sapply(1:length(dif), function(x) {
    noodors <- nrow(dif[[x]])
    sum(dif[[x]])/(noodors^2 - noodors)
  })
  difno.avg <- sapply(1:length(dif.no), function(x) {
    noodors <- nrow(dif.no[[x]])
    #cat(sum(dif.no[[x]]),(noodors^2 - noodors),'\n')
    sum(dif.no[[x]])/(noodors^2 - noodors)
  })
  overlap.avg <- sapply(1:length(overlap), function(x) {
    noodors <- nrow(overlap[[x]])
    sum(overlap[[x]])/(noodors^2)
  })
  #make pretty results and output by setting names
  names(res.compare) <- c('reliable','unreliable','mixed')
  names(overlap.avg) <- c('reliable','unreliable','mixed')
  names(dif.avg) <- c('reliable','unreliable','mixed','silent')
  names(difno.avg) <- c('reliable','unreliable','mixed','silent')
  list(res.compare,overlap.avg,dif.avg,difno.avg)
}


#put all three CampMatrix measures together, and if needed write it to file
#...:the arguments are the same as getCampDataSigList, its just passed on as it is
getCampMatScores <- function(...,foldname='.',op=1){
  params <- list(...)
  dat.ls <- do.call(getCampDataSigList,params)
  res1 <- getCampDataRel(dat.ls)
  res2 <- getCampResponseSigCor(dat.ls)
  res3 <- computeCampOverlap(dat.ls)
  cat(res1,'\n')
  cat(res2,'\n')
  cat(mean(res3[[1]]),mean(res3[[2]]),'\n')
  print(res3[[2]])
  outputfile <- paste(foldname,'/tseries_nocells.csv',sep = '')
  cat(foldname,'\n')
  write.table(res1,file = outputfile,col.names = T)
  outputfile <- paste(foldname,'/tseries_correlation.csv',sep = '')
  write.table(res2,file = outputfile,col.names = T)
  outputfile <- paste(foldname,'/tseries_reliable_overlap.csv',sep = '')
  write.table(res3[[1]],file = outputfile,col.names = T,row.names = T)
  outputfile <- paste(foldname,'/tseries_nonreliable_overlap.csv',sep = '')
  write.table(res3[[2]],file = outputfile,col.names = T,row.names = T)
  outputfile <- paste(foldname,'/tseries_nonreliable_and_reliable_overlap.csv',sep = '')
  write.table(res3[[3]],file = outputfile,col.names = T,row.names = T)
  
  list(res1,res2,res3)
}

#now, take this function, and also calculate the z-score and max value, and analyze the results on that basis
#results: # of significant responses per tiral, # of reliable cells, correlation between #reliable responses vs signal levels
#and the probability overlap/cloud matrix for reliable vs non reliable cells, also # of reliable and nonreliable cells
#params: no of relibale trials, z-score or mean

#todo :  
#adjust the z-score and mean for those cells that have -ve nos. Best way is to adjust the mean to be the threshold.

#function to traverse the directoires

#check if there is odor habituation for the unreliable cells.


#built on top of DriTraverse, checks for a condition. if fulfilled does an operation in that directory. Should be run from root level of the data dir
#if it is a file search, does the operation in that directory. 
#if it is a directory search, does the opertation inside the directory
#this function traverses through the current and all sub-directories recursively. In two modes
#mode: 1- the search pattern is applied to the directory names and only lists out the directories that match it
#2: the search pattern applies to files only and lists out all the files that match the pattern
#3: separate search patterns for directories and files
#searchstr: if mode is 3, this is a 2 elements list (dirsearch,filesearch)
#targetfn: name of the target function to call
#targetpars: the parameters for this function
#op: 1 - list out from this directory, 2 - list full pathnames
CampDirTraverse <- function(foldname='.',searchstr=c('.*','.*'),level=0,mode=1,targetfn=c(),targetpars=list(),op=1){
  #cat('In dir ',foldname,'with mode ',mode,'and search',searchstr,'\n')
  dirsearch <- switch(mode,searchstr[1],'.*',searchstr[1]) #need the swtich for mode 2
  #searches for a string 4 characters long that are not alphabets or digits 
  filesearch <- switch(mode,'^[^\\d\\w]{4}$',searchstr[2],searchstr[2]) 
  
  #list all the files in this dir that match the search string 
  files <- list.files(path=foldname,pattern = '*')
  filetypes <- file.info(paste(foldname,'/',files,sep = ''))["isdir"] # get the file types
  if (dim(filetypes)[1]==0) return(c('NA')) #if no files, return NA
  
  #process dirs and files, and apply filters separately
  dirs <- files[which(filetypes == T)] #test if dir
  filesonly <- files[which(filetypes != T)]
  dirfilter <- grepl(pattern = dirsearch,x = dirs,ignore.case = T)
  filefilter <- grepl(pattern = filesearch,x = filesonly,ignore.case = T)
  filespath <- switch(op,paste(foldname,filesonly,sep = '/'),filesonly)
  filespath <- filespath[filefilter] #get rid of empty filenames
  #cat('In',foldname,length(filespath),dirs,filespath,'\n')
  if (length(filespath)==0) cat('\nNo files that match the pattern in ',foldname)
  else {#setworking directory, execute and then set it back
    currentdir <- getwd()
    setwd(foldname)
    cat('In foldname ',foldname,' and current dir ',currentdir,getwd(),'\n')
    getCampMatScores()
    setwd(currentdir)
  }
  #else do.call(targetfn,list(foldname=foldname))
  
  if (length(dirs) > 0) { #if there are directories
    dirs <- paste(foldname,'/',dirs,sep = '') #add the directory name
    #cat('\nres',str(dirs))
    res <- sapply(dirs, CampDirTraverse,searchstr,level+1,mode,targetfn,targetpars,op) #get files in each dir recursively
    res1 <- unlist(res) #make the lists into vectors
    names(res1) <- NULL #dont need the names
    #the results are now dirs and files from this directory as well as its sub-directories
    res <- switch(mode,dirs[dirfilter],filespath,c(dirs[dirfilter],filespath)) 
    return(c(res,res1))#return (c(res1,filesonly[filefilter]))
  }
  #do we need to return this stuff?
  switch(mode,dirs[dirfilter],filespath,c(dirs[dirfilter],filespath)) 

}


#built on top of DriTraverse, checks for a condition. if fulfilled does an operation in that directory. Should be run from root level of the data dir
#this particular directory traverse operation gets useful discrimination data from all directories
#if it is a file search, does the operation in that directory. 
#if it is a directory search, does the opertation inside the directory
#this function traverses through the current and all sub-directories recursively. In two modes
#mode: 1- the search pattern is applied to the directory names and only lists out the directories that match it
#2: the search pattern applies to files only and lists out all the files that match the pattern
#3: separate search patterns for directories and files
#searchstr: if mode is 3, this is a 2 elements list (dirsearch,filesearch)
#targetfn: name of the target function to call
#targetpars: the parameters for this function
#op: 1 - list out from this directory, 2 - list full pathnames
campAllDirDiscScores <- function(foldname='.',searchstr=c('.*','.*'),level=0,mode=1,targetfn=c(),targetpars=list(),op=1){
  res <- NULL #initialize the return
  #cat('In dir ',foldname,'with mode ',mode,'and search',searchstr,'\n')
  dirsearch <- switch(mode,searchstr[1],'.*',searchstr[1]) #need the swtich for mode 2
  #searches for a string 4 characters long that are not alphabets or digits 
  filesearch <- switch(mode,'^[^\\d\\w]{4}$',searchstr[2],searchstr[2]) 
  
  #list all the files in this dir that match the search string 
  files <- list.files(path=foldname,pattern = '*')
  filetypes <- file.info(paste(foldname,'/',files,sep = ''))["isdir"] # get the file types
  if (dim(filetypes)[1]==0) return(c()) #if no files, return NA
  
  #process dirs and files, and apply filters separately
  dirs <- files[which(filetypes == T)] #test if dir
  filesonly <- files[which(filetypes != T)]
  dirfilter <- grepl(pattern = dirsearch,x = dirs,ignore.case = T)
  filefilter <- grepl(pattern = filesearch,x = filesonly,ignore.case = T)
  filespath <- switch(op,paste(foldname,filesonly,sep = '/'),filesonly)
  filespath <- filespath[filefilter] #get rid of empty filenames
  #cat('In',foldname,length(filespath),dirs,filespath,'\n')
  if (length(filespath)==0) cat('\nNo files that match the pattern in ',foldname)
  else {#setworking directory, execute and then set it back
    currentdir <- getwd()
    setwd(foldname)
    cat('\nIn foldname ',foldname,' and current dir ',currentdir,getwd(),names(targetpars),'\n')
    #print(campGetTrialDetails())
    #res <- c(campGetTrialDetails(),'elem')
    #res <- c(do.call(campCompareSatNorm,targetpars),'elem')
    res <- c(do.call(targetfn,targetpars),'elem')
    #cat(str(res))
    setwd(currentdir)
  }
  #else do.call(targetfn,list(foldname=foldname))
  
  if (length(dirs) > 0) { #if there are directories
    res.empty <- c()
    cat('\n dirs are ',dirs)
    #condition 1: the subdirectory is empty and the return is empty
    if(identical(dirs,'NA') ) return(NULL)# nothing so just return null
    dirs <- paste(foldname,'/',dirs,sep = '') #add the directory name
    res <- lapply(dirs, campAllDirDiscScores,searchstr,level+1,mode,targetfn,targetpars,op) #get files in each dir recursively
    names(res) <- dirs #preserive names
    #condition 2: valid and invalid subdirectories together. Take out the invalid ones(ones with 0 elements) 
    cat('\nfoldname: ',foldname)
    res <- removeListEmpty(res)
    #condition 3: if there are elements in this directory, and sub-directories with multiple elements
    #bring them all to the same level
    if(length(res) > 0){#multiple elements, check if they are elemental
      #cat('\nchecking elem list in ',foldname,'\n',str(res))
      elem.lst <- sapply(res, function(x) campIsElem(x))
      cat('\nelem lst is ',elem.lst,' ELEMENTS')
      nonelem.lst <- which(elem.lst==F) #the non elements
      tmp.lst <- res[which(elem.lst==T)]
      if(length(nonelem.lst)>0) for (i in nonelem.lst) {
        tmp.lst <- c(tmp.lst,res[[i]])
      }
      #cat('\nnonelem.lst:',nonelem.lst,str(tmp.lst))
      res <- tmp.lst
    }
    #check here if any of these are empty, if they are delete them from the list
    cat('\nin ',foldname,' looking at ',dirs,' with dir no ',length(dirs),' and dir returning res length',length(res))
    return(c(res))#return (c(res1,filesonly[filefilter]))
    #str(c(tmp1[c(1,2)],tmp1[[3]][c(1,2)]))
  }
  #do we need to return this stuff?
  switch(mode,dirs[dirfilter],filespath,c(dirs[dirfilter],filespath)) 
  cat('\nElem returning',length(res),'\n')
  res
}


#generic improvement on the previous function
#built on top of DriTraverse, checks for a condition. if fulfilled does an operation in that directory. Should be run from root level of the data dir
#this particular directory traverse operation gets useful discrimination data from all directories
#if it is a file search, does the operation in that directory. 
#if it is a directory search, does the opertation inside the directory
#this function traverses through the current and all sub-directories recursively. In two modes
#mode: 1- the search pattern is applied to the directory names and only lists out the directories that match it
#2: the search pattern applies to files only and lists out all the files that match the pattern
#3: separate search patterns for directories and files
#searchstr: if mode is 3, this is a 2 elements list (dirsearch,filesearch)
#targetfn: name of the target function to call
#targetpars: the parameters for this function
#op: 1 - list out from this directory, 2 - list full pathnames
campAllDirDiscScoresGen <- function(foldname='.',searchstr=c('.*','.*'),level=0,mode=1,targetfn=c(),targetpars=list(),op=1){
  res <- NULL #initialize the return
  #cat('In dir ',foldname,'with mode ',mode,'and search',searchstr,'\n')
  dirsearch <- switch(mode,searchstr[1],'.*',searchstr[1]) #need the swtich for mode 2
  #searches for a string 4 characters long that are not alphabets or digits 
  filesearch <- switch(mode,'^[^\\d\\w]{4}$',searchstr[2],searchstr[2]) 
  
  #list all the files in this dir that match the search string 
  files <- list.files(path=foldname,pattern = '*')
  filetypes <- file.info(paste(foldname,'/',files,sep = ''))["isdir"] # get the file types
  if (dim(filetypes)[1]==0) return(c()) #if no files, return NA
  
  #process dirs and files, and apply filters separately
  dirs <- files[which(filetypes == T)] #test if dir
  filesonly <- files[which(filetypes != T)]
  dirfilter <- grepl(pattern = dirsearch,x = dirs,ignore.case = T)
  filefilter <- grepl(pattern = filesearch,x = filesonly,ignore.case = T)
  filespath <- switch(op,paste(foldname,filesonly,sep = '/'),filesonly)
  filespath <- filespath[filefilter] #get rid of empty filenames
  #cat('In',foldname,length(filespath),dirs,filespath,'\n')
  if (length(filespath)==0) cat('\nNo files that match the pattern in ',foldname)
  else {#setworking directory, execute and then set it back
    currentdir <- getwd()
    setwd(foldname)
    cat('\nIn foldname ',foldname,' and current dir ',currentdir,getwd(),names(targetpars),'\n')
    #print(campGetTrialDetails())
    #res <- c(campGetTrialDetails(),'elem')
    res <- c(do.call(targetfn,targetpars),'elem')
    #cat(str(res))
    setwd(currentdir)
  }
  #else do.call(targetfn,list(foldname=foldname))
  
  if (length(dirs) > 0) { #if there are directories
    res.empty <- c()
    cat('\n dirs are ',dirs)
    #condition 1: the subdirectory is empty and the return is empty
    if(identical(dirs,'NA') ) return(NULL)# nothing so just return null
    dirs <- paste(foldname,'/',dirs,sep = '') #add the directory name
    res <- lapply(dirs, campAllDirDiscScoresGen,searchstr,level+1,mode,targetfn,targetpars,op) #get files in each dir recursively
    names(res) <- dirs #preserive names
    #condition 2: valid and invalid subdirectories together. Take out the invalid ones(ones with 0 elements) 
    cat('\nfoldname: ',foldname)
    res <- removeListEmpty(res)
    #if(length(res)==1) {#if there is only one element, then we again need to bring the element up a level
    #  res <- res[[1]] 
    #  #cat('\nreturning one element ',length(res))
    #  return(res)
    #}
    #condition 3: if there are elements in this directory, and sub-directories with multiple elements
    #bring them all to the same level
    if(length(res) > 0){#multiple elements, check if they are elemental
      #cat('\nchecking elem list in ',foldname,'\n',str(res))
      elem.lst <- sapply(res, function(x) campIsElem(x))
      cat('\nelem lst is ',elem.lst)
      nonelem.lst <- which(elem.lst==F) #the non elements
      tmp.lst <- res[which(elem.lst==T)]
      if(length(nonelem.lst)>0) for (i in nonelem.lst) {
        tmp.lst <- c(tmp.lst,res[[i]])
      }
      #cat('\nnonelem.lst:',nonelem.lst,str(tmp.lst))
      res <- tmp.lst
    }
    #check here if any of these are empty, if they are delete them from the list
    cat('\nin ',foldname,' looking at ',dirs,' with dir no ',length(dirs),' and dir returning res length',length(res))
    return(c(res))#return (c(res1,filesonly[filefilter]))
    #str(c(tmp1[c(1,2)],tmp1[[3]][c(1,2)]))
  }
  #do we need to return this stuff?
  switch(mode,dirs[dirfilter],filespath,c(dirs[dirfilter],filespath)) 
  cat('\nElem returning',length(res),'\n')
  res
}



#given a list, checks if this an elemental element of the list. Basic check is that 
#the word 'elem' is contained within it. A bit kludgy. Fix later. Liam
campIsElem <- function(elem,op=1){
  #if(is.list(elem)) return(F) #if it is a list, it ain't elemental
  cat('\ncampIsElem')
  if(length(elem)>0) {
    res <- sapply(elem, function(x){
      identical(x,'elem')
    })
    #if(length(elem)==1) res <- identical(elem,)
    if (length(which(res==T))>0) return(T)
    else return(F)
  }
  F
}


#this function gets the epth of thee lements of the list where some of the elements are themselves nested lists
getListDepth <- function(lst,level=1,op=1){
  if( is.list(lst) ) res <- sapply(lst, function(x){
    flattenList(x,level = level+1) + 1
  })  
  else return(1) #bottom node, so return 1 and then count upwards 
  cat('level\t',level,';',res)
  if(level > 1) max(res)
  else res-1
}

#thi function gets all the scores for the overlap files
summarizeOverlapScores <- function(dir='.',searchstr=c('*tseries*','*reliabl*'),op=1){
  tmp <- DirTraverse(searchstr = searchstr,mode = 2) #gets all the overlap files
  #now separate them into lists with each list containing a directories overlap files
  res <- getFilesSameDir(tmp)
  #now,get the results
  res.lst <- lapply(res, function(x){
    #get the nonreliable files
    tmp <- unlist(gregexpr('nonreliable',x))
    #get the ssum of the three matrices
    #cat('\n',x)
    matsums <- sapply(x, function(y){
      #cat('\n',y)
      matsum <- sum(read.table(y))
    })
    nonrel <- which(tmp>=1) #ones for which search resturned a postive no
    rel <- which(tmp<1)
    c(matsums[rel],sum(matsums[nonrel]))
  })
  #converst this to a dF
  res.df <- transposeDF(convertNestedListsDF(res.lst))
  row.names(res.df) <- names(res.lst)
  colnames(res.df) <- c('reliable','nonreliable')
  res.df
}

#thi function gets the  scores for the number of core,non-core cells
summarizeCellScores <- function(dir='.',searchstr=c('*tseries*','*nocells*'),op=1){
  tmp <- DirTraverse(searchstr = searchstr,mode = 2) #gets all the overlap files
  #now separate them into lists with each list containing a directories overlap files
  res <- getFilesSameDir(tmp)
  #now,get the results
  res.lst <- lapply(res, function(x){
    dat <- read.table(x)
    unlist(dat)
  })
  #converst this to a dF
  res.df <- transposeDF(convertNestedListsDF(res.lst))
  row.names(res.df) <- names(res.lst)
  colnames(res.df) <- c('significant','reliable','nonreliable')
  res.df
}


#given the results of the CalcCampDataSigList, this will group by odors, and calculate the major significant cells and so forth.
#analysis: 1 - get the response values, 2- get the rank values, 
#3, get the number of significant cells per trial per odor, 4 - get the number of common cells across trials 
#op=1
analyzeCampSigDataOdors <- function(data,analysis=1,op=1){
  #first group by the names
  odor.grp <- split(data,names(data))
  #cat(str(odor.grp))
  #now do the analysis for each individual group
  res.ls <- lapply(odor.grp, function(x){
    #cat(str(x))
    res <- analyzeCampSigDataTrials(data = x,analysis = analysis,op = op)
  })
  res.ls
}

#****************this might be the old stuff and code ************************
#analyzes an individual group which is a list of trials for the same odor
analyzeCampSigDataTrials <-function(data,analysis=1,op=1){
  if(analysis==1){#get the response values as a DF
    res <- lapply(data, function(x){
      trialval.vec <- x[[2]]
    })
    res <- convertNestedListsDF(res)
    row.names(res) <- 1:getLength(res)
  }
  if(analysis==2){#get cell ranks as a DF
    res <- lapply(data, function(x){
      trialrank.vec <- x[[1]]
    })
    res <- convertNestedListsDF(res)
    row.names(res) <- 1:getLength(res)
  }
  if (analysis==3){
    res <- sapply(data, function(x){
      trialsig.vec <- x[[1]]
      length(trialsig.vec[trialsig.vec>0]) #gets the number of significant cells per trial
    })
  }
  if (analysis==4){
    adata <- data
    names(adata) <- sprintf('t%s',as.character(1:length(data)))
    #cat(names(adata),'\n')
    res <- sapply(adata, function(x){#gets a list of significant cells for each odor
      trialsig.vec <- x[[1]]
      which(trialsig.vec>0) #get all the positions where it is >0      
    })
    res <- unlist(res,use.names = F)
    freq.table <- table(res)
    #now get the list of all numbers that occur once,twice,...
    res <- lapply(1:length(data), function(x){
      as.numeric(names(freq.table[freq.table==x]))
    })
    names(res) <- 1:length(data)
  }
  
  res
}

#function to analyze the values, either response or rank, of the cells
#data: the data from analyzeCampSigDataOdors
#resptype: values, 1 or ranking, 2
#stats: different stats to calculated, 1-mean of non-zero cells, 2 - highest value, 3- lowest value 
getCampOdorStats <-function(data,stats=1,op=1){
  dat.ls <- analyzeCampSigDataOdors(data,analysis = 1) #values
  datr.ls <- analyzeCampSigDataOdors(data,analysis = 2) #ranks
  res.ls <- lapply(1:length(dat.ls), function(x){
    #the stats to get
    res1 <- dat.ls[[x]]
    res1[res1>0] <- 1
    res1 <- apply(res1,1,sum) #no of significant responses for each cell
    res2 <- table(res1) #gfrequency of each cell's response
    res3 <- length(res1[res1>0]) #number of responsive cells across trials
    if (stats==1){ #avg of cellular resp
      avg <- dat.ls[[x]] #get the cellular resp DF
      res4 <- apply(avg,1,function(y) mean(y[y>0])) #avg. non-zero vals
      res4[is.na(res4)] <- 0 
      res4 <- round(res4,2) #round to 2 decimals
      avgr <- datr.ls[[x]]
      res5 <- apply(avgr,1,function(y) mean(y[y>0]))
      res5[is.na(res5)] <- 0
      res5 <- round(res5,2)
    }
    nosigtrials <- ceiling(length(dat.ls[[x]])/2)
    res6 <- mean(sapply(dat.ls[[x]], function(y) length(y[y>0]))) #avg. no of sig. cells per trial
    res7 <- length(res1[res1>=nosigtrials]) #counts no. of cells with reliable responses
    res8 <- which(res1>=nosigtrials) # gets those cells
    names(res8) <- NULL
    res9 <- length(res1[res1 < nosigtrials & res1!=0]) #no. of non-reliable cells with sig responses
    res10 <- which(res1 < nosigtrials & res1 != 0)
    res <- list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)
    names(res) <- c('#sig. resp.','freq. sig','# sig. cells','avg. resp','avg. rank','#sig. / trial','# rel.','rel. cells','# non. rel.','non rel. cells')
    res
  })
  names(res.ls) <- names(dat.ls) #as res.ls is analysis of dat.ls. so names should be same
  res.ls
}

#this function consolidates all the earlier functions. Given a directory, will get the responses, calculate significant cells as well as their
#their responses per trial and across odors
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for (header,data) files
#alpha is used for calculating the SDs above mean where the significance should be fixed
#meanmax: mean or max signal to be used for calculating significance thresholds 1 - max, 2 - mean
#sigval: the signal value to be returned is based on mean =1  or max=2
#rank: basis for ranking significant responses. 1 - max sig -threshold, 2 - z score, 3 - max signal signal
#analysis: 1, get the number of significant cells per trial per odor, 2 - get the number of common cells across trials, 3 - get the response values
#4- get the rank values
#resptype: values, 2 or ranking, 1
#stats: different stats to calculated, 1-mean of non-zero cells, 2 - highest value, 3- lowest value 
#op
computeCampStats <-function(foldname='.',fpatterns=c('tseries.*head.*csv','tseries.*data[0-9]+.*csv'),alpha=0.01,
                            rank=2,maxmean=1,sigval=1,analysis=1,stats=1,resptype=1,op=1){
  #get the significance data, the ranking and signal for each trial of every odor
  data.ls <- calcCampDataSigList(foldname = foldname,fpatterns = fpatterns,alpha = alpha,rank = rank,maxmean = maxmean,sigval=sigval)
  #analysis of this data
  res.anal <- analyzeCampSigDataOdors(data = data.ls,analysis = analysis)
  res.stats <- getCampOdorStats(data.ls,stats = stats)
  #measuring jacard similarities
  #getCampOverlapScores(res.stats)
  list(data.ls,res.stats)
}

#the significance itself is from Campbell's measurements
computeCampStatsSig <-function(foldname='.',fpatterns=c('tseries.*head.*csv','tseries.*sigpvals.*csv'),alpha=0.01,
                            rank=2,maxmean=1,sigval=1,analysis=2,stats=1,resptype=1,op=1){
  #get the significance data, the ranking and signal for each trial of every odor
  data.ls <- getCampDataSigList(foldname = foldname,fpatterns = fpatterns,alpha = alpha)
  #analysis of this data
  res.anal <- analyzeCampSigDataOdors(data = data.ls,analysis = analysis)
  notrials <- length(res.anal[[1]])
  nocells <- getLength(res.anal[[1]])
  res.stats <- getCampOdorStats(data.ls,stats = stats)
  res1 <- sapply(res.anal, function(x) {
    tmp <- apply(x, 1, sum)
    #cat(tmp,'\n')
    length(which(tmp>=(notrials/2) ))
  })
  res2 <- sapply(res.stats, function(x) x[[7]])
  
  list(notrials,nocells,mean(res1)/notrials,mean(res1)/notrials*(100/nocells),res1,res2)
}

#***********************old stuff probably ends here************************************

#Given the frequnecy list of cell's significant responses and list of significant cells, calculates the overlap score
#freq: is the frequency list for each cell over the 6 trials
#cell.list: is the list of significant cells for every odor across the 6 trials
#probno: whether the result shoudl be a probability or numbers
#op= compare core to core or non-core to noncore, 2 - compare non-core to core and core to non-core
getCampOverlapScores <- function(freq,cell.lst,cell2.lst=c(),notrials=6,probno=1,op=1){
  #cat(names(cell.lst),length(cell.lst),notrials,names(notrials),'end\n')
  odors <- names(cell.lst) #the names of cell.lst and notrials should coincide.
  #cat('\nidentical',identical(odors,names(cell.lst)))
  res <- sapply(1:length(cell.lst), function(x){#traverse the list of sig. cells for each odor, e.g. 7 odor list
    #working logic: cell.lst gives the cells that are in the core or non-core set being compared for option 1
    #freq.lst: its purpose is obvious from the name. given two odors, take the interesection of their cell.lst and then
    #compute probabbilities based on freq.lst. For option 2, take intersection of cell.lst and cell2.lst
    tmp <- sapply(1:length(cell.lst), function(y){#do it again
      #cat(' and ',y,' ::')
      over.set <- switch(op,intersect(cell.lst[[x]],cell.lst[[y]]),intersect(cell.lst[[x]],cell2.lst[[y]])) 
      #get the cells common to both odors, and calcaulte the prob. of them firing for both
      if (length(over.set)>0) over <- sapply(over.set, function(z){
        prob <- freq[[x]][z]*freq[[y]][z]/(notrials[odors[x]]*notrials[odors[y]])#(notrials^2) #liam, should this be hardcoded here, only if the no of trial is always 6
      })
      else over <- c(0)
      #cat('overlap sets:',str(over),'\t')
      switch(probno,sum(over),length(over.set)) #sum all the probabilities. That's no of common cells
    })
  })
  #cat(names(freq))
  if (length(names(freq)) > 1) {#multuple cols and rows
    colnames(res) <- names(freq)
    row.names(res) <- names(freq)
  }
  res
}

#Given the frequnecy list of cell's significant responses and list of significant cells, calculates the difference scores
#which is the difference in the probabilities.
#freq: is the frequency list for each cell over the 6 trials
#cell.list: is the list of significant cells for every odor across the 6 trials
#compare: 1: get the overlap, 2: get the difference in probability
#op= compare core to core or non-core to noncore, 2 - compare non-core to core and core to non-core, 3 - comapre core and non-core with cells with no response
#4 = 3 and 4
#probno: whether the result shoudl be a probability or numbers
getCampDiffScores <- function(freq,cell.lst,cell2.lst=c(),notrials=6,compare=1,probno=1,op=1){
  allcells <- 1:length(freq[[1]]) #no of cells in this odor set list
  res <- sapply(1:length(cell.lst), function(x){#traverse the list of sig. cells for each odor, e.g. 7 odor list
    #working logic: cell.lst gives the cells that are in the core or non-core set being compared for option 1
    #freq.lst: its purpose is obvious from the name. given two odors, take the interesection of their cell.lst and then
    #option 2: is a little different for diff than overlap
    #compute probabbilities based on freq.lst. For option 2, take intersection of cell.lst and cell2.lst
    tmp <- sapply(1:length(cell.lst), function(y){#do it again
      if(op==1) over.set <- campGetCellsPair(odor1 = x,odor2 = y,rel1.lst = cell.lst,rel2.lst = cell.lst,celltype = op,op=1)
      else over.set <- campGetCellsPair(odor1 = x,odor2 = y,rel1.lst = cell.lst,rel2.lst = cell2.lst,celltype = op,op=1)
      if (length(over.set)>0) over <- sapply(over.set, function(z){
        prob <- abs(freq[[x]][z]-freq[[y]][z])/(notrials) #liam, should this be hardcoded here, only if the no of trial is always 6
      })
      else over <- c(0)
      switch(probno,sum(over),length(over.set)) #sum all the probabilities. That's no of common cells
    })
  })
  #cat(names(freq))
  if (length(names(freq)) > 1) {#multuple cols and rows
    colnames(res) <- names(freq)
    row.names(res) <- names(freq)
  }
  res
}

#result: the number of cells that are considreed for the difference and the amount of difference are highly correlated
#is this inversely correlatd with overlap measures.


#gets the cells btween two paits of odors that are 
#celltype=1, reliable or unreliable in both, 2, reliable in one and not the other, 
#3 - sig. responses for one but silent for the other, 4 - options 3 and 4
#odor1: 1st odor no, odor2: second odor no
#rel1.lst: list of reliable cells for all odors
#rel2.lst: list of unreliable cells for all odors, or reliable cells if its celltype=1
#op: get probabilitiy differences, 2 - get numbers of cells in each group
campGetCellsPair <- function(odor1,odor2,rel1.lst,rel2.lst,celltype=1,op=1){
  if(celltype==1){
    over.set <- intersect(rel1.lst[[odor1]],rel2.lst[[odor2]])
    #get the cells common to both odors, and calcaulte the prob. of them firing for both
  }
  if(celltype==2 || celltype==4) {
    # types of comparisons of 2 odors: b/w cells that are reliable in one and not the other, and 
    s1 <- union(intersect(rel1.lst[[odor1]],rel2.lst[[odor2]]),intersect(rel2.lst[[odor1]],rel1.lst[[odor2]]))
    over.set <- s1
  }
  if(celltype==3 || celltype==4){    #cells that have some significant response and silent otherwise.
    #all cells that are in 1
    o1 <- union(rel1.lst[[odor1]],rel2.lst[[odor1]])
    #all cells that are in 2
    o2 <- union(rel1.lst[[odor2]],rel2.lst[[odor2]])
    over.set <- union(setdiff(o1,o2),setdiff(o2,o1))
    #cat('\ndiff',setdiff(o1,o2),'\t',setdiff(o2,o1),'\nover: ',over.set)
    #cat('\ndiff',o1,'\n',o2)
  } 
  if(celltype==4) {
    over.set <- union(s1,union(setdiff(o1,o2),setdiff(o2,o1)))
    #cat('\n',s1,'\tdiff',setdiff(o1,o2),'\t',setdiff(o2,o1),'\nover: ',over.set)
    
  }
  switch(op,over.set,length(over.set)) 
}


#Given the frequnecy list of cell's significant responses and list of significant cells, calculates the difference scores
#which is the difference in the probabilities.
#freq: is the frequency list for each cell over the 6 trials
#cell.list: is the list of significant cells for every odor across the 6 trials
#compare: 1: get the overlap, 2: get the difference in probability
#op= compare core to core or non-core to noncore, 2 - compare non-core to core and core to non-core, 3 - comapre core and non-core with cells with no response
getCampDiffScores.old <- function(freq,cell.lst,cell2.lst=c(),notrials=6,compare=1,op=1){
  #cat(str(freq))
  allcells <- 1:length(freq[[1]]) #no of cells in this odor set list
  cat('allcells',allcells,'\n')
  res <- sapply(1:length(cell.lst), function(x){#traverse the list of sig. cells for each odor, e.g. 7 odor list
    #working logic: cell.lst gives the cells that are in the core or non-core set being compared for option 1
    #freq.lst: its purpose is obvious from the name. given two odors, take the interesection of their cell.lst and then
    #option 2: is a little different for diff than overlap
    #compute probabbilities based on freq.lst. For option 2, take intersection of cell.lst and cell2.lst
    tmp <- sapply(1:length(cell.lst), function(y){#do it again
      #cat(' and ',y,' ::')
      if(op==1){
        over.set <- intersect(cell.lst[[x]],cell.lst[[y]])
        #get the cells common to both odors, and calcaulte the prob. of them firing for both
      }
      else {
        # 2 types of comparisons of 2 odors: b/w cells that are reliable in one and not the other, and 
        #cells that have some significant response and silent otherwise.
        s1 <- union(intersect(cell.lst[[x]],cell2.lst[[y]]),intersect(cell2.lst[[x]],cell.lst[[y]]))
        #all cells that are in2
        o2 <- union(cell.lst[[y]],cell2.lst[[y]])
        #all cells that are in 1
        o1 <- union(cell.lst[[x]],cell2.lst[[x]])
        o1complement <- union(o1,setdiff(allcells,o2)) #all the cells that are in 1 but not 2
        o2complement <- union(o2,setdiff(allcells,o1)) #all the cells that are in 2 but not 1
        over.set <- union(s1,union(setdiff(o1,o2),setdiff(o2,o1)))
        tmp.over <- sapply(allcells, function(z){
          tst <- 0
          if((freq[[x]][z] < 4 && freq[[y]][z] >= 4) || (freq[[x]][z] >= 4 && freq[[y]][z] < 4) || (freq[[x]][z] == 0 && freq[[y]][z] > 0) 
             || (freq[[x]][z] > 0 && freq[[y]][z] == 0)) tst <- z
          #cat('\t',z,freq[[x]][z],freq[[y]][z],tst)
          tst
        })
        #tmp.over and over.set are two different ways of doing the same thing. I don think we need o1complement
      }  
      #cat(sort(over.set),'\n',getAboveThresh(tmp.over),'\n')
      #cat(x,y,length(over.set),length(tmp.over),'\t')
      if (length(over.set)>0) over <- sapply(getAboveThresh(tmp.over), function(z){
        prob <- abs(freq[[x]][z]-freq[[y]][z])/(notrials) #liam, should this be hardcoded here, only if the no of trial is always 6
      })
      else over <- c(0)
      #cat('\n',freq[[x]][over.set],'\n',freq[[y]][over.set],'\n',over)
      sum(over) #sum all the probabilities. That's no of common cells
    })
  })
  #cat(names(freq))
  if (length(names(freq)) > 1) {#multuple cols and rows
    colnames(res) <- names(freq)
    row.names(res) <- names(freq)
  }
  res
}

#given two vectors calculates Jacard similarity
#op=1, jacard, 2 - overlap, 3 - dice
overlapScore<-function(vec1,vec2,op=1){
  switch(op,length(intersect(vec1,vec2))/length(union(vec1,vec2)),
         length(intersect(vec1,vec2))/min(length(vec1),length(vec2)))
}

#this computes the lifetime sparseness of the response data
#dat.df: the data across all trials or events. The rows are cells, while each column is the cell's
#response for each trial or event
#op: type of sparseness, 1 - your measure of norming using the max, 2 - Carandini lifetimes sparseness,
#-lifetime kurtosis
computeLifeSparseness <- function(dat.df,op=1){
  n <- nrow(dat.df) #no of neurons or rows
  if(op==1){#using max to norm
    res <- cleanNAVec(sapply(1:n, function(i) sum(dat.df[i,])/max(dat.df[i,]) ) )
    res <- (1/res)^2
    res[res==Inf] <- 0
    return(res)
  }
  if(op==2){#using carandini lifetime sparseness
    res <- cleanNAVec(sapply(1:n, function(i) (sum(dat.df[i,])^2)/sum(dat.df[i,]^2) ) )
    #cat('\ncLS',res)
    res <- (1 - res/n)/(1 - 1/n)
    return(res)
  }
  if(op==3){#lifetime kurtosis
    m <- ncol(dat.df)
    n.mn <- sapply(1:n,function(i) mean(unlist(dat.df[i,])) )
    n.sd <- sapply(1:n,function(i) sd(unlist(dat.df[i,])) )
    #cat('\ncLS',n,'sds',n.mn,n.sd)
    res <- sapply(1:n, function(i) {
      #cat('\n',unlist(dat.df[i,]),':',str(unlist(dat.df[i,]) ) )
      tmp <- (1/m) * sum(((unlist(dat.df[i,]) - n.mn[i])/n.sd[i])^4) - 3  
    }) 
    res <- cleanNAVec(res)
    return(res)
  }
}

#this function calculates the correlation between avg. response levels and no of significant responses for different
#combinations of rank, i.e. calculting the signal level, and maxmean whether we use the signal max or mean
calcCampParamsStats <- function(op=1){
  seqrank <- c(1,2,3)
  seqmaxmean <- c(1,2)
  res.ls <- lapply(seqmaxmean, function(x){
    res <- sapply(seqrank, function(y){
      tmp <- computeCampStats(rank = y,maxmean = x)
      tmp1 <- sapply(1:7, function(z) cor(tmp[[2]][[z]][[4]],tmp[[2]][[z]][[1]]))
      mean(tmp1)
    })
  })
  res.ls
}


#write functions that consolidates all the tasks so we can compare the different methods for obtaining significance
#functions to compare the cell overlaps for different odors, between the top 5% and the rest. 
#Also have to do this for responses between trials. what is the average amount of overlap


#write functions to get the corresponding values for each of the cells across the different trials, as well as their ranking (maybe average values to sstart with?)
# or a data frame where the rows are cells, and cols these values for each trial or odor trial
#also a data frame of cells and the number of trials for which they had significant responses
#you can do correlations to determine the similarity in odor responses

#function that gets all the non-zero position rankings and checks if the highest rank exceeds the number of ranked positions
#op= 1 - returns the difference, 2-returns T if rank positions and highest rank match; 3 - number of ranks over number of positions
checkRankSize <-function(vect,op=3){
  vec <- vect[vect>0]
  res <- switch(op,max(vec)-length(vec),length(vec) < max(vec),length(vec[vec>length(vec)]))
  res
}

#a function that takes the non-zero ranking positions, and re-ranks them so that if the list is n-long the ranking goes from 1 to n instead of 
#1 to n+x. e.g., 1  3 2 3 12 12 121 21 1 to 1  3 0 3 12 0 121 0 1
#op=1, vector interspersed with 0s, 2-compact vector
compactRanks <- function(vect,op=1){
  tmp <- vect[vect>0]
  #just re-rank them 
  res <- getMatRank(tmp,decrease = F)   
  if(op==2) {#get the non-zero posns, and populate them with the new ranks
    posns <- which(vect>0)
    tmp <- rep(0,length(vect))
    tmp[posns] <- res
    res <- tmp
  }
  res
}

#this is a function of all the functions to call to obtain results
campbell_analysis_functions<-function(){
  #correlations between sig. responses and rank or responses of cells to odors
  res4 <- computeCampStats(rank = 2,maxmean = 2)
  sapply(1:7, function(x) cor(res4[[2]][[x]][[4]],res4[[2]][[x]][[1]]))
  # 0.8113442 0.8445220 0.8476637 0.8616899 0.7703134 0.8131225 0.9077846
  #odor3_rank_vs_nosig_max_thresh
  fploteq(res[[2]][[3]][[5]],res[[2]][[3]][[1]],xlabel = 'avg. rank',ylabel = '# sig. responses',ticknoy = 1,ticknox = 1)
  #odor3_resp_vs_nosig_max_thresh
  fploteq(res[[2]][[3]][[4]],res[[2]][[3]][[1]],xlabel = 'avg. resp.',ylabel = '# sig. responses',ticknoy = 1,ticknox = 5)
  #odor3_resp_vs_nosig_mean_sig
  fploteq(res4[[2]][[3]][[4]],res4[[2]][[3]][[1]],xlabel = 'avg. resp.',ylabel = '# sig. responses',ticknoy = 1,ticknox = 5)
  
  #getting the best possible stats
  tmp <- calcCampParamsStats()
  #0.5577379 0.6850918 0.6798632 max, with sig-thresh,zscore,sig
  #0.7171980 0.8366343 0.8202276 mean, with sig-thresh,zscore,sig
  
  
  #determining overlap scores
  res <- computeCampStats(rank = 3,maxmean = 2)
  tmp <- lapply(res[[2]], function(x) x[[8]])
  getCampOverlapScores(lapply(res[[2]], function(x) x[[1]]),lapply(res[[2]], function(x) x[[8]]))
  
  #plotting the overlap scores
  tst2 <- read.table(file = 'plots/csv_tser_1458/tseries_nonreliable_overlap.csv')
  tst1 <- read.table(file = 'plots/csv_tser_1458/tseries_reliable_overlap.csv')
  heatmap(as.matrix(tst1),scale='none',Rowv = NA,Colv = NA,labRow = '',labCol = '',col=paste("gray",1:100,sep=""),symm = T)
  heatmap(as.matrix(tst2),scale='none',Rowv = NA,Colv = NA,labRow = '',labCol = '',col=paste("gray",51:100,sep=""),symm = T)
  image(as.matrix(tst2),col = grey(seq(.83,1,.001)),xlab = '',ylab = "",frame=F,axes=F)
  image(as.matrix(tst1),col = grey(seq(0,1,.001)),xlab = '',ylab = "",frame=F,axes=F)
  
  #plotting response vs no. of sig cells
  tmp <- getCampDataSigList()
  tmp3 <- getCampResponseSigCor(tmp)
  fploteq.DF(tmp3[[1]],xlabel='no. trials',ylabel='mean response',ticknoy=5,ticknox=2)
  
  #for getting the stats on the Matrix data in any particular directory, do

  res2 <- computeCampOverlap(res)
  
  #command for traversing through directorues
  CampDirTraverse(searchstr = c('*tseries*','*sigpval*'),mode = 2,targetfn = getCampMatScores)
  
  #computing the overlap and difference scores
  tmp <- getCampDataSigList() #get the data
  tmp2 <- campClassifyCells(tmp)
  tmp3 <- getCampDiffScores(tmp2[[1]],tmp2[[2]],tmp2[[3]],notrials = 6,probno = 1,op=3) #diff for cells that are sign for odor1, silent for odor2
  tmp4 <- getCampOverlapScores(tmp2[[1]],tmp2[[2]],tmp2[[3]],notrials = 6,probno = 1,op=1) #overalp of reliable cells
  diag(tmp4) <- 0 #have to take the diagonal out of it, maybe not actually since it is relevant for overlap but not difference
  cor(as.vector(tmp3),as.vector(tmp4)) #their correlation
  #comparing differences and the number of responding cells for case where sign for odor1, silent for odor2
  tmp3 <- getCampDiffScores(tmp2[[1]],tmp2[[2]],tmp2[[3]],notrials = 6,probno = 1,op=3)
  tmp4 <- getCampDiffScores(tmp2[[1]],tmp2[[2]],tmp2[[3]],notrials = 6,probno = 2,op=3)
  cor(as.vector(tmp3),as.vector(tmp4))
  #differences between cells that are unreliable 
  tmp4 <- getCampDiffScores(tmp2[[1]],tmp2[[2]],tmp2[[3]],notrials = 6,probno = 2,op=2)
  tmp3 <- getCampDiffScores(tmp2[[1]],tmp2[[2]],tmp2[[3]],notrials = 6,probno = 1,op=2)
  cor(as.vector(tmp3),as.vector(tmp4))
  #result: as the number of cells being considered increses the amount of difference between the odor pairs also increases. Thus, if there is more 
  #difference between two odors, they are likelier to have more diffrences in firing, too
  
  
  
}


#this function characterizes the reliability distributions of unreliable cells
#it gets all the cells that respond unrel number of times in the first notrials
#alldata: the return of getCampDataSigList
#unrel,notrials: the frequency of response per notrials trials, e.g, 1 in 3 trials would be unrel=1, notrials = 3
campGetUnrelDist <- function(alldata,odors=c(),unrel=1,notrials=3,op=1){
  #seeing if cells with a singl e response in the first 3 trials are more biased to respond again depending on the trial number
  #get all the responses across all trials for the specified odors or all odors
  alltrials <- campHabituation(alldata,op=4)
  if(length(odors) > 0) alltrials <- alltrials[odors]
  alltrials.df <- joinUnevenDFs(alltrials)
  #determine no of trials and fix the names of all the rows
  total.trials <- ncol(alltrials.df)
  rownames(alltrials.df) <- 1:nrow(alltrials.df) #the rowname is the row no
  #get the unreliable cell rows in for the first notrials trials
  unrel.cells <- getRowsDFNonzeros(alltrials.df[,1:notrials],n=unrel,op=2) #get all cells with 'unrel' response(s) exactly in the first 'notrials' trials
  unrel.rows <- alltrials.df[rownames(unrel.cells),] #get the rows for these unrel cells specifically 
  #now, get all those cells from here that fire again in the rest of the trials, and the frequency of their firing
  relrows <- c(list(unrel.rows),sapply((unrel+1):(unrel+total.trials-notrials),function(i) getRowsDFNonzeros(unrel.rows,n=i,op=2) ) )
  lenrows <- length(relrows)
  names(relrows) <- 1:lenrows
  #only consider those rows that have cells in them
  nonempty.rows <- which(sapply(1:lenrows,function(i) length(relrows[[i]])) > 0)
  relrows <- relrows[nonempty.rows]
  lenrows <- length(relrows)
  #now, make a cumulative list of rownames,
  relrows.cum <- list(rownames(relrows[[lenrows]]) )
  freqnames <- c(1)
  #make the cumulative list
  for(i in seq_wrap(from=1,to=length(relrows)-1 )){
    if(length(relrows[[i]]) > 0){ #only add those frequencies that have cells
      freqnames <- c(freqnames,i) #the frequency that is currently being processed
      relrows.tmp <- c(rownames(relrows[[lenrows-i]]),unlist(relrows.cum[[length(relrows.cum)]])) 
      #cat('\ncum',i,lenrows-i,length(relrows.tmp),length(rownames(relrows[[lenrows-i]])) )
      relrows.cum <- c(relrows.cum,list(relrows.tmp) )
    }
  }
  names(relrows.cum) <- freqnames
  #separate out the reliability rows in the cumulative lists, and order them based on the first resposive trial
  relrows.cum <- rev(relrows.cum)
  res <- lapply(1:(length(relrows.cum)-1),function(i) {
    relnames <- setdiff(relrows.cum[[i]],relrows.cum[[i+1]])
    #cat('\nnams',i,length(relnames),':',relnames)
    reltrials <- alltrials.df[orderMatRowsPrecedence(alltrials.df[relnames,])[,3],]
  } )
  if(isDataType(relrows.cum[[length(relrows.cum)]])==1) res <- c(res,list(alltrials.df[relrows.cum[[length(relrows.cum)]],]) )
  else res <- c(res,list(alltrials.df[orderMatRowsPrecedence(alltrials.df[relrows.cum[[length(relrows.cum)]],])[,3],]))
  res
}


#generates stats of how often you will see certai frequencies in a population
#This is for figuring out if the frequency of firing 1/2 or 1/3 etc. that is seen in a population is feasible
#cellfreq: how often the cell responds for unrel.trials # of trials
#rangefreq: the number of frequencies you should consider, the numerator
#cellfreq: the number of occurences out of unrel.trials
#unrel.trials: the number of trials for calculating the base probability
#probavg: how you should average probabilities: 1 - a log based base prob, 2 - an additive based base prob
#op=1, generate the population composition, 2 - requires observed, and compare with expected and give the probability
#weight: how all the reliabilitues should be weighted, could be a distribution or some other type of function
#freqop: whetehre you are going to go (1) 1/2, 1/3... or (2) 2,3,...
campGenRelStats <- function(nocells,rangefreq=c(3:12),unrel.trials=3,notrials=6,obs=c(),cellfreq=1,probavg=1, 
                            weight=function(x) x,freqop=1,op=1){
  #claculate the poulation composition, need for op=1 and 2
  orig.prob <- cellfreq/unrel.trials
  #first get the predictions for each cell frequency for how many cells in the next 'x' trials will be silent, fire once, twice, ...
  switch(probavg,cell.freqs <- lapply(rangefreq,function(i) {
    theory.prob <- switch(freqop,cellfreq/i,i/cellfreq)
    prob <- sapply(0:(notrials-unrel.trials),function(j) dbinom(x = j,size = (notrials-unrel.trials),prob = theory.prob) * nocells )
    #cat('\nprob',cellfreq/i,roundVecToSum(prob,total = nocells),'\t',prob)
    roundVecToSum(prob,total = nocells)
  }),
  cell.freqs <- lapply(seq((orig.prob-0.5),(orig.prob+0.1),0.05),function(i) {
    prob <- sapply(0:(notrials-unrel.trials),function(j) dbinom(x = j,size = (notrials-unrel.trials),prob = i) * nocells )
    #cat('\nprob',round(i,2),roundVecToSum(prob,total = nocells),'\t',prob)
    roundVecToSum(prob,total = nocells)
  }) )
  #cat('\ncGNr',orig.prob,str(cell.freqs),rangefreq)
  #put it in a data frame format
  freqs.df <- transposeDF(convertNestedListsDF(cell.freqs))
  #get the avg
  avg <- sapply(1:nrow(freqs.df),function(i) {
    theory.prob <- switch(freqop,cellfreq/rangefreq[i],rangefreq[i]/cellfreq)
    #cat('\ncGNR',i,theory.prob)
    prob.wt <- unlist(freqs.df[i,])*weight(theory.prob)
  })
  #cat('\ncGNR',str(t(avg)),as.matrix(t(avg)) )
  #cat(str(avg))
  freqs.df <- rbind(freqs.df,roundVecToSum(apply(freqs.df,2,mean),nocells) )
  #make sure the names are in the format of trial frequency 1/x
  rnames <- switch(probavg,sapply(rangefreq,function(i) paste('',as.character(i),sep = '')),
                   sapply(seq((orig.prob-0.5),(orig.prob+0.1),0.05),function(i) as.character(round(i,2)) ) )
  #cat(str(freqs.df),str(rnames))
  rownames(freqs.df) <- c(rnames,'avg.')
  colnames(freqs.df) <- 0:(notrials-unrel.trials)
  

  if(op==2){#calculate the probability of the observed population given the various base probabilities
    #use the chi square test  
    probs <- sapply(1:nrow(freqs.df),function(i) {
      tmp <- unlist(freqs.df[i,])
      kisq <- cleanNAVec( (tmp - obs)^2/tmp )
      #cat('\n',i,': ',sum(kisq))
      significance <- 1 - pchisq(sum(kisq),df=length(obs)-1)
    }) 
    names(probs) <- rnames
  }
  switch(op,freqs.df,probs)
  #t(avg)
}

#for dataset will cpmpute the overlap score for all odor pairs by reliability or another metric speciified in op
#op=1, per odor, 2- per cell 
#overop: the overlap condition or op for campOVerlapAllAnalog
#trialop: option for campGetTrialDetails; 1 - get it from the directory, 2 - get it from the data
#thresh: the threshold to filter similar or dissimilar odor correlations
#threshop: 0 - no threshold, 1 - above threshold, 2- below threshold, 3 - between thresh1 and thresh2
campComputeOverlap <- function(datasig,skip=c(),overop=9,trialop=2,thresh=0.5,threshop=1,op=1){
  #reliability vs overlap, do it for all cells per odorpair, and then all odor-pairs
  #overlap.allcells <- campOverlapAllAnalog(data = datasig,op = overop)
  overlap.allcells <- campOverlapAllAnalog(datasig,skip = skip,op=overop,trialop = trialop) #op=9, overlap all cells
  #cat('\ncCO',str(overlap.allcells))
  overlap.odorpairs <- lapply(overlap.allcells,function(x) sumDF(x))
  over.allodors <- joinListDFs(overlap.odorpairs)
  over.mean <- meanDF(over.allodors)
  over.byrel <- split(over.allodors,over.allodors[,1])
  #return a list of overlap where each list element is a vector of odor overlap for that level of reliability
  res <- lapply(over.byrel,function(x) x[,2])
  res <- campThreshByCorr(switch(op,res,overlap.allcells),datasig = datasig,thresh = thresh,threshop = threshop)
  res
}

#thresholds the results of res.ls by applying thresh to the correlations between odors and thresholding
#datasig:return of getCampDataSigList()
#res.ls: the list or vector that is to be thresholded. Assume, that the elements have names in the form odor1,odor2
#thresh: the threshold to filter similar or dissimilar odor correlations
#threshop: 0 - no threshold, 1 - above threshold, 2- below threshold, 3 - between thresh1 and thresh2
campThreshByCorr <- function(res.ls,datasig,thresh,threshop,op=1){
  res <- res.ls
  if(thresh==0) return(res)  
  res.cor <- computeListMatCorr(alldata = datasig,matop = 2,op=1)#similarity on averaged vectors, and pearson
  res.cor.vec <- unlist(flattenDFtoList(res.cor))
  if(threshop==1) res <- res.ls[names(which(res.cor.vec > thresh))]
  if(threshop==2) res <- res.ls[names(which(res.cor.vec < thresh))]
  if(threshop==3) res <- res.ls[names(which(res.cor.vec > thresh[1] & res.cor.vec < thresh[2]))]
  res
}

#computes the correlation between reliable and unreliable cells for different odor-pairs,
#filters them by odor similarity
#skip: skip odor pairs with this odor
#op: 1- above a thresh. similarity, 2 = below, 3 - in between
campRUCellCorr <- function(dat.lst,thresh=0.5,skip=c(''),corop=1,op=1){
  #get the correlation between odor pairs for all cells, rel. and unrel. cells as a vector
  allcells <- computeListMatCorr(dat.lst,celltype = 4,op = corop)
  relcells <- computeListMatCorr(dat.lst,celltype = 2,op=corop)
  unrelcells <- computeListMatCorr(dat.lst,celltype = 3,op=corop)
  allcells.vec <- unlist(flattenDFtoList(allcells))
  relcells.vec <- unlist(flattenDFtoList(relcells))
  unrelcells.vec <- unlist(flattenDFtoList(unrelcells))
  #now do the filtering
  relcells.res <- switch(op,relcells.vec[which(allcells.vec>thresh)],relcells.vec[which(allcells.vec<thresh)],
                         relcells.vec[which(allcells.vec>thresh[1] & allcells.vec<thresh[2])])
  unrelcells.res <- switch(op,unrelcells.vec[which(allcells.vec>thresh)],unrelcells.vec[which(allcells.vec<thresh)],
                           unrelcells.vec[which(allcells.vec>thresh[1] & allcells.vec<thresh[2])])
  #take out skip odor-pairs
  relcells.res <- relcells[getSkipStrIndex(names(relcells.res),skip=skip)]
  unrelcells.res <- unrelcells[getSkipStrIndex(names(unrelcells.res),skip=skip)]
  #return as a list for plotting with fstripchartvecs
  res <- list(unrelcells.res,relcells.res)
  names(res) <- c('unrel.','rel.')
  res
}


hallemfunctions <- function(){
  #converting the hallem names to rob;s names
  #algo: get rob's names, search for a substring in hallem's names,
  #then make sure it is the same one and then subsitute it
  #idea: write a function to name match which will normalize the differences
  #between - and blank,i.e., first-second == first second
  tmp1 <- unique(tmp[,1])
  tmp1 <- as.character(tmp1)
  which(grepl(pattern = 'octano',x = hallem[,1],fixed = T))
  hallem1[105,1]
  hallem1[105,1] <- tmp1[4]
  
  #gtting robnames goign
  hallem1 <- hallem
  hallem2 <- read.csv('../../../camptest/hallem_robnames.csv',header = T)
  hallem1[,1] <- hallem2[,2]
  hallem1 <- ConvertDfCols(hallem1,1,op=2)
  
  #this gets the similarity between oddor pairs based on hallem, and the corresponding 
  #overlap at each reliability level for each odor pair
  tmp <- getCampDataSigList()
  res <- campFreqCellsAnalog(tmp)
  res3 <- campOverlapAllAnalog(tmp,nolevels = 2)
  res4 <- campGetHallemCorr(hallem1,getSkipStrIndex(names(res[[2]]),skip = 'paraff',op=2))  
  res5 <- campPrepareAnalog(res3)
  res6 <- campMatchHallemRob(res5,res4) #return with first col the hallem cor, and the others showing the overlap for each KC reliability measure for each odor pair
  fploteq(x=res6[,1],y=res6[,-c(1,8:13)],xlabel = 'similarity(hallem)',ylabel = 'overlap(rob)',ticknox = 2,ticknoy = 5)
  fploteq(x=res6[,1],y=res6[,-c(1,2:7)],xlabel = 'similarity(hallem)',ylabel = 'overlap(rob)',ticknox = 2)
  
  #splitting the data into reliable and unreliable cells
  res7 <- campMatchHallemRob(res5,res4,op=2) 
  fploteq(x=res7[,1],y=res7[,-c(1)],xlabel = 'similarity(hallem)',ylabel = 'overlap(rob)',ticknox = 2,ticknoy = 5)
  
  #getting a characterization of overlap and difference and plotting some examples of each
  res6[which(res6[,1]<0.1),] #get odor pairs below a certain threshold
  res <- campFreqCellsAnalog(tmp) #getting the frequency of cells for each odor
  temp1 <- campGetTrialDetails()[[1]]
  #get the overlap between two similar odor pairs
  temp2 <- campOverlapAllPairAnalog(odor2 = '2-heptanone',odor1 = 'pentyl-acetate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=7)
  fploteq(temp2[,1],temp2[,2],ylabel = 'overlap',xlabel = 'reliability',ticknox = 1,ticknoy = 2)
  temp2 <- campOverlapAllPairAnalog(odor1 = '2-heptanone',odor2 = 'pentyl-acetate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=7)
  fploteq(temp2[,1],temp2[,2],ylabel = 'overlap',xlabel = 'reliability',ticknox = 1,ticknoy = 2)
  #get the difference now
  temp2 <- campOverlapAllPairAnalog(odor2 = '2-heptanone',odor1 = 'pentyl-acetate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=8)
  fploteq(temp2[,1],temp2[,2],ylabel = 'overlap',xlabel = 'reliability',ticknoy = 20,ticknox = 1,fixx = c(1,4),fixy = c(-.2,.4))
  temp2 <- campOverlapAllPairAnalog(odor1 = '2-heptanone',odor2 = 'pentyl-acetate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=8)
  fploteq(temp2[,1],temp2[,2],ylabel = 'useful D.',xlabel = 'reliability',ticknoy = 30,ticknox = 1,fixy = c(-.1,1.1))
  #dissimilar odors, difference
  temp2 <- campOverlapAllPairAnalog(odor1 = 'diethyl-succinate',odor2 = 'ethyl-lactate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=8)
  temp2 <- campOverlapAllPairAnalog(odor2 = 'diethyl-succinate',odor1 = 'ethyl-lactate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=8)
  fploteq(temp2[,1],temp2[,2],ylabel = 'useful D.',xlabel = 'reliability',ticknoy = 20,ticknox = 1)
  #dissimiar odors, overlap
  temp2 <- campOverlapAllPairAnalog(odor1 = 'diethyl-succinate',odor2 = 'ethyl-lactate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=7)
  temp2 <- campOverlapAllPairAnalog(odor2 = 'diethyl-succinate',odor1 = 'ethyl-lactate',freq = res[[2]],notrials = temp1,nolevels = 1,relscore = 2,op=7)
  #same plot for the 4 temp2 plots above
  fploteq(temp2[,1],temp2[,2],ylabel = 'overlap',xlabel = 'reliability',ticknoy = 20,ticknox = 1)
  #result: sometimes the reliability cells do not match against overlap and useful D. because while overlap will be 0 if one of the odors has a 0 response, it
  #would be positive/negative for useful D.
  
  #gets the summary of reliable, significant, and unreliable cells for this data set
  temp3 <- getCampDataRel(tmp,op=2)
  fstripchartvecs(as.list(temp3),markersize = 0.8,ylabel = '% cells',semthick = 2,pairplot = 0,tickno = 10,fixy = c(1,60))
  
  #computing for the 0-mean hallem data set
  hallem3 <- normalizeMat(hallem1[,-1],scalemean = 100,changemean = T,op=2)
  hallem2[,-1] <- hallem3
  res41 <- campGetHallemCorr(hallem2,getSkipStrIndex(names(res[[2]]),skip = 'paraff',op=2))
  res71 <- campMatchHallemRob(res5,res41,op=2)
  fploteq(x=res71[,1],y=res71[,-c(1)],xlabel = 'similarity(hallem)',ylabel = 'overlap(rob)',ticknox = 2,ticknoy = 5)
  campHallRobGetStats(res71,relno = 2)
  
  #cpmparing frequencies
  names(tst) <- names(res[[2]])
  tst[which(as.vector(apply(tst,1,prod)) > 0),]
  
  #significant correlations
  tmp <- getCampDataSigList()
  tst <- getGroupedListsDf(tmp,op=2)
  tst1 <- computeListMatCorr(tst)
  mean(sapply(1:6,function(x) cor(unlist(res4[x,]),unlist(tst2[c(x),]))))
  fploteq(unlist(flattenDFtoList(res4)),unlist(flattenDFtoList(tst1[-6,-6])),xlabel = 'hallem',ylabel = 'rob',ticknox = 2)
  cor(unlist(flattenDFtoList(res4)),unlist(flattenDFtoList(tst1[-6,-6]))) #for ideal data set, correlation is 0.65
  
  #doing it for differences
  res32 <- campOverlapAllAnalog(tmp,nolevels = 2,op = 2)
  res52 <- campPrepareAnalog(res32)
  res72 <- campMatchHallemRob(res52,res4,op=2)
  fploteq(x=res72[,1],y=res72[,-c(1)],xlabel = 'similarity(hallem)',ylabel = 'diff (rob)',ticknox = 2,ticknoy = 5)
  campHallRobGetStats(res72,relno = 2) #get correlations: unreliable: -0.75, reliable: -0.54
  
  #doing it with probability differences
  #to provide results of aggregated reliability across cells
  res33 <- campOverlapAllAnalog(tmp,nolevels = 2,relscore = 2,op = 3)
  #results of every cell
  res33 <- campOverlapAllAnalog(tmp,nolevels = 2,relscore = 2,op = 4)
  #results of correlations 
  tmp <- getCampDataSigList()
  tst <- getGroupedListsDf(tmp,op=2)
  tst1 <- computeListMatCorr(tst)
  cbind(unlist(flattenDFtoList(res4))[-c(6,7,14,21,28,35)],unlist(flattenDFtoList(tst1[-6,-6]))[-c(6,7,14,21,28,35)])
  fploteq(res33$`diethyl-succinate,pentyl-acetate`[,1],res33$`diethyl-succinate,pentyl-acetate`[,2],xlabel = 'KC reliability',ylabel = 'useful (disc.)')
  fploteq(res33$`ethyl-lactate,ethyloctanoate`[,1],res33$`ethyl-lactate,ethyloctanoate`[,2],xlabel = 'KC reliability',ylabel = 'useful (disc.)')
  #get all correlations with one command
  tst2 <- computeListMatCorr(getGroupedListsDf(getCampDataSigList(),matop=2))
  tst3 <- campMatchHalRobCor(hallem2,tst2) #get the closest correlation pairs
  res34 <- campOverlapAllAnalog(getCampDataSigList(),nolevels = 2,relscore = 2,op = 4)#and get all overlaps and compare pairs from tst3
  
  #making a data frame of odors and their cell differences, and the means of the cells
  res.110109 <- campMatchHalRobCor(hallem2,skip = c('empty','humulene'),op=2)
  #do this for a whole bunch of experiments
  #now, put the data together for all these experiments
  res1 <- campProcessOdorSets(list(res.0904,res.0904.2,res.110108,res.110109,res.110109.2))
  res2 <- t(sapply(res1,function(x) c(x[1,1],x[1,2]))) #shows the useful difference for all cells with reliability 1
  fploteq(res2[,1],res2[,2],xlabel = 'overlap',ylabel = 'difference',ticknox = 2,ticknoy = 2) #cor: -0.37
  res3 <- t(sapply(res1,function(x) c(x[ncol(x),1],x[ncol(x),2]))) #shows the useful difference for all cells with reliability 2
  fploteq(res3[,1],res3[,2],xlabel = 'overlap',ylabel = 'difference',ticknox = 2,ticknoy = 8) #cor: -0.39
  res4 <- t(sapply(res1,function(x) c(x[nrow(x),1],x[nrow(x),2]))) #shows the useful difference for all cells with the most reliability
  fploteq(res4[,1],res4[,2],xlabel = 'overlap',ylabel = 'difference',ticknox = 2,ticknoy = 2) # correlation: -0.76
  #result: the difference increases with reliability for odorpairs that are dissimilar. When similarity is low, the differences between
  #unreliable cells are not that great and oscillate between just above and below zero and its low. But, for reliable cells, the other
  #odor activation is most likely silent or unreliable, so the difference will be higher.
  
  
  #the odor pairs that we tested: similar, dissimilar, and middle
  fploteq(res33$`pentanal,hexanal`[,1],res33$`pentanal,hexanal`[,2],xlabel = 'KC reliability',ylabel = 'useful (disc.)',ticknox = 5,ticknoy = 5)
  fploteq(res33$`hexanal,benzaldehyde`[,1],res33$`hexanal,benzaldehyde`[,2],xlabel = 'KC reliability',ylabel = 'useful (disc.)',ticknox = 5,ticknoy = 5)
  fploteq(res33$`pentyl-acetate,6-methyl-5-hepten-2-one`[,1],res33$`pentyl-acetate,6-methyl-5-hepten-2-one`[,2],xlabel = 'KC reliability',ylabel = 'useful (disc.)',ticknox = 5,ticknoy = 5)
  
  #synapse strengthening and weakening with saturation.
  fplotFn(x=c(0,15),eq = function(x) x^1/(x^1 + 1^1),xlabel = 'trials',ylabel = 'synapse strength',ticknox = 3,ticknoy = 2)
  fplotFn(x=c(0,15),eq = function(x) 1/(1 + (x/1)),xlabel = 'trials',ylabel = 'synapse strength',ticknox = 3,ticknoy = 2)
  
  #comparing odors pairs with and without saturation
  res2 <- campCompareSatNorm(hallem2)
  fploteq(res2[[1]][,1],res2[[1]][,c(4,7)],xlabel = 'similarity',ylabel = 'useful D.',ticknox = 2,ticknoy = 4) #unreliable cells
  fploteq(res2[[1]][,1],res2[[1]][,c(5,8)],xlabel = 'similarity',ylabel = 'useful D.',ticknox = 3,ticknoy = 2) #reliable cells
  fploteq(res2[[1]][,1],res2[[1]][,c(5)]-res2[[1]][,c(8)],xlabel = 'similarity',ylabel = 'useful D.',ticknox = 3,ticknoy = 1) #change from unreliable to reliable
  fploteq(res2[[1]][,2],res2[[1]][,c(3)]-res2[[1]][,c(6)],xlabel = 'similarity',ylabel = 'slope change',ticknox = 3,ticknoy = 2) #slope change from unreliable to reliable  
  fitlinearlaw(res2[[1]][,2],res2[[1]][,c(3)]-res2[[1]][,c(6)])
  
  #comparing saturation and normal, reliable vs unreliable, dissimilar vs similar
  res2 <- campCompareSatNorm(hallem2)
  res3 <- res2[[1]][which(res2[[1]][,2] < 0.15),]
  res4 <- res2[[1]][which(res2[[1]][,2] > 0.3),]
  fstripchartvecs(list(res4[,5],res4[,8]),markersize = 0.8,ylabel = 'useful D.',semthick = 2) #reliable
  fstripchartvecs(list(res4[,4],res4[,7]),markersize = 0.8,ylabel = 'useful D.',semthick = 2) #unreliable 
  #comparing hallem and rob correlations for this set
  res2 <- campCompareSatNorm(hallem2) #correlation hallem, rob : 0.497
  fploteq(res2[[1]][,1],res2[[1]][,2],ticknox = 10,ticknoy = 5,rndfact = 1,fixx = c(0,.3),fixy = c(0,.3),xlabel = 'hallem',ylabel = 'rob')
  res2 <- campCompareSatNorm(hallem2,matchlevel = 0.8) #correlation hallem, rob : 0.816
  fploteq(res2[[1]][,1],res2[[1]][,2],ticknox = 2,ticknoy = 5,rndfact = 1,xlabel = 'hallem',ylabel = 'rob')
  
  #generic function for traversing all the directories, used here to get all the dataset sdetails.
  tst <- campAllDirDiscScoresGen(searchstr = c('*tseries*','*sigpval*'),mode = 2,targetfn = campGetTrialDetails,targetpars = list())
  
  #set the working directory first
  setwd("/nadata/cnl/data/shyam/campbell/data_dir/camptest/data_dir")
  #traversing through the directories and getting the campCompareSatNOrm data, #gotox
  tmp2 <- campAllDirDiscScores(searchstr = c('*tseries*','*sigpval*'),mode = 2,targetfn = campCompareSatNorm,targetpars = list(hallem2,halrobmatch=1,matchlevel=0.2))
  res.alldirsflySN <- tmp2
  tmp3 <- tmp2[which(sapply(tmp2, function(x) is.null(x[[1]]))==F)] #just the data for the valid directories
  tmp4 <- lapply(tmp3,'[[', 1)
  #tmp5 <- convertNestedListsDF(tmp4)
  tmp5 <- joinListDFs(lstdfs = tmp4)
  res <- campProcessDirSets(tmp3,op=1) #get the first results as a data frame
  res1 <- campProcessDirSets(tmp3,op=2)
  fstripchartvecs(res1[[1]][1:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 0)#dissimilar
  fstripchartvecs(res1[[2]][1:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 0)#midsimmilar
  fstripchartvecs(res1[[3]][1:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 0)#similar
  fstripchartvecs(res1[[1]][1:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,methodstr = 'overplot')#dissimilar
  fstripchartvecs(res1[[3]][1:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,methodstr = 'overplot')#similar
  fstripchartvecs(res1[[2]][1:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,methodstr = 'overplot')#midsimilar
  #plot similarity of hallem vs rob, corr, 0.86
  fploteq(tmp5[,1],tmp5[,2],ticknox = 5,ticknoy = 2,rndfact = 10,fixx = c(-.2,.8),fixy = c(-.1,.7),xlabel = 'hallem',ylabel = 'rob')
  #plotting slope change for the 3 groups
  res2 <- campProcessDirSets(tmp3,op=3)
  fstripchartvecs(res2[[1]],markersize = 0.8,ylabel = 'slope change.',semthick = 2,tickno=5)
  fstripchartvecs(res2[[2]],markersize = 0.8,ylabel = 'slope change.',semthick = 2,pairplot = 1,methodstr = 'overplot') # all three groups separately
  
  
  #testing reliability vs overlap
  res5 <- campComputeOverlap(tmp)
  fstripchartvecs(lapply(res5,function(x) x[,2]),markersize = 0.8,ylabel = 'overlap',semthick = 2)
  
  #do odor similarity, and then get cells
  res <- campMatchHalRobCor(hallem2,optype = 9,matchlevel = .7)
    
  
  #testing how no of trials influences no of unreliable cells and showing that 
  res <- campProcessDirSets(tmp3,op=1) #get the first results as a data frame
  plotLinearFit(res$trials,res$ur.cells,graphparams = list(xlabel='# trials',ylabel='# unrel. cells'))
  tst <- meanDF(cbind.data.frame(res$trials,res$ur.cells))
  #plotLinearFit(tst[,1],tst[,2],graphparams = list(xlabel='# trials',ylabel='# unrel. cells',ticknoy=5))
  plotLinearFit(tst[,1],tst[,2],graphparams = list(xlabel='# trials',ylabel='# unrel. cells',ticknoy=10,fixy = c(0,60)))
  
  plotLinearFit(res$trials,res$ursum.nor,graphparams = list(xlabel='# trials',ylabel='# unrel. sum'))
  tst <- meanDF(cbind.data.frame(res$trials,res$ursum.nor))
  #both no of unreliable cells and their overlap increases with no of trials
  #plotLinearFit(tst[,1],tst[,2],graphparams = list(xlabel='# trials',ylabel='# unrel. sum',ticknoy=5))
  plotLinearFit(tst[,1],tst[,2],graphparams = list(xlabel='# trials',ylabel='# unrel. sum',ticknoy=10,fixy=c(-1,4)))
  plotLinearFit(res$trials,res$urel.act,graphparams = list(xlabel='# trials',ylabel='# unrel. overlap',ticknoy = 5))
  #result: as the nymber of trials increases, the number of unreliable cells increases, but their overall contribution stays the same
  #because the mean unreliable cell contribution goes down because they are more likely to overlap with more cells.
  
  
  #cutoff calculations and change in the cutoff 
  tst <- campCalcCutoffRange(res,ptsrange = 30,cutrange = seq(0,10,0.5),cutop = 2)
  #change in cutoff from normal to saturated, similar 
  fploteq(tst[,1],tst[,c(7)]-tst[,c(4)],ticknox = 5,markersize = 0.8,xlabel = 'cutoff',ylabel = 'inc. odors (dis)',ticknoy = 1)
  
  #cutoff reliable vs all cells
  tst1 <- campCalcCutoffRange(res,ptsrange = 30,cutrange = seq(0,10,0.5),cutop = 1,op=1)
  c(green='reliable cells, normal',blue='reliable cells, saturated',black='all cells, normal',red='all cells, saturated')
  fploteq(tst1[,1],tst1[,2:5],markersize = 0.8,xlabel = 'cutoff',ylabel = '# odors pairs',ticknox = 3)
  tst1 <- campCalcCutoffRange(res,ptsrange = 30,cutrange = seq(0,7,0.5),cutop = 1,op=1)
  fploteq(tst1[,1],100*(tst1[,c(3)]-tst1[,c(5)])/tst1[,c(5)],markersize = 0.8,xlabel = 'cutoff',ylabel = '% improvement',ticknox = 3,ticknoy = 1)
  
  
  #habituation
  #no of times each trial is the highest response : flies
  tst1 <- campHabituation(tmp,op=1)
  tst2 <- getDensDfStats(transposeDF(convertNestedListsDF(tst1)),start = 1)
  fploteq(1:6,unlist(tst2[1,]),ticknox = 1,fixx = c(1,6),stddev = unlist(tst2[3,]),fixy = c(0,16),ticknoy = 2,xlabel = 'trial #',ylabel = '# peaks',markersize = 1.2)
  #mouse
  tst1 <- campHabituation(tst,op=1)
  
  #looking whether response go down after first significant response
  tst1 <- campHabituation(tmp,op=2)
  tst3 <- getDensDfStats(transposeDF(convertNestedListsDF(tst1)),start = 1)
  fploteq(1:2,unlist(tst3[1,]),ticknox = 10,fixx = c(1,2),stddev = unlist(tst3[3,]),fixy = c(0,40),ticknoy = 2,xlabel = 'first trial peak',ylabel = '# peaks')
  #look at the first and second significant response
  tst <- campHabituation(tmp,op=3)
  tst1 <- joinUnevenDFs(tst) #join all of the odors
  tst2 <- tst1[which(tst1[,1]>0),] #only get the ones with at least one sig response
  tst3 <- tst2[which(tst2[,2]>0),] #get ones with both significant responses
  fstripchartvecs(list(tst3[,1],tst3[,2]),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 1,methodstr = 'overplot')
  fstripchartvecs(list(tst4[,1],tst4[,2]),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 0,methodstr = 'jitter') #without lines goign across
  tst4 <- tst3[which(tst3[,1]<0.5),] #blow up responses from 0 to 0.5
  fstripchartvecs(list(tst4[,1],tst4[,2]),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 1,methodstr = 'overplot')
  
  #plot by number of sigresponses
  tst <- campHabituation(tmp,op=3)
  tst1 <- joinUnevenDFs(tst) #join all of the odors
  tst2 <- getRowsDFNonzeros(tst1,n=4,op=3)
  fstripchartvecs(as.list(data.frame(tst2)),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot')
  
  tmp1 <- joinUnevenDFs(campHabituation(getCampDataSigList(),op=3))
  tst3 <- getRowsDFNonzerosAll(tmp1,op=3) #tst3 gives a list where each element is a df of 1,2,3...etc sig. responses
  #e.g., plotting 3 sig responses
  fstripchartvecs(as.list(data.frame(tst3[[3]])),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot')
  
  #looking at all responses across odors
  res <- campHabituation(getCampDataSigList(),op=4)
  res2 <- joinUnevenDFs(res)
  rownames(res2) <- 1:nrow(res2)
  res1 <- getRowsDFNonzeros(res2[,1:4],n=1,op=2) #get no of cells with 1 response in first 4 trials
  nrow(getRowsDFNonzeros(res2[rownames(res1),],n=2,op=1)) #of these show cells with at least 2 responses
  
  #seeing if cells with a singl e response in the first 3 trials are more biased to respond again depending on the trial number
  res <- campHabituation(getCampDataSigList(),op=4)
  res2 <- joinUnevenDFs(res)
  rownames(res2) <- 1:nrow(res2)
  res1 <- getRowsDFNonzeros(res2[,1:3],n=1,op=2) #get all cells with one response exactly in the first 3 trials
  res3 <- res2[rownames(res1),] #get these cells specifically 
  res4 <- res3[which(res3[,1]>0),] #number that respond in first trial
  nrow(res4) #whats this number
  nrow(getRowsDFNonzeros(res4,n=2,op=1)) #how many respond again, repeat by changing 1 to 2 and then 3
  res5 <- getRowsDFNonzeros(res3,n=1,cols = 2,op=5) #get the responses for the cells with one response in the first 3 trials
  fstripchartvecs(as.list(data.frame(res5.1[which(res5.1[,2]==0),])),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot') #single freq
  fstripchartvecs(as.list(data.frame(res5.1[which(res5.1[,2]>0),])),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot') #more than one
  #comparing respones between the first and second significant responses
  temp <- getRowsDFNonzeros(res2,n=2,op=3)
  names(temp) <- c('T1','T2')
  fstripchartvecs(as.list(data.frame(temp)),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot') #pairplot
  fstripchartvecs(as.list(data.frame(temp)),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 0) #no connecting lines
  
  #comparing peak responses across all trials, only significant responses
  temp1 <- as.list(data.frame(res2))
  names(temp1) <- 1:lenght(temp1)
  fstripchartvecs(as.list(data.frame(res2)),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 0,fixy = c(0,0.4)) #truncate at 0.4
  #comparing peak responses across all trials, only significant responses
  temp <- lapply(1:ncol(res2), function(x) getAboveThresh(res2[,x]))
  names(temp) <- 1:length(temp)
  fstripchartvecs(temp,markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 0)
  
  #no of significant responses per trial
  temp2 <- sapply(1:ncol(res2), function(x) length(which(res2[,x]>0)) )
  fploteq(1:length(temp2),temp2,xlabel = 'trial #',ylabel = '# sig. responses',fixx = c(1,6))
  
  #testing if no of trials vs resp is correlated, avg. cells for each set of trials
  temp <- campTrialNoVsResp(getCampDataSigList())
  temp1 <- t(temp1)
  temp2 <- temp1[which(temp1[,2]!=0),] #get the trials that have that many cells
  cor(temp2[,1],temp2[,2])
  NlsFitCDF(temp2[,2],dist =2) #good fit for a gamma
  #if you want to do all cells for all trials
  temp <- campTrialNoVsResp(getCampDataSigList(),op=2)
  temp2 <- temp[which(temp[,2]!=0),]
  #no of responses per trial odorwise
  tst <- campHabituation(tmp,op=1)
  tst2 <- getDensDfStats(transposeDF(convertNestedListsDF(tst)),start = 1)
  #plotting reliability vs response levels
  temp2 <- getDfSplit(temp)
  fstripchartvecs(lapply(temp2,function(x) x[,2]),markersize = 0.8,ylabel = 'avg. response',semthick = 2,pairplot = 0,tickno = 4)
  
  #gives the number of cells you would expect to respond if their frequency was 1 in 3, 4, ...
  fploteq(3:12,floor(142-sapply(3:12, function(x) dbinom(x=c(0),size = 3,prob = 1/x) * 142)),xlabel = 'frequency of response',ylabel = '# responsive cells')
  
  #getting the euclidean statistics that Glenn had requested
  tmp1 <- campGetTrialDistances(tmp,op=3,rettype = 2)
  fstripchartvecs(as.list(tmp1),markersize = 0.8,ylabel = 'trial distance',semthick = 2,pairplot = 2,tickno = 4,methodstr = 'overplot')
  tmp1 <- campGetTrialDistances(tmp,op=1,rettype = 2)
  fstripchartvecs(as.list(tmp1),markersize = 0.8,ylabel = 'origin distance',semthick = 2,pairplot = 2,tickno = 1,methodstr = 'overplot')
  
  
}

campOverlapAllPairAnalog.old <- function(odor1,odor2,freq,resp,notrials,nolevels,cells=1,
                                         relscore=1,op=1){
  vec1 <- freq[[odor1]]
  vec2 <- freq[[odor2]]
  if(op==11 || op==13 || op==14){#the vector is reliability * avg. response rate
    vec1 <- vec1*resp[[odor1]]
    vec2 <- vec2*resp[[odor2]]
  }
  if(isDataType(cells)==4){#if the cells list exists, get only the selected cells
    allcells <- union(cells[[1]],cells[[2]])
    #cat('\ncOAPA',allcells)
    vec1 <- vec1[allcells]
    vec2 <- vec2[allcells]
  }
  #cat(vec1[1:20],'\n',vec2[1:20],'\n')
  #calculate the reliability score
  if(relscore==1) rel <- sapply(1:length(vec1), function(x){
    floor(((vec1[x]+vec2[x])/2)*10*nolevels)/(10*nolevels)
  })
  else rel <- sapply(1:length(vec1), function(x) vec1[x])
  #cat('\ncampoverlao',notrials)
  #make sure you pass only those notrials for the two pairs of odors
  over <- campOverlapFunction(vec1,vec2,notrials[c(odor1,odor2)],cells=cells,op=op)
  #cat('\nover return: ',str(over))
  if(op==13 || op ==14) {
    res.df <- over #population level computation, only one number
  }
  else res.df <- cbind(rel,over) #make the results into a matrix
  #cat('\nover1',str(res.df) )
  #and sort the matrix by relliability, and sum all the overlaps with the same reliability
  switch(op,sumDF(getThreshDF(res.df,op=1)),sumDF(getThreshDF(res.df,op=1)),
         sumDF(getThreshDF(res.df,op=1)),getThreshDF(res.df,op=1),getThreshDF(res.df,op=1),getThreshDF(res.df,op=1),
         res.df[which(res.df[,2]>0),],res.df[which(res.df[,2]!=0),],getThreshDF(res.df,op=1),getThreshDF(res.df,op=1),
         getThreshDF(res.df,op=1),getThreshDF(res.df,thresh = -1,op=1),res.df,res.df)
}


#see defintions below
campOverlapAllAnalog.old <- function(data,nolevels=2,skip=c('paraffin'),relscore=2,op=1){
  res.sig <- getGroupedListsDf(data,op=1) #gets the significant response Dfs
  notrials <- campGetTrialDetails()[[1]] #no of trials, assumes that all odors have same # of Trials
  #have to calculate the freq. list for all the odors, which is a single list
  freq.lst <- lapply(res.sig, function(x){
    #first the freq cell list
    freq.cells <- apply(x, 1, sum)
  })
  odors <- getSkipStrIndex(names(freq.lst),skip = skip,op = 2) #get odors without skip
  pairs <- convertMatVec(length(odors)) #get the pairing of odors 
  res.ls <- lapply(pairs, function(x){#calculate overlap/difference for all odorpairs
    campOverlapAllPairAnalog(odor1=odors[x[1]],odor2=odors[x[2]],freq=freq.lst,
                             notrials=notrials,nolevels=nolevels,relscore = relscore,op=op)
  })
  #just get the odor pair names
  res.names <- sapply(pairs, function(x){
    paste(odors[x[1]],odors[x[2]],sep = ',')
  })
  names(res.ls) <- res.names
  res.ls
}



campGetOverlapTopNPop.old <- function(data.lst,thresh.win=c(0,1),topn=0,deciles=c(),symop=1,overop=9,sameop=1,op=1){
  sig.topn <- campGetTopNCells(data.lst = data.lst,topn = topn,deciles = deciles,op=2) #get the topN neurons for all odors
  #get the overlap for all cells of each odor pair
  #overlap.allcells <- campComputeOverlap(datasig = data.lst,overop = overop,op=2)
  overlap.allcells <- campOverlapAllAnalog(data = data.lst,skip = c(),op = overop,trialop = 2,cells = sig.topn)
  #cat('\n',str(overlap.allcells),'\n',str(sig.topn))
  #Algo:what we have to do is for each odor pair, designate the first of the pair as the primary odor, then
  #get the topn neurons for that odor, then go back to the overlap list entry and only get those cells and their overlap.
  #now, get the odor pairs whose correlations are within the thresh.win, and get the overlaps of the cells only for these pairs
  odor.corr <- computeListMatCorr(alldata = data.lst,celltype=4) #correlations of all odor pairs
  odor.corr <- flattenDFtoList(odor.corr)
  odorpairs <- names(odor.corr)[which(odor.corr>=thresh.win[1] & odor.corr<=thresh.win[2])] #names of odor pairs that satisfy the criterion
  overlap.target <- overlap.allcells[odorpairs] #the overlap numbers for the targetted population
  # res <- lapply(names(overlap.target),function(x){
  #   #choose the appropriate top N cells from one or both odors, based on symop
  #   odors <- strsplit(x=x,split = ',')[[1]]
  #   if(overop==13 || overop == 14 || overop == 15) {
  #     odornames <- union(names(sig.topn[odors[1]][[1]]),names(sig.topn[odors[2]][[1]]))
  #     cbind.data.frame(length(odornames),overlap.target[[x]])
  #   }
  #   else{
  #     odornames <- names(sig.topn[odors[1]][[1]])
  #     sumres <- overlap.target[[x]][odornames,] #get the particular cells specific in sig.topn by using their names
  #     #cat('\n',x,mean(sumres[,1]),sum(sumres[,2]) )
  #     cbind.data.frame(mean(sumres[,1]),sum(sumres[,2]))
  #   }
  # })
  #cat('\n',str(res),str(overlap.target))
  res <- overlap.target
  #res <- overlap.target
  names(res) <- names(overlap.target)
  if(sameop==2){#take out similar odor pairs
    res <- res[getIndexSameSubStr(names(res),op=2)] #op=2 gets strings with dissimilar substrings
  } 
  #cat('\nCGOY',str(res))
  res.df <- joinListDFs(res)
  #list(names(sig.topn),names(overlap.allcells),overlap.allcells[odorpairs],res,res.df )
  res <- switch(op,res,res.df,res.df[,2])
}


#given a data set will generate a lst of dataframes, where each data frame is just the values of the 
#reliable or unreliable cells
#op=1, significant cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells, 5 - specified cells
#headerop: how we get header information, 1 = look through the folder, 2 - look through alldata
campGenRUCellMat.old <- function(alldata,headerop=1,op=1){
  dat.ls <- getGroupedListsDf(alldata,op=2) #get response values
  trials <- campGetTrialDetails(alldata = alldata,op=2)[[1]]
  if(op==4) return(dat.ls)
  #get reliable and unreliable cells
  cell.class <- campClassifyCells(alldata,op=2) #defailt: get it through alldata
  if(headerop==1) cell.class <- campClassifyCells(alldata) 
  #all sign cells = rel + unrel cells
  cell.class[[1]] <- lapply(1:length(cell.class[[1]]),function(i) c(cell.class[[2]][[i]],cell.class[[3]][[i]]) )
  veclen <- length(alldata[[1]][[1]]) # no of cells
  nocols <- sapply(1:length(dat.ls),function(i) ncol(dat.ls[[i]])) #no of trials for every odor
  #cat('\n',sapply(1:length(dat.ls),function(i) ncol(dat.ls[[i]])))
  #go through the response matrix for each odor, and get only either the reliable or unreliable cells
  cell.mat <- lapply(1:length(cell.class[[op]]),function(i) {
    vec <- genVecOnes(veclen,cell.class[[op]][[i]]) #get a vector of ones in posns of (un)reliable cells
    mat <- genVecCopies(vec,nocols[i]) #make a matrix
    #cat('\n',str(mat),str(dat.ls[[i]]))
    as.data.frame(filterWithMat(dat.ls[[i]],mat) ) # a matrix of 
  })
  names(cell.mat) <- names(dat.ls)
  cell.mat
}

campGenRUCellMat1.old <- function(alldata,headerop=1,retype=1,op=1){
  dat.ls <- getGroupedListsDf(alldata,op=2) #get response values
  rel.ls <- getGroupedListsDf(alldata,op=1) #get reliability values
  trials <- campGetTrialDetails(alldata = alldata,op=2)[[1]]
  #get reliability and response rates
  res <- campGetReliabResp(data.lst = alldata) #first get the data frames of reliability and response rates
  #now, get the posns for the options specified
  posns.lst <- switch(op,
                      lapply(names(dat.ls),function(i){#=1, significant cells
                        posns <- which(res[[1]][i]>0)
                      }),
                      lapply(names(dat.ls),function(i){#2 - reliable cells
                        posns <- which(res[[1]][i]>(trials[i]/2))
                      }),
                      lapply(names(dat.ls),function(i){#3 - unreliable cells
                        posns <- which(res[[1]][i]<=(trials[i]/2) & res[[1]][i]>0)
                      }),
                      lapply(names(dat.ls),function(i){#4 - all cells
                        posns <- 1:nrow(dat.ls[[i]])
                      }),
                      lapply(names(dat.ls),function(i){#5 - specified cells
                        posns <- cells[[i]]
                      }))
  names(posns.lst) <- names(dat.ls) #name results appropriately
  cells <- lapply(names(posns.lst),function(i){#now pick all the cells for each trial of every odor
    trial.cells <- lapply(1:trials[i],function(j) {
      tmp <- dat.ls[[i]][,j]*rel.ls[[i]][,j]
      resp <- tmp[posns.lst[[i]] ]
      #cat('\n',str(tmp),str(res[[2]][i]) )
      #resp <- setMatNonzero(unlist(res[[2]][i])[posns.lst[[i]] ])*tmp #have to filter out the non-significant cells
      #if(i=='2-heptanone') cat('\n',tmp,'\n',resp,'\n',setMatNonzero(unlist(res[[2]][i])[posns.lst[[i]] ]))
      names(resp) <- posns.lst[[i]] #the names of each vector carry that position
      resp
    })
    trial.cells <- switch(retype,convertNestedListsDF(trial.cells),trial.cells)
  })
  names(cells) <- names(dat.ls)
  cells
}

#function that gets the correlation between a list of matrices
#dat.ls: the list of matrices or data frames
#sel=1, correlation between every vector column across the two matrices and 
#the number is their average, sel=2, correlation between selected vectors in each matrix,
#selected col specified in by 'col'
#op: 1 - all sig. cells, 2 - reliable cells, 3 - unreliable cells, 4 - all cells
computeListMatCorr.old <- function(alldata,sel=1,col=1,op=1){
  #dat.ls <- getGroupedListsDf(alldata,op=2)
  dat.ls <- campGenRUCellMat.old(alldata,op=op) #select the kind of cells we want to get correlations for
  if(length(dat.ls)<2) return(NULL) #if there is only one or no odors, there is nothing to compare it with
  res.ls <- sapply(1:length(dat.ls), function(x){
    #cat('\n x',x)
    res <- sapply(1:length(dat.ls), function(y){
      #cat('\n',x,y,str(dat.ls[[x]]),str(dat.ls[[y]]))
      if(sum(dat.ls[[x]])==0 || sum(dat.ls[[y]])==0) allcor <- c(0,0)
      else  allcor <- switch(sel,cor(dat.ls[[x]],dat.ls[[y]]),
                             cor(dat.ls[[x]][,col],dat.ls[[y]][,col]))
      #cat('\tcLMC',x,';',y,':',allcor,str(allcor) )
      mean(allcor)
    })
  })
  #cat('\ncLMC',str(res.ls),str(dat.ls))
  rownames(res.ls) <- names(dat.ls)
  colnames(res.ls) <- names(dat.ls)
  res.ls
}

#compares (the statistics of) the overlap between different percentile ranges
#poutput: the number of overlapping cells that don't include similar elements
#matop:
#params:
#retype:
#op:1 - the no of overlapping cells, 2 - dot product of overlapping cells, 3 - relative dot product, i.e., dot product/norm(a)*norm(b)
# 4 - norm product for the percentile, 5 - the non-overlap between the two odors, 6 - percentage overlap, 7 - percent nonoverlap
# 8 - continuous overlap/non-overlap
#trialop is 2 for op=2,3,4.
#trialop: 1 get the result by computing it across all trial combinations, and averaging, leave out same-trial pairs
#trialop: 2 get the result by computing result for averaged trial response
campCompPercentRanges.old <- function(alldata,cells=c(),matop=2,params=c(10,10),trialop=1,op=1){
  dat.ls <- campGenRUCellMat(alldata,op=5,matop = 1,params = params,retype = 3) #select the kind of cells we want to get correlations for
  #cat('\n',str(dat.ls))
  if(length(dat.ls)<2) return(NULL) #if there is only one or no odors, there is nothing to compare it with
  notrials <- campGetTrialDetails(alldata = alldata,op=2)[[1]] #get no of trials for each odor
  stim.norms <- sapply(names(dat.ls[[2]][[2]]),function(x) lpnorm(unlist(dat.ls[[2]][[2]][x])) )
  #iterate through odor pairs, get the number of cells in common, and also get their dot. product.
  res.ls <- sapply(names(dat.ls[[1]]),function(i){
    res.in <- sapply(names(dat.ls[[1]]),function(j){
      #get number of cells in common
      stim1 <- (rownames(dat.ls[[1]][[i]])) #cells for stim 1
      stim2 <- rownames(dat.ls[[1]][[j]]) #cells for stim2
      total <- union(stim1,stim2)
      overlap <- intersect(stim1,stim2) #cells in common
      nonoverlap <- setdiff(union(stim1,stim2),intersect(stim1,stim2)) #cells that are not in common
      dotprod <- sum(unlist(dat.ls[[1]][[i]][overlap,])*unlist(dat.ls[[1]][[j]][overlap,]))
      #cat('\n',i,j,':',unlist(dat.ls[[1]][[i]][overlap,]),',',unlist(dat.ls[[1]][[j]][overlap,]),':',overlap,dotprod,'\n')
      normprod <- lpnorm(unlist(dat.ls[[1]][[i]][stim1,]))*lpnorm(unlist(dat.ls[[1]][[j]][stim2,]))
      comm.over <- sum(min(dat.ls[[1]][[i]][overlap,],dat.ls[[1]][[j]][overlap,]) )
      noncomm.over <- sum(abs(dat.ls[[1]][[i]][overlap,]-dat.ls[[1]][[j]][overlap,]) )
      if(length(nonoverlap)==0 || i==j) {
        comm.over <- 1
        noncommon <- 1
      }
      else  noncommon <- sum(dat.ls[[1]][[i]][setdiff(stim1,overlap),]) + sum(dat.ls[[1]][[j]][setdiff(stim2,overlap),])
      #cat('\n',i,j,comm.over,noncomm.over,noncommon,setdiff(stim1,overlap),':',setdiff(stim2,overlap) )
      switch(op,length(overlap),dotprod,dotprod/as.numeric(stim.norms[i]*stim.norms[j]),normprod,length(nonoverlap),
             length(overlap)/length(total),length(nonoverlap)/length(total),comm.over/(noncomm.over+noncommon))
    })
  })
  res.no <- sum(res.ls[upper.tri(res.ls,diag = F)])
  list(res.ls,res.no)
}

campComparePercentile.old <- function(dat1.df,dat2.df,dat1.avg,dat2.avg,trialop=1,op=1){
  #get number of cells in common
  stim1 <- rownames(dat1.df) #cells for stim 1
  stim2 <- rownames(dat2.df) #cells for stim2
  total <- length(union(stim1,stim2)) #total number of cells from both stimuli
  overlap <- intersect(stim1,stim2) #cells in common
  #cat('\toverlap',length(overlap),'total:',total,'nonoverlap',length(union(setdiff(stim1,overlap),setdiff(stim2,overlap))) )
  nonoverlap <- setdiff(union(stim1,stim2),intersect(stim1,stim2)) #cells that are not in common
  if(op <= 9 ) choice <- switch(op,1,2,3,4,2,3,4,5,6)
  else choice <- op
  #all these have to be done over all trials, and then averaged or done on the average of all trials
  if(trialop==1 && !(op %in% c(2:4)) ){#second condn: for op=2:4, we look at trial averaged response only
    res.ls <- sapply(1:ncol(dat1.df),function(i){
      res.en <- sapply(1:ncol(dat2.df),function(j){
        #cat('\n',str(dat1.df[setdiff(stim1,overlap),i]),':',stim1,':',overlap,':',setdiff(stim1,overlap),'\n',choice )
        if(choice==3 || choice ==4){
          total <- length(union(names(dat1.df)[which(dat1.df[,i]>0)],names(dat2.df)[which(dat2.df[,j]>0)]))
        }
        tmp <- switch(choice,
                      countNonZeroes(dat1.df[overlap,i]*dat2.df[overlap,j]), #gets the overlap of cell between 2 trials       
                      #the cells of stim1 not overlapping and stim 2 not overlapping
                      countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i]),
                      countNonZeroes(dat1.df[overlap,i]*dat2.df[overlap,j])/total, #gets the overlap of cell between 2 trials       
                      (countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i]))/total
        )
        if(choice==5){#continuous overlap
          stim1 <- rownames(dat1.df)[which(dat1.df[,i]>0)]
          stim2 <- rownames(dat2.df)[which(dat2.df[,j]>0)]
          common <- intersect(stim1,stim2)
          common.over <- sum(unlist(sapply(common,function(s) min(dat1.df[s,i],dat2.df[s,j]) ) ) )
          common.nonover <- sum(abs(dat1.df[common,i]-dat2.df[common,j]) )
          noncommon <- sum(dat1.df[setdiff(stim1,stim2),i]) + sum(dat2.df[setdiff(stim2,stim1),j])
          if(identical(dat1.df,dat2.df) ) tmp <- 1
          else tmp <- common.over/(common.nonover+noncommon)  
        }
        if(choice==6) {
          tmp <- countNonZeroes(dat1.df[overlap,i]*dat2.df[overlap,j])/(countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i]))
          #cat('\t',countNonZeroes(dat1.df[setdiff(stim1,overlap),i]),countNonZeroes(dat2.df[setdiff(stim2,overlap),i]) )
        }
        if(choice==100) {#the test condition
          tmp <- total
        }
        if(choice==10) {#total
          tmp <-  length(union(names(dat1.df)[which(dat1.df[,i]>0)],names(dat2.df)[which(dat2.df[,j]>0)])) 
        }
        tmp
      })
      #cat('\n',(countNonZeroes(dat1.df[setdiff(stim1,overlap),i]) + countNonZeroes(dat2.df[setdiff(stim2,overlap),i])) )
      res.en #overlap of dat1 col i with dat2 all cols
    })
    res.mn <- mean(cleanNAVec(res.ls))
  }
  if(trialop==2){#we are looking at average response comparisons
    res.mn <- switch(choice,
                     length(overlap),
                     length(c(setdiff(stim1,overlap),setdiff(stim2,overlap))),
                     length(overlap)/total,
                     length(c(setdiff(stim1,overlap),setdiff(stim2,overlap)))/total
    )
    if(choice==5){
      common <- intersect(stim1,stim2)
      #cat('\ncommon',common,':',stim1,':',stim2,str(dat1.avg) )
      common.over <- sum(unlist(sapply(common,function(s) min(dat1.avg[s,1],dat2.avg[s,1]) ) ) )
      common.nonover <- sum(abs(dat1.avg[common,1]-dat2.avg[common,1]) )
      noncommon <- sum(dat1.avg[setdiff(stim1,stim2),1]) + sum(dat2.avg[setdiff(stim2,stim1),1])
      if(identical(dat1.df,dat2.df) ) res.mn <- 1
      else res.mn <- common.over/(common.nonover+noncommon) 
      #cat('\ncommon',common,':',common.over,':',common.nonover,':',noncommon )
    }
    if(choice==6) 
      #cat('\n',length(c(setdiff(stim1,overlap),setdiff(stim2,overlap))) )
      res.mn <- length(overlap)/length(c(setdiff(stim1,overlap),setdiff(stim2,overlap)))
    if(choice==100) {#the test condition
      tmp <- total
    }
  }
  #cat('\tres:',res.mn)
  res.mn
}

#matches the Hallem and Rob correlations and gets the useful difference for those odor pairs that are within 0.1 correlations
#according to both Hallem and Rob
#hallem: the hallem data set
#rob: correlations between Rob's odors
#odors: the odor names for doing correlations
#optype: 1 - overlap; 2, difference, 3 - probability difference, 
#hallrobmatch: how to balance the two sets: 1 - only odors present in both, no similar odor pairs 
#2 - only odors present in both, with similar odor pairs,
#3 -rob odors, similar odor pairs ok, 4 - rob odors, similar odor pairs ok
#4 - probability difference, all cells, 5 - probability difference, all cells w/ saturation
#also in all cases i
#op=1, the sum of reliability cells of the same type, 2 - mean of all reliability cells of the same type, 
#todo: stuff to weed out empty and mismatched and same odor pairs
campMatchHalRobCor.old <- function(hallem,skip=c('paraff','empty'),matchlevel = 0.1,halrobmatch=1,optype=4,op=1){
  #getall the data, and then rob's odor correlations
  alldata <- getCampDataSigList() #gets the original data
  notrials <- campGetTrialDetails()[[1]] #get no of trials
  #get the Rob correlations
  rob <- computeListMatCorr(getGroupedListsDf(alldata,matop=2))
  #get the hallem correlations for the odors from rob
  hallemdf <- campGetHallemCorr(hallem,rownames(rob))
  #print(hallemdf)
  halrob.match <- 0 #there is no match
  if(length(hallemdf) != 0) {#there is overlap in odors between hallem and rob  
    odors <- intersect(rownames(rob),rownames(hallemdf))
    #since this is a comparison, these are the only odors we are interested in, will also take out similar odor-pairs
    res <- pickMatchedVals(hallemdf,rob[odors,odors],prec = matchlevel)
    if (length(res)>0) {#there is a match so we are good.
      res <- t(res) #transpose so that the odors are in rows
      colnames(res) <- c('match','hallem','rob')
      halrob.match <- 1 #there is a match!
    }
  }  
  if (halrob.match==0){#there is no overlap, so do correlation for only rob's odors
    res <- cbind.data.frame(unlist(flattenDFtoList(rob))) #make it into a DF
    rownames(res) <- names(flattenDFtoList(rob))
    colnames(res) <- c('rob')
  }
  overlap.cols <- campOverlapAnalogRel(alldata = alldata,notrials=notrials,skip = skip,relscore = 2,optype = optype,op=op)
  if(hallrobmatch<=2){#take out similar odor pairs
    odors <- rownames(overlap.cols)[which(sapply(rownames(overlap.cols), function(x) areSubStringsSame(x)) == F)]
    odors <- getSkipStrIndex(odors,skip=skip,op=2) #odors without skip
  }
  cat(str(res))
  print(res[odors,])
  #now, we have to match Rob and Hallem according to the conditions
  #res <- cbind.data.frame(res[rownames(overlap.cols),],overlap.cols)
  res <- cbind.data.frame(res[odors,],overlap.cols[odors,])
  res
}

