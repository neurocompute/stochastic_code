
#simon_data_analysis.R: code for data analysis: Shyam Srinivasan (C), shyam@snl.salk.edu
#version 1.1 contains code for a variable number of valences, i.e., MBONs and DANs. In this the number of MBONs and DANs should match.
#this file contains functions for analyzing simons data

#have to first construct the equivalent of the getCampDataSigList function. Once that is done, we are done, we can just 
#call all the campbell analysis functions after that

#this functio is the equivalent of getCampDataSigList: basically reads the data in, and then for each ell for each odor determines the response
#value and computes if it's over background.
#returns a a list where each element contains the following vectors, the sigcells, the response values, and basemean,base sd, and z-scores
#foldname: the current folder or the folder that contains the Matlab generated data files
#fpatterns: the patterns for the sigpvals file
#alpha is used for calculating the SDs above mean where the significance should be fixed
#meanmax and sigval and ranking too are not needed here
#op: types of return functions
#dirop: whether the foldname is op=1, '.', 2 : relative, or 3: absolute e.g., /home/shyam....
simGetDataSigList <- function(foldname='.',fpatterns=c('.*series.*data*'),alpha=0.01,odorname=T,op=1,
                              dirop=1){
  datafile <- list.files(path = foldname,pattern = fpatterns[1])
  if(foldname != '.' || dirop==2) {#we are not in the directory, so save the current one and set the path
    currentdir <- getwd()
    newdir <- paste(currentdir,'/',foldname,'/',sep = '') #go to the actual directory
    setwd(newdir)
  }
  dat.df <- read.csv(file = datafile)
  #get each trial separately, the last col. contains the trial number
  trials.lst <- getDfSplit(dat.df,factorcol = ncol(dat.df))
  #rearrange so that the lists are in the right order,i.e., 1 2 3 instead of 1 10 11
  trials.lst <- trials.lst[sortStringNos(names(trials.lst))]
  #check if the names are properly assigned
  #now, for each dataframe get the sig. cells, response values, basemean,base sd, sig. as z-score
  #cat('\nsimGEt',length(trials.lst))
  trialdata.ls <- lapply(trials.lst, function(i) {
    #cat('\n',str(i))
    simGetTrialDataStats(trial.df = i,stats = 1,alpha = alpha,op=4)
  }) 
  #make sure the that the names are appropriate and return
  odornames <- simReadTrialOdorMap()  
  names(trialdata.ls) <- odornames[,ifelse(odorname,3,2)] #the third col contains the names, 2nd contains numbers
  if(foldname != '.' || dirop==2)  setwd(currentdir) #set the directory back
  trialdata.ls
}




#this function takes in a data framme of a trial and returns the following as a list
# sig. cells, response values, basemean,base sd, sig. as z-score
#stats: 1 - mean, 2 - max
#op: 1- return signal, 2 - return base_mean, base_sd, thresh.base 3 - signal, base.mean, base.sd, thresh.base
#4: sig. cells, response values, basemean,base sd, sig. as z-score
simGetTrialDataStats <- function(trial.df=c(),foldname='.',fpatterns=c('.*series.*head.*csv'),stats=1,alpha=0.01,op=1){
  headerfile <- list.files(path = foldname,pattern = fpatterns[1])
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile,header = T,colClasses = "character")
  #cat('\nsGTDS',str(headerfile) )
  frame.size <- as.numeric(headerfile$fp[1])
  #cat('\nhere')
  odor.onset <- floor(as.numeric(headerfile$stimlat[1])*frame.size )
  back.frames <- c(1:(odor.onset-1) ) #calculate back and odor frames
  odor.frames <- odor.onset - 1 + c(1:ceiling(frame.size*(as.numeric(headerfile$stimdur[1]) + as.numeric(headerfile$offset[1]) + as.numeric(headerfile$extratime[1])) ) )
  #cat(back.frames,'\n',odor.frames,'\n')
  #each row is one time frame so we iterate through the cols for each cell
  backmean <- sapply(2:(ncol(trial.df)-2), function(i) mean(trial.df[back.frames,i]))
  backsd <- sapply(2:(ncol(trial.df)-2), function(i) sd(trial.df[back.frames,i]))
  nosds <- qnorm(1-alpha) 
  thresh.base <- backmean + nosds*backsd #the threshold for baseline
  #cat('\nthresh',thresh.base,backmean,backsd)
  resp <- switch(stats,sapply(2:(ncol(trial.df)-2), function(i) mean(trial.df[odor.frames,i])),
                 sapply(2:(ncol(trial.df)-2), function(i) max(trial.df[odor.frames,i])) )
  sig.cells <- setAboveThreshOne(resp - thresh.base,op=2)
  zscore <- (resp - backmean)/backsd
  #cat('\nthresh',thresh.base,'\n',resp,'\n',resp-thresh.base,'\n',(resp-thresh.base)*sig.cells,'\n')
  #mimic the returns of camGetTrialDataStats
  resp <- resp - thresh.base #have to subtract the threshold from the signal
  res <- switch (op,resp,list(backmean,backsd,thresh.base),list(resp,backmean,backsd,thresh.base),
                 list(sig.cells,resp,backmean,backsd,zscore))
  res
}

#this function given a data set, the return of getCAmpDataSigList, will give you two lists
#list(all significant cell response means, all significant cells base mean)
simGetDataSigBase <- function(dat,op=1){
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
    sig.zscore <- dat[[x]][[5]]
    #cat('\n',x,str(dat[[x]]))
    signals <- sig.signal[which(sig.resp==1)]
    basemeans <- sig.base[which(sig.resp==1)]
    sig.cells <- sig.resp[which(sig.resp==1)]
    sig.zscore <- sig.zscore[which(sig.resp==1)]
    list(signals+basemeans,basemeans,sig.cells,sig.zscore)
  })
  res.ls
}



#given a list of the odor responses, the return of simGEtDataSigList, will get the odor correlations
#op=1, correlation between every vector column across the two matrices and 
#the number is their average, op=2, correlation between selected vectors in each matrix
simGetOdorCorrelations <- function(data.lst,op=1){
  data.grp <- getGroupedListsDf(data.lst = data.lst,op=2) #all the response vals
  data.sig <- getGroupedListsDf(data.lst = data.lst,op=1) #all the sig cells
  #filter one data set with another
  sig.resp <- lapply(1:length(data.grp), function(i){
    res <- data.grp[[i]] * data.sig[[i]]
    names(res) <- names(data.grp[[i]])
    res
  })
  names(sig.resp) <- names(data.grp)
  #cat(str(sig.resp))
  res <- computeListMatCorr(sig.resp,op=op)
  res <- makeDfFromList(flattenDFtoList(res))
  res
}


#function that gets a mapping in terms of data frame, where the first column is the trial index
#and the second col is the odor number
#input: could be the name of a file or the dataframe of the data read from the csv file
#op:1 - vector, 2 - df
simGetTrialOdorMap <-function(foldname='.',fpatterns=c('series*_*data*'),odorname=T,op=1){
  if(is.character(foldname)){ #if input is a string it's a file name
    datafile <- list.files(path = foldname,pattern = fpatterns[1])
    dat.df <- read.csv(file = datafile)
  }
  else cat('\nCant get odor map')
  #get each trial separately, the last col. contains the trial number
  trials.lst <- getDfSplit(dat.df,factorcol = ncol(dat.df))
  #rearrange so that the lists are in the right order,i.e., 1 2 3 instead of 1 10 11
  trials.lst <- trials.lst[sortStringNos(names(trials.lst))]
  #rearrange the names of the list so that they are odor bassed, skip for now, can do later
  odormap <- sapply(trials.lst, function(i) simGetTrialOdor(i,op=1))
  #get the odor names and put them in too
  odors <- read.csv('odors.csv',header = F)
  odor.names <- sapply(1:length(trials.lst), function(i) odors[odormap[i],1])
  odormap <- cbind.data.frame(1:length(trials.lst),odormap,odor.names)
  names(odormap) <- c('trial index','odor#','name')
  switch(op,odormap[,2],odormap)
}


#function that reads the odor mapping as in the trial number to odor number mapping
#and the second col is the odor number
#input: could be the name of a file or the dataframe of the data read from the csv file
#odorname: T - instead of numbers, substitute the odor name, F - leave the numbers
#op:1 - vector, 2 - df
simReadTrialOdorMap <-function(foldname='.',fpatterns=c('odorindex.*header*','series*_*data*'),odorname=T,op=1){
  headerfile <- list.files(path = foldname,pattern = fpatterns[1]) #read the header file
  #cat('\nsimReadT',headerfile)
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile,header = T,colClasses = "character")
  else headerfile <- sismGetTrialOdorMap() #if the odor map file does not exist, create it
  headerfile <- ConvertDfCols(headerfile,cols = c(1),op=1) #set the 1st column to numericsim
  if(odorname==T && ncol(headerfile)==2){ #if we want odor names and it is not alreaddy part of the file
    odors <- read.csv('odors.csv',header = F)
    odor.names <- sapply(1:nrow(headerfile), function(i) odors[headerfile[i,2],1])
    headerfile <- insertColDf(headerfile,newcol = odor.names,posn = 3)
  }
  headerfile
}


#this function gets the odor index from the trial data frame
#op: 1 - odor number, 2 - mapping c(trial index,odor index)
simGetTrialOdor <- function(trial.df,op=1){
  #the last col is the trial index and the second last column is odor no
  switch(op,trial.df[1,ncol(trial.df)-1],c(trial.df[1,ncol(trial.df)],trial.df[1,ncol(trial.df)-1]) )
}

#makes a header in the form of the Campbell header
#fpattern: contains the details for file with details of simon's data
#output: the name of the output header file
#op=1
simMakeHeader <- function(foldname='.',fpatterns=c('*params.*csv','series*_*data*'),output='series_header.csv',op=1){
  #get the values from the detail file
  paramsfile <- list.files(path = foldname,pattern = fpatterns[1])
  #cat('\nsimMakeH',paramsfile,str(paramsfile))
  if(length(paramsfile) == 1) paramsfile <- read.csv(file = paramsfile,header = F)
  #read in the odor trial data
  odornames <- simReadTrialOdorMap(foldname = foldname)  
  len <- nrow(odornames)
  stimlat <- rep(paramsfile[2,2],len) #odor onset=stimlat
  stimdur <- rep(paramsfile[3,2],len) #1s odor pulse
  extratime <- rep(0.5,len)
  offset.odor <- rep(1,len) #odor offset, 1s
  fp <- rep(paramsfile[1,2],len)  #frame rate, no of frames/sec
  #cat(str(odornames),len,'\n',stimlat,'\n',stimdur,'\n',extratime,'\n',offset.odor,'\n',fp)
  header <- cbind.data.frame(odornames[,2],stimlat,stimdur,extratime,offset.odor,fp)
  names(header) <- c('odor','stimlat','stimdur','extratime','offset','fp')
  #now get the prefix, e.g., mouse164, and write the file name
  datafile <- list.files(path = foldname,pattern = fpatterns[2]) #get the prefix to be attached
  output_name <- getPattern(string = datafile,pattern = '(m[a-z]+[0-9]+[_]+)+',postpattern = fpatterns[2])
  write.csv(x = header,file = paste(output_name,output,sep = ''),row.names = F)
  T
}

#makes a trial odor map file. Basically a column of mappings of odor index number to odor number
simMakeOdorMapFile <- function(foldname='.',fpatterns=c('series*.*data*'),output='trial_odorindex_header.csv',op=1){
  #get the odor map, and write it to a file so that you can save time
  odormap <- simGetTrialOdorMap(foldname = foldname,fpatterns = fpatterns,op=2)
  #read the datafile name so that you can get the mouse prefix
  datafile <- list.files(path = foldname,pattern = fpatterns[1])
  output_name <- getPattern(string = datafile,pattern = '(m[a-z]+[0-9]+[_]+)+',postpattern = fpatterns[1])
  write.csv(odormap,file = paste(output_name,output,sep = ''),row.names = F)
  T
}

#makes a trial details file
#op=1: write the names as odor names, 2 - write them as numbers
simMakeTrialDetailsFile <- function(foldname='.',fpatterns=c('series*.*data*'),output='trial_details.csv',op=1){
  #get the odor map, and write it to a file so that you can save time
  odormap <- simReadTrialOdorMap(foldname = foldname)
  #read the datafile name so that you can get the mouse prefix
  datafile <- list.files(path = foldname,pattern = fpatterns[1])
  output_name <- getPattern(string = datafile,pattern = '(m[a-z]+[0-9]+[_]+)+',postpattern = fpatterns[1])
  data.df <- read.csv(datafile,header = F)
  nocells <- ncol(data.df)-3 #the cells are in the columns
  odor.freq <- table(odormap[,2])
  odornames <- makeMapDfCols(odormap,cola = 2,colb = 3)
  trialdetails <- cbind.data.frame(odornames[names(odor.freq)],as.vector(odor.freq),rep(nocells,length(odor.freq)))
  names(trialdetails) <- c('odors','freq','nocells')
  write.csv(trialdetails,file = paste(output_name,output,sep = ''),row.names = F)
  T
}
  
#function that reads the odor mapping as in the trial number to odor number mapping
#and the second col is the odor number
#input: could be the name of a file or the dataframe of the data read from the csv file
#op:1 - vector, 2 - df
simReadTrialDetailsFile <-function(foldname='.',fpatterns=c('.*trial_details.csv','series*_*data*'),op=1){
  headerfile <- list.files(path = foldname,pattern = fpatterns[1]) #read the header file
  if(length(headerfile) == 1) headerfile <- read.csv(file = headerfile,header = T,colClasses = "character")
  else {
    headerfile <- simMakeTrialDetailsFile() #if the trial detail file does not exist, create it
    headerfile <- list.files(path = foldname,pattern = fpatterns[1]) #read the header file
    headerfile <- read.csv(file = headerfile,header = T,colClasses = "character")
  }
  headerfile <- ConvertDfCols(headerfile,cols = c(2,3),op=1) #set columns 2 and 3 to numeric
  headerfile
}


#generic function: function to compute the effect of alpha, the significance, on the number of reliable and unreliable cells
#foldname: the folder that holds the dataset
#op = 1 get the number of reliable and unreliable cells.
# 2: return the classified data, the getCampDAtaSigList structures
# 3: list of reliable and unreliable cells
# alpharange: the range or values of alphas for which the function is evaluated
simGetAlphaEffectRUCells <- function(foldname=".",alpharange=c(0.05,0.01,0.005,0.001,0.0005),op=1){
  #cat('\ncampGetAlpha')
  res.ls <- lapply(alpharange,function(x){
    #cat('\nalpha: ',x)
    dat.sig <- simGetDataSigList(foldname = foldname,alpha = x )
    trial.class <- campClassifyTrialCells(dat.sig)
    cells <- joinListDFs(trial.class[[4]])/length(dat.sig[[1]][[1]]) * 100
    #print(cells)
    res <- switch(op,apply(cells,2,mean), #reliable and unreliable cells
                  dat.sig,                #the data set
                  list(trial.class[[2]],trial.class[[3]]))    
  })
  names(res.ls) <- alpharange #set the names to alpha value (significance) that is being served
  res.ls
}

#function that calculates the UD for similar and dissimilar odors for a range of alphas, i.e., significance/
# alpharange: the range or values of alphas for which the function is evaluated
# foldname: the folder that holds the dataset
# cols: the cols that you are interested in. empty is get everything
simComputeUDAlphas <- function(foldname='.',alpharange=c(0.05,0.01,0.005,0.001,0.0005),cols=c(),op=1){
  #get the different campDatasiglist with different alphas
  datsets.ls <- simGetAlphaEffectRUCells(foldname = foldname,alpharange = alpharange,op=2)
  res.ls <- lapply(datsets.ls, function(x){
    res <- campComputeUD(x)
    if (length(cols)==0) res.df <- convertNestedListsDF(res)
    else res.df <- convertNestedListsDF(lapply(res,'[',cols))
    res.df
  })
  #cat('\n',str(res.ls))
  #res.ls
  res.df <- joinListDFs(res.ls)
  res.df  
}

#getting data for Simon. Basically, the list of reliable and unreliable cells and the plane number
#output: cell no, reliable or unreliable or silent, plane no, trial no(s)
#foldername: of the place where the data resides
simGenRUCellsPlaneData <- function(foldname='.',fpatterns=c('.*series.*data*'),alpha=0.01,op=1){
  datafile <- list.files(path = foldname,pattern = fpatterns[1])
  dat.df <- read.csv(file = datafile)
  #get the sig data list and then the classified trials
  data.sig <- simGetDataSigList(foldname = foldname,fpatterns = fpatterns,alpha = alpha)
  trial.data <- campClassifyTrialCells(data.sig)
  nocells <- length(trial.data[[1]][[1]]) 
  #get the plane number
  plane.nos <- colnames(dat.df)[2:(ncol(dat.df)-2)]
  #lets go through all the cells and mark them accordingly by looking in the rows of the dat.df to 
  trial.cells <- as.data.frame(trial.data[[1]])
  colnames(trial.cells) <- names(trial.data[[1]])
  res.df <- cbind(plane.nos,trial.cells)
  #now, take out the silent cells
  cell.activity <- apply(trial.cells, 1, sum)
  silent.cells <- which(cell.activity == 0)
  cat('\nsilent cells',silent.cells,'\n')
  #cbind(res.df,cell.activity)
  res.df <- res.df[-silent.cells,]
  res.df
}


#all of the analysis calls and functions
simFunctions <- function(){
  #run these files to do a quick setup
  #1. first copy the odors.csv and odor_params files and then change thhe data file name to contain series_data0.csv
  #2. run the commands in the following order
  simMakeOdorMapFile()
  simMakeHeader()
  simMakeTrialDetailsFile()
  
  tst <- simGetDataSigList(alpha = 0.05) #gets the responses in this dir. similar to campDataSigList
  #correlation mean response vs trials
  tmp3 <- getCampResponseSigCor(tst,op=2) #gets the number of trials, and mean response of each cell, op =2 as a df
  tmp4 <- meanDF(tmp3[[1]]) #get the averages
  #plot and get correlations, cor ~ 0.96
  fploteq(tmp4[,1],tmp4[,2],xlabel='no. trials',ylabel='mean response',ticknoy=5,ticknox=2)
  cor(tmp4[,1],tmp4[,2]) #cor ~ 0.96
  cor(tmp3[[1]][,1],tmp3[[1]][,2]) #cor ~ 0.56 , as a group the correlation is lower
  fploteq(tmp3[[2]][,1],tmp3[[2]][,2],xlabel='no. trials',ylabel='mean response',ticknoy = 2,rndfact = c(1,5))
  #all correlations, mean is 0.8
  sapply(tmp3, function(i) cor(i[,1],i[,2])) 
   
  #reliable and unreliable ccells
  tst1 <- simGetDataSigList(alpha = 0.05)
  temp1 <- campClassifyCells(tst1)
  #just plot all the cells
  fstripchartvecs(temp1[[2]],markersize = 0.8,ylabel = 'cells',semthick = 0,tickno=5,markersizemean = 0)
  #result: strange that there are no cells numbered 1 to 50, what's going on?
  #reliable and unreliable cells
  fstripchartvecs(list(res2,res3),markersize = 0.8,ylabel = '% cells',semthick = 1,tickno=2,markersizemean = 0,fixy = c(0,0.6),pairplot = 2,methodstr = 'overplot')
  
  
  #overlap
  res3 <- campOverlapAllAnalog(tst1,nolevels = 1) #get the overlap, summed for each reliability level
  res4 <- joinUnevenDFs(res3) # join the scores for all odor pairs
  res5 <- getDfSplit(res4,factorcol = 1) #split by reliability level and then plot
  fstripchartvecs(res5,markersize = 0.8,ylabel = 'overlap',semthick = 0,tickno=3,rndfact = 5)
  
  
  #correlations between the MB responses
  tst <- simGetDataSigList(alpha = 0.01)
  temp2 <- computeListMatCorr(getGroupedListsDf(tst),op=1)
  temp2 <- simGetOdorCorrelations(tst)
  write.csv(as.data.frame(temp2),file = '/home/shyam/extras/sync/measurements/olfaction/simon_analysis/odor_correlations.csv')
  #for each odor find the odor with the most correlation
  sapply(1:10, function(i) which(temp2[i,]==max(temp2[i,])))
  res33 <- campOverlapAllAnalog(tst,op = 9)
  res33 <- campAnalogStats(res33,op=5)
  temp3 <- unlist(flattenDFtoList(temp2))
  res34 <- insertColDf(res33,newcol = temp3[rownames(res33)],colname = 'cor.')
  #plot first as pair plot and then jsut plot
  fstripchartvecs(list(res34[,4],res34[,5]),markersize = 0.8,ylabel = 'overlap',semthick = 1,tickno=2,markersizemean = 1,pairplot = 2,methodstr = 'overplot')
  fstripchartvecs(list(res34[,4],res34[,5]),markersize = 0.8,ylabel = 'overlap',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  #overlap of similar and dissimilar odors
  #dissimilar
  fstripchartvecs(list(res34[which(res34[,1]<0.15),4],res34[which(res34[,1]<0.15),5]),markersize = 0.8,ylabel = 'overlap',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  #similar
  fstripchartvecs(list(res34[which(res34[,1]>0.5),4],res34[which(res34[,1]>0.5),5]),markersize = 0.8,ylabel = 'overlap',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  #midsimilar
  fstripchartvecs(list(res34[which(res34[,1]>0.15) && which(res34[,1]<0.5),4],res34[which(res34[,1]>0.15) && which(res34[,1]<0.5),5]),markersize = 0.8,ylabel = 'overlap',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
    
  #useful discrimination
  setwd("/nadata/cnl/data/shyam/simon/mouse164")
  res33 <- campOverlapAllAnalog(tst,op = 4)
  res33 <- campAnalogStats(res33,op=5)
  temp3 <- unlist(flattenDFtoList(temp2))
  res34 <- insertColDf(res33,newcol = temp3[rownames(res33)],colname = 'cor.')
  #plot discrimination differences all odors
  fstripchartvecs(list(res34[,4],res34[,5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 2,methodstr = 'overplot')
  fstripchartvecs(list(res34[,4],res34[,5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  #discrimination, similar, dissimilar, and mid-similar
  fstripchartvecs(list(res34[which(res34[,1]<0.15),4],res34[which(res34[,1]<0.15),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  fstripchartvecs(list(res34[which(res34[,1]>0.5),4],res34[which(res34[,1]>0.5),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  fstripchartvecs(list(res34[which(res34[,1]>0.15) && which(res34[,1]<0.5),4],res34[which(res34[,1]>0.15) && which(res34[,1]<0.5),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  
  #with saturation
  res35 <- campOverlapAllAnalog(tst,op = 5)
  res35 <- campAnalogStats(res35,op=5)
  res36 <- insertColDf(res35,newcol = temp3[rownames(res35)],colname = 'cor.')
  fstripchartvecs(list(res36[,4],res36[,5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 2,methodstr = 'overplot')
  fstripchartvecs(list(res36[,4],res36[,5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  fstripchartvecs(list(res36[which(res36[,1]<0.15),4],res36[which(res36[,1]<0.15),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  fstripchartvecs(list(res36[which(res36[,1]>0.5),4],res36[which(res36[,1]>0.5),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  fstripchartvecs(list(res36[which(res36[,1]>0.15) && which(res36[,1]<0.5),4],res36[which(res36[,1]>0.15) && which(res36[,1]<0.5),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 0)
  
  
  #figure in paper
  #comparison of saturated vs non-saturated
  fstripchartvecs(list(res34[,4],res36[,4],res34[,5],res36[,5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  #means: normal - 2.2 vs 1.04, saturated - 1.86 vs 7.6
  getStats(res36[,4])
  #similar
  fstripchartvecs(list(res34[which(res34[,1]<0.15),4],res36[which(res36[,1]<0.15),5],res34[which(res34[,1]<0.15),4],res36[which(res36[,1]<0.15),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  #dissimilar
  fstripchartvecs(list(res34[which(res34[,1]>0.5),4],res36[which(res36[,1]>0.5),4],res34[which(res34[,1]>0.5),5],res36[which(res36[,1]>0.5),5]),markersize = 0.8,ylabel = 'useful D.',semthick = 1,tickno=2,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  
  #result: there is always a reduction in useful discrimination for the reliable cells with saturation because if a reliable cell overlaps with an unreliable 
  #cell, then with saturation while the reliable cell's contributions do not go up, the unreliable cell's contributions increase with time.
  
  #function that calculates the odor habituation responses over a number of trials
  #op=1, for every odor, for each cell across all trials, get the trial where it hits its peak
  #op=2, now, across odor see if its always the first trial that is the peak response
  #op=3, across odors, gets the significant responses for each cell, if not all cells, padded out
  #op-4, just get all the responses for each odor across trials as df
  
  #habituation
  #no of times each trial is the highest response : flies
  tst1 <- campHabituation(tst,op=1)
  tst2 <- getDensDfStats(transposeDF(convertNestedListsDF(tst1)),start = 1)
  fploteq(1:8,unlist(tst2[1,]),ticknox = 1,fixx = c(1,8),stddev = unlist(tst2[3,]),ticknoy = 1,xlabel = 'trial #',ylabel = '# peaks',markersize = 1.2,rndfact = c(1,1))
  #looking whether response go down after first significant response
  tst1 <- campHabituation(tst,op=2)
  tst3 <- getDensDfStats(transposeDF(convertNestedListsDF(tst1)),start = 1)
  fploteq(1:2,unlist(tst3[1,]),fixx = c(1,2),stddev = unlist(tst3[3,]),ticknoy = 2,xlabel = 'first trial peak',ylabel = '# peaks')
  #look at the first and second significant response
  tst1 <- campHabituation(tst,op=3)
  tst2 <- joinUnevenDFs(tst1) #join all of the odors
  tst3 <- tst2[which(tst2[,1]>0),] #only get the ones with at least one sig response
  tst4 <- tst3[which(tst3[,2]>0),] #get ones with both significant responses
  fstripchartvecs(list(tst4[,1],tst4[,2]),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 1,methodstr = 'overplot')
  fstripchartvecs(list(tst4[,1],tst4[,2]),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 0,methodstr = 'jitter')
  tst5 <- tst4[which(tst4[,1]<0.5),] #blow up responses from 0 to 0.5
  fstripchartvecs(list(tst5[,1],tst5[,2]),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 1,methodstr = 'overplot')
  
  #plot by number of sigresponses
  tst1 <- campHabituation(tst,op=3)
  tst2 <- joinUnevenDFs(tst1) #join all of the odors
  tst3 <- getRowsDFNonzeros(tst2,n=4,op=3)
  fstripchartvecs(as.list(data.frame(tst3)),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot')
  
  tmp1 <- joinUnevenDFs(campHabituation(tst,op=3))
  tst3 <- getRowsDFNonzerosAll(tmp1,op=3) #tst3 gives a list where each element is a df of 1,2,3...etc sig. responses
  #e.g., plotting 3 sig responses
  fstripchartvecs(as.list(data.frame(tst3[[3]])),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 2,methodstr = 'overplot')
  fstripchartvecs(as.list(data.frame(tst3[[3]])),markersize = 0.8,ylabel = 'response',semthick = 2,pairplot = 0,methodstr = 'overplot')
  
  #looking at all responses across odors
  res <- campHabituation(getCampDataSigList(),op=4)
  res2 <- joinUnevenDFs(res)
  rownames(res2) <- 1:nrow(res2)
  res1 <- getRowsDFNonzeros(res2[,1:4],n=1,op=2) #get no of cells with 1 response in first 4 trials
  nrow(getRowsDFNonzeros(res2[rownames(res1),],n=2,op=1)) #of these show cells with at least 2 responses
  
  
  
  
  #getting the euclidean statistics that Glenn had requested
  tmp1 <- campGetTrialDistances(tst,op=3,rettype = 2)
  fstripchartvecs(as.list(tmp1),markersize = 0.8,ylabel = 'trial distance',semthick = 2,pairplot = 2,tickno = 4,methodstr = 'overplot')
  tmp1 <- campGetTrialDistances(tst,op=1,rettype = 2)
  fstripchartvecs(as.list(tmp1),markersize = 0.8,ylabel = 'origin distance',semthick = 2,pairplot = 2,tickno = 5,methodstr = 'overplot')
  
  
  #revision: finding out 
  
  
}

#function that compares the contributes of reliable and unreliable cells with and without saturation
#op=1, summed data, 2 - mean data
simCompareSatNorm <- function(cor=c(),skip = c('empty','paraff'),op=1){
  alldata <- getCampDataSigList() #get the data
  #get the correlations for hallem and rob
  #get non-saturated data
  res33 <- campOverlapAllAnalog(alldata,skip = skip,op = 4)#op=4, normal
  res33 <- campAnalogStats(res33,op=5) #get all the stats
  res33 <- res33[odors,]
  res35 <- campOverlapAllAnalog(alldata,skip = skip,op = 6)#op=6, which cells are active for both odors
  res35 <- campAnalogStats(res35,op=6) #get the cell overlap stats
  res35 <- res35[odors,]
  #saturated data
  res34 <- campOverlapAllAnalog(alldata,skip = skip,op = 5)#op=4, saturated 
  res34 <- campAnalogStats(res34,op=5) #get all the stats, slope, difference ...
  res34 <- res34[odors,]
  #put it all togetehr, odor similarity(hallem and rob), slope, unreliable and reliable cell contributtions: normal first then saturated
  #cat(str(res33),str(rescor))
  res <- cbind(rescor[odors,c(2,3)],res33[odors,2:4],res34[odors,2:4],res33[odors,c(7,8)],rescor[odors,ncol(rescor)],res35[odors,c(1,2)])
  names(res)[3:13] <- c('slope.nor','ursum.nor','rsum.nor','slope.sat','ursum.sat','rsum.sat','ur.cells','r.cells','trials','urel.act','rel.act') 
  list(res,res33,res34,rescor,rescor.sat)
}
