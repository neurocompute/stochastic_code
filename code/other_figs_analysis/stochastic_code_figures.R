#This file contains code that generates the figures in the paper Srinivasan et al., 2023. Each figure is a separate function.
#The naming convention is as follows. stoch_fig_x_a would be main figure X, panel a. stoch_fig_sx_a would be supplement figure X, panel a

#Datasets: before executing the code make sure that the datasets are present, and note their directories. The main dataset for fly is 
#fly-main_dataset and for mouse is 164. When we use extra datasets, we will mention them explicitly, otherwise, the datasets used 
#are the main ones. The fly main set is loaded into the variable tmp and the mouse one into variable tst.

#loading the fly main and mouse datasets: make it a global variable.
#tmp <- getCampDataSigList(foldname = 'path/to/dir', alpha = 0.05)
#tst <- simGetDataSigList(foldname = 'path/to/dir', alpha = 0.05)

#shows how to load the datasets into variables tmp and tst used throughout the file
#as above, make the commands within these datasets or variables that hold them into global variables, so you can use them
#Before executing any of the commands for the any of the figures below, execute the loading of datasets commands first.
stoch_code_loading_data <-function(){
  #fly main dataset
  tmp <- getCampDataSigList(foldname = 'path/to/dir', alpha = 0.05)
  #mouse main dataset
  tst <- simGetDataSigList(foldname = 'path/to/dir', alpha = 0.05)
}

#Figure 1: fly and mouse statistics of responses
stoch_code_fig1 <- function(){
  #flies
  #plotting the # of significant, reliable and unreliable cells for every odor
  #getting the number of significant responses, reliable and unreliable cells that fire per odor
  #1c. plots  rel. and unrel. for each odor
  temp3 <- getCampDataRel(tmp,op=2)
  fstripchartvecs(as.list(temp3[,2:3]),markersize = 0.8,ylabel = '% cells',semthick = 2,pairplot = 0,tickno = 4,fixy = c(1,60))
  #the nos
  apply(temp3,2,mean)
  apply(temp3,2,sd)/sqrt(nrow(temp3))
  
  #1d. getting the number of reliable and unreliable cell responses per trial
  temp <- campClassifyTrialCells(tmp)
  fstripchartvecs(as.list(joinListDFs(temp[[4]])),markersize = 0.8,ylabel = '% cells',semthick = 2,pairplot = 0,tickno = 5,fixy = c(1,40))
  #getting the numbers
  apply(joinListDFs(temp[[4]]),2,mean) #6.52381 8.97619 
  apply(joinListDFs(temp[[4]]),2,sd)/sqrt(7) #2.127042 2.553249 
 
  #1e: plotting reli. vs cor. 
  res1 <- campGetReliabStats(tmp,respop = 2)
  fploteq(res1[[4]][[1]][,1],res1[[4]][[1]][,2],stddev = res1[[4]][[1]][,3],markersize = 0.8,xlabel = 'reliability',ylabel = 'mean response',ticknoy = 5,rndfact=c(1,5),fixx = c(1,6),fixy=c(0,1),stdthick = 15,eq = function(x) .02*(exp(0.6*x)) )

  #mouse  
  #1f: plots  rel. and unrel. per odor
  temp3 <- getCampDataRel(tst,op=2)
  temp3 <- temp3/length(tst[[1]][[1]]) *100
  fstripchartvecs(as.list(temp3[,2:3]),markersize = 0.8,ylabel = '% cells',semthick = 2,pairplot = 0,tickno = 3,fixy = c(0,60))
  
  #1g: rel and unrel cells per trial
  temp <- campClassifyTrialCells(tst)
  temp1 <- joinListDFs(temp[[4]])/length(tst[[1]][[1]]) * 100
  fstripchartvecs(as.list(temp1),markersize = 0.8,ylabel = '% cells',semthick = 2,pairplot = 0,tickno = 10,fixy = c(1,40))

  #1h: trials vs avvg. response
  res3 <- campGetReliabStats(tst,respop = 2)
  fploteq(res3[[4]][[1]][,1],res3[[4]][[1]][,2]/100,stddev = res3[[4]][[1]][,3]/100,markersize = 0.8,xlabel = 'reliability',ylabel = 'mean response',ticknoy = 5,rndfact = c(1,5),fixy=c(0,1),fixx = c(1,8),stdthick = 15,eq = function(x)  + 0.01*(exp(0.5*x)) )
  
  #writing to a file
  write.xlsx(temp3[,2:3],file='fig1.xlsx',sheetName = 'fig1c',append = F)
  write.xlsx(joinListDFs(temp[[4]]),file='fig1.xlsx',sheetName = 'fig1d',append = T)
  write.xlsx(as.data.frame(res1[[4]][[1]]),file='fig1.xlsx',sheetName = 'fig1e',append = T)
  write.xlsx(temp3[,2:3],file='fig1.xlsx',sheetName = 'fig1f',append = T)
  write.xlsx(temp1,file='fig1.xlsx',sheetName = 'fig1g',append = T)
  write.xlsx(as.data.frame(res3[[4]][[1]]),file='fig1.xlsx',sheetName = 'fig1h',append = T)
  
}

#Fig. 2: showing the cumulative distribution of the response, overlap
stoch_code_fig2 <- function(){
  
  #cumulative distributions of reliability
  #flies: 2B
  temp <- campClassifyTrialCells(tmp,trials = 6,op = 2)
  temp2 <- as.vector(table(unlist(temp[[1]])) )
  NlsFitCDF(unlist(sapply(1:6,function(i) rep(i/6,temp2[i+1]))),dist = 2,graphparams = list(ticknox=2,ticknoy=2,rndfact=c(1,1),xlabel='reliability',ylabel='cum. freq.',fixx=c(0,1)) )
  
  #mouse: 2E
  temp3 <- campClassifyTrialCells(tst,trials = 8,op = 2)
  temp4 <- as.vector(table(unlist(temp3[[1]])) )
  NlsFitCDF(unlist(sapply(1:8,function(i) rep(i/8,temp4[i+1]))),dist = 2,graphparams = list(ticknox=2,ticknoy=2,rndfact=c(1,1),xlabel='reliability',ylabel='cum. freq.',fixx=c(0,1)))
  
  #cumulative distribution of responses
  #flies, fig. 2c
  temp1 <- campGenRUCellMat(tmp,op=1)
  temp3 <- unlist(temp1)
  #best gamma fit
  NlsFitCDF(temp3[which(temp3>0)],dist = 2,graphparams = list(ticknox=4,ticknoy=5,rndfact=c(1,5),xlabel='response rate',ylabel='cum. freq'),fittype = 4,params = c(1.15,.32))
  
  #mouse, fig. 2f
  res1 <- campGenRUCellMat(tst,op=1)
  res2 <- unlist(res1)
  NlsFitCDF(res2[which(res2>0)]/100,dist = 2,graphparams = list(ticknox=4,ticknoy=2.5,rndfact=c(1,5),xlabel='response rate',ylabel='cum. freq') )

  #cumulative distribution of the overlap
  #flies, Fig. 2d
  res5 <- campComputeOverlap(tmp,op=2)
  res4 <- joinListDFs(res5)
  NlsFitCDF(res4[,2],dist = 2,graphparams = list(ticknox=1,ticknoy=1,rndfact=c(2,2),xlabel='overlap',ylabel='cum. freq.'))

  #mouse, Fig. 2g
  res7 <- campComputeOverlap(tst,op=2)
  res8 <- joinListDFs(res7)
  NlsFitCDF(res8[,2],dist = 2,graphparams = list(ticknox=1,ticknoy=1,rndfact=c(2,2),xlabel='overlap',ylabel='cum. freq.'))
  # 0.1881611 0.3521271

  #do each plot after the correponding code.
  res1 <- RelCumFreq(unlist(sapply(1:6,function(i) rep(i/6,temp2[i+1]))))
  write.xlsx(res1,file='fig2.xlsx',sheetName = 'fig2b',append = F)
  res2 <- RelCumFreq(temp3[which(temp3>0)])
  write.xlsx(res2,file='fig2.xlsx',sheetName = 'fig2c',append = T)
  res3 <- RelCumFreq(res4[,2])
  write.xlsx(res3,file='fig2.xlsx',sheetName = 'fig2d',append = T)
  res4 <- RelCumFreq(unlist(sapply(1:8,function(i) rep(i/8,temp4[i+1]))))
  write.xlsx(res4,file='fig2.xlsx',sheetName = 'fig2e',append = T)
  res5 <- RelCumFreq(res2[which(res2>0)]/100)
  write.xlsx(res5,file='fig2.xlsx',sheetName = 'fig2f',append = T)
  res6 <- RelCumFreq(res8[,2])
  write.xlsx(res6,file='fig2.xlsx',sheetName = 'fig2g',append = T)

}

#figure 1: examining correlations between fly OSNs and KCs
#to do: generate hallem 2
#showing the second part of Fig. 4: overlap vs reliability for similar and dissimilar odors
stoch_code_fig4 <-function(){
  #for all the odor sets
  #for each odor set, go to the corresponding unzipped directory, and run this function
  #as an aillustration, for the set main, unzip the main fly dataset and in that diretory run the 3 commands below.
  ressig.main <- campMatchHalRobCor(hallem2,celltype = 1,matchlevel = 1)
  resrel.main <- campMatchHalRobCor(hallem2,celltype = 2,matchlevel = 1)
  resunrel.main <- campMatchHalRobCor(hallem2,celltype = 3,matchlevel = 1)
  
  #run these threes commands in the 0904 set directory, and so on
  ressig.0904 <- campMatchHalRobCor(hallem2,celltype = 1,matchlevel = 1)
  resrel.0904 <- campMatchHalRobCor(hallem2,celltype = 2,matchlevel = 1)
  resunrel.0904 <- campMatchHalRobCor(hallem2,celltype = 3,matchlevel = 1)
  
  ressig.110109.2 <- campMatchHalRobCor(hallem2,celltype = 1,matchlevel = 1)
  resrel.110109.2 <- campMatchHalRobCor(hallem2,celltype = 2,matchlevel = 1)
  resunrel.110109.2 <- campMatchHalRobCor(hallem2,celltype = 3,matchlevel = 1)
  
  ressig.110109 <- campMatchHalRobCor(hallem2,celltype = 1,matchlevel = 1)
  resrel.110109 <- campMatchHalRobCor(hallem2,celltype = 2,matchlevel = 1)
  resunrel.110109 <- campMatchHalRobCor(hallem2,celltype = 3,matchlevel = 1)
  
  ressig.110108 <- campMatchHalRobCor(hallem2,celltype = 1,matchlevel = 1)
  resrel.110108 <- campMatchHalRobCor(hallem2,celltype = 2,matchlevel = 1)
  resunrel.110108 <- campMatchHalRobCor(hallem2,celltype = 3,matchlevel = 1)
  
  ressig <- joinUnevenDFs(list(ressig.main,ressig.0904,ressig.110108,ressig.110109,ressig.110109.2))
  
  resrel <- joinUnevenDFs(list(resrel.main,resrel.0904,resrel.110108,resrel.110109,resrel.110109.2))
  
  resunrel <- joinUnevenDFs(list(resunrel.main,resunrel.0904,resunrel.110108,resunrel.110109,resunrel.110109.2))
  #The main figure 4a-c
  #fitting only those correlations that are within 0.2 of each other in Rob and Hallem. 
  temp <- which((abs(ressig[,2]-ressig[,3]) <= 0.2)) #all odor pair indices
  temp1 <- ressig[temp,]
  #r2 = 0.89, filename:correlation_antenna_vs_mb_allcells
  fploteq(temp1[,2],temp1[,3],ticknox = 4,ticknoy = 4,rndfact = c(10,10),eq = function(x) 0.047 + 0.78*x,xlabel = 'OSNs',ylabel = 'KCs',fixx=c(-.3,0.9))
  write.xlsx(as.data.frame(temp1),file='fig4.xlsx',sheetName = 'fig4a',append = F)
  
  temp1 <- resrel[temp,]
  #r2 = 0.8, filename:correlation_antenna_vs_mb_reliable
  fploteq(temp1[,2],temp1[,3],ticknox = 4,ticknoy = 4,rndfact = c(10,10),eq = function(x) 0.068 + 0.89*x,xlabel = 'OSNs',ylabel = 'KCs',fixy=c(-0.3,0.9),fixx=c(-.3,0.9))
  write.xlsx(as.data.frame(temp1[,2:3]),file='fig4.xlsx',sheetName = 'fig4b',append = T)
  
  temp1 <- resunrel[temp,]
  #r2 = 0.2, filename:correlation_antenna_vs_mb_unreliable
  fploteq(temp1[,2],temp1[,3],ticknox = 4,ticknoy = 4,rndfact = c(1,10),eq = function(x) 0.039 + 0.05*x,xlabel = 'OSNs',ylabel = 'unreliable KCs',fixy=c(-0.3,0.9))
  write.xlsx(as.data.frame(temp1[,2:3]),file='fig4.xlsx',sheetName = 'fig4c',append = T)
  
  #flies: 4de
  res5 <- campGetOverlapCells(tmp,thresh.win = c(.35,1),op=3)
  fstripchartvecs(res5,markersize = 0.8,ylabel = 'overlap',semthick = 2,flipaxis = 2,stretchx = 1,angle = F,xlabel = 'reliability') #flipped axis
  #overlap_vs_reliability_flies
  res6 <- campGetOverlapCells(tmp,thresh.win = c(0,.15))
  fstripchartvecs(res6,markersize = 0.8,ylabel = 'overlap',semthick = 2,flipaxis = 2,stretchx = 1,angle = F,xlabel = 'reliability')

  res7 <- convertUnequalListToDf(res5)
  res8 <- convertUnequalListToDf(res6)
  write.xlsx(res7,file='fig4.xlsx',sheetName = 'fig4d',append = T)
  write.xlsx(res8,file='fig4.xlsx',sheetName = 'fig4e',append = T)
  
  #mouse, 4fg
  res5 <- campGetOverlapCells(tst,thresh.win = c(.3,1))
  fstripchartvecs(res5,markersize = 0.8,ylabel = 'overlap',semthick = 2)  
  #flipped axis
  fstripchartvecs(res5,markersize = 0.8,ylabel = 'overlap',semthick = 2,flipaxis = 2,stretchx = 1,angle = F,xlabel = 'reliability')
  res6 <- campGetOverlapCells(tst,thresh.win = c(0,.15))
  fstripchartvecs(res6,markersize = 0.8,ylabel = 'overlap',semthick = 2,tickno = 10,fixy = c(0,1))
  fstripchartvecs(res6,markersize = 0.8,ylabel = 'overlap',semthick = 2,flipaxis = 2,stretchx = 1,angle = F,xlabel = 'reliability',tickno = 4,rndfact = 5)
  
  res7 <- convertUnequalListToDf(res5)
  res8 <- convertUnequalListToDf(res6)
  write.xlsx(res7,file='fig4.xlsx',sheetName = 'fig4f',append = T)
  write.xlsx(res8,file='fig4.xlsx',sheetName = 'fig4g',append = T)
  
}


#the contribution of reliable and unreliable cells to odor discrimination and sparsity
stoch_code_fig5 <- function(){
  
  #odor similarity for each percentile of cell responses
  #flies; 5A
  #cosine_similarity_sparse_coding_figure_fly
  temp4 <- campGetOverlapTopNPopSeq(data.lst = tmp,thresh.win = c(0.5,1),topnseq = seq(25,100,25),sameop = 2,overop=14,op=3,deciles = c(),respop = 2)
  temp3 <- campGetOverlapTopNPopSeq(data.lst = tmp,thresh.win = c(0,.15),topnseq = seq(25,100,25),sameop = 2,overop=14,op=3,deciles = c(),respop = 2)
  fstripchartvecs(list(temp4,temp3),markersize = 0.8,ylabel = 'cosine similarity',semthick = 2,markersizemean = 1.2,fixy = c(0,0.9),tickno = 3)
  
  res3 <- convertListToDF(temp3)
  res4 <- convertListToDF(temp4)
  write.xlsx(res4,file='fig5.xlsx',sheetName = 'fig5a1',append = F)
  write.xlsx(res3,file='fig5.xlsx',sheetName = 'fig5a2',append = T)
  
  #mouse; 5C
  #cosine_similarity_sparse_coding_figure_mouse
  temp4 <- campGetOverlapTopNPopSeq(data.lst = tst,thresh.win = c(0.5,1),topnseq = seq(25,100,25),sameop = 2,overop=14,op=3,deciles = c(),respop = 2)
  temp3 <- campGetOverlapTopNPopSeq(data.lst = tst,thresh.win = c(0,.15),topnseq = seq(25,100,25),sameop = 2,overop=14,op=3,deciles = c(),respop = 2)
  fstripchartvecs(list(temp4,temp3),markersize = 0.8,ylabel = 'cosine similarity',semthick = 2,markersizemean = 1.2,tickno = 3,fixy = c(0,.9))
  
  res3 <- convertListToDF(temp3)
  res4 <- convertListToDF(temp4)
  write.xlsx(res3,file='fig5.xlsx',sheetName = 'fig5c2',append = T)
  write.xlsx(res4,file='fig5.xlsx',sheetName = 'fig5c1',append = T)
  
  
  #comparing number of top and bottom cells in each reliability category  
  #flies: 5B
  res2 <- campCompTopRelCells(tmp,topn = 25,deciles = 25,op=1)
  res3 <- campCompTopRelCells(tmp,topn = 100,deciles = 25,op=1)
  fstripchartvecs(list(res2,res3),markersizemean = 1.2,col = T,flipaxis = 1,semthick = 1,angle = F,xlabel = 'reliability',ylabel = '# cells')
  
  res4 <- convertUnequalListToDf(res2)
  res5 <- convertUnequalListToDf(res3)
  write.xlsx(res4,file='fig5.xlsx',sheetName = 'fig5b1',append = T)
  write.xlsx(res5,file='fig5.xlsx',sheetName = 'fig5b2',append = T)
  
  
  #mouse: 5D
  res2 <- campCompTopRelCells(tst,topn = 10,deciles = 10,op=1,respop = 2)
  res3 <- campCompTopRelCells(tst,topn = 100,deciles = 10,op=1,respop = 2)
  fstripchartvecs(list(res2,res3),markersizemean = 1.2,col = T,flipaxis = 1,semthick = 1,angle = F,xlabel = 'reliability',ylabel = '# cells',markersize = 0.8,tickno = 5)
  
  res4 <- convertUnequalListToDf(res2)
  res5 <- convertUnequalListToDf(res3)
  write.xlsx(res4,file='fig5.xlsx',sheetName = 'fig5d1',append = T)
  write.xlsx(res5,file='fig5.xlsx',sheetName = 'fig5d2',append = T)
  
  
  #classifier results for applying LDA to the data
  #flies and mice; classification accuracy: 5E
  #the data; accuracy
  temp1 <- rbind.data.frame(c(0.6,0.2,0.4),c(0.57,0.28,0.5),c(0.4,0.1,0.5),c(0.64,0.21,0.36),c(0.55,0.2,0.4),c(0.64,0.14,0.42))
  #f: flies, m:mouse
  rownames(temp1) <- c('lda.m','lda.f','knn.m','knn.f','svm.m','svm.f')
  colnames(temp1) <- c('Rel.','Unrel.','all')
  stackedHorzBarPlot(temp1,horz = T,spaces = 0.2,ticknox = 5,fontsize = 0.9,sepwidth = 1)
  #flies
  HorzBarPlot(temp2[2,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  #mouse
  HorzBarPlot(temp1[1,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  
  write.xlsx(temp1[1:2,],file='fig5.xlsx',sheetName = 'fig5e',append = T)
  
  #Fig. 5F; auc
  temp2 <- rbind.data.frame(c(0.8,0.59,0.64),c(0.72,0.72,0.83),c(0.74,0.39,0.59),c(0.71,0.36,0.65),c(0.88,0.7,0.84),c(0.93,0.66,0.82))
  rownames(temp2) <- c('lda.m','lda.f','knn.m','knn.f','svm.m','svm.f')
  colnames(temp2) <- c('Rel.','Unrel.','all')
  #flies
  HorzBarPlot(temp2[2,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  #mouse
  HorzBarPlot(temp2[1,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  
  write.xlsx(temp2[1:2,],file='fig5.xlsx',sheetName = 'fig5f',append = T)
  
}


#plotting the data for Useful discrimination with normal and extended training
stoch_code_fig6 <- function(){
  
  #flies:6AB
  #set the working directory first. Refer to the figure and all the datasets that were used. Then unzip all of them and go to the 
  #parent directory and execute these commands
  setwd("/path/to/dir")
  #traversing through the directories and getting the campCompareSatNOrm data, #gotox
  tmp2 <- campAllDirDiscScores(searchstr = c('*tseries*','*sigpval*'),mode = 2,targetfn = campCompareSatNorm,targetpars = list(hallem2,halrobmatch=1,matchlevel=0.2))
  #res.alldirsflySN <- tmp2# contains all the results of all directories <- tmp2
  tmp3 <- tmp2[which(sapply(tmp2, function(x) is.null(x[[1]]))==F)] #just the data for the valid directories
  res1 <- campProcessDirSets(tmp3,op=2)
  
  #fly dissimilar
  #reliable
  fstripchartvecs(res1[[1]][3:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,ticklabs = c('N','E'),methodstr = 'overplot',fixy = c(0,20),tickno = 2.5)
  #unreliable
  fstripchartvecs(res1[[1]][1:2],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,ticklabs = c('N','E'),methodstr = 'overplot',fixy = c(0,20),tickno = 5)
  sapply(res1[[1]][1:4],mean) #for stats
  
  test5 <- convertListToDF(res1[[1]][1:4])
  test5 <- convertListToDF(convertListNos2Str(res1[[1]][1:4]))
  write.xlsx(test5,file='fig6.xlsx',sheetName = 'fig6a',append = F)
  write.csv(test5,file = 'fig6a')
  
  #fly similar
  #reliable
  fstripchartvecs(res1[[3]][3:4],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,ticklabs = c('N','E'),methodstr = 'overplot',fixy = c(0,20),tickno = 12.5)
  #unreliable
  fstripchartvecs(res1[[3]][1:2],markersize = 0.8,ylabel = 'useful D.',semthick = 2,pairplot = 1,ticklabs = c('N','E'),methodstr = 'overplot',fixy = c(0,20),tickno = 6.25)
  sapply(res1[[3]][1:4],mean) #for stats
  
  test6 <- convertListToDF(res1[[3]][1:4])
  write.xlsx(test6,file='fig6.xlsx',sheetName = 'fig6b',append = T)
  write.csv(test6,file = 'fig6b')
  
  #mouse: 6CD  
  temp2 <- simGetOdorCorrelations(tst)
  temp3 <- unlist(flattenDFtoList(temp2))
  temp4 <- unlist(flattenDFtoList(computeListMatCorr(tst,matop = 1)))
  res33 <- campAnalogStats(campOverlapAllAnalog(tst,op = 4),dat.ls = tst,op=5)
  res34 <- insertColDf(res33,newcol = temp3[rownames(res33)],colname = 'cor.')
  res35 <- campAnalogStats(campOverlapAllAnalog(tst,op = 5),dat.ls = tst,op=5)
  res36 <- insertColDf(res35,newcol = temp3[rownames(res35)],colname = 'cor.')
  
  #choose similar odors
  res8 <- getIndexSameSubStr(rownames(res34),op=2)
  res2 <- intersect(which(res34[,1]<0.15),res8)  #dissimilar
  res3 <- intersect(which(res34[,1]>0.5),res8)   #similar
  #dissimilar unreliable;dissimilar_reliable_ud_mouse;dissimilar_unreliable_ud_mouse
  fstripchartvecs(list(res34[res2,4],res36[res2,4]),markersize = 0.8,ylabel = 'useful D.',ticklabs = c('N','E'),fixy = c(0,20),semthick = 1,tickno=2.5,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  #dissimilar reliable
  fstripchartvecs(list(res34[res2,5],res36[res2,5]),markersize = 0.8,ylabel = 'useful D.',ticklabs = c('N','E'),fixy = c(0,20),semthick = 1,tickno=5,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  
  test1 <- convertListToDF(list(res34[res2,4],res36[res2,4],res34[res2,5],res36[res2,5]))
  write.xlsx(test1,file='fig6.xlsx',sheetName = 'fig6c',append = T)
  
  #similar unreliable;similar_reliable_ud_mouse;similar_unreliable_ud_mouse
  fstripchartvecs(list(res34[res3,4],res36[res3,4]),markersize = 0.8,ylabel = 'useful D.',ticklabs = c('N','E'),fixy = c(0,20),semthick = 1,tickno=5,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  #similar rliable 
  fstripchartvecs(list(res34[res3,5],res36[res3,5]),markersize = 0.8,ylabel = 'useful D.',ticklabs = c('N','E'),fixy = c(0,20),semthick = 1,tickno=12.5,markersizemean = 1,pairplot = 1,methodstr = 'overplot')
  
  test2 <- convertListToDF(list(res34[res3,4],res36[res3,4],res34[res3,5],res36[res3,5]))
  write.xlsx(test2,file='fig6.xlsx',sheetName = 'fig6d',append = T)
  write.csv(test2,file = 'fig6d')
  
  #fly: 6E
  res4 <- campCompUDAllSingleRel(dat.ls = tmp,op=1)
  res5 <- joinListDFs(res4)
  res6 <- sapply(1:6,function(i) sum(res5[c(2*i,2*i -1),4]) - sum(res5[c(2*i,2*i -1),3]))
  #reliability_vs_disc_increase_fly
  fploteq(1:6,res6,xlabel = 'reliability',ylabel = 'disc. change',fixx = c(1,6),ticknoy = 4,eq = function(x) 0*x)
  
  test3 <- cbind.data.frame(1:6,res6)
  write.xlsx(test3,file='fig6.xlsx',sheetName = 'fig6e',append = T)
  write.csv(test3,file = 'fig6e')
  
  #mouse: 6F
  tst7 <- campCompUDAllSingleRel(dat.ls = tst,op=1)
  tst8 <- joinListDFs(tst7)
  #magnitude difference
  tst9 <- sapply(1:8,function(i) sum(tst8[c(2*i,2*i -1),4]) - sum(tst8[c(2*i,2*i -1),3]))
  #reliability_vs_disc_increase_mouse
  fploteq(1:8,tst9,xlabel = 'reliability',ylabel = 'disc. increase',fixx = c(1,8),ticknoy = 4,eq = function(x) 0*x)
  
  test4 <- cbind.data.frame(1:8,tst9)
  write.xlsx(test4,file='fig6.xlsx',sheetName = 'fig6f',append = T)
  write.csv(test4,file = 'fig6f')
  
}

#this generates the code for the analysis of the fly data in the methods figure 8
stoch_code_fig8 <- function(){
  
  #Fig. 8F
  res1 <- campFreqCellsAnalog(tmp)
  #heptanone
  test <- read.csv('tseries_data3.csv',stringsAsFactors = F,header = F);test3 <- as.matrix(test)
  #example_plot_ethyloctanoate_1st_trial
  fploteq(1:40,transposeDF(as.data.frame(test3[c(18,108,110,124),1:40])),markersize = 0.4,fixy = c(-.25,0.95),ticknoy = 2.5, ticknox = 5, xlabel = 'time (frames)',ylabel = 'Ca^2+',rndfact = c(1,2))
  write.xlsx(as.data.frame(test3),file='fig8.xlsx',sheetName = 'fig8f',append = F)
  write.csv(transposeDF(as.data.frame(test3[c(18,108,110,124),1:40])),file='fig8f.csv')
  
  #Fig. 8G
  res3 <- campGetDataSigBase(tmp) 
  res4 <- lapply(res3,'[[',1)
  res5 <- lapply(res3,'[[',2)
  temp <- RelCumFreqList(list(unlist(res5),unlist(res4)))
  fploteq(temp,markersize = 1,ticknox = 2,rndfact = c(5,1),ticknoy = 5,xlabel = 'response',ylabel = 'frequency')
  write.xlsx(as.data.frame(temp[[1]]),file='fig8.xlsx',sheetName = 'fig8g1',append = T)
  write.xlsx(as.data.frame(temp[[2]]),file='fig8.xlsx',sheetName = 'fig8g2',append = T)
  write.csv(as.data.frame(temp[[1]]),file='fig8g1.csv')
  write.csv(as.data.frame(temp[[2]]),file='fig8g2.csv')
  
  #Fig. 8H; base_sign_significant_cells_vs_response
  temp1 <- list(unlist(res5),unlist(res4))
  names(temp1) <- c('base','sig.')
  fstripchartvecs(temp1,sem = 3,ylabel = 'response',markersize = 0.8,markersizemean = 1.2)
  temp2 <- convertListToDF(temp1)
  write.csv(temp2,file='fig8h.csv')
  
  #Fig. 8I
  temp2 <- which(unlist(res5) > min(unlist(res4)))
  temp3 <- list(temp1[[1]][temp2],temp1[[2]][temp2])
  names(temp3) <- c('base','sig')  
  fstripchartvecs(temp3,sem = 3,ylabel = 'response',markersize = 0.8,markersizemean = 1.2,pairplot = 1,tickno = 2.5)
  temp4 <- convertListToDF(temp3)
  write.csv(temp4,file='fig8i.csv')
  
}

#main modeling figure 3
#this function contains the code to generate the modeling results of main Figure 3 exploring the effect of noise 
#in various components on circuit responses
stoch_code_fig3 <- function(){
  
  #Simulations to see the effect of noise.
  circuit1 <- createFlyCircuit(glomno = 50,kcno = 100,noodors = 6,aplCut = 10)
  
  #Fig. 3C: glomnoise
  res.glomnoise.ls <- exploreNoiseSeqRuns(circuit1,noruns = 5,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1))
  temp <- processNoiseSeqRunsResults(res.glomnoise.ls)
  fploteq(seq(.3,1,0.1),apply(temp,1,mean)[-(1:3)],stddev = apply(temp,1,sd)[-(1:3)],markersize = 0.8,stdthick = 2,stdepsilon = .25,ticknox = 2,ticknoy = 4,fixy = c(0,4),xlabel = 'noise',ylabel = 'rel/unrel cells',hline = 0.72)
  
  #Fig. 3D: kcapl
  res.kcaplnoise.ls <- exploreNoiseSeqRuns(circuit1,noruns = 5,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=3)
  temp <- processNoiseSeqRunsResults(res.kcaplnoise.ls)
  fploteq(seq(.3,1,0.1),apply(temp,1,mean)[-(1:3)],stddev = apply(temp,1,sd)[-(1:3)]/sqrt(5),markersize = 1,stdthick = 2,stdepsilon = .25,ticknox = 2,ticknoy = 8,xlabel = 'noise',ylabel = 'rel/unrel cells')
  
  #Fig. 3C: glomerular noise
  res.glomnoise.ls <- exploreNoiseSeqRuns(circuit1,noruns = 5,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1))
  temp <- processNoiseSeqRunsResults(res.glomnoise.ls)
  fploteq(seq(.3,1,0.1),apply(temp,1,mean)[-(1:3)],stddev = apply(temp,1,sd)[-(1:3)],markersize = 0.8,stdthick = 2,stdepsilon = .25,ticknox = 2,ticknoy = 4,fixy = c(0,4),xlabel = 'noise',ylabel = 'rel/unrel cells',hline = 0.72)
  
  #Fig. 3E
  test5 <- cirGetParamsTargetAllApls(temp5,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 )
  res5 <- cirAnalGoodParams(test5,op=2,par = 5,noiseop = 22)
  NlsFitCDF(res5[[5]][,6],graphparams = list(xlabel='pn noise',ylabel='freq.',ticknoy=2,ticknox=5),dist = 5)
  
  #Fig. 3F
  tst1 <- seq(0.25,0.425,0.025)
  tst2 <- seq(0,.1,.05)
  tst3 <- c(list(tst1),repList(tst2,3))
  temp12 <- exploreNoiseParams(noruns = 10,notrials = 6,noiseset = 2,noiserange = tst3,aplrange = c(5:9),op=22,classop = 2,kcapl_conn_par=tst4,aplkc_conn_par=list(8,c(20,3,0.4,0.45)));
  res2 <- cirAnalGoodParams(cirGetParamsTargetAllApls(temp12,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 ),op=2,par = 1)
  NlsFitCDF(res2[[5]][,6],graphparams = list(xlabel='apl noise',ylabel='freq.',ticknoy=2))
  
  #Fig. 3G
  temp15 <- exploreNoiseParams(noruns = 10,notrials = 6,noiseset = 2,noiserange = list(seq(0.25,0.35,0.05),seq(0,0.5,0.1),seq(0,0.5,0.1),seq(0,0.2,0.1),seq(0.1,0.5,0.1)),aplrange = c(6:11),classop = 2,kcapl_conn_par=list(7,c(20,2,0.4,0.45)),aplkc_conn_par=list(8,c(20,2,0.4,0.45)),thresh = list(0,c(1,0)),op=30)
  res15 <- cirAnalGoodParams(cirGetParamsTargetAllApls(temp15,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 ),op=2,par = 5,noiseop = 30)
  NlsFitCDF(res15[[5]][,2],graphparams = list(xlabel='kc noise',ylabel='freq.',ticknoy=2),dist = 1)
  
  #Fig. 3H
  temp5 <- exploreNoiseParams(noruns = 8,notrials = 6,noiseset = 2,noiserange = list(seq(0,0.5,.025),seq(0,0.5,0.1),seq(0,0.4,.05),seq(0,0.4,0.05)),aplrange = c(5:12),op=22,classop = 2)
  test5 <- cirGetParamsTargetAllApls(temp5,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 )
  res5 <- cirAnalGoodParams(test5,op=2,par = 5,noiseop = 22)
  fploteq(res5[[5]][,3],res5[[5]][,5],xlabel = 'apl noise',ylabel = 'pn-kc noise')
  
  #Fig. S7A: glommbnoise
  res.glommbnoise.ls <- exploreNoiseSeqRuns(circuit1,noruns = 5,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=2)
  temp <- processNoiseSeqRunsResults(res.glommbnoise.ls)
  fploteq(seq(.3,1,0.1),apply(temp,1,mean)[-(1:3)],stddev = apply(temp,1,sd)[-(1:3)],markersize = 0.8,stdthick = 2,stdepsilon = .25,ticknox = 2,ticknoy = 8,fixy = c(0,4),xlabel = 'noise',ylabel = 'rel/unrel cells',hline = 0.72)
  
  #Fig. S7B: apl->kc
  circuit3 <- createFlyCircuit(glomno = 50,kcno = 100,noodors = 6,aplCut = 6)
  res.aplkcnoise3.ls <- exploreNoiseSeqRuns(circuit3,noruns=2,notrials = 6,noiseset = 2,noiserange = list(seq(0,1,.025)),op=4)
  temp3 <- processNoiseSeqRunsResults(res.aplkcnoise3.ls,seqop = 2,op=102)
  fploteq(seq(0.2,1,0.025),unlist(temp3[[1]][1,9:41]),stddev = unlist(temp3[[2]][1,9:41]),markersize = 0.8,stdthick = 4,stdepsilon = .1,ticknox = 1,ticknoy = 4,rndfact = c(2,2),xlabel = 'noise',ylabel = 'rel/unrel cells',fixx = c(0.2,1),fixy = c(0,4),hline = 0.72) 
  
  #Fig. S7C: pn,pn-mb
  res.pnpnmbnoise.ls <- exploreNoiseSeqRuns(circuit1,noruns = 5,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=10)
  temp <- processNoiseSeqRunsResults(res.pnpnmbnoise.ls)
  fploteq(seq(.3,1,0.1),apply(temp,1,mean)[-(1:3)],stddev = apply(temp,1,sd)[-(1:3)],markersize = 0.8,stdthick = 2,stdepsilon = .25,ticknox = 2,ticknoy = 8,fixy = c(0,4),xlabel = 'noise',ylabel = 'rel/unrel cells',hline = 0.72)
  
}



#figure S1: for these plots, first run the main figure code to get res1 and res2, and then run this code
stoch_code_figs1 <- function(){

  #fly plots: A
  res1 <- campGetReliabStats(tmp,respop = 2)
  fploteq(res1[[4]][[1]][,1],res1[[4]][[1]][,2],stddev = res1[[4]][[1]][,3],markersize = 0.8,xlabel = 'reliability',ylabel = 'mean response',fixx = c(1,6),stdthick = 15,eq = function(x) .02*(exp(0.6*x)),logs = c(F,T) )
  write.xlsx(as.data.frame(res1[[4]][[1]]),file='figs1.xlsx',sheetName = 'figs1a',append = F)
  
  res3 <- campGetReliabStats(tst,respop = 2)
  #mouse plots: B
  fploteq(res3[[4]][[1]][,1],res3[[4]][[1]][,2]/100,stddev = res3[[4]][[1]][,3]/100,markersize = 0.8,xlabel = 'reliability',ylabel = 'mean response',fixx = c(1,8),stdthick = 15,eq = function(x)  + 0.01*(exp(0.5*x)),logs = c(F,T) )
  write.xlsx(as.data.frame(res3[[4]][[1]]),file='figs1.xlsx',sheetName = 'figs1b',append = T)
  write.csv(as.data.frame(res3[[4]][[1]]),file='figs1b.csv')
  
  #flies: C
  temp3 <- campGetReliabResp(tmp)
  temp3 <- campGetReliabResp(tmp,respop = 2) #expected value average response
  temp4 <- matchDataFrames(temp3[[1]],temp3[[2]])
  temp5 <- listDF(temp4)
  fstripchartvecs(temp5[getSortPosns(names(temp5) )][2:7],ylabel = 'avg. response',tickno = 4,xlabel = 'reliability')
  write.xlsx(convertUnequalListToDf(temp5[getSortPosns(names(temp5) )][2:7]),file='figs1.xlsx',sheetName = 'figs1c',append = T)
  
  #mouse: D
  temp <- campGetReliabResp(tst)
  temp <- campGetReliabResp(tst,respop = 2) #expected value average response
  temp1 <- matchDataFrames(temp[[1]],temp[[2]])
  temp2 <- listDF(temp1)
  fstripchartvecs(temp2[getSortPosns(names(temp2) )][2:9],ylabel = 'avg. response',tickno = 4,xlabel = 'reliability')
  write.xlsx(convertUnequalListToDf(temp2[getSortPosns(names(temp2) )][2:9]),file='figs1.xlsx',sheetName = 'figs1d',append = T)
  
  
}


#contains code to generate Fig. S2 showing CV, fano factor, and SD for flies and mice
stoch_code_figS2 <- function(){
  #signal vs sd for flies, S2C
  temp <- campSigNoiseComp(tmp,op=1,result = 2)
  temp1 <- joinUnevenDFs(temp)
  fploteq(temp1[,1],temp1[,2],ticknox = 5,ticknoy = 3,xlabel = 'signal',ylabel = 'SD',eq = function(x) 0.04 + 0.297*x )
  write.xlsx(as.data.frame(temp1),file='figs2.xlsx',sheetName = 'figs2c',append = F)
  
  # signal vs sd for mouse; S2D
  temp <- campSigNoiseComp(tst,op=1,result = 2)
  temp1 <- joinUnevenDFs(temp)
  fploteq(temp1[,1],temp1[,2],ticknox = 5,ticknoy = 3,xlabel = 'signal',ylabel = 'SD',eq = function(x) 3.19+0.52*x )
  write.xlsx(as.data.frame(temp1),file='figs2.xlsx',sheetName = 'figs2d',append = T)
  write.csv(as.data.frame(temp1),file='figs2d.csv')
  

  #fly, S2A, CV, mean= 0.4, ff, 0.22
  res5 <- campSigNoiseFF(tmp,op=2)
  fstripchartvecs(res5,markersize = 0.4,semthick = 0.6,tickno = 5,markersizemean = 0.8)
  write.xlsx(convertUnequalListToDf(res5),file='figs2.xlsx',sheetName = 'figs2a1',append = T)
  
  res5 <- campSigNoiseFF(tmp,op=3)
  fstripchartvecs(res5,markersize = 0.4,semthick = 0.6,tickno = 5,markersizemean = 0.8)
  write.xlsx(convertUnequalListToDf(res5),file='figs2.xlsx',sheetName = 'figs2a2',append = T)
  
  #mouse, S2B cv and ff by freq. CV, mean= 0.76, ff, 0.71
  res4 <- campSigNoiseFF(tst,op=2)
  fstripchartvecs(res4,markersize = 0.4,semthick = 0.6,tickno = 5,markersizemean = 0.8)
  write.xlsx(convertUnequalListToDf(res4),file='figs2.xlsx',sheetName = 'figs2b1',append = T)
  
  res4 <- campSigNoiseFF(tst,op=3)
  fstripchartvecs(res4,markersize = 0.4,semthick = 0.6,tickno = 5,markersizemean = 0.8)
  write.xlsx(convertUnequalListToDf(res4),file='figs2.xlsx',sheetName = 'figs2b2',append = T)
  
}


#plotting the correaltions between OSNS and KCs, ML/Classifier algorithms applied to data
stoch_code_figS3 <- function(){
  
  #the figures for the supplement Fig. S3a-c
  #first run the code in the function stoch_fig3 to populate the ressig, resrel, and resunrel structures.
  fploteq(ressig[,2],ressig[,3],eq = function(x) 0.16 + 0.32*x,ticknox = 4,ticknoy = 5,fixy=c(-.1,.8),fixx = c(-.3,0.9),rndfact = c(.2,2),xlabel = 'hallem',ylabel = 'rob' )
  write.xlsx(as.data.frame(ressig[,2:3]),file='figs3.xlsx',sheetName = 'figs3a',append = T)
  
  fploteq(resrel[,2],resrel[,3],eq = function(x) 0.17 + 0.44*x,ticknox = 4,ticknoy = 5,fixy=c(-.1,.8),fixx = c(-.3,0.9),rndfact = c(.2,2),xlabel = 'hallem',ylabel = 'rob' )
  write.xlsx(as.data.frame(resrel[,2:3]),file='figs3.xlsx',sheetName = 'figs3b',append = T)

  fploteq(resunrel[,2],resunrel[,3],eq = function(x) 0.04 + 0.04*x,ticknox = 4,ticknoy = 10,fixy=c(-.1,0.2,0.5,.8),fixx = c(-.3,0.1,0.5,0.9),rndfact = c(.4,1),xlabel = 'OSNs',ylabel = 'KCs',fixop = 1 )
  write.xlsx(as.data.frame(resunrel[,2:3]),file='figs3.xlsx',sheetName = 'figs3c',append = T)
  
  #classifier figures: G,H
  #classifier results for applying kNN and SVM to the data
  #flies and mice; classification accuracy: G
  #the data; accuracy
  temp1 <- rbind.data.frame(c(0.6,0.2,0.4),c(0.57,0.28,0.5),c(0.4,0.1,0.5),c(0.64,0.21,0.36),c(0.55,0.2,0.4),c(0.64,0.14,0.42))
  #f: flies, m:mouse
  rownames(temp1) <- c('lda.m','lda.f','knn.m','knn.f','svm.m','svm.f')
  colnames(temp1) <- c('Rel.','Unrel.','all')
  stackedHorzBarPlot(temp1,horz = T,spaces = 0.2,ticknox = 5,fontsize = 0.9,sepwidth = 1)
  #flies
  HorzBarPlot(temp1[4,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  HorzBarPlot(temp1[6,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  #mouse
  HorzBarPlot(temp1[3,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  HorzBarPlot(temp1[5,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  write.xlsx(temp1[3:6,c(2,1,3)],file='figs3.xlsx',sheetName = 'figs3g',append = T)
  
  #Fig. S3H; auc
  temp2 <- rbind.data.frame(c(0.8,0.59,0.64),c(0.72,0.72,0.83),c(0.74,0.39,0.59),c(0.71,0.36,0.65),c(0.88,0.7,0.84),c(0.93,0.66,0.82))
  rownames(temp2) <- c('lda.m','lda.f','knn.m','knn.f','svm.m','svm.f')
  colnames(temp2) <- c('Rel.','Unrel.','all')
  #flies
  HorzBarPlot(temp2[4,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  HorzBarPlot(temp2[6,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  #mouse
  HorzBarPlot(temp2[3,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  HorzBarPlot(temp2[5,c(2,1,3)],horz = F,spaces = 0.1,ticknox = 1,sepwidth = 1,color = 'c',sem = F,op=3,outputop = 0)
  write.xlsx(temp2[3:6,c(2,1,3)],file='figs3.xlsx',sheetName = 'figs3h',append = T)
  
  #Fig. S3D 
  #computing the classifier performance of accuracy and AUC versus correlation: Figure S3D: figure classifier AUC validation. Figure XX
  res6 <- campDetCorVsClass(tst,classifier = 10)
  cor(res6[[5]][lower.tri(res6[[5]],diag = T)],res6[[3]][[6]][lower.tri(res6[[3]][[6]],diag = T)])
  
  #comparing how odors that are similar (correlation) have low AUC and vice versa to show AUC's ability to distinguish odors.
  tst1 <- res6[[5]][lower.tri(res6[[5]],diag = T)] # corr
  tst2 <- res6[[3]][[6]][lower.tri(res6[[3]][[6]],diag = T)] #auc
  tst3 <- which(tst2 < 0.25) #get all the AUCs less thatn 0.25, i.e., similar ones
  #D: similar
  fstripchartvecs(list(tst1[tst3],tst2[tst3]),tickno = 2.5,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  write.xlsx(cbind.data.frame(tst1[tst3],tst2[tst3]),file='figs3.xlsx',sheetName = 'figs3d',append = T)
  
  
  tst3 <- which(tst2 > 0.75)
  #D: dissimilar
  fstripchartvecs(list(tst1[tst3],tst2[tst3]),tickno = 2,fixy = c(0,1),ylabel = 'measure',pairplot = 1,methodstr = 'overplot',semthick = 1.2)
  write.xlsx(cbind.data.frame(tst1[tst3],tst2[tst3]),file='figs3.xlsx',sheetName = 'figs3d2',append = T)
  
  #Fig. S3E,F
  #correlations just within KCs
  #flies: allcells_vs_rucells_flies
  temp1 <- computeMatCor(computeListMatCorr(tmp,celltype = 4),computeListMatCorr(tmp,celltype = 3),op=3)
  temp2 <- computeMatCor(computeListMatCorr(tmp,celltype = 4),computeListMatCorr(tmp,celltype = 2),op=3)
  fploteq(list(temp1,temp2),eq = list(function(x) -0.08 + 0.27*x,function(x) -0.5 + 1.74*x),fixy = c(-0.5,1),ticknox = 5,fixx = c(0,1),rndfact = c(5,5),xlabel = 'all cells',ylabel = 'R/U cells')
  test1 <- cbind.data.frame(temp1[[1]],temp1[[2]],temp2[[2]])
  names(test1) <- c('all cells','U','R')
  write.xlsx(test1,file='figs3.xlsx',sheetName = 'figs3e',append = T)
  write.csv(test1,file='figs3e.csv')
  
  #mice: allcells_vs_rucells_mice
  temp3 <- computeMatCor(computeListMatCorr(tst,celltype = 4),computeListMatCorr(tst,celltype = 3),op=3)
  temp4 <- computeMatCor(computeListMatCorr(tst,celltype = 4),computeListMatCorr(tst,celltype = 2),op=3)
  fploteq(list(temp3,temp4),eq = list(function(x) -0.01 + 0.1*x,function(x) -0.22 + 1.61*x),fixy = c(-0.5,1),ticknox = 5,fixx = c(0,1),rndfact = c(5,5),xlabel = 'all cells',ylabel = 'R/U cells')
  test2 <- cbind.data.frame(temp3[[1]],temp3[[2]],temp4[[2]])
  names(test2) <- c('all cells','U','R')
  write.xlsx(test2,file='figs3.xlsx',sheetName = 'figs3f',append = T)
  write.csv(test2,file='figs3f.csv')
  
}

#Fig. S4: showing that overlap increases with reliability
stoch_code_figS4 <- function(){
  #showing that overlap increases with reliability: Suppl Fig. Reliability vs overlap
  #flies: A
  res5 <- campComputeOverlap(tmp,threshop = 0)
  fstripchartvecs(res5,markersize = 0.8,ylabel = 'overlap',semthick = 2)
  write.xlsx(convertUnequalListToDf(res5),file='figs4.xlsx',sheetName = 'figs4a',append = T)
  
  #mouse: B
  res6 <- campComputeOverlap(tst,threshop = 0)
  fstripchartvecs(res6,markersize = 0.8,ylabel = 'overlap',semthick = 2)
  write.xlsx(convertUnequalListToDf(res6),file='figs4.xlsx',sheetName = 'figs4b',append = T)
  
  #calculates the number of cells in each reliability class
  #flies: D
  test2 <- countNoRelClass(tmp,op=2)
  #nocells_vs_reliability_fly
  fploteq(1:test2[[3]],test2[[1]],stddev = test2[[2]],xlabel = 'reliability',ylabel = '# cells',fixx = c(1,6),ticknoy = 2,rndfact = c(1,5),stdthick = 4)
  res2 <- cbind.data.frame(1:test2[[3]],test2[1:2])
  names(res2) <- c('reliability','# cells','sem')
  write.xlsx(res2,file='figs4.xlsx',sheetName = 'figs4d',append = T)
  
  #mouse: E
  test3 <- countNoRelClass(tst,op=2)
  #nocells_vs_reliability_mouse
  fploteq(1:test3[[3]],test3[[1]],stddev = test3[[2]],xlabel = 'reliability',ylabel = '# cells',fixx = c(1,8),ticknoy = 2,rndfact = c(1,5),stdthick = 4)
  res2 <- cbind.data.frame(1:test3[[3]],test3[1:2])
  names(res2) <- c('reliability','# cells','sem')
  write.xlsx(res2,file='figs4.xlsx',sheetName = 'figs4e',append = T)
  
  #flies: F
  test7 <- campCompRelOverlap(tmp)
  fploteq(1:6,test7[[3]][1:6],rndfact = c(1,0.25),xlabel = 'reliability',ylabel = 'bias',fixx = c(1:6),ticknoy = 2,fixy = c(0,10) )
  res2 <- cbind.data.frame(1:6,test7[[3]][1:6])
  names(res2) <- c('reliability','bias')
  write.xlsx(res2,file='figs4.xlsx',sheetName = 'figs4f',append = T)
  
  #mouse: G
  test7 <- campCompRelOverlap(tst)
  fploteq(1:8,test7[[3]][1:8],rndfact = c(1,0.25),xlabel = 'reliability',ylabel = 'bias',fixx = c(1:8),ticknoy = 4,fixy = c(0,80) )
  res2 <- cbind.data.frame(1:8,test7[[3]][1:8])
  names(res2) <- c('reliability','bias')
  write.xlsx(res2,file='figs4.xlsx',sheetName = 'figs4g',append = T)
  
}

#S5 shows the frequency plots of the overlaps
stoch_code_figS5 <-function(){
  #fly
  res5 <- campComputeOverlap(tmp,op=2,thresh = 0.5,threshop = 1)
  res4 <- joinListDFs(res5)
  res6 <- campComputeOverlap(tmp,op=2,thresh = 0.15,threshop = 2)
  res7 <- joinListDFs(res6)
  #together, same plot, S5a
  NlsFitCDFList(list(res4[,2],res7[,2]),dist = 2,graphparams = list(ticknox=1,ticknoy=1,rndfact=c(5,2),xlabel='overlap',ylabel='cum. freq.',fixx=c(0,1),markersize=0.8))
  #S5b
  fplotFn(eq = list(function(x) dgamma(x,shape = 0.15,scale = 0.44),function(x) dgamma(x,shape = 0.21,scale = 0.67)),ticknox = 0.2,tickunit = 1,ticknoy = 2,xlabel = 'overlap',ylabel = 'probability')
  tst4 <- RelCumFreqList(list(res4[,2],res7[,2]))
  write.xlsx(as.data.frame(tst4[[1]]),file='figs5.xlsx',sheetName = 'figs5a1',append = T)
  write.xlsx(as.data.frame(tst4[[2]]),file='figs5.xlsx',sheetName = 'figs5a2',append = T)
  
  #mouse
  res1 <- campComputeOverlap(tst,op=2,thresh = 0.5,threshop = 1)
  res2 <- joinListDFs(res1)
  res3 <- campComputeOverlap(tst,op=2,thresh = 0.15,threshop = 2)
  res4 <- joinListDFs(res3)
  
  #S5d
  NlsFitCDFList(list(res2[,2],res4[,2]),dist = 2,graphparams = list(ticknox=1,ticknoy=1,rndfact=c(5,2),xlabel='overlap',ylabel='cum. freq.',fixx=c(0,1),markersize=0.8))
  #S5e
  fplotFn(eq = list(function(x) dgamma(x,shape = 0.24,scale = 0.09),function(x) dgamma(x,shape = 0.22,scale = 0.45)),ticknox = 0.2,tickunit = 1,ticknoy = 2,xlabel = 'overlap',ylabel = 'probability')
  tst4 <- RelCumFreqList(list(res2[,2],res4[,2]))
  write.xlsx(as.data.frame(tst4[[1]]),file='figs5.xlsx',sheetName = 'figs5d1',append = T)
  write.xlsx(as.data.frame(tst4[[2]]),file='figs5.xlsx',sheetName = 'figs5d2',append = T)
  
    
  #looking at fly and mouse correlations distributions, S5cf
  #S5c; fly  
  temp <- computeListMatCorr(tmp,matop = 2,op=1)
  NlsFitCDF(temp[lower.tri(temp)],dist = 5,graphparams = list(ticknox=2,ticknoy=4,fixx=c(-0.2,0.6),rndfact=c(10,5),xlabel='correlations',ylabel='cum. freq'))
  write.xlsx(as.data.frame(RelCumFreq(temp[lower.tri(temp)])),file='figs5.xlsx',sheetName = 'figs5c',append = T)
  
  
  #S5f; mouse  
  temp1 <- computeListMatCorr(tst,matop = 2,op=1)
  NlsFitCDF(temp1[lower.tri(temp1)],dist = 5,graphparams = list(ticknox=2,ticknoy=4,fixx=c(-0.2,0.6),rndfact=c(10,5),xlabel='correlations',ylabel='cum. freq'))
  write.xlsx(as.data.frame(RelCumFreq(temp1[lower.tri(temp1)])),file='figs5.xlsx',sheetName = 'figs5f',append = T)
  
}

#showing that the reliability of cells is on a spectrum
stoch_code_figS6 <- function(){
  
  #showing that unreliable cells are on a spectrum; flies: A
  temp <- campGetUnrelDist(tmp) #will get the unreliable cells, where each list element gives the amount of times the
  #cell fires again. All ordered by the first time the cell fires in the first 3 trials
  temp2 <- joinUnevenDFs(temp[1:3])
  unrel.fly <- rbind(temp[[4]],joinUnevenDFs(temp[3:1]))
  drawMatGridPlot(unrel.fly,grayscale = T)
  write.xlsx(as.data.frame(unrel.fly),file='figs6.xlsx',sheetName = 'figs6a',append = T)
  
  #mouse: B
  temp <- campGetUnrelDist(tst)
  temp1 <- rbind(joinUnevenDFs(temp[1:3]),temp[[4]])
  #reverse:
  temp3 <- rev(temp)
  unrel.mouse <- rbind(joinUnevenDFs(temp3))
  write.xlsx(as.data.frame(unrel.mouse),file='figs6.xlsx',sheetName = 'figs6b',append = T)
  
  #16 trials flies: C
  #setwd("/nadata/cnl/data/shyam/campbell/data_dir/camptest/data_dir/101007/tseries-10062010-1245-001dir")
  #res.101007_1 <- getCampDataSigList()
  tmp1 <- res.101007_1
  temp <- campGetUnrelDist(tmp1)
  temp1 <- rbind(joinUnevenDFs(temp[1:3]),temp[[4]])
  #"reverse"
  temp3 <- rev(temp)
  temp1 <- rbind(joinUnevenDFs(temp3[2:4]))
  unrel.fly16 <- rbind.data.frame(temp3[[1]],temp1)
  drawMatGridPlot(temp1,grayscale = T)
  write.xlsx(as.data.frame(unrel.fly16),file='figs6.xlsx',sheetName = 'figs6c',append = T)
  write.csv(unrel.fly16,file='figs6c.csv')
  
  #theoretical explanations
  res1 <- campGenRelStats(nocells = 142,rangefreq = c(3:16),obs = c(98,33,10,1))
  res1 <- rbind(res1,c(98,33,10,1))
  #rownames(res1) <- c(rownames(res1)[1:15],'data')
  rownames(res1) <- c(paste('p=1/',rownames(res1)[1:14],sep = ''),rownames(res1)[15],'data')
  stackedHorzBarPlot(res1,horz = T,fontsize = 0.6,ticknox = 1,rndfact = 5,sepwidth = 1)
  write.xlsx(res1,file='figs6.xlsx',sheetName = 'figs6d',append = T)
  write.csv(res1,file='figs6d.csv')
  
  res2 <- campGenRelStats(nocells = 30,rangefreq = c(3:24),notrials = 16,obs = c(16,10,3,1,0,0,0,0,0,0,0,0,0,0),op=1)
  res2 <- rbind(res2,c(16,10,3,1,0,0,0,0,0,0,0,0,0,0))
  rownames(res2) <- c(paste('p=1/',rownames(res2)[1:22],sep = ''),'avg','data')
  #rownames(res1) <- c(rownames(res1)[1:23],'data')
  stackedHorzBarPlot(res2,horz = T,fontsize = 0.6,ticknox = 2,rndfact = 5,sepwidth = 1)
  write.csv(res2,file='figs6f.csv')
  
  
  res3 <- campGenRelStats(nocells = 566,rangefreq = c(3:16),notrials = 8,obs = c(314,142,66,24,13,7))
  res3 <- rbind(res3,c(314,142,66,24,13,7))
  rownames(res3) <- c(paste('p=1/',rownames(res3)[1:14],sep = ''),'avg','data')
  stackedHorzBarPlot(res3,horz = T,fontsize = 0.6,ticknox = 4,rndfact = 5,sepwidth = 1)
  write.csv(res3,file='figs6e.csv')

}

#compares the correlations between pairs of odors with the cosine similarities of the top 25, 26-50, 51-75, 76-100 percentile of cells
stoch_code_figs8 <- function(){
  #Fig. S8
  #fly:corr_cosine_different_topk_fly
  res3 <- campOverSimRankData(alldata = tmp)
  res3 <- campOverSimRankData(alldata = tmp,respop = 2) #avg. as E(X)
  res4 <- getFitEqnList(res3)
  fploteq(x=res3,ticknox = 3,ticknoy = 2,xlabel = 'correlations',ylabel = 'cosine similarity',eq = res4)
  test3 <- joinListDFs(res3)
  write.xlsx(test3,file='figs8.xlsx',sheetName = 'figs8a',append = T)
  
  #mouse:corr_cosine_different_topk_mouse
  res3 <- campOverSimRankData(alldata = tst)
  res3 <- campOverSimRankData(alldata = tst,respop = 2)  #avg. as E(X)
  res4 <- getFitEqnList(res3)
  fploteq(x=res3,ticknox = 3,ticknoy = 2,xlabel = 'correlations',ylabel = 'cosine similarity',eq = res4)
  write.xlsx(joinListDFs(res3),file='figs8.xlsx',sheetName = 'figs8b',append = T)
  
}

#gets the data for the effect of significance level on discrimination analysis
stoch_code_figS9 <- function(){
  
  #Fig. S9A; reliable vs unreliable cells changing significance threhsold
  #the first command has to be carried out in the directory with 
  setwd("/path/to/dir")
  res <- campGetAlphaEffectRUCells(alpharange = c(0.1,0.05,0.01,0.005,0.001,0.0005),op=3)
  res1 <- campCompareAlphaRUcells(res,compst = 2,compend = 3,op=2)
  res2 <- campCompareAlphaGroups(res)
  res5 <- convertNestedListsDF(res2)
  #reliable_vs_unreliable_cellslost_change_alpha
  fploteq(unlist(res5[2,]),unlist(res5[3,]),ticknox = 5,ticknoy = 5,xlabel = 'reliable',ylabel = 'unreliable',eq=function(x) 2 + 5.44*x )
  write.xlsx(res5[2:3,],file='figs9.xlsx',sheetName = 'figs9a',append = T)
  
  
  #Fig. S9B,C
  res5 <- campComputeUDAlphas(cols = c(3,4),alpharange=c(0.1,0.05,0.01,0.005,0.001,0.0005),meas.sim = c(0.15,0.4))
  res6 <- sapply(seq(1,11,2), function(i) c(sum(res5[c(i,i+1),1]),sum(res5[c(i,i+1),2])) )
  res7 <- sapply(seq(1,11,2), function(i) c(sum(res5[c(i,i+1),3]),sum(res5[c(i,i+1),4])) )
  res8 <- sapply(seq(1,11,2), function(i) sum(res5[c(i:(i+1)),3])/sum(res5[c(i:(i+1)),1]))
  #Fig. S9Bdissim_normal_vs_sim_sat_changing_alpha
  fploteq(c(0.1,0.05,0.01,0.005,0.001,0.0005),res6[1,]/res7[2,],ticknox = 2,fixy = c(0,2),ticknoy = 1,logs = c(T,F),xlabel = 'significance',ylabel = 'relative discrimination')
  test1 <- cbind.data.frame(c(0.1,0.05,0.01,0.005,0.001,0.0005),res6[1,]/res7[2,])
  names(test1) <- c('alpha','relative discrimination')
  write.xlsx(test1,file='figs9.xlsx',sheetName = 'figs9b',append = T)
  write.csv(test1,file='figs9b.csv')
  
  #Fig. S9C
  fploteq(c(0.1,0.05,0.01,0.005,0.001,0.0005),res8,ticknox = 2,fixy = c(0,1),ticknoy = 2,logs = c(T,F),xlabel = 'significance',ylabel = 'sim./dissim UD')
  test2 <- cbind.data.frame(c(0.1,0.05,0.01,0.005,0.001,0.0005),res8)
  names(test2) <- c('alpha','sim/dissim UD')
  write.csv(test1,file='figs9c.csv')
  
  #mouse:
  setwd("/path/to/dir")
  #fig D
  res <- simGetAlphaEffectRUCells(alpharange = c(0.1,0.05,0.01,0.005,0.001,0.0005))
  res1 <- transposeDF(as.data.frame(res))
  #Fig. S9D
  fploteq(res1[,1],res1[,2],xlabel = 'reliable',ylabel = 'unreliable',ticknoy = 2,eq = function(x) 0.3 + 2.2*x,ticknox = 5)
  write.csv(res1,file='figs9d.csv')
  
  res5 <- simComputeUDAlphas(cols = c(3,4),alpharange = c(0.1,0.05,0.01,0.005,0.001,0.0005))
  res7 <- sapply(seq(1,11,2), function(i) c(sum(res5[c(i,i+1),3]),sum(res5[c(i,i+1),4])) )
  #Fig. S9E
  fploteq(c(0.1,0.05,0.01,0.005,0.001,0.0005),res7[2,]/res7[1,],ticknox = 2,fixy = c(0,3),ticknoy = 1,logs = c(T,F),xlabel = 'significance',ylabel = 'relative discrimination')
  test2 <- cbind.data.frame(c(0.1,0.05,0.01,0.005,0.001,0.0005),res7[2,]/res7[1,])
  names(test2) <- c('alpha','relative discrimination')
  write.csv(test2,file='figs9e.csv')
  
  #Fig. S9F
  res9 <- sapply(seq(1,11,2), function(i) sum(res5[c(i:(i+1)),1])/sum(res5[c(i:(i+1)),4]))
  fploteq(c(0.1,0.05,0.01,0.005,0.001,0.0005),res9,ticknox = 2,fixy = c(0,2),ticknoy = 2,logs = c(T,F),xlabel = 'significance',ylabel = 'relative discrimination')
  test2 <- cbind.data.frame(c(0.1,0.05,0.01,0.005,0.001,0.0005),res9)
  names(test2) <- c('alpha','sim/dissim UD')
  write.csv(test2,file='figs9f.csv')
  
}


#supplementary modeling figure supporting the conclusions of the main Figure, Fig. 3
stoch_code_figS10 <- function(){
  
  #Fig. S10A1
  temp12 <- exploreNoiseParams(noruns = 10,notrials = 6,noiseset = 2,noiserange = tst3,aplrange = c(5:9),op=22,classop = 2,kcapl_conn_par=tst4,aplkc_conn_par=list(8,c(20,3,0.4,0.45)));
  res2 <- cirAnalGoodParams(cirGetParamsTargetAllApls(temp12,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 ),op=2,par = 1)
  NlsFitCDF(res2[[5]][,2],graphparams = list(xlabel='apl-threshold noise',ylabel='freq.',ticknoy=2,ticknox=1),dist = 1)
  
  #Fig. S10A2 
  temp5 <- exploreNoiseParams(noruns = 8,notrials = 6,noiseset = 2,noiserange = list(seq(0,0.5,.025),seq(0,0.5,0.1),seq(0,0.4,.05),seq(0,0.4,0.05)),aplrange = c(5:12),op=22,classop = 2)
  test5 <- cirGetParamsTargetAllApls(temp5,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 )
  NlsFitCDF(res5[[5]][,4],graphparams = list(xlabel='apl-kc noise',ylabel='freq.',ticknoy=2),dist = 5)
  
  #Fig. S10b
  test5 <- cirGetParamsTargetAllApls(temp5,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 )
  res5 <- cirAnalGoodParams(test5,op=2,par = 5,noiseop = 22)
  fploteq(res5[[5]][,3],res5[[5]][,5],xlabel = 'apl noise',ylabel = 'pn-kc noise')
  
  #Fig. S10C
  temp14 <- exploreNoiseParams(noruns = 10,notrials = 6,noiseset = 2,noiserange = list(seq(0.2,0.45,0.025),seq(0,0.1,0.05),seq(0,0.1,0.1),seq(0,0.1,0.05),seq(0,0.5,0.1)),aplrange = c(6:10),classop = 2,kcapl_conn_par=list(7,c(20,2,0.4,0.45)),aplkc_conn_par=list(8,c(20,2,0.4,0.45)),op=30)
  res14 <- cirAnalGoodParams(cirGetParamsTargetAllApls(temp14,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 ),op=2,par = 1,noiseop = 30)
  HorzBarPlot(100-res14[[2]][2:6],color = 's',sem = F,op=2,ticknox = 1.5)
  
  #Fig. S10D, left and right
  temp6 <- exploreNoiseParams(noruns = 8,notrials = 6,noiseset = 2,noiserange = list(seq(0.3,0.5,0.01),seq(0,0.1,0.05),seq(0,0.1,0.05),seq(0,0.1,0.05)),aplrange = c(5:9),op=22,classop = 2)
  test6 <- cirGetParamsTargetAllApls(temp6,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 )
  tst6 <- cirAnalGoodParams(test6,op=2,par=1)
  #left
  NlsFitCDF(tst6[[5]][,3],graphparams = list(xlabel='apl noise',ylabel='freq.',ticknoy=2))
  
  temp7 <- exploreNoiseParams(noruns = 8,notrials = 6,noiseset = 2,noiserange = list(seq(0.3,0.5,0.025),seq(0,0.1,0.05),seq(0,0.1,0.05),seq(0,0.1,0.05)),aplrange = c(5:9),op=22,classop = 2,kcapl_conn_par=list(7,c(10,2,0.4,0.45)),aplkc_conn_par=list(8,c(20,2,0.4,0.45)))
  test7 <- cirGetParamsTargetAllApls(temp7,targetpars = list(c(0.72,5.26,7.23,6.1,29)),op=1,parno = 0 )
  res1 <- cirAnalGoodParams(test7,op=2,par = 1)
  #right
  NlsFitCDF(res1[[5]][,3],graphparams = list(xlabel='apl noise',ylabel='freq.',ticknoy=2))
  
  #Fig. S10 exploting # synapses vs % top,   #top_match_vs_synapses_3_distributions
  res3 <- sapply(seq(2,40,4),function(i) cirExploreWTAparams(params = c(i,3,.2,0.12),op=3) )
  res4 <- sapply(seq(2,40,4),function(i) cirExploreWTAparams(params = c(i,1,1.15,0.12),op=3) )
  res5 <- sapply(seq(2,40,4),function(i) cirExploreWTAparams(params = c(i,2,0.4,0.5),op=3) )
  fploteq(seq(2,40,4),cbind.data.frame(res3[1,],res4[1,],res5[1,]),ticknox = 2,ticknoy = 2,xlabel = '# synapses',ylabel = '% top')

}



#Fig. S11: testing some examples from parameter exploration, and showing the effect of number of cells on the results  
stoch_code_figS11 <- function(){
  #examples for paper
  test3 <- createFlyCircuit(glomno = 50,kcno = 2000,aplCut = 9,kcapl_conn_par=list(7,c(10,2,0.4,0.45)),aplkc_conn_par=list(8,c(10,2,0.4,0.45)),noodors = 6,classop = 2,thresh = list(1,c(2,10)))
  test4 <- exploreNoise(test3,op=22,noiseparams = list(list(1,c(0,0.35),1),list(1,c(0,0.1),1),list(1,c(0,0.1),1),list(1,c(0,0),1)),notrials = 6)
  test5 <- cirEvalSimResp(test4,meas.sim = c(0.15,0.3));test5[1:2];
  #Fig. S11A-D: effect of # KCs
  test.examples <- list(test3,test4) #kc = 150
  test5 <- cirEvalSimResp(test4,meas.sim = c(0.15,0.3));test5[1:2];test.examples <- c(test.examples,list(test3,test4)) #kc=250
  test5 <- cirEvalSimResp(test4,meas.sim = c(0.15,0.3));test5[1:2];test.examples <- c(test.examples,list(test3,test4)) #kc=500
  test5 <- cirEvalSimResp(test4,meas.sim = c(0.15,0.3));test5[1:2];test.examples <- c(test.examples,list(test3,test4)) #kc=1000
  test5 <- cirEvalSimResp(test4,meas.sim = c(0.15,0.3));test5[1:2];test.examples <- c(test.examples,list(test3,test4)) #kc=2000
  #plotting them, 150, 500, 2000
  tmp1 <- campGetFitVecsData(test.examples[[2]])
  tmp2 <- campGetFitVecsData(test.examples[[6]])
  tmp3 <- campGetFitVecsData(test.examples[[10]])
  #cumulative_overlap_150_500_2000
  NlsFitCDFList(list(tmp1[[3]],tmp2[[3]],tmp3[[3]]),dist = 2,colorop = 2,graphparams = list(ticknox=2,ticknoy=10,rndfact=c(1,5),ylabel='frequency',xlabel='overlap'))
  #cumulative_overlap_150_500_2000
  NlsFitCDFList(list(tmp1[[2]],tmp2[[2]],tmp3[[2]]),dist = 2,colorop = 2,graphparams = list(ticknox=1,ticknoy=1,rndfact=c(5,5),ylabel='frequency',xlabel='reliability'))
  #cumulative_response_150_500_2000
  NlsFitCDFList(list(tmp1[[1]],tmp2[[1]],tmp3[[1]]),dist = 2,colorop = 2,graphparams = list(ticknox=2,ticknoy=5,rndfact=c(5,5),ylabel='frequency',xlabel='response levels',markersize=0.5))
  
  
  #other scharacteristics
  tmp2 <- lapply(seq(2,10,2),function(i) cirEvalSimResp(test.examples[[i]],meas.sim = c(0.15,0.3))[[4]])
  tmp3 <- lapply(1:length(tmp2), function(i) tmp2[[i]]/c(0.72,6.1,29,4.2,4))
  tmp4 <- transposeList(tmp3)
  #Fig. A; changing_kcs_150_to_2000
  fstripchartvecs(tmp4,markersize = 0.8,sem = 0,tickno = 5,fixy = c(0,1.5))
  #result: check later if the if the reliability goes up with more cells.

  #Fig. E: good; bad examples:
  test1 <- createFlyCircuit(glomno = 50,kcno = 150,aplCut = 8,kcapl_conn_par=list(7,c(20,2,0.4,0.45)),aplkc_conn_par=list(8,c(20,2,0.4,0.45)),noodors = 6,classop = 2)  
  test2 <- exploreNoise(test1,op=30,noiseparams = list(list(1,c(0,0.24),1),list(1,c(0,0.1),1),list(1,c(0,0.1),1),list(1,c(0,0),1),list(1,c(0,0.1),1)),notrials = 6); res1 <- cirEvalSimResp(test2,meas.sim = c(0.15,0.3));res1[[4]];
  test1.examples <- list(test1,test2); #perfect parameters
  test1.examples <- c(test1.examples,list(test1,test2)) # very little noise, urel.sum: 0.41, 1.51
  test1.examples <- c(test1.examples,list(test1,test2)) # too much noise rel.sum: 4, 1.71
  #tmp 1.2, 5.32  5.32 3.28; 1 to 4 & 2.5 to 1.76 

  res2 <- lapply(seq(2,6,2),function(i) cirEvalSimResp(test1.examples[[i]],meas.sim = c(0.15,0.3))[[4]])
  res3 <- lapply(1:length(res2), function(i) res2[[i]]/c(0.72,6.1,29,4.2,4))
  res4 <- transposeList(res3)
  #good_bad_examples_comparing_properties
  fstripchartvecs(res4[1:4],markersize = 0.8,sem = 0,tickno = 5,fixy = c(0,5))
  
}


#simulations for generating the modeling plots 3, 
stoch_code_fig3S6 <- function() {
  res <- createFlyCircuit(glomno = 50,kcno = 100,noodors = 6,aplCut = 10)
  res1 <- exploreNoise(res,notrials = 6,noiseparams = list(1,c(0,.6),1),noiseset = 2,resop = 1,op=1)

  #Figs. 3,S6
  #glom noise
  res.glomnoise <- exploreNoiseSeq(res,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1))
  #pmmb
  res.glommbnoise <- exploreNoiseSeq(res,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=2)
  #kcapl
  res.kcaplnoise <- exploreNoiseSeq(res,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=3)
  #apl=kc
  res.aplkcnoise <- exploreNoiseSeq(res,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=4)
  #pn,pn-mb
  res.pnpnmbnoise <- exploreNoiseSeq(res,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=10)
  #pn,pn-mb, and apl-kc
  res.pn_pnmb_aplkcnoise <- exploreNoiseSeq(res,notrials = 6,noiseset = 2,noiserange = seq(0,1,.1),op=10)
  
  
}

#gets the dataset statistics for all datasets in Tables S3 and 4
stoch_code_tableS34 <- function(){
  
  #getting the basic data for all datasets
  #mouse
  res2 <- campGetReliabStatsLists(list(tst,mouse.163.data,mouse.164.data,mouse.7.data,mouse.8.data,mouse.9.data),respop = 2,op=2)
  names(res2) <- c('tst','mouse.163.data','mouse.164.data','mouse.7.data','mouse.8.data','mouse.9.data')
  transposeDF(convertNestedListsDF(res2))[-1,-(1:3)][c(2,1,3,4,5),]
  xtable::xtable(transposeDF(convertNestedListsDF(res2))[-1,-c(1:3,9)][c(2,1,3,4,5),])
  #fly
  res3 <- campGetReliabStatsLists(list(tmp,fly.01052011.001,fly.01052011.001.110108,fly.01052011.001.110109,fly.01052011.001.1101092,fly.09042009.001,fly.10092009.012,fly.01052011.002),respop = 2,op=2)
  names(res3) <- c('tmp','fly.01052011.001','fly.01052011.001.110108','fly.01052011.001.110109','fly.01052011.001.1101092','fly.09042009.001','fly.10092009.012','fly.01052011.002')
  transposeDF(convertNestedListsDF(res3))[-1,-(1:3)][c(6,1:5,7),]
  xtable::xtable(transposeDF(convertNestedListsDF(res3))[-1,-c(1:3,9)][c(6,1:5,7),])
  
}