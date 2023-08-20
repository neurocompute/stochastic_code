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

#main modeling figure 3, and supplementary modeling figure S7
#this function contains the code to generate the modeling results of main Figure 3 exploring the effect of noise 
#in various components on circuit responses
stoch_code_fig3 <- function(){
  
  #to run the code in this function, you will need to source two files, separately. For Figures 3C-E and S7A-C, please run
  #circuits0.R, and then the code here. This file contains the first model with a single synapse from APL back to KCs 
  #For all other figures, source circuits.R: this file contains the code the second model that explores the effect of multiple synapses
  #from APL to every KC.
  
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
  
  #For analysis of the data from the Lin et al. paper presented in Figure S7, use the following code on the Lin Data.
  #The data can be obtained by writing to the authors of Lin et al., 2014
  #Fig. S7H,I
  tmp3 <- getLinMatList()
  #I: plot all cumulative lists
  tmp4 <- RelCumFreqList(computeLinMatCum(tmp3[[2]])) #get the odors as RelcumLists
  #fit to average
  NlsFitCDF(unlist(computeLinMatCum(tmp3[[2]])),dist=2)
  tmpfn <- function(x) pgamma(x,shape = 6.11,scale = 13.98) #assign it to a fuction
  #plot that functio now
  fploteq(tmp4,eq = tmpfn,plottype = 8,markersize = .75,xlabel = 'odor responses',ylabel = 'cum freq',ticknox = 3,ticknoy = 5)
  #I: get one from the list and plot it. e.g., tmp4[[1]] will plot the first list.
  fploteq(tmp4[[1]],eq = tmpfn,plottype = 8,markersize = .75,xlabel = 'odor responses',ylabel = 'cum freq',ticknox = 3,ticknoy = 5)  
  
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


