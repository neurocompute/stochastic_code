
#circuits.R: code for implementing perceptrons and other neuronal networks: Shyam Srinivasan (C), shyam@snl.salk.edu
#version 1.1 contains code for a variable number of valences, i.e., MBONs and DANs. In this the number of MBONs and DANs should match.

#this file contains the classes that will be needed to evaluate the fly olfactory circuit and maybe even neural networks
#nomenclature: NameClass (all uppercase) for classes, nameFn (first letter lowercase) for functions, 
#namevar (all lowercase) for variables.


#setup constants here as a list: Neuron, connection matrix, and FlyMBg object
#Anytime these are changed, they have to be done so in three places: here, init, and set/get
#Specific rules for inserting. The numbers have to be contiguous because they serve as indices to the vector 
#of names, so you have to increase the number of everyone after your insertee.
#const.N$id gives you the structure's id inside const.FM for example res$apl$id = const.FM$apl
const.N <- list(neurons=1,neurondesc=2,dist=3,basenoise=4,neurontype=5,params=6,posns=7,noise=8,id=9,thresh=10)
const.CM <- list(connmat=1,type=2,source=3,target=4,params=5,size=6,syndist=7,noise=8,
                 learningrule=9,gain=10,dopno=11,dopparams=12,dopdist=13,doptype=14,dopconn=15,curtime=16,
                 synactfn=17)
const.FM <- list(glom=1,mbkc=2,mbkc_net=3,apl=4,dans=5,mbons=6,
                 glommb_conn=7,kcapl_conn=8,aplkc_conn=9,kcmbon_conn=10,noodors=11,valences=12,trainingval=13,
                 teststim=14,tatparams=15,origstim=16)
 
#teststim: contains the stimuli number that have been added for testing how well the system has learnt, as well as the ones that were 
#used for training. teststim = c(vector of odors to be used for testing). Initially both vectors are the same
#trainingval: the training valences for all the odors, i.e. odor 1 to valence 1, 2 to 2 and so on; = list(c(training valences),c(odors to be used for training)))
#tatparams: the parameters that describe how training takes place and testing takes place, for now traintype:1-alternate, 2 - sequence, 3 - random
#testtype: 1 - is valence learnt, 2 - pairwise discrimination. Overall c(traintype,testtype)
#origstim: is the original stimuli that does not contain any mixes, just the pure odorants that are likely to be randomly chosen
#noise: the amount of noise generated in this neuron or connection network of the form c(type of noise distribution,params,
#additive or multiplicative)

#Noise: how to do we deal with noise. The basic procedure has 3 steps. 
#1. you set the noise parameters with the function, setNoiseParams. This is the same for neurons and 
#connection matrices. It goes in, and assigns the noise params to the noise property. The op in 
#setnoiseparams chooses the components to which you must add noise. 
#2. Explore noise: In explorenoise, you call setnoiseparams, and then call the addnoise with the same 
#op specifying the structure that you want to add noise to. Then, you compute the rest of the function
#3. You will have to also make changes to the exploreNoiseSeq function
#4. Make changes to computeNoNoiseComp
#For connection matrices you don't have to do anything else. For neurons, you will have to add noise
#within the procedure 
#Procedure for neurons: everything the same until addnoise. Over there, add noise for connection matrices.
#For neurons, noise will have to be added in the computeFiring.neurons procedure.  



#number: number of neurons of this type
#neurontype: type of neuron 1- excitatory, 2 - inhibitory, 3 - modulatory
#neurondesc: descroption of neuron
#posns: the position of the neurons. 0 or single number means no posnal information, a vector or matrix indicates 
#posn in 1 or 2 dimensional space
#basenoise: set to true or False depending on whether this neuron should have noise or not 
#params: of the distribution, e.g., for dist=1, params=x means sd=x
#dist: distribution of firing rates for this neuron, 0 - constant number, 1- gaussian, 2 - exponential, 3 - gamma, 4 - binary
#5 uniform distribution
#10 - user specified for all neurons, which is given in params
#params: parameters for the firing rates based on distribution, for dist=0, one number
#for dist=10, a vector of firing rates and length(params) = number
#noStimuli: no of different stimuli that these set of neurons code, e.g., for 30 odors, this would be 30. default is 1
#noise: the amount of noise generated in this neuron or connection network of the form 
#c(type of noise, distribution,params,noistype, i.e., additive or multiplicative). new change c() to list()
#list(type e.g.,1-gaussoion,(mean,sd),0-multiplicative or 1 additive)
#additive or multiplicative), this is the kind of noise that changes with each trial, so initialize here, but change 
#only in computeFiring
#extraparams: extra params if needed, when specified do the distribution from this.
#op=1, 
#id: specifies the identity of the neuron, like const.FM$apl (which is 4)
#thresh: the threshold for the neuron to fire. c(type,no) assumes the same threshold for all neurons, but, can change later.
#here type: 1 - absolute thresh, i.e., sig - no, 2 - percentile thresh, no is between (0,100) i.e., sig - thresh, 
#where thresh is bottom no percentile,e.g., if no=5, no. that marks the end of the 5th percentile. 
newNeurons <- function(number=1,basenoise=F,neurontype=1,neurondesc='excitatory',posns=0,dist=1,params=c(10,1),exparams=c(),
                       noise=list(0,c(0,0),0),noStimuli=1,id=1,thresh=c(1,0),op=1){
  #generate the 
  neurons <- genNeurons(number = number,basenoise = basenoise,dist = dist,params = params,noStimuli = noStimuli,noise=noise)
  #posns, figure this one out when we actually use it. Might be tricky like figuring out the rostral caudal gradient
  if(length(posns)>1){#code for putting positional information
    
  }
  names(neurons) <- 1:noStimuli #name all the stimuli on the lists
  cells <- list(neurons,neurondesc,dist,basenoise,neurontype,params,posns,noise,id,thresh)
  names(cells) <- c('neurons','neurondesc','dist','basenoise','neurontype','params','posns','noise','id','thresh')
  res <- structure(cells,class = 'Neurons') #make it a class and return it
  #print(cells)
  #cat('\n new neurons',res$id,res$thresh,'done',res[[9]],res[[10]],'\n')
  res
}

#given a neuron object, will threshold all of them.
#neuron: an object of he class of Neurons
thresholdNeuron <- function(neuron,op=1){
  stim <- getData(neuron,const.N$neurons)
  thresh <- getData(neuron,const.N$thresh)
  if(thresh[2]==0) return(neuron) #0-threhsold, so nothing to do
  threshneuron <- setData(neuron,threshStimuli(stim,thresh),const.N$neurons)
  threshneuron
}

#given a list of neurons and threshold returns a thresholded list of neurons 
#thresh: the threshold for the neuron to fire. c(type,no) assumes the same threshold for all neurons, but, can change later.
#here type: 1 - absolute thresh, i.e., sig - no, 2 - percentile thresh, no is between (0,100) i.e., sig - thresh, 
#where thresh is bottom no percentile,e.g., if no=5, no. that marks the end of the 5th percentile. 
threshStimuli <- function(stimlist,thresh,op=1){
  if(thresh[2]==0) return(stimlist) #0-threhsold, so nothing to do
  #set the threshold, for the mean multiply by two since mean is presumably half
  #cat('\n sum stimuli',sum(unlist(stimlist)))
  thresh.value <- switch(thresh[1],thresh[2],getBottomNPercentileThresh(unlist(stimlist),bottomn = thresh[2],op=1) )
  res.stim <- lapply(stimlist, function(x){
    res <- threshValTarget(x - thresh.value) #threshold the neurons    
    res
  })
  #cat('\nthreshval',thresh.value,'\t',thresh)
  res.stim
}


#function that generates neurons according to the parameters specified. 
#The return is a bunch of neuronal firing rates as vectors 
#Better this way, as it can then be used 
#by various functions
#dist: 0 - all neurons are the same, 10 - neurons specified by params, 1 - gaussian, 2 - exponential
# 3 - gamma, 4 - poisson, 5 - uniform
# params: the parameters of the distribution
genNeurons <- function(number,basenoise=F,dist=1,params=c(10,1),noise=list(0,c(0,0),0),noStimuli=1){
  #the firing rate depends on the value of dist
  if(dist==0) {#all neurons are the same number
    if (length(params)==1) neurons <- lapply(1:noStimuli, function(i) neurons <- rep(params,number))
    else neurons <- lapply(1:noStimuli, function(i) neurons <- rep(params[i],number))
  }
  if(dist==10) {
    if (length(params)==1) neurons <- lapply(1:noStimuli, function(i) neurons <- params)
    else neurons <- lapply(1:noStimuli, function(i) neurons <- params)
  }
  if(dist>0 && dist<=5){#firing rate is based on a distribution
    neurons <- lapply(1:noStimuli, function(i) 
      switch(dist,rnorm(number,mean = params[1],sd=params[2]),
             rexp(number,rate = params[1]),
             rgamma(number,shape = params[1],scale=params[2]),
             rpois(number,lambda = params[1]),
             runif(n = number,min = params[1],max = params[2])) )
  }
  #now add noise
  if(basenoise){ #add noise to the neurons
    neurons <- lapply(1:noStimuli, function(i){
      temp <- neurons[i] *(1 + switch(noise[1],rnorm(length(neurons[i]),noise[[1]][1],noise[[2]][2]),
                                      runif(length(neurons[i]),min = noise[[1]][1],max = noise[[1]][2])))
      #make sure none of the firing rates are -ve
      neurons <- setVecThresh(temp,thresh = 0)
    })
    noise <- noise #this is the kind of noise that changes with each trial, so initialize here, but change only in computeFiring
  }
  neurons
  
}

#generates connection matrix, default Gamma, whose params are specified by params
#the matrix is size rows x cols
#source: the structure from which the connections are coming in, e.g., OB
#target: the structure to which the connections are going, e.g., PCx
#type: 1 - distributed, 2 - topographic
#nosyn - average number of contacts that an input neuron makes onto the output neuron population
#syndist - the synaptic distribution; 1- agaussian matrix, 2 -gamma distribution, 3 - compound poisson Gamma distribution
# 4 - uniform distribution, 5 - constant, number in params.
# 6 - each source and target neuron make multiple synapses, params= c(avg. no of synapses by each source neuron onto a target,connection distribution)
#if compound poisson gamma params = (poisson mean,shape,scale)
#7 - multiple synapse from each source neuron. and these are represented as a vector
#8 - the matrix connection, the number of synapses is given by the vector
#noise: the amount of noise generated in this neuron or connection network of the form c(type of noise distribution,params,
#additive or multiplicative), this is the kind of noise that changes with each trial, so initialize here, but change only in computeFiring
#params - params of the distribution
#origin - a string, the source population of neurons
#target - a string, target population of neurons
#dop=c() - no dopamine, not a learning matrix, default parameters:list(2,5,c(2,0.9,1.1),1)
#list(no of dopamine neurons, dopamine network type, params for the network connections,dopamine type=reward,punish)
# 1 for reward, -1 for punish
#learningrule: the rule for updating the strength of the synapses.(alpha,beta,ec50,T:sat fn,F:linear function) additive decrease s-alpha, multiplicative inc: s*beta 
#gain = 1, the gain of the feedforard or feedback connection. Baked into the connection matrix, but you can 
#get the original strength by 1/gain multiplication
newConnMatrix <- function(size=c(1,1),type=1,syndist=3,params=c(1,1.15,.16),dop=c(),
                          source='source',target='target',noise=list(0,c(0,0),0),learningRule=c(1/38,2,0.2,T),
                          gain=1,op=1){
  #kcmbon_conn <- newConnMatrix(mbno,length(getData(mbons)[[const.N$neurons]]),
  #                             dop=length(getData(dans)[[const.N$neurons]]),type = 5,params = c(16,2,1)) 
  source <- size[1] #source and target, source forms the cols
  target <- size[2]
  #  glommb_conn <- newConnMatrix(c(glomno,kcno),syndist = 10,params=glommb_conn,dop=c())
  if(type==1){#distributed circuit
    #cat('\nnewconn',size)
    connMat <- genConnMat(rows = size[2],cols = size[1],params = params,type = syndist,noise = noise)
  }
  if(type==2){#topographic circuit
    #calculate the number of targets that a source neuron contacts
    source_target <- size[2]/size[1]
    #generate matrix, with rows as the target and cols are the source
    #lets just do the source as rows first and then transpose in the end
    connMat.ls <- lapply(1:size[1], function(x){
      #figure out which neurons this should connect to
      #basically, put zeroes, ones in the topographic places, and then zeroes again 
      if (source_target > 1) target_neurons <- c(rep(0,(x-1)*source_target),rep(1,source_target),rep(0,(size[2]/source_target-x)*source_target))
      else {
        #cat(x,source_target,ceiling(x*source_target),x*source_target,'\n')
        target_neurons <- c(rep(0,size[2]))
        target_neurons[ceiling(x*source_target)] <- 1
      }
      target_neurons
    })
    connMat <- convertNestedListsDF(connMat.ls)
  }
  if(length(dop) > 0){#this circuit does learning, so set up those params too
    #needs to be normalized if we are going to use saturating functions, add the factor so that we dont hve a 1 as that screws up the saturation calculation
    connMat <- connMat/(max(connMat)*1.010101) #normalize the matrix so that the maximum synapse strength is set to 1
    dopNo <- dop[[1]] #no of dop neurons
    dopDist <- dop[[2]] #the synaptic distribution
    dopParams <- dop[[3]] #the params for the dop connection, same as the params of the connection matrix
    dopType <- dop[[4]] #the kind of dopamine, reward or punishment, 1 for reward, -1 for punish
    #generate another connection matrix, in this case the source is the dop neurons, and the targets are the synapses
    #of the kcs, so the KCs themselves
    dopConn <- genConnMat(rows = size[1],cols = dopNo,params = dopParams,type = dopDist,noise = noise) #should noise be independent?
    curTime <- 0 #this is the start so the time counter is 0
    synactfn <- function(x) x #the synapse activiation function
    connMatrix <- list(connMat,type,source,target,params,size,syndist,noise,learningRule,gain,dopNo,dopParams,dopDist,dopType,dopConn,curTime,synactfn)
    names(connMatrix) <- c('connMat','type','source','target','params','size','syndist','noise','learningRule','gain',
                           'dopNo','dopParams','dopDist','dopType','dopConn','curTime','synActFn')
    res <- structure(connMatrix, class = c('Learning','ConnMatrix')) #it is a learning connection matrix
  }
  else {#not a learning matrix
    connMatrix <- list(connMat,type,source,target,params,size,syndist,noise,learningRule,gain)
    names(connMatrix) <- c('connMat','type','source','target','params','size','syndist','noise','learningRule',
                           'gain')
    if(syndist == 7) res <- structure(connMatrix, class= c('ConnMatrix','MultSyn')) #multiple synapses bw neurons
    else res <- structure(connMatrix, class= c('ConnMatrix','SingleSyn')) #single synapse between neurons
  }
  res #now, return the connection matrix object
}

#setup generic methods to change the values inside the object


#function calculates or recalculates the synaptic weight normalizinng function 
#self: the circuit structure
#mbkc_neurons: the firing patterns of the mbkc for different stimuli
#valence: reward(approach) or punishment(avoid) 1 or -1
#aplcut: the apl percentage cutoff, requited in order to calculate, if 0, mbkc.neurons is actually mbkc_net so you dont need aplcut
#op=1,  
normSynWts <- function(self,mbkc.neurons,valence=1,aplCut=0,op=1){
  #gets the high kc values for each odor and take their mean.
  if (aplCut>0) high.mbkc <- mean( sapply(1:length(mbkc.neurons), function(i) max(mbkc.neurons[[i]] - getTopNPercentileThresh(mbkc.neurons[[i]],topn = aplCut,op=2))) )
  else high.mbkc <- mean( sapply(1:length(mbkc.neurons), function(i) max(mbkc.neurons[[i]])) )
  #get the kcmbon connection matrix and DANS for this particular valence
  kcmbon_conn.lst <- getData(self,const.FM$kcmbon_conn)
  kcmbon_conn <- kcmbon_conn.lst[[valence]] #the learning matrix for this valence
  dans <- getData(self,const.FM$dans)[[valence]] #the dopamine neurons 
  #now, do the adjustment function for the multiplicative weight update
  dop_kc.nw <- getData(kcmbon_conn,const.CM$dopconn)
  dop_kc.nw <- apply(dop_kc.nw,2,mean) #mean synaptic strengths from each dopamine neuron
  dop.neurons <- getData(dans,const.N$neurons) #the dopamine neurons firing
  total.dop <- sapply(1:length(dop.neurons[[1]]), function(i) dop.neurons[[1]][i] * dop_kc.nw[i]) 
  high.mbkc <- sum(total.dop) * high.mbkc #roughly the highest update value possible 
  beta <- getData(kcmbon_conn,const.CM$learningrule)[2] #get beta 
  #function to normalize synaptic weighting: the idea is that when x is 0, i.e., the kc is almost silent, the multiplicative factor is 1
  #and when x is close to highest possible value of high.mbkc, it is 1/beta. Gives the equation below. returns 1 for x=0, and beta for x=high
  #the thresh equation ensures that once the x is above high.bkc(since it is a mean), the function evaluates to beta: 
  syn <- function(x)  threshUpperVal(((beta-1)/(high.mbkc))*x + 1,val = beta) 
  #set the function
  kcmbon_conn <- setData(kcmbon_conn,val=syn,index=const.CM$synactfn) #first set the learning matrix
  #kcmbon_conn.lst <- getData(self,const.FM$kcmbon_conn) #now update this matrix in the list of matrices
  kcmbon_conn.lst[[valence]] <- kcmbon_conn
  circuit.mb <- setData(self,val=kcmbon_conn.lst,index=const.FM$kcmbon_conn) #update the list of matrices in the circuit structure
  #cat('\nnormsynwts',circuit.mb$kcmbon_conn.avo$synActFn(.2),circuit.mb$kcmbon_conn.app$synActFn(.2),valence)
  circuit.mb
}




#this computes the change in the synapses and neurons for the various objects with training
trainingUpdate <-function(self,...){
  UseMethod("trainingUpdate",self)
}


#this computes the change in the synapses and neurons for the various objects with training
#the input is one odor specified by odorno
#only learning. no calculation of mbon rates which should be done in the testing phase
#timestep: no of time steps, i.e., trials to iterate through
#odorno: the odor number on which training is done. 
#valence: the valence number that we would like to train.
trainingUpdate.FlyMBg <-function(self,odorno=1,valence=1,timestep=1){
  pt <- proc.time()
  temp <- computeFiring(self = self,odor=odorno) #calculate the firing rate for all the neurons
  #now, update the connection strengths of the kc_mbon_conn matrices based on mbkc_net firing
  #first get the learning matrix, mbon, and dans associated with a particular valence
  pt1 <- proc.time()
  mbon <- getData(self,const.FM$mbons)[[valence]]
  dans <- getData(self,const.FM$dans)[[valence]]
  kcmbon_conn.lst <- getData(self,const.FM$kcmbon_conn)
  kcmbon_conn <- kcmbon_conn.lst[[valence]]
  
  mbkc_net <- getData(temp,const.FM$mbkc_net) #get the mb_kc firing
  dan <- getData(self,const.FM$dans)[[valence]]
  pt2 <- proc.time()
  #update the learning matrix according to KC firing
  #cat('\ntrainingUPDATe valence',unlist(kcmbon_conn$connMat),';\t',unlist(mbkc_net$neurons[[odorno]]),';\t',unlist(getData(dan,const.N$neurons)[[1]]) )
  kcmbon_conn <- trainingUpdate(kcmbon_conn,neuron.input = mbkc_net$neurons[[odorno]],
                                dop.input=getData(dan,const.N$neurons)[[1]],currentTime=timestep)
  kcmbon_conn.lst[[valence]] <- kcmbon_conn #update the learning matrix
  #cat('\ntrainingUPDATe 2: valence',valence,' : odor',odorno,': ',unlist(kcmbon_conn$connMat) )
  temp <- setData(temp,val = kcmbon_conn.lst,const.FM$kcmbon_conn) #set the mb firing rate
  pt3 <- proc.time()
  temp
}



#self: this connection matrix that needs a training update
#neuron.input: the kcs or source neurons and their firing rates. The amount by which the synapses are strengthened is a function of dop input and source input
#dop.input: the firing rate of the input dopamine neurons. This determines the change in synaptic weights 
#traintime: the current time, so this will update the curtime and also use this to calcculate the degradation in synaptic weights
#sat: T - whether we use a saturating function for update, F - we use a linear proportional function
trainingUpdate.ConnMatrix <- function(self,neuron.input,dop.input,currentTime,op=1){
  #algo: if KC is active, its synapses undergo multipliicative learning, and inactive KCs undergo addditive changes
  #for the active synapses, factor the multiplicative factor by KC input and dop input,
  #get the conn matrix params: synapse strengths, and learning rules and dop input
  pt <- proc.time()
  
  cmat.flat <- getData(self,const.CM$connmat,op = 1)
  learningRule <- getData(self,const.CM$learningrule) #get the learning rules
  synActFn <- getData(self,const.CM$synactfn) #get the synapse acitvation fn for factoring in the dop update
  dopConn <- getData(self,const.CM$dopconn) #the dopamine to neurons connections, rows - neurons, cols - dopamine
  cmat.dim <- dim(cmat.flat) #the dimension of the matrix for you to put it back
  dim(cmat.flat) <- NULL #make it into a vector, column wise
  pt1 <- proc.time()
  dop.strength <- dop.input %*% t(dopConn) #calculates the dopamine input strength for every KC
  #calculate the kC and dopamine input factors that are going to increase the weight
  inpSize <- length(neuron.input) # no of input neurons
  update.vals <- sapply(1:inpSize, function(i){
    #for each kc neuron, go through each dopamine input and calculate the strenth increase based on
    # kc.firing*dop.firing*dop.kc.synapse.strength
    kc <- neuron.input[i] * dop.strength[i] 
    kc
  }) #returning a vector of dop. strength inputs for each KC
  pt2 <- proc.time() 
  active_pos <- which(update.vals != 0) #active synapses, put in 
  active_pos <- unlist(lapply(1:(length(cmat.flat)/inpSize),function(i) active_pos + ((i-1)*inpSize)  ) )
  inactive_pos <- setdiff(1:length(cmat.flat),active_pos) #inactive synapses
  #inactive_pos <- unlist(lapply(1:(length(cmat.flat)/inpSize),function(i) inactive_pos + (i-1)*(inpSize) ) )
  pt3 <- proc.time() 
  #now, go through the connection matrix and update the synapse on this basis
  posns <- rep(1:inpSize,length(cmat.flat)/inpSize) #all the positions 
  #cat('\n posns',posns,'\t',length(cmat.flat),'\t',inpSize,'\naactuvee',active_pos,'\n',inactive_pos)
  if(learningRule[4]==T){#synapses change according to a saturating function
    #do active and inactive synapses separately
    cmat.active <- sapply(active_pos, function(i) {#active synapses
      #generate the appropriate positions, so, for an input size of 5, 3 remains 3 and 6 becomes 1
      #val <- updateSynapseSat(synapse = cmat.flat[i],update = update.vals[pos],learningRule = learningRule,activeSyn = T,synActFn = synActFn)
      x0 <- invSatFn(cmat.flat[i],ec50 = learningRule[3],n = 1)
      x1 <- x0 / (synActFn(update.vals[posns[i]]) ) #multiplicative learning
      val <- satFn(x1,ec50 = learningRule[3],n = 1)
    })
    cmat.inactive <- sapply(inactive_pos, function(i) {#inactive synapses
      #generate the appropriate positions, so, for an input size of 5, 3 remains 3 and 6 becomes 1
      #val <- updateSynapseSat(synapse = cmat.flat[i],learningRule = learningRule,activeSyn = F,synActFn = synActFn)
      x0 <- invSatFn(cmat.flat[i],ec50 = learningRule[3],n = 1)
      x1 <- x0 + learningRule[1] #multiplicative learning
      val <- satFn(x1,ec50 = learningRule[3],n = 1)
      #if (update.vals[pos]==0) cat(',learning1',learningRule[1],'update',cmat.flat[i]+learningRule[1],'cmat',cmat.flat[i])
      #else cat('\n',i,' learning2',learningRule[2],'update val',update.vals[pos],'update',cmat.flat[i] * learningRule[2] * update.vals[pos],'cmat',cmat.flat[i])
      #cat(' val',val)
    })
    pt4 <- proc.time()
    cmat <- cmat.flat #now populate the cmat vector with noth sets of results
    cmat[active_pos] <- unlist(cmat.active) #unlist needed in case it is 0-length list
    cmat[inactive_pos] <- unlist(cmat.inactive)
    #testing if any of these are not a number which later causes a problem in computeFiring: liam
  }
  else {#synapses change according to a linear function
    cmat <- sapply(1:length(cmat.flat), function(i) {
      #generate the appropriate positions, so, for an input size of 5, 3 remains 3 and 6 becomes 1
      if(update.vals[posns[i]] == 0) val <- updateSynapseLin(synapse = cmat.flat[i],learningRule = learningRule,activeSyn = F,synActFn = synActFn)
      else val <- updateSynapseLin(synapse = cmat.flat[i],update = update.vals[posns[i]],learningRule = learningRule,activeSyn = T,synActFn = synActFn)
      val
    })
    cmat <- threshUpperVal(threshValTarget(cmat,val=0,target = 0.000001),val=1) #control the lower and uppper ends of the weights
  }
  dim(cmat) <- cmat.dim #assemble it into a matrix again
  conn <- self #cant modify self as it is a funnction parameter, so use conn to update
  #pt4 <- proc.time()
  curTime <- getData(conn,const.CM$curtime)
  if(length(curTime) <= 1) curTime <- list(c(0,0,0),c(0,0,0),c(0,0,0),c(0,0,0))
  #cat('\n time',(pt4-pt)[1:3],(pt4-pt3)[1:3],(pt4-pt2)[1:3],'change synapse',(pt3-pt2)[1:3],'update',(pt2-pt1)[1:3])
  curTime <- list(curTime[[4]] + (pt4-pt3)[1:3],curTime[[3]] + (pt3-pt2)[1:3],curTime[[2]] + (pt2-pt1)[1:3],curTime[[1]] + (pt1-pt)[1:3])
  
  conn <- setData(conn,val=cmat,index=const.CM$connmat) #now, set this data and return 
  conn <- setData(conn,val=curTime,index=const.CM$curtime) #set the current time
  conn #the new updated connection matrix
}


#this function updates synapses that are active, i.e., the KC is activated, the nonlinear Hill Fn update
#learning rule specifies the multiplicative factor
#synapse: the current value of the synapse
#update: factor by which it should be changed, accounts for KC and dopamine firing, 0 if it is an inactive synapse
#learning rule: the multiplicative learning rule, n = fill coefficient
#activeSyn: this is an active synapse
#sat: T - whether we use a saturating function for update, F - we use a linear proportional function
updateActiveSynapseSat <- function(synapse,update=0,learningRule,n=1,activeSyn=T,synActFn=c(),op=1){
  #the synapse value varies according to a saturation function, so we move it on this curve according to beta
  #find out the correponding x posn of synapse, then add the learnin induced multiplicative increase, and calculate the new y posn
  x0 <- invSatFn(synapse,ec50 = learningRule[3],n = n)
  x1 <- x0 / (synActFn(update) ) #multiplicative learning
  y1 <- satFn(x1,ec50 = learningRule[3],n = n)
  #cat(' synapse',synapse,' x0',x0,' x1',x1,' y1',y1)
  y1
}

updateInactiveSynapseSat <- function(synapse,update=0,learningRule,n=1,activeSyn=T,synActFn=c(),op=1){
  #the synapse value varies according to a saturation function, so we move it on this curve according to beta
  #find out the correponding x posn of synapse, then add the learnin induced multiplicative increase, and calculate the new y posn
  x0 <- invSatFn(synapse,ec50 = learningRule[3],n = n)
  x1 <- x0 + learningRule[1] #learning or additive
  y1 <- satFn(x1,ec50 = learningRule[3],n = n)
  #cat(' synapse',synapse,' x0',x0,' x1',x1,' y1',y1)
  y1
}


#this function updates synapses that are active, i.e., the KC is activated, the linear update
#learning rule specifies the multiplicative factor
#synapse: the current value of the synapse
#update: factor by which it should be changed, accounts for KC and dopamine firing, 0 if it is an inactive synapse
#learning rule: the multiplicative learning rule, n = fill coefficient
#activeSyn: this is an active synapse
#sat: T - whether we use a saturating function for update, F - we use a linear proportional function
updateSynapseLin <- function(synapse,update=0,learningRule,n=1,activeSyn=T,synActFn=c(),op=1){
  # no saturation so, just follows the function, y = x
  x0 <- synapse
  x1 <- ifelse(activeSyn,x0 / (synActFn(update) ),x0 + learningRule[1]) #multiplicative learning or additive
  x1
}

#this function sets up a new class object for the glomerular-KC-MBON-DAN circuit
#glom: no of glomeruli
#kcno: no of kcs
#danno: no of dans per type, we have two types of DANs
#mbonno: no of mbons per compartment, we have two types of MBONs
#nooodors: total number of odors.
#novalences: either the number of reward punsihment type behaviors, i.e., a single number or the valences like c(1,2,3)
#tatparams: the parameters that describe how training takes place and testing takes place, for now traintype:1-alternate, 2 - sequence, 3 - random
#testtype: 1 - is valence learnt, 2 - pairwise discrimination. Overall c(traintype,testtype)
#glommb_conn:the glomeruli to KC connection based on chucks study. 3 param sets: 
#1=chooses glomeruli for eacch KC,2= hypergeometric distribution, 3=synaptic distribution, gamma 
#a vector of valences = c('app','avo') or c('sugar','wateer','shock','heat')
#params: use this to pass other params, for example, for glomeruli, params=list(1,dist,dist.params)
#learningrule: the rule for updating the strength of the synapses.(alpha,beta,ec50,T:sat fn,F:linear function) additive decrease s-alpha, multiplicative inc: s*beta 
#thresh: set the threshold for neurons list(neuronop,threshparams). neuronop: 1 - mbkc, 2 - mbkc, 3 - mbkc and mbkcnet 
newGlomMBcircuit <- function(glomno=50,kcno=2000,danno=5,mbonno=1,noodors=4,params=c(),learningRule=c(1/38,2,0.2,T),
                             glommb_conn=list(c(8,0.715),c(0.07),c(2,4,4)),kcmbon_conn= list(c(13.5,1.15,0.16),list(danno,5,c(2,0.9,1.1),1)),
                             novalences=2,tatparams=c(1,1),thresh=c(1,0),op=1){
  #neurons with multiple entries: glom, mb, mbkc_net, and mbons, apl
  #single entry: dans
  #all connection structures are also single.
  #cat('\nglommbcircuit')
  glom <- newNeurons(number = glomno,dist = 2,noStimuli = noodors,id = const.FM$glom) #projection neuron firing

  glommb_conn <- newConnMatrix(c(glomno,kcno),syndist = 10,params=glommb_conn,dop=c())
  mbkc <- newNeurons(number = kcno,dist = 0,params = c(0),noStimuli = noodors,id = const.FM$mbkc,thresh = thresh)
  mbkc_net <- newNeurons(number = kcno,dist = 0,params = c(0),noStimuli = noodors,id=const.FM$mbkc_net)
  kcapl_conn <- newConnMatrix(c(kcno,1),syndist = 6,dop=c(),params = 1)
  apl <- newNeurons(number = 1,dist = 0,params = c(0),noStimuli = noodors,neurondesc = 'inhibitory',id=const.FM$apl)
  aplkc_conn <- newConnMatrix(c(1,kcno),syndist = 6,dop=c(),params = 1)
  
  #setup the valence n/w
  if(length(novalences)>1) valences <- novalences
  else valences <- 1:novalences
  # # dans = # valencs
  dans <- lapply(1:novalences, function(i) newNeurons(number = danno,dist = 0,params = c(1),id = const.FM$dans) )
  names(dans) <- valences
  # mbons = # valences, each mbons also contains the mbon firing rates given a particular odor.
  mbons <- lapply(1:novalences, function(i) newNeurons(number = mbonno,dist = 0,params = c(0),noStimuli = noodors,id = const.FM$mbons) )
  names(mbons) <- valences
  
  #now you have to put in the connection matrix and learning rules for the MB-MBON connections
  #13.5 synapses from every KC to a DAN; one average 2 ssynapses from every DAN to KCs
  kcmbon_conn <- lapply(1:novalences, function(i) newConnMatrix(c(kcno,mbonno),params = c(13.5,1.15,0.16),
                                                                dop = list(danno,5,c(2,0.9,1.1),1) ) )
  names(kcmbon_conn) <- valences
  
  #initialize the training vals, default is no of valences and alternating, and origstim and testing are the same too
  origstim <- cbind.data.frame(genSequence(length(valences),size = noodors,op = tatparams[1]),1:noodors,1:noodors) 
  colnames(origstim) <- c('val','odor#','mix')
  #add glomeruli for training val and teststim
  allglom <- addData(glom,stimul=glom)
  trainingval <- origstim
  trainingval[,2] <- 1:nrow(trainingval) + origstim[nrow(origstim),2]


  allglom <- addData(allglom,stimul=glom)
  teststim <- origstim
  teststim[,2] <- 1:nrow(teststim) + trainingval[nrow(trainingval),2]
  
  glomMbCircuit <- list(allglom,mbkc,mbkc_net,apl,dans,mbons,
                        glommb_conn,kcapl_conn,aplkc_conn,kcmbon_conn,noodors,valences,trainingval,teststim,tatparams,origstim)
  names(glomMbCircuit) <- c('glom','mbkc','mbkc_net','apl','dans','mbons',
                            'glommb_conn','kcapl_conn','aplkc_conn','kcmbon_conn','noodors','valences','trainingval',
                            'teststim','tatparams','origstim')
  res <- structure(glomMbCircuit,class = c('FlyMBg','Circuit','GLOMMB')) #make it a class and return it
  res
}

#modification on the previous circuit where theere are a whole bunch of connections between 
#ka-apl and apl-kc instead of just one. according to what Mehrab and Glenn wanted
#this function sets up a new class object for the glomerular-KC-MBON-DAN circuit
#glom: no of glomeruli
#kcno: no of kcs
#danno: no of dans per type, we have two types of DANs
#mbonno: no of mbons per compartment, we have two types of MBONs
#nooodors: total number of odors.
#novalences: either the number of reward punsihment type behaviors, i.e., a single number or the valences like c(1,2,3)
#tatparams: the parameters that describe how training takes place and testing takes place, for now traintype:1-alternate, 2 - sequence, 3 - random
#testtype: 1 - is valence learnt, 2 - pairwise discrimination. Overall c(traintype,testtype)
#glommb_conn:the glomeruli to KC connection based on chucks study. 3 param sets: 
#1=chooses glomeruli for eacch KC,2= hypergeometric distribution, 3=synaptic distribution, gamma 
#a vector of valences = c('app','avo') or c('sugar','wateer','shock','heat')
#params: use this to pass other params, for example, for glomeruli, params=list(1,dist,dist.params)
#kcapl_conn and apl_kc_conn: list(synapse distribution type,params of synapse distribution: c(no of synapse,distribution type
#, distrinbution params)), only 1 synapse -> c(6,c(1)) or list(7,c(6,1,1.15,0.12))
#third option you can also just specify the connection synapses: list(dist.type:8, 
#the synapse dist: paramsc(1,1.15,0.12),connection vector e.g., c(6,7,8,4,3)).  
#use this option when you want to generate correlated connection patterns
#learningrule: the rule for updating the strength of the synapses.(alpha,beta,ec50,T:sat fn,F:linear function) additive decrease s-alpha, multiplicative inc: s*beta 
#thresh: set the threshold for neurons list(neuronop,threshparams). neuronop: 0 - no threshold, 1 - mbkc, 2 - mbkc_net, 3 - mbkc and mbkcnet 
newPNKCAPLcircuit <- function(glomno=50,kcno=2000,danno=5,mbonno=1,noodors=4,params=c(),learningRule=c(1/38,2,0.2,T),
                              kcapl_conn_par=list(7,c(6,1,1.15,0.12)),aplkc_conn_par=list(8,c(6,1,1.15,0.12)),
                              glommb_conn=list(c(8,0.715),c(0.07),c(2,4,4)),kcmbon_conn= list(c(13.5,1.15,0.16),list(danno,5,c(2,0.9,1.1),1)),
                              novalences=2,tatparams=c(1,1),thresh=list(0,c(1,0)),op=1){
  #neurons with multiple entries: glom, mb, mbkc_net, and mbons, apl
  #single entry: dans
  #all connection structures are also single.
  glom <- newNeurons(number = glomno,dist = 2,noStimuli = noodors,id=const.FM$glom) #projection neuron firing
  glommb_conn <- newConnMatrix(c(glomno,kcno),syndist = 10,params=glommb_conn,dop=c())
  
  mbkc <- newNeurons(number = kcno,dist = 0,params = c(0),noStimuli = noodors,neurondesc = 'excitatory',id=const.FM$mbkc)
  
  apl <- newNeurons(number = 1,dist = 0,params = c(0),noStimuli = noodors,neurondesc = 'inhibitory',id=const.FM$apl)
  kcapl_conn <- newConnMatrix(c(kcno,1),syndist = kcapl_conn_par[[1]],dop=c(),params = kcapl_conn_par[[2]])
  mbkc_net <- newNeurons(number = kcno,dist = 0,params = c(0),noStimuli = noodors,neurondesc = 'exc. wtaed',id=const.FM$mbkc_net)
  #specify the connection params for the feedback
  syn.cor.vec <- getConnSynapseNo(kcapl_conn$connMat) #gets the vector of the number of synapses from each KC to aPL
  
  conn.cor <- convertMatInteger(genCorrVecsSpec(syn.cor.vec[[1]],n=2))[,2] #generate the correlated vector
  if(min(conn.cor)<0) conn.cor <- conn.cor + abs(min(conn.cor)) #remove all -ve values
  conn.cor <- floor((conn.cor/mean(conn.cor)) * aplkc_conn_par[[2]][1]) #normalize to have the mean specified by aplkc_conn params
  aplkc_conn_params <- c(aplkc_conn_par,list(conn.cor))
  #cat('\napl-kc params:',str(aplkc_conn_params),conn.cor )
  aplkc_conn <- newConnMatrix(c(1,kcno),syndist = aplkc_conn_params[[1]],dop=c(),params = aplkc_conn_params[-1])

  #setup the valence n/w
  if(length(novalences)>1) valences <- novalences
  else valences <- 1:novalences
  # # dans = # valencs
  dans <- lapply(1:novalences, function(i) newNeurons(number = danno,dist = 0,params = c(1),id=const.FM$dans) )
  names(dans) <- valences
  # mbons = # valences, each mbons also contains the mbon firing rates given a particular odor.
  mbons <- lapply(1:novalences, function(i) newNeurons(number = mbonno,dist = 0,params = c(0),noStimuli = noodors,id=const.FM$mbons) )
  names(mbons) <- valences
  
  #now you have to put in the connection matrix and learning rules for the MB-MBON connections
  #13.5 synapses from every KC to a DAN; one average 2 ssynapses from every DAN to KCs
  kcmbon_conn <- lapply(1:novalences, function(i) newConnMatrix(c(kcno,mbonno),params = c(13.5,1.15,0.16),
                                                                dop = list(danno,5,c(2,0.9,1.1),1) ) )
  names(kcmbon_conn) <- valences
  
  #initialize the training vals, default is no of valences and alternating, and origstim and testing are the same too
  origstim <- cbind.data.frame(genSequence(length(valences),size = noodors,op = tatparams[1]),1:noodors,1:noodors) 
  colnames(origstim) <- c('val','odor#','mix')
  #add glomeruli for training val and teststim
  allglom <- addData(glom,stimul=glom)
  trainingval <- origstim
  trainingval[,2] <- 1:nrow(trainingval) + origstim[nrow(origstim),2]
  
  allglom <- addData(allglom,stimul=glom)
  teststim <- origstim
  teststim[,2] <- 1:nrow(teststim) + trainingval[nrow(trainingval),2]
  
  PNKCAPLcircuit <- list(allglom,mbkc,mbkc_net,apl,dans,mbons,
                        glommb_conn,kcapl_conn,aplkc_conn,kcmbon_conn,noodors,valences,trainingval,teststim,tatparams,origstim)
  names(PNKCAPLcircuit) <- c('glom','mbkc','mbkc_net','apl','dans','mbons',
                            'glommb_conn','kcapl_conn','aplkc_conn','kcmbon_conn','noodors','valences','trainingval',
                            'teststim','tatparams','origstim')
  res <- structure(PNKCAPLcircuit,class = c('FlyMBg','Circuit','PNKCAPL')) #make it a class and return it

  #should we initialize it here, calcaulte APL and get the thing going, Maybe?
  res
}


#this function sets up a new class object for the OB-PCx circuit
#the latest function is in the file ob_pcx_circuit.R
#glom: no of glomeruli
#pcx: no of pcx neurons
#nooodors: total number of odors.
#glommb_conn:the glomeruli to KC connection based on chucks study. 3 param sets: 
#1=chooses glomeruli for eacch KC,2= hypergeometric distribution, 3=synaptic distribution, gamma 
#a vector of valences = c('app','avo') or c('sugar','wateer','shock','heat')
#params: use this to pass other params, for example
#pcneurons: ratios of the different PC neurons in the following order: pc1i,pc2i,pc2e,pc3e.
#gcneurons: number of gcs
#pcpc, pcob, obpc, and pcpc-prams: pc2i2e,pc1i2e,pc2e3e,pc2e2i,pc2i3e : parameters for the connection networks
#learningrule: the rule for updating the strength of the synapses.(alpha,beta,ec50,T:sat fn,F:linear function) additive decrease s-alpha, multiplicative inc: s*beta 
newOBPcxcircuit.old <- function(glomno=200,pcno=2800,noodors=4,obpcparams=c(),pcpcparams=c(),pcobparams=c(),obobparams=c(),op=1){
  #neurons in thi circuit: glom, pc1i,pc2i,pc2e,pc3e, and gcs
  #next interation should also include obei (ob inhibitory neurons in EPL), obpi (periglomerular inhibitory neurons)
  #networks: glompc1i,glompc2e,pc1ipc2e,pc2epc3e,pc2epc2i,pc2ipc3e,glompc3e
  glom <- newNeurons(number = glomno,dist = 2,noStimuli = noodors) #projection neuron firing
  cat('\nError here')
  #5/56, 25/56 etc gives the number of each type of cell
  pc1i <- newNeurons(number = pcno*(5/56),dist = 0,noStimuli = noodors,params = c(0),neurontype = 'inh')
  pc2e <- newNeurons(number = pcno*(25/56),dist = 0,noStimuli = noodors,params = c(0),neurontype = 'ext')
  cat('\n1Error here')
  pc2i <- newNeurons(number = pcno*(8/56),dist = 0,noStimuli = noodors,params = c(0),neurontype = 'inh')
  pc3e <- newNeurons(number = pcno*(20/56),dist = 0,noStimuli = noodors,params = c(0))
  #number of granule cells, which is 100x of glom, gc1 are neonatal cells, gc2 are adult newborn gcs
  gc1frac <- 0.8 #fraction of neonatal gcs
  cat('\n2Error here')
  gc1 <- newNeurons(number = glomno*100*(gc1frac),dist = 2,noStimuli = noodors,params = c(0)) #gcs firing
  gc2 <- newNeurons(number = glomno*100*(1-gc1frac),dist = 2,noStimuli = noodors,params = c(0)) #gcs firing
  
  cat('\n3Error here')
  
  #now,do the networks: from ob to pc
  glompc1i_conn <- newConnMatrix(size = c(glom,pc1i),type = 2,params = obpcparams,syndist = 3) #type 3 n/w
  cat('\n3.1Error here')
  glompc2e_conn <- newConnMatrix(size = c(glom,pc2e),type = 2,params = obpcparams,syndist = 3) #type 3 n/w
  glompc3e_conn <- newConnMatrix(size = c(glom,pc3e),type = 2,params = obpcparams,syndist = 3) #type 3 n/w
  
  cat('\n4Error here')
  #within pc
  # pc1i2e_conn <- newConnMatrix(size = c(pc1i,pc2e),type = 2,params = pc1i2e,syndist = 20) #type=20,spatial conn matrix
  pc2e2i_conn <- newConnMatrix(size = c(pc1i,pc2e),type = 2,params = pc2e2i,syndist = 20) #type=20,spatial conn matrix
  pc2i2e_conn <- newConnMatrix(size = c(pc1i,pc2e),type = 2,params = pc2i2e,syndist = 20) #type=20,spatial conn matrix
  # pc2e3e_conn <- newConnMatrix(size = c(pc1i,pc2e),type = 2,params = pc2e3e,syndist = 20) #type=20,spatial conn matrix
  # pc2i3e_conn <- newConnMatrix(size = c(pc1i,pc2e),type = 2,params = pc2i3e,syndist = 20) #type=20,spatial conn matrix
  
  cat('\n5Error here')
  
  #from pc to gcs in ob
  # pc2eglom_gc1 <- newConnMatrix(size = c(pc2e,gc1),type = 2,params = pc2egc1,syndist = 20) #type=20,spatial conn matrix
  # pc2eglom_gc2 <- newConnMatrix(size = c(pc2e,gc2),type = 2,params = pc2egc2,syndist = 20) #type=20,spatial conn matrix
  # pc3eglom_gc1 <- newConnMatrix(size = c(pc3e,gc1),type = 2,params = pc3egc1,syndist = 20) #type=20,spatial conn matrix
  # pc3eglom_gc2 <- newConnMatrix(size = c(pc3e,gc2),type = 2,params = pc3egc2,syndist = 20) #type=20,spatial conn matrix
  cat('\nError here')
  
  
  
}


#initialization of the object after the construction is done
init <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("init",self)
  }
}

#this function is needed to set the parameters or the connection matrix between apl and kcs so that, 
#and setup the KC-MBON networks, basically normalize the update function
#you get the required percentage of KCs firing
#aplCut: apl value so that only 10 % of KCs fire
init.FlyMBg <- function(self,aplCut=5,op=1){
  #the difference between init and compute firing is two fold: this sets up APL gain, that doesn't
  #And, there is no noise in these calculations.
  mbkc <- computeFiring(getData(self,const.FM$glom),getData(self,const.FM$glommb_conn),
                        dest=getData(self,const.FM$mbkc)) #get the KC firing rates
  circuit.mb <- setData(self,index=const.FM$mbkc,val=mbkc)
  mbkc.neurons <- getData(mbkc,const.N$neurons) #just the neurons to calculate the apl cutoff
  #this calculates the APL gain needed to achive aplCUt and updates, apl, aplkc_conn, and mbkc_net
  circuit.mb <- calcAPLConn(circuit.mb,aplCut = aplCut,op=1)
  #now, that you have the MB KC firing going, do the MBONs
  kcmbon_conn.lst <- getData(circuit.mb,const.FM$kcmbon_conn)
  mbon.lst <- getData(circuit.mb,const.FM$mbons)
  for(i in 1:length(kcmbon_conn.lst) ){
    mbon.lst[[i]] <- computeFiring(getData(circuit.mb,const.FM$mbkc_net)
                                   ,kcmbon_conn.lst[[i]],mbon.lst[[i]],op=1)
  }
  circuit.mb <- setData(circuit.mb,val = mbon.lst,index=const.FM$mbons) #set the mbon firing rates
  #valences next, update the synaptric weights
  valences <- getData(self,const.FM$valences)
  for(i in 1:length(valences)){
    circuit.mb <- normSynWts(circuit.mb,getData(getData(circuit.mb,const.FM$mbkc_net),const.N$neurons)
                             ,valence = valences[i],aplCut = aplCut)
  }
  #also have to run or do learning for a couple of training odors to get the network primed for doin stuff
  #either a subset of existing odors, or just random odors picked from the glomerular distribution
  circuit.mb
}

#todo: generate sample odors for doing apl stuff by just shuffling the existing mbkcs, do apl, then
#reset to the earlier stuff


#function that calculates, based on the apl cutoff, the connection strength from KCs to APL and APL
#back to KCs. Unilike the earlier function, here, we will modify the existing connections, and not make 
#new connection matrices.
#aplcut: the number of kcs that should be firing
#self: the Flymb structure that needs to be modified
#op=1, the results as a new fly circuit, 2 : just the results
calcAPLConn <- function(self,aplCut=5,op=1){
  #procedure: Compute the firing rates of APL. 
  #Then, for every stimulus, compute the feedback constant that gives you aplCut number of neurons
  #the averag of all of them is the APL feedback gain.
  apl <- getData(self,const.FM$apl)
  aplkc_conn <- getData(self,const.FM$aplkc_conn)
  mbkc <- getData(self,const.FM$mbkc)
  kcapl_conn <- getData(self,const.FM$kcapl_conn)
  #now calculate the number that apl should be scaled by to be in the top aplcut percentile
  res <- updateFbConnGain(mbkc = mbkc,kcapl_conn = kcapl_conn,aplkc_conn = aplkc_conn,apl=apl,
                          aplCut = aplCut,op=5)
  #update everything
  res.mb <- setData(self,index=const.FM$apl,val=res[[3]])
  res.mb <- setData(res.mb,index=const.FM$aplkc_conn,val=res[[2]])
  mbkc_net <- setData(getData(self,const.FM$mbkc_net),index=const.N$neurons,val=getData(res[[4]],const.N$neurons) )
  res.mb <- setData(res.mb,index=const.FM$mbkc_net,val=mbkc_net)
  switch(op,res.mb,res)
}


#given the list of mb vectors, and the kc->apl connection matrix, will update the apl->Kc conn
#matrix with the gain that will, on average, give aplCut # neurons firing
#mbvec.lst: the list of kc responses for a bunch of odors
#kcapl_conn: the connection matrix going from kcs to apl
#aplkc_conn: the connection matrix going  back from apl to kcs
#apl: is the apl or intermediate target structure
#aplCut: the aplcut # neurons you want
#op=1: the average gain, 2: updated connection matrix, 3 - apl 4. updated kc responses, 
#5: updated gain, updated apl, kc responses, and conn.nmatric
updateFbConnGain <- function(mbkc,kcapl_conn,aplkc_conn,apl=c(),aplCut=5,shufflen=400,op=1){
  #now calculate the number that apl should be scaled by to be in the top 1-aplcut percentile
  apl.new <- computeFiring(mbkc,kcapl_conn,dest=apl)
  #generate shuffled mbkc odors so that we can narrow down on APL gain.
  mbkc.shuffle.neurons <- genShuffledResponses(getData(mbkc,const.N$neurons),n=shufflen,op=1)
  #mbkc.shuffle.neurons <- getData(mbkc,const.N$neurons)
  mbkc.shuffle <- setData(mbkc,index=const.N$neurons,val=mbkc.shuffle.neurons)
  apl.shuffle <- computeFiring(mbkc.shuffle,kcapl_conn,dest=apl)
  apl.shuffle.neurons <- getData(apl.shuffle,const.N$neurons)
  syn.wt <- sapply(1:length(aplkc_conn$connMat), function(j) sum(aplkc_conn$connMat[[j]]) )
  #algorithm: claculate the APL gain value range for each odor stimulus that gives you
  #aplCut # mbkc neurons firing. Then, take the average of all the intervals. That's your gain.
  apl.gains <- lapply(1:length(mbkc.shuffle.neurons), function(i){
    calcAPLGain(mbkc.shuffle.neurons[[i]],apl = apl.shuffle.neurons[[i]],syn.vec = syn.wt,
                aplCut = aplCut)
  })
  apl.gain.ints <- lapply(apl.gains,'[[',2) #apl.gains is list of c(#mbkc>0, gain interval)
  #get mean interval vals that are in the most intervals, and take the mean if more than one
  #but, first, check extreme condition if apl.gains is null
  if((length(apl.gain.ints) == 0) || is.null(unlist(apl.gain.ints)) ) {#get the mean of the closest candidates
    allgains <- combineListMat(lapply(apl.gains,'[[',1)) #gets all the apl.gains
    candidates <- getClosestCandidates(allgains[2,],
                                       target = ceiling(length(mbkc.shuffle.neurons[[1]])*(aplCut)/100),op=2)
    apl.gain <- mean(allgains[1,candidates])
  }
  else apl.gain <- mean(getIntervalsCommon(apl.gain.ints)[,1]) 
  
  #update the connection matrix from apl back to Kcs with the apl gain
  aplkc_conn.new <- setData(aplkc_conn,index=const.CM$connmat,
                            val=lapply(aplkc_conn$connMat, function(x) x * apl.gain))
  
  #calculate the feedback and then subtract from existing mbkc
  mbkc.fb <- computeFiring(apl.new,aplkc_conn.new) #no dest. here, because this is just feedback
  mbkc.new <- twoSetOperations(mbkc,mbkc.fb,op=1) #subtrracts the second from the first
  switch(op,apl.gain,aplkc_conn.new,apl.new,mbkc.new,list(apl.gain,aplkc_conn.new,apl.new,mbkc.new) )
}

#function that calculates the gain given aplcut, neuron firing rates, apl rate, and synatpic weights
#algorithm: go thetough the range of apl gain calcaulted for each kc, and see which ones gives you 
#the desired target number of cells. Then, 
calcAPLGain <- function(mb.vec,apl,syn.vec,aplCut,step.size=0.05,op=1){
  #have to search
  gain.range <- cleanNAVec((mb.vec/syn.vec)/apl) #gets rid of INf
  zero.synapses <- countZeroes(syn.vec)
  #keep iterating in either direction until you get close to answer and it doesn't change.
  #determine the cutoff number of neurons you want
  target <- ceiling(length(mb.vec) * (aplCut)/100)
  delta <- (max(gain.range) - min(gain.range)) * step.size
  #variables to go through the loop
  res.no <- c(); res.gain <- c(); gain.guess <- min(gain.range)
  active.mb <- length(mb.vec) #keeps track of all the mb > 0
  #gain=0, all cells active, iterates from the min. gain to the gain that shuts off all cells. target 
  #has to be in the middle
  while ((active.mb > 0) && (active.mb>zero.synapses) && (apl > 0)) { #if apl < 0, infinite loop
    mb.val <- mb.vec - (apl * syn.vec * gain.guess)
    res.gain <- c(res.gain,gain.guess)
    active.mb <- length(which(mb.val>0))
    res.no <- c(res.no,active.mb )
    gain.guess <- gain.guess + delta
  }
  #the gain and the number active at this gain value
  #find the aplCut number and the interval for all values that produce this no
  #select all the vals, then interval is c(min,max)
  canditates <- which(res.no==target)
  if(length(canditates)>0) tar.inter <- c(min(res.gain[canditates]),max(res.gain[canditates])) 
  else {#if nothing matches, just don't return anything
    tar.inter <- c()
  }
  #cat('\ntar:',target,'tar. interval',tar.inter,'\n',res.no,'\n',res.gain)
  list(rbind(res.gain,res.no),tar.inter)
}


#function that calculates the gain given aplcut, neuron firing rates, apl rate, and synatpic weights
#algorithm: go thetough the range of apl gain calcaulted for each kc, and see which ones gives you 
#the desired target number of cells. Then, 
calcAPLGain.old <- function(mb.vec,apl,syn.vec,aplCut,step.size=0.05,op=1){
  #have to search
  thresh <- getTopNPercentileThresh(mb.vec,topn = 100-aplCut+1)
  gain.range <- cleanNAVec((mb.vec/syn.vec)/apl) #gets rid of INf
  #if the numbers are too close to each other, increase the range, basically within one fold
  if(isVecInRange(gain.range,range = 1)) gain.range <- c(gain.range,max(gain.range)*2)
  #keep iterating in either direction until you get close to answer and it doesn't change.
  #determine the cutoff number of neurons you want
  cat('\nmin max range:',min(gain.range),max(gain.range))
  target <- ceiling(length(mb.vec) * (aplCut)/100)
  delta <- (max(gain.range) - min(gain.range)) * step.size
  res.change <- 1
  gain.guess <- max(gain.range)
  res.no <- c()
  res.gain <- c()
  #cat('\n',gain.guess,'\t',gain.range,cleanNAVec(gain.guess))
  #iteratesfrom maximum to minimum and stores the gain in res.gain and the number of +ve kcs in res.no
  while (gain.guess > min(gain.range)) {
    mb.val <- mb.vec - (apl * syn.vec * gain.guess)
    res.gain <- c(res.gain,gain.guess)
    res.no <- c(res.no,length(which(mb.val>0)) )
    gain.guess <- gain.guess - delta
  }
  #the gain and the number active at this gain value
  #find the aplCut number and the interval for all values that produce this no
  #select all the vals, then interval is c(min,max)
  canditates <- which(res.no==target)
  if(length(canditates)>0) tar.inter <- c(min(res.gain[canditates]),max(res.gain[canditates])) 
  else {#if nothing matches, just don't return anything
    tar.inter <- c()
  }
  #cat('\ntar:',target,'tar. interval',tar.inter,'\n',res.no,'\n',res.gain)
  list(rbind(res.gain,res.no),tar.inter)
}


#function that gets the sparsity of the stimuli
#stim:1, the set to anaylysze: 1 - all odors, 2 - origm. 3 - training, 4 - test
#op:1 type of match
cirGetSparsity <- function(circuit,stim=1,op=1){
  stimuli <- switch(stim,1:getData(circuit,const.FM$noodors),getData(circuit,const.FM$origstim)[,2],
                    getData(circuit,const.FM$trainingval)[,2],getData(circuit,const.FM$teststim[,2]) )
  #before <- getData(getData(circuit,const.FM$mbkc),const.N$neurons)
  after <- getData(getData(circuit,const.FM$mbkc_net),const.N$neurons)
  #cat('\nstim',stimuli)
  vecsize <- length(after[[1]])
  res <- sapply(stimuli, function(i){
    countNonZeroes(after[[i]])/vecsize
  })
  mean(res) #return the average sparsity
}

#function that calculates the order mismatch for the n odors, because the WTA network
#stim:1, the set to anaylysze: 1 - all odors, 2 - origm. 3 - training, 4 - test
#op:1 type of match
cirStimulusWTAMatch <- function(circuit,stim=1,op=1){
  stimuli <- switch(stim,1:getData(circuit,const.FM$noodors),getData(circuit,const.FM$origstim)[,2],
                    getData(circuit,const.FM$trainingval)[,2],getData(circuit,const.FM$teststim[,2]) )
  before <- getData(getData(circuit,const.FM$mbkc),const.N$neurons)
  after <- getData(getData(circuit,const.FM$mbkc_net),const.N$neurons)
  #cat('\nstim',stimuli)
  res <- sapply(stimuli, function(i){
    calcVecsMatch(before[[i]],after[[i]],op=3)
  })
  sparsity <- cirGetSparsity(circuit)
  res <- mean(res)/sparsity
  res
}

#tests the effect of changing the connection params on how many of the top cells loose their posisiotn
#because of the WTA.
#op: params to manipulate: 1: no of synapses, 2 - change distribution & distribution params
# 3 : synapses and distribution
#params: the params
cirExploreWTAparams <- function(noruns=10,glomno = 50,kcno = 100,noodors = 6,aplCut = 5,classop = 2,
                                kcapl_conn_par=list(7,c(10,3,0.3,0.4)),
                                aplkc_conn_par=list(8,c(10,3,0.3,0.4)),params,op=1){
  aplkc_conn <- aplkc_conn_par;kcapl_conn <- kcapl_conn_par;
  if(op==1) {#synapses
    aplkc_conn[[2]][1] <- params;kcapl_conn[[2]][1] <- params
  }
  if(op==2) {#distribution type
    aplkc_conn[[2]][-1] <- params;kcapl_conn[[2]][-1] <- params
  }
  if(op==3) {#synpass & distribution type
    aplkc_conn[[2]] <- params;kcapl_conn[[2]] <- params
  }
  #cat('\n',str(aplkc_conn))
  res <- sapply(1:noruns, function(i){
    circuit <- createFlyCircuit(glomno = glomno,kcno = kcno,noodors = noodors,aplCut = aplCut,classop = classop,
                                kcapl_conn_par=kcapl_conn,aplkc_conn_par=aplkc_conn)
    cirStimulusWTAMatch(circuit,stim = 2)
  })
  c(mean(res),sd(res)/sqrt(noruns))
}



#liam: the next three are test and helper functions
#calculates the APL effect on mb cells, using a toy model of distributions 
#gainfactor: the amount by which the gain should be adjusted
#aplcut: the top % of Kcs that should respond
calcAPLConnTest <- function(vec,svec1,svec2,gainfactor=0.5,aplcut=20,op=3){
  #calcaultes the apl firing rates, and based on that the gain, which is 
  #calculated at that number that reduces everything below the top 20 to 0 
  apl <- as.vector(vec %*% svec1)
  gain <- getTopNPercentileThresh(vec,20)/apl
  inh.vec <- apl * svec2 * gain * gainfactor
  res <- vec - inh.vec
  #arrange the vectors by descending order of firing rate of the original mb signals
  res.mat <- getMatSortPosns(cbind(vec,res,inh.vec),decrease = T,op=op)
  res.mat
}


#function that calculates the differences in position of the neurons in the rbnaking order because of APL WTA
#veclst: the odor mb vectots as a list
#svec1 and 2: the kc-apl and apl-kc synapses
#top: the ranking order change for the top neurons
testCalcACT <- function(veclst,svec1,svec2,top=10,op=1){
  res.lst <- sapply(veclst, function(x){
    tmp <- calcAPLConnTest(x,svec1,svec2)
    #arrange the matrix by the mbkc_net response vector
    tmp1 <- getMatSortPosns(tmp,col=2,decrease = T,op=2)
    cat('\n',tmp1[,1])
    res <- sapply(1:top, function(i){
      (tmp1[i]-i) # +ve indicates that the mbkc net rank has moved up
      # if(tmp1[i,1] != i) res1 <- 1
      # else res1 <- 0
    })
    #sum(res)
  })
  res.lst
}

testGenRanCorVec <- function(vec,n=4,rho=0.7,op=1){
  cormat <- genCorrVecs(vec,n=n,rho = rho)
  for (i in 1:n) {
    tmp <- cormat[,i]
    maxn <- max(tmp)
    tmp1 <- ceiling(10 * tmp/maxn)
    #cat('\n',tmp,':',tmp1)
    cormat[,i] <- tmp1
  }
  cormat
}

# commands to run to test
# res2 <- testGenRanCorVec(runif(20,.5,.58))
# testCalcACT(res5$mbkc$neurons[1:4],res2[,1],res2[,3])


# REsults: tentatively there are two constraints here. You need to have about an equal number of synapses
# between apl and each of the KCS. If not, you can get a wide variation in synaptic strenghts, and then you can
# have the case where the highest firing neurons might not actually be the highest firing neurons after
# apl inhibition. This scenario is also the reason for having some correlation between number of synapses.
# If they are not, you can get wide variation again. It could also be the reason why you need a wide variety
# of synapses. To average out imbalances.
# Basically, noise from trial to trial in APL firing can give you noise at the threshold edges.
# Yo ualso need this to maintain a proper threshold. Otherwise, there would be wide fluctuations. 


#updates the structyre in question
updateObj <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("updateObj",self)
  }
}


#updates the conection matrix object
#paramtype: tells you which param is to be updated, specified by const.CM
#gain: updates the gain, based on params=gain
#params: the new params
updateObj.ConnMatrix <- function(self,paramtype,params,op=1){
  #for now just do them using an if and can switch when we have mroe
  if(paramtype==const.CM$gain){
    #the original synaptic weights w are obtained by w*1/gain
    #for the update just continue with w * gain and newgain = gain *params
    newgain <- getData(self,const.CM$gain) * params
    connmat <- getData(self,const.CM$connmat)
    if(isDataType(connmat)==2){#it's a matrix
      newconnmat <- newgain * connmat  
    }
    else {#it's a list of connection matrices
      newconnmat <- lapply(connmat, function(x){
        newmat <- x * params
      })
      names(newconnmat) <- names(connmat)
    }
    #update the mtrix with new values and return
    updatemat <- setData(self,val=newgain,index=const.CM$gain)
    updatemat <- setData(updatemat,val=newconnmat,index=const.CM$connmat)
    return(updatemat)
  }
  
  return(self) #nothing changed. just return the existing matrix
}


#sets the parameters for the inhibitory neurons at whatever mean rate we need to get the required sparsity
setSparsity <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("setSparsity",self)
  }
}

#sets the sparsity of the MB KC population firing
#aplcut: the sparsity level that we need
setSparsity.FlyMBg <- function(self,aplCut=5,op=1){
  #it;s the same as init now, since the apl calculation will have to be redone and
  # then everything downsteam would change, too  
  circuit.mb <- init(self,aplCut=aplCut)
  circuit.mb# <- setData(self,val=kcapl_conn,const.FM$kcapl_conn)
}


#sets the number of valences downstream of KCs, and then the MBON-KC_DAN n/w
setValences <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("setValences",self)
  }
}

# newGlomMBcircuit <- function(glomno=50,kcno=2000,danno=5,mbonno=1,noodors=4,params=c(),learningRule=c(1/38,2,0.2,T),
#                              glommb_conn=list(c(8,0.715),c(0.07),c(2,4,4)),novalences=2,tatparams=c(1,1),op=1){

#sets the number of valences downstream of KCs, and then the MBON-KC_DAN n/w
#nodans: tells you the number of dans per type of DAN
setValences.FlyMBg <- function(self,nombons=2,nodans=0,op=1){
  #get no of odors
  glom <- getData(self,const.FM$glom)
  noodors <- size(glom,op=1)
  #get no of kcs
  kcs <- getData(self,const.FM$mbkc_net)
  kcno <- size(kcs,op=2)
  #get no of dans per type of DAN
  if(nodans == 0)  {
    dans <- getData(self,const.FM$dans)
    danno <- size(dans[[1]],op=2) #no of dans per dan type
  }
  else danno <- nodans
  #get no of mbons per type of mbon
  mbons <- getData(self,const.FM$mbons)
  mbonno <- size(mbons[[1]],op=2) #no of dans per dan type
  #cat('\nSetValences',noodors,kcno,danno,mbonno)
  valences <- 1:nombons
  # # dans = # valencs
  dans <- lapply(valences, function(i) newNeurons(number = danno,dist = 0,params = c(1)) )
  names(dans) <- valences
  # mbons = # valences, each mbons also contains the mbon firing rates given a particular odor.
  mbons <- lapply(valences, function(i) newNeurons(number = nombons,dist = 0,params = c(0),noStimuli = noodors) )
  names(mbons) <- valences
  #now you have to put in the connection matrix and learning rules for the MB-MBON connections
  #13.5 synapses from every KC to a DAN; one average 2 ssynapses from every DAN to KCs
  kcmbon_conn <- lapply(valences, function(i) newConnMatrix(c(kcno,mbonno),params = c(13.5,1.15,.16),dop = list(danno,5,c(2,0.9,1.1),1) ) )
  names(kcmbon_conn) <- valences
  circuit.mb <- setData(self,val=dans,const.FM$dans) #set dans. mbons. and kcmbon_nw
  circuit.mb <- setData(circuit.mb,val=mbons,const.FM$mbons)
  circuit.mb <- setData(circuit.mb,val=kcmbon_conn,const.FM$kcmbon_conn)
  circuit.mb <- computeFiring(circuit.mb,op=1)
  #get the function for the normalization function for adjusting weights up and down
  mbkc_net <- getData(self,const.FM$mbkc_net) #just the neurons
  mbkc_net.neurons <- getData(mbkc_net,const.N$neurons) 
  for(i in valences){
    circuit.mb <- normSynWts(circuit.mb,mbkc_net.neurons,valence = valences[i],aplCut = 0)
  }
  circuit.mb
  
}  

#modifies, the aplCut, nombons or the no of original stimuli of the circuit
modifyFlyMBg <- function(circuit,aplCut=0,nombons=0,nostim=0,
                         tatparams=c(1,1),testset=1,mixparams=list(list(0.2),c(1)),trainset=1,noiseset=1,combnop=2,op=1){
  if(aplCut==0) newcircuit <- circuit #dont change the aplCut for the circuit
  else newcircuit <- setSparsity(circuit,aplCut = aplCut)
  
  if(nombons>0) {#change the number of MBONS and KCMBON n/w
    newcircuit <- setValences(newcircuit,nombons=nombons)
  }
  
  if(nostim>0){#set no of orig. stim
    #change all stimuli
    newcircuit <- setStimuli(circuit,nostimuli=nostim,stimtype=1)
    #cat('\nsize',size(newcircuit$glom,op=1))
    #change the testing and training stimuli too, since they have to connected to original stimuli
    newcircuit <- prepFlyCircuit(circuit = newcircuit, tatparams = tatparams,testset = testset,trainset = trainset,mixparams = mixparams,
                                 combnop = combnop)
  }
  newcircuit  
}

#a copy constructor: the elements that could stay the same are the PN-KC network, and the KC-MBON synapses initially
#elements that could change are the aplCut, stimuli as well as the number of stimuli, the number of MBONs, and the KCMBON n/w
#op changes the following =1, straightforward copy, 2 - training stimuli, 3 - test stimuli, 4 - both, 5 - set aplcut, 6 - MBONs and KCMBON nw, 
#7, change no of orig. stimuli,
#settype: depeneding on op, specifies the appropriate param, If op = 2, do trainset, if op = 3, follow testset, if op = 4, settype = c(trainset,testset),
#if op = 5, aplCut, if op = 6 no of MBONS, 7 - if op = 7 no of orig. stim, testset, trainset 
#trainset: 1 - original stimuli, 2 - mixed stimuli, 3 - orginal mixed with random background odor stimuli, 4 - mixed stimuli which also includes a background odor
#testset: 1 - same as training, 2 - mixed stimuli, 3 - mixed with noise, 4 - mixed stimuli which also includes a background odor
#params is the set of params: params : list(mixparams,combno)
#mixparams: list(list(mix),mixop);mix-the proportion of mixing between stimulus pairs: list(c(.1,.2)) e.g., mixing in the proportion .1:.9 and .2:.8
#mixop: how they should be mixed, by substituting neuron firing (1) or weighted averaging (2)
#combnop= options for choosing mixes of odors; 1, all possible combinations of vectors from each type of label, 2 - combinations of vectors from each type of label only once, e.g., 
copyConFlyMBg <- function(circuit,settype=1,params=list(list(list(0.2),c(1)),2),op=1){
  #get all the attributes of this circuit size-wise
  newcircuit <- circuit
  if(op==1) return(circuit)
  #initialize the input arguments, basically set to change nothing, 0 implies nothing to change
  trainset <- 1
  testset <- 0
  aplCut <- 0
  nombons <- 0
  mixparams <- params[[1]]
  combnop <- params[[2]]
  nostim <- 0
  #modfiy the input defaults based on op choice
  if(op==2){# change training stimuli, stimuli type specified by settype
    trainset <- settype
  }
  if(op==3){# change testing stimuli, stimuli type specified by settype
    testset <- settype
  }
  if(op==4){# change training and testing stimuli, stimuli type specified by settype
    trainset <- settype[1]
    testset <- settype[2]
  }
  if(op==5){#set aplCut
    aplCut <- settype
  }
  if(op==6){#change the MBON and KCMBON nw
    nombons <- settype
  }
  if(op==7){#change no of original stimuli
    nostim <- settype[1]
    trainset <- settype[2]
    testset <- settype[3]
  }
  newcircuit <- modifyFlyMBg(circuit = circuit,aplCut = aplCut,nombons = nombons,nostim = nostim)
  newcircuit <- prepFlyCircuit(circuit = circuit,testset = testset,trainset = trainset,mixparams = mixparams,
                               combnop = combnop)
}

#function that creates the fly circuit with the parameters you want, initializes it, and then you call the training update or testobj function
#noodors: determines the size of the original stimuli
#tatparams= c(trainparams,testparams); order in which they are trained; trainparams: 1 - alternate, 2 - sequential, 3 - random
#testparams: 1 - is valence training confirmed, 2 - pairwise discrimination
#noiseset: specifies the odors to which noise shuold be added. 0 all odors, 1- test set, 2 - training set
#trainset: 1 - nothing to do, keep what you have, 2 - mixed stimuli, 3 - orginal mixed with random background odor stimuli, 4 - mixed stimuli which also includes a background odor
#testset: 0 - nothing to do keep whatever you have, 1 - same as training, 2 - mixed stimuli, 3 - mixed with noise, 4 - mixed stimuli which also includes a background odor
#mixparams: list(list(mix),mixop);mix-the proportion of mixing between stimulus pairs: list(c(.1,.2)) e.g., mixing in the proportion .1:.9 and .2:.8
#mixop: how they should be mixed, by substituting neuron firing (1) or weighted averaging (2)
#combnop=1, all possible combinations of vectors from each type of label, 2 - combinations of vectors from each type of label only once, e.g., 
#if there are 4 vectors 1,2,3,4 if type 1 if 1,2 are pair then 1,4 cannot be a pair
#classop: specifies the circuit that is created. 1 - newGlomMBcircuit, 2 - newPNKCAPLcircuit
#op:1 return new circuit, 2 - list(newcircuit,oldcircuit)
#  kcmbon_conn <- lapply(1:novalences, function(i) newConnMatrix(c(kcno,mbonno),params = c(13.5,1.15,0.16),
#dop = list(danno,5,c(2,0.9,1.1),1) ) )
#thresh: set the threshold for neurons list(neuronop,threshparams). neuronop: 0 - none, 1 - mbkc, 2 - mbkc_net, 3 - mbkc and mbkcnet 
createFlyCircuit <- function(glomno = 5,kcno = 10,danno = 5,mbonno = 1,noodors=4,novalences=2,aplCut=5,learningRule=c(1/38,2,0.2,T),
                             kcapl_conn_par=list(7,c(6,1,1.15,0.12)),aplkc_conn_par=list(8,c(6,1,1.15,0.12)),
                             glommb_conn=list(c(8,0.715),c(0.07),c(2,4,4)),kcmbon_conn= list(c(13.5,1.15,0.16),list(danno,5,c(2,0.9,1.1),1)),
                             tatparams=c(1,1),noiseparams=list(c(1,0,0.1),c(1,0,0.1)),noiseop=0,testset=1,mixparams=list(list(0.2),c(1)),
                             trainset=1,noiseset=1,combnop=2,classop=1,thresh=list(0,c(1,0)),op=1){
  #choose the type of FLy circuit you want
  # circuit.type.create <- switch(classop,newGlomMBcircuit,newPNKCAPLcircuit) 
  # circuit <- circuit.type.create(noodors = noodors,glomno = glomno,kcno = kcno,danno = danno,mbonno = mbonno,novalences = novalences,
  #                             glommb_conn = glommb_conn,tatparams = tatparams,kcapl_conn_par=kcapl_conn_par,
  #                             aplkc_conn_par=aplkc_conn_par,thresh=thresh)
  circuit <- switch(classop,
                    newGlomMBcircuit(noodors = noodors,glomno = glomno,kcno = kcno,danno = danno,mbonno = mbonno,novalences = novalences,
                                     glommb_conn = glommb_conn,tatparams = tatparams,kcmbon_conn = kcmbon_conn,thresh=thresh),
                    newPNKCAPLcircuit(noodors = noodors,glomno = glomno,kcno = kcno,danno = danno,mbonno = mbonno,novalences = novalences,
                                      glommb_conn = glommb_conn,tatparams = tatparams,kcapl_conn_par=kcapl_conn_par,
                                      aplkc_conn_par=aplkc_conn_par,kcmbon_conn = kcmbon_conn,thresh=thresh)
                    ) 

  #initialize the cirucit if you need to.
  newcircuit <- init(self = circuit,aplCut = aplCut)
  newcircuit <- prepFlyCircuit(circuit = newcircuit,tatparams = tatparams,testset = testset,trainset = trainset,mixparams = mixparams,
                               combnop = combnop)
  # circuit <- init(circuit,aplCut = aplCut)
  # newcircuit <- prepTrainingSet(circuit = circuit,nostim = noodors,mix = mixparams,trainset = trainset, combnop = combnop)
  # newcircuit <- prepTestingSet(circuit = newcircuit,nostim = noodors,mix = mixparams,testset = testset,combnop = combnop)
  # newcircuit <- prepTestingSet(circuit = newcircuit,noiseparams = noiseparams,noiseop = noiseop,mixparams = mixparams,
  #                              noiseset=noiseset,trainset=1,testset = testset)
  switch(op,newcircuit,list(newcircuit,circuit))
}

#creates the testing and training stimuli for this fly circuit based on the input params.
#circuit: the circuit whose training or testing regimen is to be changed
#aplCut: 0, dont do anything since we are likely modifying an existing circuit so leave as is, if >0, change it
#noiseset: adding noise to the training circuit
#nombons: 0, dont do anything
prepFlyCircuit <- function(circuit,tatparams=c(1,1),testset=1,mixparams=list(list(0.2),c(1)),trainset=1,noiseset=1,combnop=2,
                           op=1){
  #set the training and test stimuli
  newcircuit <- prepTrainingSet(circuit = circuit,mix = mixparams,trainset = trainset, combnop = combnop)
  newcircuit <- prepTestingSet(circuit = newcircuit,mix = mixparams,testset = testset,combnop = combnop)
}

#function which prepares the training set for the circuit. There are three main options all to do with choosing 
#stimuli. They can either be pure stimuli, mix stimuli. Noiise does not make sense, so leave out stimuli with noise added for now.
#trainset: 1 - original stimuli, 2 - mixed stimuli, 3 - orginal mixed with random background odor stimuli, 4 - the original stimuli
#testset: 1 - same as training, 2 - mixed stimuli, 3 - mixed with noise, 4 - mixed stimuli which also includes a background odor
#mixparams: list(list(mix),mixop);mix-the proportion of mixing between stimulus pairs: list(c(.1,.2)) e.g., mixing in the proportion .1:.9 and .2:.8, 
#important consideration: if there are multuple mixes that should be tried, they should be their own list item: list(mix1,mix2,...) not list(c(mix1,mix2,...))
#nostim= no of original stimuli
#mix: if trainset = 2, gives the number of mixes possible with the stimuli
#combnop=1, all possible combinations of vectors from each type of label, 2 - combinations of vectors from each type of label only once, e.g., 
#if there are 4 vectors 1,2,3,4 if type 1,2,1,2, if 1,2 are pair then 1,4 cannot be a pair
prepTrainingSet <- function(circuit,mix=list(list(0.2),c(1)),trainset=1,combnop=2,op=1){
  #cat('\nTrainingset',trainset,size(circuit$glom,op=1))
  if(trainset == 1) return(circuit) #nothing to do, we already have the stimuli we need
  if(trainset == 2){#2: mixed stimuli
    newcircuit <- makeMixStimuli(circuit = circuit,mix = mix,stimtype = 2,op=1,combnop=combnop)
  }
  if(trainset == 3){#noise added stimuli
    newcircuit <- makeMixStimuli(circuit = circuit,mix = mix,stimtype = 2,op=2,combnop=combnop)
  }
  if(trainset==4){#same as original stimuli set
    origstim <- getData(circuit,const.FM$origstim)
    glom <- getData(getData(circuit,const.FM$glom),const.N$neurons) #the neuron stimuli list
    newcircuit <- setStimuli(circuit,stimuli=glom[origstim[,2]],valences=origstim[,c(1,3)],stimtype=3)
  }
  newcircuit
}


#function which prepares the testing stimuli of the circuit.
#testset: 0, do nothing 1 - same as training, 2 - mixed stimuli, 3 - mixed with background odor, 4 - the original stimuli?
#2 - different from training combination of training stimuli specified by the mixture proporations
#mixparams: list(list(mix),mixop);mix-the proportion of mixing between stimulus pairs: list(c(.1,.2)) e.g., mixing in the proportion .1:.9 and .2:.8, 
#important consideration: if there are multuple mixes that should be tried, they should be their own list item: list(mix1,mix2,...) not list(c(mix1,mix2,...))
#mixop: how they should be mixed, by substituting neuron firing (1) or weighted averaging (2)
#noiseset: specifies the odors to which noise shuold be added. 0 all odors, 1- test set, 2 - training set
#combnop=1, all possible combinations of vectors from each type of label, 2 - combinations of vectors from each type of label only once, e.g., 
#if there are 4 vectors 1,2,3,4 if type 1,2,1,2, if 1,2 are pair then 1,4 cannot be a pair
prepTestingSet <- function(circuit,mix=list(list(0.2),c(1)),testset=1,combnop=2,op=1){
  #cat('\nTsetset',testset,size(circuit$glom,op=1))
  if(testset == 0) return(circuit)
  if(testset==1){#same as training set
    trainingval <- getData(circuit,const.FM$trainingval)
    glom <- getData(getData(circuit,const.FM$glom),const.N$neurons) #the neuron stimuli list
    newcircuit <- setStimuli(circuit,stimuli=glom[trainingval[,2]],valences=trainingval[,c(1,3)],stimtype=3)
  }
  if(testset == 2){#2: mixed stimuli, will be different from trainset
    newcircuit <- makeMixStimuli(circuit = circuit,mix = mix,stimtype = 3,op=1,combnop=combnop)
  }
  if(testset == 3){#noise added stimuli
    newcircuit <- makeMixStimuli(circuit = circuit,mix = mix,stimtype = 3,op=2,combnop=combnop)
  }
  if(testset==4){#same as original stimuli set
    origstim <- getData(circuit,const.FM$origstim)
    glom <- getData(getData(circuit,const.FM$glom),const.N$neurons) #the neuron stimuli list
    newcircuit <- setStimuli(circuit,stimuli=glom[origstim[,2]],valences=origstim[,c(1,3)],stimtype=2)
  }
  # if(noiseop>0){#add noise to testing set
  #   newcircuit <- addNoise(newcircuit,neuronnoise=noiseparams[[1]],connnoise=noiseparams[[2]],noiseset=noiseset,op=noiseop)
  # }
  newcircuit
}


#given pure odors, will generate the required mixed stimuli
#op: 1 - mixthe simtuli, 2 - mix the given stimuli at the ratio indicated by mix with a random background odor
#stimtype: specifies whether it is the training or test set: 1 - train set, 2- test set
#back: random background odor either a vector of odors or two params specifying the distribution of the odor parameters
#odors: the odors that form the basis of the mix odor
#stimtype: the stimulus type to which these stimuli should be added. Default is testing stimuli
#1 - original stimuli, 2 - training stimuli, 3 - testing stimuli, default:3
makeMixStimuli <- function(circuit,mix=list(list(0.2),1),back=c(),stimtype=2,combnop=2,op=1){
  #generate the mix data, always from original stimuli
  origstim <- getData(circuit,const.FM$origstim)
  stim <- getData(getData(circuit,const.FM$glom),const.N$neurons)[origstim[,2]]
  if(op==1){#mix stimuli with each other
    stimuli.mix <- genVectorCombos(veclst = stim,combno = length(mix[[1]][[1]])+1,mix = mix[[1]],op=combnop,veclabels = origstim[,1])
    stim.df <- getValenceMixed(stimnames = names(stimuli.mix),origstim = origstim)
  }
  if(op==2){#mix stimuli with random background
    params <- getData(getData(circuit,const.FM$glom),const.N$params)
    randomodor <- newNeurons(number = length(stim[[1]]),dist = 2,noStimuli = 1,params = params)
    #stimuli.mix <- genNoiseMixtures(stimuli = stim,mix = mix[[1]],op=mix[[2]],randodor = getData(randomodor,const.N$neurons)[[1]] )
    #cat('\n',getData(randomodor,const.N$neurons)[[1]])
    stimuli.mix <- genVectorCombos(veclst = stim,combno = length(mix[[1]][[1]])+1,mix = mix[[1]],op=combnop,
                                   fixedvec = getData(randomodor,const.N$neurons)[[1]],veclabels = origstim[,1] )
    stim.df <- getValenceMixed(stimnames = names(stimuli.mix),origstim = origstim)
    #cat('\nmMS,',str(stimuli.mix),str(stim.df))
  }
  if(op==3){#mix the stimuli with each other and random background
    params <- getData(getData(circuit,const.FM$glom),const.N$params)
    randomodor <- newNeurons(number = length(stim[[1]]),dist = 2,noStimuli = 1,params = params)
    stimuli.mix <- genVectorCombos(veclst = stim,combno = length(mix[[1]][[1]])+1,mix = mix[[1]],op=combnop,
                                   fixedvec = getData(randomodor,const.N$neurons)[[1]],veclabels = origstim[,1] )
    stim.df <- getValenceMixed(stimnames = names(stimuli.mix),origstim = origstim)
  }
  newcircuit <- setStimuli(circuit,stimuli=stimuli.mix,valences=stim.df[,c(1,3)],stimtype=stimtype)

}

#outputs the odors in the form of getCampDataSigList. So, signfiicnat response vector, odor response vector
#no mean, background, and z-score.
#The input is actually a bunch of circuits, the output of exploreNoise
#op: 1 - all odors, 2 - origstim, 3 - training values, 2 - test stim
cirGetFlyStimOut <- function(circuit,op=1){
  mbkc.neurons <- getData(getData(circuit,const.FM$mbkc_net),const.N$neurons)
  #now, filter for the responses specified by op
  stimuli <- switch(op,1:getData(circuit,const.FM$noodors),getData(circuit,const.FM$origstim)[,2],
                    getData(circuit,const.FM$trainingval)[,2],getData(circuit,const.FM$teststim[,2]) )
  mbkc.neurons <- mbkc.neurons[stimuli]   
  sig.neurons <- lapply(1:length(mbkc.neurons), function(i) {
    sig <- setAboveThreshOne(mbkc.neurons[[i]])
  })
  siglist <- lapply(1:length(mbkc.neurons), function(i) list(sig.neurons[[i]],mbkc.neurons[[i]]))
  names(siglist) <- stimuli
  siglist  
}

#gets the valences for mixed stimuli
#stimnames: the mixed stimuli names
getValenceMixed <- function(stimnames,origstim,op=1){
  newstim <- sapply(1:length(stimnames), function(i){
    odor <- as.numeric(strsplit(stimnames[i],',')[[1]][1]) #odor name
    posn <- which(origstim[,2] == (odor) )
    val <- ifelse(length(posn)>0,origstim[posn,1],0) #if the random odor is greater, valence is 0 
    c(val,i,stimnames[i]) #return valence, odor #, and stimulus name containing mix proportions
  })
  newstim <- as.data.frame(t(newstim))
  names(newstim) <- c('val','odor#','mix')
  newstim <- ConvertDfCols(newstim,cols = c(1,2),op=1)
  newstim <- ConvertDfCols(newstim,cols = c(3),op=2)
  newstim
}

#this function trains the circuit for a particular odor or set of odors, or a user-specified odor
#generic method to add stuff to thee objects 
trainObj <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("trainObj",self)
  }
}

#function that trains the circuit for a particular odor or set of odors. Only training. no computations of odor rates
#odors = the odor nos to train for the specified number of trials
#trials: the no of trials to train for
#noiseparams: the noise parameters for training, specifies the components to which noise should be added
#1: noise to PNs, 2 - noise to 
#if a vector (v1,v2,...) then v1, v2, v3 for odor 1,2,3...
#op=1
trainObj.FlyMBg <- function(self,odors=c(),trials=12,op=1){
  trainedcircuit <- self #initialize trained circuit
  trainingvals <- getData(self,const.FM$trainingval) #get the training valences and odor nos
  valen <- trainingvals[,1]
  odornos <- trainingvals[,2]
  #go through all the odor nos specified by odornos and trials   
  #we wil have to us a for loop to go hrough all the odors
  for(i in seq_len(length(odornos)) ){ # go through each odor, seq_len to avoid backward sequences
    for(j in seq_len(trials)){# now, for each odor, go through all of its trials
      #train the circuit for this odor with this valence
      #cat('\nodor ',i,'trial ',j,'valence ',valen[i],';')
      trainedcircuit <- trainingUpdate(trainedcircuit,odorno = odornos[i],valence = valen[i])
    }
  }
  trainedcircuit
}
  

#this function tests the circuit for a particular odor, or a user-specified odor
#generic method to add stuff to thee objects 
testObj <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("testObj",self)
  }
}

#computes the mbon firing rates for the test odors. 
#the animals discrimination ability is the relative firring rates of these two mbons.
#odor: test the network for these odors, in the form of a list of stimuli
#valences: the valences for the input odors, they have to be named.
#if a vector (v1,v2,...) then v1, v2, v3 for odor 1,2,3...
#op=1
testObj.FlyMBg <- function(self,odors=c(),valences=c(),op=1){
  testcircuit <- self
  if(length(odors)>0) {#set teststim to point to this data
    if(length(names(valences)) > 0 ) vnames <- names(valences)
    else vnames <- 1:length(valences)
    newcircuit <- setStimuli(circuit,stimuli=odors,valences=cbind(valences,vnames),stimtype=3)
  } 
  else  newcircuit <- self #get the odor nos for testing
  #go through the test odors and compute their firiing rates
  newcircuit <- computeFiring(newcircuit,op=1) #will comput the firing rate for the mbons for odor in this new n/w
  
}

#checks the value of the specified circuit components
#comp: component to be checked. 1 - mbkc_net, 2 - kcmbon_conn, 3 - mmbon
#subcomp: if neurons specify the stimuli no, if kcmbon_conn, specify the valence that should be checked, if mbon,subcomp is odor no
#op=1 - individual components, 2 - all components in df format 
checkCircuitComp <- function(circuit,comp=2,subcomp=1,op=1){
  temp <- switch(comp,getData(circuit,const.FM$mbkc_net),getData(circuit,op=3,valence=subcomp),
                 getData(circuit,const.FM$mbons))
  #now to get the subcomponent
  if (comp==1){#neurons
    res <- getData(temp,op = 3,stimno=subcomp)
  }
  if (comp==2){#kcmbon_conn
    res <- temp
  }
  if(comp==3){
    res <- lapply(temp, function(i) getData(i,op=3,stimno=subcomp))
  }
  if(op==2 || c('FlyMBg') %in% class(circuit)){
    res.mbkcnet <- getData(getData(circuit,const.FM$mbkc_net),const.N$neurons)
    res.mbkcnet <- transposeDF(convertNestedListsDF(res.mbkcnet))
    rownames(res.mbkcnet) <- 1:nrow(res.mbkcnet)
    res.mbon <- lapply(getData(circuit,const.FM$mbons), function(i) getData(i,const.N$neurons))
    res.mbon <- transposeDF(convertNestedListsDF(res.mbon))
    row.names(res.mbon) <- 1:nrow(res.mbon)
    colnames(res.mbon) <- 1:length(getData(circuit,const.FM$valences))
    res.kcmbon <- lapply(getData(circuit,const.FM$kcmbon_conn), function(i) getData(i,const.CM$connmat))
    res.kcmbon <- transposeDF(convertNestedListsDF(res.kcmbon))
    colnames(res.kcmbon) <- 1:ncol(res.kcmbon)
    rownames(res.kcmbon) <- 1:length(getData(circuit,const.FM$valences))
    
  }
  list(res.mbkcnet,res.mbon,res.kcmbon)
}


#functions to test capacity and other things.
#nostim: number of original stimuli
#params: for the different neurons and connection matrices,
#testop: thr strategy to be used for testing efficacy of training, 1- testing if the valence that is trained stays
#2 - pairwise discrimination between odors that were trained for different valences
#op=1, alernatiing valences, 2 - consecutive valences like (1,1,1,1,2,2,2,2), 3 - valences divided randonly
#noisparams: c(1,0,.1) default gaussion noise 10%, 1 stand for gaussian noise
#noiseop = the noise option, default 0 - no noise, 1 - noise in glomeruli, added to the test stimuli
#trainset: 1 - original stimuli, 2 - mixed stimuli, 3 - orginal mixed with random background odor stimuli, 4 - mixed stimuli which also includes a background odor
#testset: 1 - same as training, 2 - mixed stimuli, 3 - mixed with noise, 4 - mixed stimuli which also includes a background odor
#mixparams: list(list(mix),mixop);mix-the proportion of mixing between stimulus pairs: list(c(.1,.2)) e.g., mixing in the proportion .1:.9 and .2:.8, 
#important consideration: if there are multuple mixes that should be tried, they should be their own list item: list(mix1,mix2,...) not list(c(mix1,mix2,...))
#mixop: how they should be mixed, by substituting neuron firing (1) or weighted averaging (2)
#noiseset: specifies the odors to which noise shuold be added. 0 all odors, 1- test set, 2 - training set
trainAndTest <- function(tcircuit=c(),nostim=4,notrials=6,novalence=2,kcno=10,glomno=5,params=c(),aplCut=5,tatparams=c(1,1),trainset=2,
                         noiseparams=list(c(1,0,0.1),c(1,0,0.1)),noiseop=0,testset=1,mixparams=list(list(0.2),1),noiseset=1,op=1){
  pt <- proc.time()
  #trains and tests a system of the specified size
  if (length(tcircuit) > 0 ) {
    circuit <- tcircuit
    nostims <- size(getData(circuit,const.FM$glom))
    novalences <- length(getData(circuit,const.FM$valences))
    circuit.nolearn <- computeFiring(tcircuit,op=1) #circuit without the learning or noise. mbons computed too
    circuit <- addNoise(circuit,op=noiseop,noiseneuron=noiseparams[[1]],connnoise=noiseparams[[2]],noiseset=1) #add noise to the testing circuit
  }
  else {
    nostims <- nostim
    novalences <- novalence
    circuit <- createFlyCircuit(noodors = nostim,aplCut = aplCut,kcno = kcno,novalences = novalences,tatparams = c(op,1),glomno = glomno,
                                trainset = trainset,noiseparams = noiseparams,noiseop = noiseop,testset=testset,mixparams = mixparams,noiseset=noiseset)
    circuit.nolearn <- computeFiring(circuit,op=1) #circuit without the learning or noise. mbons computed too
    #circuit <- addNoise(circuit,op=noiseop,noiseneuron=noiseparams[[1]],connnoise=noiseparams[[2]]) #add noise to the circuit
  }
  pt1 <- proc.time()
  #train and test
  circuit <- trainObj(circuit,trials=notrials) #the learning part
  circuit <- testObj(circuit) #computes the mbon outputs with this new matrix
  pt2 <- proc.time()
    
  #get the output of the mbons and see if they match with the positive and negative valences
  anals.res <- analyzeTestObj(circuit = circuit)
  pt3 <- proc.time()
  # cat('\nTrain and Test Time - analysis: ',(pt3-pt2)[1:3],'training and compute time: ',(pt2-pt1)[1:3],'setup time:',(pt1-pt)[1:3],
  #     '\tnotrials',notrials)
  c(anals.res,list(circuit.nolearn)) #return accuracy df, new n/w, accuracy, old n/w
}

#functions to test capacity and other things. Same as earlier, excspt for noise and mixtures
#size: of the number of responses you want to test
#params: for the different neurons and connection matrices,
#testop: thr strategy to be used for testing efficacy of training, 1- testing if the valence that is trained stays
#2 - pairwise discrimination between odors that were trained for different valences
#op=1, alernatiing valences, 2 - consecutive valences like (1,1,1,1,2,2,2,2), 3 - valences divided randonly
#noisparams: c(1,0,.1) default gaussion noise 10%
#noiseop = the noise option, default 0 - no noise, 1 - noise in glomeruli, added to the test stimuli
#testset: how should the testing stimuli be specified: 0- same as training, 1 - same as training but separate
#2 - different from training combination of training stimuli specified by the mixture proporations
#does not apply to input circuits
#mixparams: list(mix,mixop);mix-the proportion of mixing between stimulus pairs: c(.1,.2) e.g., mixing in the proportion .1:.9 and .2:.8
#mixop: how they should be mixed, by substituting neuron firing (1) or weighted averaging (2)
#noiseset: specifies the odors to which noise shuold be added. 0 all odors, 1- test set, 2 - training set
trainAndTestNoiseMix <- function(tcircuit=c(),nostim=4,notrials=6,novalence=2,kcno=10,glomno=5,params=c(),aplCut=5,tatparams=c(1,1),
                         noiseparams=list(c(1,0,0.1),c(1,0,0.1)),noiseop=0,testset=1,mixparams=list(c(0.1,0.2),1),noiseset=1,op=1){
  pt <- proc.time()
  #trains and tests a system of the specified size
  if (length(tcircuit) > 0 ) {
    circuit <- tcircuit
    nostims <- size(getData(circuit,const.FM$glom))
    novalences <- length(getData(circuit,const.FM$valences))
    circuit.original <- computeFiring(circuit,op=1) #circuit without the learning or noise. mbons computed too
    circuit <- prepTestingSet(circuit,noiseparams = noiseparams,noiseop = noiseop,noiseset = noiseset,mixparams = mixparams,testset = testset)
    anals.naive <- analyzeTestObj(circuit = circuit.original) #results of the naive circuit
  }
  else {
    nostims <- nostim
    novalences <- novalence
    circuitlst <- createFlyCircuit(noodors = nostim,aplCut = aplCut,kcno = kcno,novalences = novalences,tatparams = tatparams,
                                noiseparams = noiseparams,noiseop = noiseop,testset=testset,mixparams = mixparams,noiseset=noiseset,op=2)
    circuit.original <- computeFiring(circuitlst[[2]],op=1) #circuit without the learning or noise. mbons computed too
    circuit <- circuitlst[[1]]
    anals.naive <- analyzeTestObj(circuit = circuit.original) #results of the naive circuit
  }
  #for the original circuit, have a separate duplicate testing set
  
  pt1 <- proc.time()
  #train and test both circuits
  circuit <- trainObj(circuit,trials=notrials) #the learning part
  circuit.original <- trainObj(circuit.original,trials=notrials) #the learning part
  circuit <- testObj(circuit) #computes the mbon outputs with this new matrix
  circuit.original <- testObj(circuit.original) #computes the mbon outputs with this new matrix
  pt2 <- proc.time()
  
  #get the output of the mbons and see if they match with the positive and negative valences
  anals.res <- analyzeTestObj(circuit = circuit)
  anals.res.orig <- analyzeTestObj(circuit = circuit.original)
  pt3 <- proc.time()
  # cat('\nTrain and Test Time - analysis: ',(pt3-pt2)[1:3],'training and compute time: ',(pt2-pt1)[1:3],'setup time:',(pt1-pt)[1:3],
  #     '\tnotrials',notrials)
  #return accuracy df, new n/w, accuracy, valences,n/w w/o noise or mixing, accuracy df, accuracy,acccuracy df of naive n/w, accuracy
  c(anals.res,list(circuit.original),anals.res.orig[c(1,3)],anals.naive[c(1,3)]) 
}


#analyzes the results of the test of the circuit
analyzeTestObj <- function(circuit,op=1){
  tatparams <- getData(circuit,const.FM$tatparams)
  nostims <- nrow(getData(circuit,const.FM$teststim)) #no of testing stimuli
  #get the output of the mbons and see if they match with the positive and negative valences
  valences <- getData(circuit,const.FM$teststim)[,1] #the valences used for testing
  if(tatparams[2] == 1){#test if the trained valence is learnt
    res.df <- testCorrectValence(circuit = circuit,valences = valences)
    accuracy <- sum(res.df[,4])/nostims
  }
  if(tatparams[2] == 2) {#test if two odors of different valences can be distinguished
    res.df <- testPairs(circuit = circuit,valences = valences)
    #accuracy <- sum(res.df[,ncol(res.df)])/nrow(res.df)
    accuracy <- sum(res.df$cor)/nrow(res.df)
  }
  list(res.df,circuit,accuracy,valences)
}

#view the results of the TAT side by side with the original circuit, i.e., the one before learning.
#tatobj: the result of tatobj
#op= show the behavior results side by side
viewTATObj <- function(tatobj,op=1){
  #show them as side by side columns
  res1 <- cbind(analyzeTestObj(tatobj[[5]])[[1]],analyzeTestObj(tatobj[[2]])[[1]])
  res2 <- c(analyzeTestObj(tatobj[[5]])[[3]],analyzeTestObj(tatobj[[2]])[[3]],
            (analyzeTestObj(tatobj[[2]])[[3]]-analyzeTestObj(tatobj[[5]])[[3]])/analyzeTestObj(tatobj[[5]])[[3]])
  list(res1,res2)
}

#tests if the mbon results are the same as the valence that they were trained for
#op=1, the valence is the MBON with the lower firing rate, 2 - the ratio of the two MBONs, and apply threshold, 3 - absolute threshold. The MBON 
#with firing rate lower than threshold. If both are lower, choose the lower one 
testCorrectValence<- function(circuit,valences,thresh=1,op=1){
  #cat('\ntestCorrect ',valences)
  novalences <- getData(circuit,const.FM$valences)
  nostims <- nrow(getData(circuit,const.FM$teststim)) #no of stimuli
  teststim <- getData(circuit,const.FM$teststim)[,2]
  #get the output of the mbons and see if they match with the positive and negative valences
  res <- sapply(1:length(teststim),function(i){#for each odor see which valence is lowest
    mbons <- sapply(1:length(novalences),function(j) mean(getData(circuit,op=5,odorno=teststim[i],valence=j)) )
    pos <- which(min(mbons) == mbons)[1] # include [1] to take the first instance in case there are multiple mins
    correct <- ifelse(valences[i]==pos,1,0) 
    c(pos,correct,mbons)
  })
  #record number of incorrect as a fraction
  mix <- getData(circuit,const.FM$teststim)$mix
  res.df <- as.data.frame(t(rbind(teststim,valences,res) ),stringsAsFactors=F) #data frame to accomodate the 'mix' string
  colnames(res.df)[1:4] <- c('odor#','val','beh','correct')
  res.df <- insertColDf(res.df,newcol = mix,posn = ncol(res.df)+1,colname = 'mix')
  res.df
}

#test for pairwise discrimination. At the end of training, this is what we would like. To see if given 
# 2 odors, the animal will move towards one odor and away from another
#circuit: the circuit which has been fully trained
#valences: the order in which odors have been trained to different valences
#op
testPairs <- function(circuit,valences,op=1){
  novalences <- getData(circuit,const.FM$valences) #get the different valiences
  noval <- length(novalences) #no of valences
  nostims <- length(valences) #no of stimuli
  teststim <- getData(circuit,const.FM$teststim)
  res <- sapply(1:nrow(teststim)  ,function(i){#for each odor get each valence mbons firing
    #get the mean firing rate for all the mbons of that valence
    mbons <- sapply(1:noval,function(j) mean(getData(circuit,op=5,odorno=teststim[i,2],valence=j)) )
    c(mbons,teststim[i,1]) #also record the trained valence for this odor
  })
  res <- t(res) #so that mbons and valences are in the columns
  #now, with all the mbons calculated, start doing comparisons between pairs of odors with
  #opposing valences
  #first, make a list of all valences
  val.lst <- lapply(1:noval,function(i) {
    posns <- which(valences==i) #all the posns that have this valence
    otherposns <- setdiff(1:length(valences),posns) #posns of other valences
    #compare this valence with others, 
    res.comp <- lapply(posns,function(j){
      vals <- lapply(otherposns,function(k){
        #see if odor j's trained valence mbon is lower and odor k's trained valence mbon is lower
        res.mbon <- ifelse( (res[j,res[j,noval+1]] < res[j,res[k,noval+1]]) && (res[k,res[k,noval+1]]< res[k,res[j,noval+1]]) ,1,0)
        list(c(teststim[j,2],teststim[k,2],res[j,res[j,noval+1]],res[j,res[k,noval+1]],res[j,noval+1],res[k,res[j,noval+1]],
          res[k,res[k,noval+1]],res[k,noval+1],res.mbon),c(teststim$mix[j],teststim$mix[k]) )
      })
      res.vals <- transposeDF(convertNestedListsDF(lapply(vals,'[[',1) ) )
      name.vals <- lapply(vals,'[[',2)
      colnames(res.vals) <- c('od1','od2','mbon1','mbon2','val1','mbon1','mbon2','val2','cor')#'mixod1','mixod2'
      #get the correspondinng mix odor information, too
      res.vals <- insertColDf(res.vals,newcol = as.character(unlist(lapply(name.vals,'[[',1))),posn = ncol(res.vals)+1,colname = 'mixod1')
      res.vals <- insertColDf(res.vals,newcol = as.character(unlist(lapply(name.vals,'[[',2))),posn = ncol(res.vals)+1,colname = 'mixod2')
      res.vals
    })    
    names(res.comp) <- posns
    res.comp
  })
  names(val.lst) <- 1:length(novalences)
  val.lst <- flattenLists(val.lst)
  val.dfs <- joinListDFs(val.lst)
  val.dfs #list(res,val.dfs)
}

#do it a bunch of times to get averages
#no: number of runs
#copyop: specifies whether successive runs should use the same FlyMbg or new ones 1: different (default), 2 = the same one
#op=1, return the valence performance df, accuracy, new circuit, old circuit, 2 - the valence performance df, 3 - accuracy
trainAndTestLst <- function(no,tcircuit=c(),nostim=4,notrials=6,novalences=2,kcno=10,glomno=5,params=c(),aplCut=5,tatparams=c(1,1),
                            trainset=2,noiseparams=list(c(1,0,0.1),c(1,0,0.1)),noiseop=0,mixparams=list(list(0.2),1),noiseset=1,
                            testset=1,valenceop=1,copyop=1,op=1){
  #cat('\nIn train and test LST')
  if(length(tcircuit) == no){#same circuit for all the runs, don't need it now, but will for stochastic runs
    #the circuits for the runs are stored in tcircuit, so length(tcircuit) = no
  }
  res <- lapply(1:no,function(i) {
    temp <- trainAndTest(tcircuit = tcircuit,nostim = nostim,notrials = notrials,novalence = novalences,kcno = kcno,params = params,aplCut = aplCut,
                         tatparams = tatparams,trainset = trainset,noiseparams = noiseparams,noiseop = noiseop,mixparams = mixparams,
                         noiseset = noiseset,testset = testset,op=valenceop)
    temp <- switch(op,temp[c(1,3,2,5)],temp[c(1)],temp[c(3)],temp[c(2)],temp[c(5)]) #the valence performance df, accuracy, new circuit, old circuit
  })
  res
}

#check the params of trainAndTestLstSeq to make sure that they are ok
checkTATParams <- function(params,paramno=1,op=1){
  res <- params
  if(paramno==5){
    vals <- sapply(params,function(i) ifelse(i>=1 && i<= 3,T,F))
    res <- res[vals]
  }
  if(paramno==6){
    vals <- sapply(params,function(i) ifelse(i<2,T,F))
    res <- res[vals]
  }
  res
}

#do it while you vary your parameters along a sequence
#paramno= parameter over which we should iterate, 1 - no stim, 2 - no of trials, 3 - kcno, 4 - aplCut, 5 - novalences
#the sequence should be given in the usual form
#tatparams
#op=1, manipulate one parameter, 2 - 2 parameters
#trainset: 1 - original stimuli, 2 - mixed stimuli, 3 - orginal mixed with random background odor stimuli, 4 - mixed stimuli which also includes a background odor
#testset: 1 - same as training, 2 - mixed stimuli, 3 - mixed with noise, 4 - mixed stimuli which also includes a background odor
trainAndTestLstSeq <- function(no,tcircuit=c(),nostim=4,notrials=6,novalences=2,kcno=10,glomno=5,params=c(),aplCut=5,paramno=1,
                               tatparams=c(1,1),trainset=2,noiseparams=list(c(1,0,0.1),c(1,0,0.1)),noiseop=0,mixparams=list(list(0.2),1),
                               testset=1,noiseset=1,op=1){
  #get the paramno that is changing
  par <- switch(paramno,nostim,notrials,kcno,aplCut,novalences) #the changeable parameter
  par <- checkTATParams(par)
  #the other parameters that are fixed for this round, but can also changee
  otherpars <- switch(paramno,list(notrials=notrials,kcno=kcno,aplCut=aplCut,novalences=novalences),
                      list(nostim=nostim,kcno=kcno,aplCut=aplCut,novalences=novalences),
                      list(notrials=notrials,nostim=nostim,aplCut=aplCut,novalences=novalences),
                      list(notrials=notrials,kcno=kcno,nostim=nostim,novalences=novalences),
                      list(notrials=notrials,kcno=kcno,nostim=nostim) )
  res <- lapply(par,function(i){
    cat('\nParam ',paramno,'\tval',i)
    tmpparams <- c(otherpars,list(no=no,glomno=glomno,params=params,tatparams=tatparams,noiseparams=noiseparams,noiseop=noiseop,
                                  trainset=trainset,testset=testset,mixparams=mixparams,noiseset=noiseset)) #the third set are params that are fixed
    pars <- switch(paramno,list(nostim=i),list(notrials=i),list(kcno=i),list(aplCut=i),list(novalences=i) )
    allparams <- c(pars,tmpparams)
    if(copyno==1){#you want the same circuit across runs, so get the circuit first
        
    }
    tmp <- do.call(trainAndTestLst,allparams)
    #results are in the form of #the valence performance df, accuracy, new circuit, old circuit
  })
  names(res) <- par #so that we can identify it
  #get the results in the form of of the indices, sort of a transpose of lists
  res.data <- lapply(1:length(res[[1]][[1]]), function(i){
    lapply(res,function(j) {
      lapply(j,'[[',i)
    })
  })
  res.data #the valence performance df, accuracy, new circuit, old circuit
}


#analyze the data returned by trainAndTestLstSeq
#they are in the form of #the valence performance df, accuracy, new circuit, old circuit
#op=1, analyze the precent right, 2 - the  improvement in accuracy after comparing with 
#field, 1 - return df, col1= paramno, col2 = measure
analyzeTATLst <- function(dat.lst,field=1,op=1){
  #takes the percent right for the different runs and gives the statistics in a form that can be plotted
  if(op==1){
    res <- sapply(dat.lst[[2]], function(i){
      #cat(str(i))
      getStats(unlist(i))
    })
  }
  if(op==2){
    res <- sapply(1:length(dat.lst[[3]]), function(i){
      improvement <- sapply(1:length(dat.lst[[3]][[i]]), function(j) {
        new.accuracy <- analyzeTestObj(dat.lst[[3]][[i]][[j]])[[3]]
        old.accuracy <- analyzeTestObj(dat.lst[[4]][[i]][[j]])[[3]]
        #cat(new)
        improvement <- (new.accuracy-old.accuracy)/old.accuracy
      })
      getStats(improvement)
    })
    colnames(res) <- names(dat.lst[[3]])
  }
  res <- t(res)
  #cat('here',names(dat.lst[[3]]))
  res <- cbind(as.numeric(rownames(res)),res)
  res
}

#prepares the cuircuits for reruns of the trainAndTestLstSeq calls
#allpars: all the param combinations for the network, a DF
#i: the particular instance or row number from the all pars
#orig.nw: the list of original networks which are to be modified
prepTATLCircuits <- function(allpars,i,orig.nw,op=1){
  default.pars <- pars.all[1,]
  change.pars <- which( (pars.all[i,] - pars.all[1,]) != 0 ) #the params diff from the original
  if(3 %in% change.pars){#kcnos, change all
    return(c()) #nothing to keep  
  }
  new.circuits <- lapply(orig.nw,function(j) j[[4]]) #get the original networks
  if(1 %in% change.pars){#change nostim: keep n/w and IC, change stim
    new.circuits <- lapply(1:length(new.circuits),function(j) copyConFlyMBg(new.circuits[[j]],settype = 1))
  }
  if(2 %in% change.pars){#change notrials: keep n/w and IC, stimuli
    new.circuits <- new.circuits
  }
  if(4 %in% change.pars){#change nostim: keep n/w and IC, stim. only change aplCut
    new.circuits <- lapply(1:length(new.circuits),function(j) copyConFlyMBg(new.circuits[[j]],settype = 1))
  }
  if(5 %in% change.pars){#change nostim: keep n/w and IC, stim, change KC_MBON_conn
    new.circuits <- lapply(1:length(new.circuits),function(j) copyConFlyMBg(new.circuits[[j]],settype = 1))
  }
  new.circuits #return the modified circuits
}

testfn1 <- function(op=1,...){
  
  
}


#do it while you vary your 2 parameters along a sequence
#paramno= parameter over which we should iterate, 1 - no stim, 2 - no of trials, 3 - kcno, 4 - aplCut, 5 - novalences
#the sequence should be given in the usual form
#tatparams
#trainset: 1 - original stimuli, 2 - mixed stimuli, 3 - orginal mixed with random background odor stimuli, 4 - mixed stimuli which also includes a background odor
#testset: 1 - same as training, 2 - mixed stimuli, 3 - mixed with noise, 4 - mixed stimuli which also includes a background odor
#op=1, return the valence performance df, accuracy, new circuit, old circuit, 2 - the valence performance df, 3 - accuracy
trainAndTestLstsSeq <- function(no,tcircuit=c(),nostim=4,notrials=6,novalences=2,kcno=10,glomno=5,params=c(),aplCut=5,paramnos=c(1,2),
                               tatparams=c(1,1),trainset=2,noiseparams=list(c(1,0,0.1),c(1,0,0.1)),noiseop=0,mixparams=list(list(0.2),1),
                               testset=1,noiseset=1,op=1){
  #get the paramnoz that are changing
  pars.lst <- list(nostim,notrials,kcno,aplCut,novalences)
  pars <- pars.lst[paramnos]
  pars <- permuteGroups(pars,k=length(pars),op=2)
  otherparnos <- setdiff(1:length(pars.lst),paramnos)
  otherpars <- pars.lst[otherparnos]
  #generate a DF where each row is a list of params
  pars.all <- cbind.data.frame(pars,repVecMat(unlist(otherpars),n = nrow(pars)) )
  pars.all <- reorderDFCols(pars.all,c(paramnos,otherparnos) ) #reorder so the params to be in the order 1=stim,...
  #other fixed params
  fixedparams <- list(no=no,glomno=glomno,params=params,tatparams=tatparams,noiseparams=noiseparams,noiseop=noiseop,
                    trainset=trainset,testset=testset,mixparams=mixparams,noiseset=noiseset) #the third set are params that are fixed
  if(rerun==1){#you want rerun the same circuits across runs
    #run it once for the default params and get the circuts, You can use the same circuits for the different sequences
    varpars <- c(list(nostim=pars.all[1,1],notrials=pars.all[1,2],kcno=pars.all[1,3],aplCut=pars.all[1,4],novalences=pars.all[1,5]) )
    allparams <- c(varpars,fixedparams,list(op=1))
    orig.res <- do.call(trainAndTestLst,allparams)
  }
  #now run for all param instances
  res <- lapply(1:nrow(pars.all),function(i){
    #cat('\nParam vals ',unlist(pars.all[i,]))
    if(rerun==0) circuit.par <- tcircuit #either keep the default circuits for each run or generate new circuits everytime
    else circuit.par <- prepTATLCircuits(pars.all,i,orig.res)
    varpars <- c(list(nostim=pars.all[i,1],notrials=pars.all[i,2],kcno=pars.all[i,3],aplCut=pars.all[i,4],novalences=pars.all[i,5]) )
    allparams <- c(varpars,fixedparams,list(op=op,tcircuit=circuit.par))
    tmp <- do.call(trainAndTestLst,allparams)
    #results are in the form of #the valence performance df, accuracy, new circuit, old circuit
  })
  names(res) <- concatenateMatStr(cbind(repVecMat(paramnos,n = nrow(pars)),pars) )#so that we can identify it
  #for now, just get the percentages
  #get the results in the form of of the indices, sort of a transpose of lists
  cat('\nTATlst',length(res[[1]][[1]]),length(res))
  if(op==1) nores <- 4
  else nores <- 1
  res.data <- lapply(1:nores, function(i){
    #cat('\n i',i)
    lapply(res,function(j) {
      #cat('\t j',length(j) )
      lapply(j,'[[',i)
    })
  })
  res.data[c(1,2)] #the valence performance df, accuracy, new circuit, old circuit
}

analyzeTATLstsSeq <- function(tatres,op=1){
  #get the results of the accuracy runs the mean
  resvec <- sapply(tatres[[1]], function(i) mean(unlist(i)) )
  resnames <- strsplit(names(resvec),split = ',')
  nopars <- length(resnames[[1]])/2
  #cat('\nnopars',nopars)
  namevec <- sapply(resnames,function(i) as.numeric(i)[-(1:nopars)])
  namevec <- t(namevec)
  name.row <- unique(namevec[,1])
  name.col <- unique(namevec[,2])
  matval <- matrix(rep(0,length(name.row)*length(name.col)),nrow = length(name.row))
  #assign row and col names
  rownames(matval) <- sort(name.row)
  colnames(matval) <- sort(name.col)
  #cat(str(matval))
  for(i in seq_wrap(1,length(resvec))){
    #cat('\t',namevec[i,1],namevec[i,2])
    matval[as.character(namevec[i,1]),as.character(namevec[i,2])] <- resvec[i]
  }
  matval
}

#function that explores noise results for sequence of parameters over multiple runs
#All parameters except runs are the same as exploreNoiseSeq
#resop: option for the how the results output should be formatted, 1 = mbkc_net firing rates, 2 - mbkc, 3 - glomerular 
#4 - list(glomerular,mbkcnet) , 10 - all the circuits
#runs: no of times this is run
##verbose: F - nothing, T - keeps track of how far along you are in the simulations 
exploreNoiseSeqRuns <- function(self,noruns,notrials,noiserange=seq(0,0.3,0.1),noiseparams=list(1,c(0,0.1),1),
                            noiseset=1,outcomp=1,verbose=F,op=1,resop=4){
  res.ls <- lapply(1:noruns,function(x){
    if(verbose) cat('\n',x,' out of ',noruns,' runs.')
    res <- exploreNoiseSeq(self = self,notrials = notrials,noiserange = noiserange,noiseparams = noiseparams,
                           noiseset = noiseset,outcomp = outcomp,verbose=verbose,op = op,resop = resop)
  })  
  names(res.ls) <- 1:noruns
  res.ls
}


#this function explores the effect of noise and circuit parameters. So, as circuit params like the pn->kc connections change,
#what is the effect of noise. Does it improve with the number of PN-KC claws? It here is the amount of noise conveyed.
#circparams: the circuit parameters: list(#glom,#kcs,aplkcgain,#mbons,#dans,pn-kc params,kc-mbon params,dan-mbon params,learning rule)
#circop: list(levela): each element specifies the parameter from circparams that needs to be modified. They could be atomic or a 
#nested list, e.g., the multiple nested parameters that specify connection matrices, specify the param number to be changed of the form
#levela = list(c(paramno from level1,sub-paramno),...), e.g. level2 = list(c(6,1,1)) would modify glommbconn: sub-param 1,subsubparam 1
#circmod: the amount by which the circuit parameters should be changed, could be an additive percentage: so
#list(c(step %,left,right),otherparams); specifies how far to the left and right of the present value should we go, and step 
#specifies the step size. other params if needed. If only one param, apply this to all parameters, else one for every paramter.
#targetpars: the target parameters against which you should be comparing results, the default are the fly params
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 
#12 - DAN-KC-MBON connections, 13 - noise to everything, 0 no noise, 
#100 - add noise according to noise params to every structure
exploreNoiseCircuitParams<-function(noruns=2,notrials=2,noodors=6,noiserange=seq(0,0.3,0.1),noiseparams=list(1,c(0,0.1),1),
                                 circparams=list(glom=50,kcs=150,aplgain=6,mbons=2,dans=5,glommb=list(c(8,0.715),c(0.07),c(2,4,4)),
                                                 kcmbon=c(13.5,1.15,0.16),danmbon=list(5,c(2,0.9,1.1),1),learning=c(1/38,2,0.2,T)),noiseset=2,outcomp=1,
                                 op=1,resop=4,circmod=list(10,1,1),circop=list(1),chooesop=2,retop=1,
                                 targetpars=list(c(0.72,0.07),c(5.26,0.7),c(7.23,0.84),c(6.1,1.9),c(29,3.16)) ){
  # exploreNoiseCircuitSeq <- function(self,notrials,noiserange=list(seq(0,0.3,0.1)),noiseparams=list(1,c(0,0.1),1),
  #                                            noiseset=1,circuitoutcomp=1,circparams=list(c(-1,1,1)),circop=1,op=1,resop=1){
  #first change the circuit params, and then just call exploreNoiseSeq with the other params intact
  #generate the param range on both sides
  assignVarVals(lstval = changeCircuitParams(circparams = circparams,circmod = circmod,circop = circop),
                lstvar = c('par.combo','parlst')) 
  list(par.combo,parlst)
  
  res.ls <- lapply(1:length(parlst),function(i){#go through each parameter list combination
    #create the fly circuit for the runs
    circuit <- createFlyCircuit(glomno = parlst[[i]]$glom,kcno = parlst[[i]]$kcs,noodors = noodors,aplCut = parlst[[i]]$aplgain,
                                danno = parlst[[i]]$dans,mbonno = parlst[[i]]$mbons)
    runs <- exploreNoiseSeqRuns(circuit,noruns=noruns,notrials = notrials,noiseset = noiseset,noiserange = noiserange,
                                noiseparams = noiseparams,op=op,resop = resop)
    #cat('\nexploreNoisecircuitpars',str(runs))
    run.res <- cirGetAllParamsTarget(runs,targetpars = targetpars,chooseop = chooseop,op=2)
    list(runs,run.res)
    # res.ls <- lapply(aplrange,function(x){
    #   mcircuit <- modifyFlyMBg(circuit = circuit,aplCut = x) #set sparsity
    #   #now, run this circuit for noruns for different noise regimes for the specified components
    #   runs <- exploreNoiseSeqRuns(mcircuit,noruns=noruns,notrials = notrials,noiseset = noiseset,noiserange = noiserange,
    #                               op=op)
    #   run.res <- cirGetAllParamsTarget(runs,targetpars = targetpars,chooseop = chooseop,op=2)
    #   list(runs,run.res)
    # })
    #runs <- lapply(res.ls,'[[',1)
    #runs.res <- lapply(res.ls,'[[',2)
    # names(runs) <- aplrange
    # names(runs.res) <- aplrange
    # switch(retop,runs,runs.res,list(runs,runs.res))
  })  
  names(res.ls) <- convertDfToList(transposeDF(par.combo))
  res.ls
  #procedure: initialize circuit according to circ params, then generate all the combination of circuit params, and then
  #call them all inside the loop, initialize each time, and then  call explore noise. In the end collate all the results,
  #mainly make sure you get the names right.
  
  
}

#circparams: the circuit parameters: list(#glom,#kcs,aplkcgain,#mbons,#dans,pn-kc params,kc-mbon params,dan-mbon params,learning rule)
#circop: list(levela): each element specifies the parameter from circparams that needs to be modified. They could be atomic or a 
#nested list, e.g., the multiple nested parameters that specify connection matrices, specify the param number to be changed of the form
#levela = list(c(paramno from level1,sub-paramno),...), e.g. level2 = list(c(6,1,1)) would modify glommbconn: sub-param 1,subsubparam 1
#circmod: the amount by which the circuit parameters should be changed, could be an additive percentage: so
#list(c(left,right,step),otherparams); specifies how far to the left and right of the present value should we go, and step 
#specifies the step size. other params if needed. If only one param, apply this to all parameters, else one for every paramter.
changeCircuitParams <- function(circparams=list(glom=50,kcs=150,aplgain=6,mbons=2,dans=5,glommb=list(c(8,0.715),c(0.07),c(2,4,4)),
                                                kcmbon=c(13.5,1.15,0.16),danmbon=list(5,c(2,0.9,1.1),1),learning=c(1/38,2,0.2,T)),
                                op=1,circmod=list(10,1,1),circop=c(1)){
  #generate the param range on both sides
  allparams <- lapply(1:length(circop),function(i) {
    #now look at the corresponding parameter
    if(length(circop[[i]])==1)  parval <- circparams[[circop[[i]]]]
    else parval <- getValNestedList(circparams,circop[[i]])
    cat('\nChangecircuitparams',str(parval))
    param <- seq(parval*(1-(circmod[[1]]/100)*circmod[[3]]),
                 parval*(1+(circmod[[1]]/100)*circmod[[3]]),parval*circmod[[1]]/100)
    #cat('\n',i,param)
    param
  }) 
  #cat('\n',str(allparams))
  par.combo <- permAllGroups(allparams,retop = 2)
  colnames(par.combo) <- c()
  parlst <- lapply(1:nrow(par.combo),function(i){
    #go throug hthe param list, take each row, and then apply the values
    #using setvalNestedlist and the circop indices
    vals <- unlist(par.combo[i,]);itemlst <- circparams
    for(j in 1:length(vals)){
      #cat('\nchanging ',i,'item ',j,vals[j])
      itemlst <- setValNestedList(itemlst,circop[[j]],vals[j])
    }
    itemlst
  })
  list(par.combo,parlst)
}



#this function: explores the two main free paramters, noise and aplcut to generate a whole list of possibilities.
#then you have to navigate through all of them to arrive at the most optimal set of parameters. STore results in a multi-
#dimensional array? or single dimensional list
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
#targetpars: the target parameters against which you should be comparing results, the default are the fly params,
#chooseop:
#cirop: 
#retop: the retun type op, 1 - results as distance from target, 2 - results of the runs in terms of cells/trial&alltrials, 
#3 - both 1&2 
#classop: 1- glmmb, 2 - PNKCAPL
#verbose: F - nothing, T - keeps track of how far along you are in the simulations 
exploreNoiseParams<-function(noruns,notrials,noodors=6,noiserange=seq(0,0.3,0.1),noiseparams=list(1,c(0,0.1),1),
                             aplrange=seq(5,15,5),glomrange=c(50),kcrange=c(150),noiseset=2,outcomp=1,op=1,resop=1,
                             chooesop=2,retop=1,kcapl_conn_par=list(7,c(6,1,1.15,0.12)),aplkc_conn_par=list(8,c(6,1,1.15,0.12)),
                             targetpars=list(c(0.72,0.07),c(5.26,0.7),c(7.23,0.84),c(6.1,1.9),c(29,3.16)),thresh=list(0,c(1,0)),
                             verbose=F,classop=2){
  #create the circuit
  circuit <- createFlyCircuit(glomno = glomrange[1],kcno = kcrange[1],noodors = noodors,aplCut = aplrange[1],
                              classop = classop,kcapl_conn_par = kcapl_conn_par,aplkc_conn_par = aplkc_conn_par,thresh = thresh)
  if(verbose) cat('\nAPL range: ',aplrange)
  res.ls <- lapply(aplrange,function(x){
    if (verbose) cat('\n: APL gain: ',x)
    mcircuit <- modifyFlyMBg(circuit = circuit,aplCut = x) #set sparsity
    #now, run this circuit for noruns for different noise regimes for the specified components
    runs <- exploreNoiseSeqRuns(mcircuit,noruns=noruns,notrials = notrials,noiseset = noiseset,noiserange = noiserange,
                                verbose=verbose,op=op)
    run.res <- cirGetAllParamsTarget(runs,targetpars = targetpars,chooseop = chooseop,op=2)
    list(runs,run.res)
  })
  runs <- lapply(res.ls,'[[',1)
  runs.res <- lapply(res.ls,'[[',2)
  names(runs) <- aplrange
  names(runs.res) <- aplrange
  switch(retop,runs,runs.res,list(runs,runs.res))
}


#dat.ls: is the data in getCampDataSigList format
#processes the results of explorenoiseseqruns in a format conducive to plotting
#op=1, a list of two vectors: mean and sd , 2 = dataframe of whichever results are needed
#seqop: the options for the processnoiseseqresults function. 
#seqop=1, get ratio of rel to unrel cells per trial, =2, no of rel cells/trial
#seqop=3, no of unrel cells/trial, 4 - # rel cells/odor, 5 - # unrel/odor
#op=1, whatever option is specified by seqop, 100 - combination of those that return a single value (mean)
#op=101 combination of those that return a single value (SD), 102: get both mean and SD
processNoiseSeqRunsResults <- function(dat.ls,seqop=1,op=1){
  if(op==100 || op==101 || op== 102){#multiple options
      res1 <- processNoiseSeqRunsResults(dat.ls,seqop = 1,op=1)
      res2 <- processNoiseSeqRunsResults(dat.ls,seqop = 2,op=1)
      res3 <- processNoiseSeqRunsResults(dat.ls,seqop = 3,op=1)
      res4 <- processNoiseSeqRunsResults(dat.ls,seqop = 4,op=1)
      res5 <- processNoiseSeqRunsResults(dat.ls,seqop = 5,op=1)
      resm.ls <- lapply(list(res1,res2,res3,res4,res5), function(x) {
        switch(op-99,apply(x,1,mean),apply(x,1,sd)/sqrt(ncol(x)) )
      })
      if (op < 102) {
        res <- transposeDF(convertNestedListsDF(resm.ls) )
        rownames(res) <- c('ratio','rel/tr','unrel/tr','rel/od','unrel/od')
        colnames(res) <- rownames(res1)
      } else {#102: 
        #cat('\nlist',str(list(res1,res2,res3,res4,res5)))
        res.mn <- lapply(list(res1,res2,res3,res4,res5), function(x) apply(x,1,mean) )
        res.sd <- lapply(list(res1,res2,res3,res4,res5), function(x) apply(x,1,sd)/sqrt(ncol(x)) )
        #print(res.sd)
        res.mn.df <- transposeDF(convertNestedListsDF(res.mn) )
        rownames(res.mn.df) <- c('ratio','rel/tr','unrel/tr','rel/od','unrel/od')
        colnames(res.mn.df) <- rownames(res1)
        res.sd.df <- transposeDF(convertNestedListsDF(res.sd) )
        rownames(res.sd.df) <- c('ratio','rel/tr','unrel/tr','rel/od','unrel/od')
        colnames(res.sd.df) <- rownames(res1)
        res <- list(res.mn.df,res.sd.df)
      }
      return(res )
  }
  #make a matrix or datafrrame 
  res.ls <- lapply(dat.ls,function(x){
    res <- processNoiseSeqResults(x,op = seqop)
  })
  res <- convertNestedListsDF(res.ls)
  #cat('\nprocessNoiseSeqRunsResults',str(res.ls))
  #print(res)
  res
}


#function that looks through the results of processNoiseSeqRunsResults and then tells you if which params 
#get you closest to your results
#parno: the param no,i.e., the row to look at, if 0 do all params
#targetpars: the target params a list of two vectors (means,sem) or just means in which case sems <- 10% of means
#chooseop: the way of choosing fit. 1 - all overlaps, 2 - min overlao, 3 - max overlap, 
#op: 1 - return value of the difference between them
cirGetParamsTarget <- function(dat.ls,parno=1,targetpars,chooseop=1,op=3){
  allpars <- processNoiseSeqRunsResults(dat.ls = dat.ls,op = 102)
  if(length(targetpars)==1) targetsems <- sapply(targetpars[[1]], function(x) x*0.2) #sem is 10 % of signal
  else targetsems <- targetpars[[2]]
  if(parno==0) parnos <- 1:nrow(allpars[[1]])
  else parnos <- parno
  close <- sapply(1:ncol(allpars[[1]]),function(i){
    #get overlap for all params. If they are all 0, this is 0
    alloverlap <- sapply(1:length(parnos), function(j){
      overlap <- intervalOverlap(c(targetpars[[1]][[j]],targetsems[j]),c(allpars[[1]][j,i],allpars[[2]][j,i]),inputop = 2)
    })
    if(prodListVectors(alloverlap)==0) res.over <- 0
    else res.over <- lpnorm(alloverlap)
  })
  if(sum(close)==0) return(switch(op,0,0,c(0,0)) ) #no bites, so return nothing
  names(close) <- names(allpars[[1]]) #assign names
  #get the appropropriate posns according to op
  posns <- switch(op,which(close>0),min(close[which(close>0)]),max(close[which(close>0)]))
  posns <- switch(op,posns,which(posns==close),which(posns==close))
  res <- close[posns]
}

#parno: the param no,i.e., the row to look at, if 0 do all params
#targetpars: the target params a list of two vectors (means,sem) or just means in which case sems <- 10% of means
#chooseop: the way of choosing fit. 1 - all overlaps, 2 - min overlao, 3 - max overlap, 
#op: 1 - return value of the difference between them
cirGetParamsTargetAllApls <- function(dat.ls,parno=1,targetpars,chooseop=1,op=3){
  apl.grps <- 1:length(dat.ls)
  apl.names <- names(dat.ls)
  res.ls <- lapply(apl.grps, function(i){
    res <- cirGetParamsTarget(dat.ls = dat.ls[[i]],parno = parno,targetpars = targetpars,op = op)
    names(res) <- paste(apl.names[i],',',names(res),sep = '')
    res
  })
  res.ls  
}

#function that looks through the results of processNoiseSeqRunsResults and then tells you if which params 
#get you closest to your results
#dat.ls:
#targetpars: all target params in the same order as parno which is the number of rows
#def.width: the default width for determining the extent of the interval on both sides of the parameter mean
#targetpars: the target params c(mean,sem)
#chooseop: the way of choosing fit.  1 - min overlao, 2 - max overlap, 3 - difference (observed-target), 
#4 -norm ratio (observed-target)/target 
#op: relative overlap, 2 - relative difference
cirGetAllParamsTarget <- function(dat.ls,targetpars,chooseop=1,def.width=0.15,op=2){
  allpars <- processNoiseSeqRunsResults(dat.ls = dat.ls,op = 102) #get the number of params
  #cat('\ncirGetAllParamsTarget',str(allpars)) 
  res.ls <- lapply(1:ncol(allpars[[1]]),function(i){#iterate through parameter combination
    parscore <- sapply(1:nrow(allpars[[1]]),function(j){#through each of the constraints
      overlap <- intervalOverlap(targetpars[[j]],c(allpars[[1]][j,i],allpars[[2]][j,i]),inputop = 2,op = op)
      if(overlap==0) overlap<-10
      overlap
    })
    score <- sum(parscore)
  })
  names(res.ls) <- names(allpars[[1]])
  res <- unlist(res.ls)
  res
}

#function that analyzes the results of the good params. The return of the results of cirGetAllparamstarget
#dat.ls: the output of cirGetAllPars, of the form, list(c(parms1,params2...),c()..). each vector correponds
#to an APLcutoff number. The values are the euclidean distance from the desired results
#noiseop: specifies the noise option. like 22 is apl, apl-kc, pn-mb, pn
#par: the param you want to analyze: 1 - apl in op=22., if you say 0, it is apl cut-off
#op: 1 - return the analyzed resutls. list(avg. vals,params freq,apl.freq) 
#2 - list(avg. vals,params freq,apl.freq,sorted par vals)
cirAnalGoodParams <- function(dat.ls,noiseop=22,par=1,op=1){
  dat.vec <- flattenLists(dat.ls)
  nonzeroes <- which(dat.vec > 0)
  #get all the param values as cols in a df
  cat('\n',names(dat.vec),nonzeroes)
  parvals.df <- convertStringNos2Df(names(dat.vec)[nonzeroes],op=1) 
  parvals.df <- cbind.data.frame(parvals.df,dat.vec[nonzeroes])
  #cirassign gives you the labels of the parameters that were changed according to noiseop
  cnames <- c('params',cirAssignNoiseLabels(noiseop = noiseop,op=1),'dist')
  colnames(parvals.df) <- cnames
  totalno <- nrow(parvals.df) #total no of successful params
  #sort DF by noise in APL
  parvals.sort.df <- getMatSortPosns(parvals.df,col = 3,op=3)
  
  #cat(str(parvals.df),str(dat.vec))
  #average value of the params; omit the first col
  avg.val <- apply(parvals.sort.df[,-1],2,mean)
  avg.sem <- apply(parvals.sort.df[,-1],2,mean)/sqrt(totalno)
  
  #no of parameters and their frequency of having no noise
  par.freq <- sapply(2:ncol(parvals.df), function(i) length(which(parvals.sort.df[,i]==0)) )
  par.freq <- round((totalno-par.freq) * (100/totalno),0)
  par.freq.df <- rbind.data.frame(as.vector(par.freq))
  names(par.freq.df) <- colnames(parvals.sort.df)[2:ncol(parvals.df)]
  
  #frequencies of the different noise values for APL
  apl.freq <- table(parvals.sort.df[,par+2])
  apl.freq.df <- rbind.data.frame(as.vector(apl.freq))
  names(apl.freq.df) <- names(apl.freq)
  
  #the frequencies as a plot
  freq.plt <- NlsFitCDF(parvals.df[,3],graphparams = list(xlabel='apl noise',ylabel='freq.',ticknoy=2), dist = 1)
  cat('\nstatistics',freq.plt) #print out the results
  
  res <- list(avg.val,par.freq.df,apl.freq.df,totalno,parvals.sort.df,freq.plt)
  names(res) <- c('avg. par vals','params freq','apl par freq','#success','pars DF','freq.plot')
  switch(op,res[-(length(res)-1)],res)
}

#function to see how closely the dataset matches the observed results
cirEvalSimResp <- function(dat.ls,meas.sim=c(0.15,0.35),op=1){
  #do the response stats first
  resp.stats <- campGetReliabStats(dat.ls)
  resp.stats <- c(resp.stats[[1]],resp.stats[[2]][1],resp.stats[[3]][1:2])
  
  ud.stats <- campComputeUD(dat.ls,meas.sim = meas.sim)
  ud.chg <- rbind(ud.stats$dissim.sat/ud.stats$dissim.st,ud.stats$sim.sat/ud.stats$sim.st)
  #example_param_response_gamma_fit
  responses <- getAboveThresh(unlist(lapply(dat.ls,'[[',2)) )
  
  #get responses and ud characteristics
  resp.ud <- c(resp.stats[[5]][1]/resp.stats[[5]][2],resp.stats[[6]],resp.stats[[7]],ud.chg[,3])
  names(resp.ud) <- c('rel./unrel.','rel.','unrel','D. dissim.','D sim.')
    
  list(resp.stats,ud.chg,responses,resp.ud)
}


#finds the params that satisfy a certain condition
#op: the kind of contraint, 1 - minimum of const for all the pars summed
#2- none of the params are more than const
#retop: what to return. 1 - the posns where this occurs, 2 - the rows where this occurs
cirFindParamsConstraints <- function(dat.df,const=c(),pars = 1,retop=2,op=1){
  if(op==1){
    res <- apply(dat.df[,pars], 1, sum)
    posns <- which(res <= const)
  }
  if(op==2){
    posns.lst <- lapply(pars, function(i) which(dat.df[,i] <= const) )
    posns <- Reduce(intersect,posns.lst)
  }
  switch(retop,posns,dat.df[posns,],length(posns))
}

#noiseop: specifies the noise option. like 22 is apl, apl-kc, pn-mb, pn
#op: 1 - add apl-cut-off at the at the start, 2 - don't add
cirAssignNoiseLabels <- function(noiseop,op=1){
  if(noiseop==22) res <- c('apl', 'apl-kc', 'pn-mb', 'pn')
  if(noiseop==30) res <- c('apl', 'apl-kc', 'pn-mb', 'pn','kc')
  switch(op,c('apl cutoff',res),res)
}


#op: 1 - range of values from low to high noise, 2 - optimum parameter
cirOutputFormPars <- function(dat.ls,targetpars,choosop=1,def.width=0.15,retop=1,op){
  low.noise <- cirGetAllParamsTarget(dat.ls = dat.ls,targetpars = targetpars,choosop = 1,retop = retop)
  high.noise <- cirGetAllParamsTarget(dat.ls = dat.ls,targetpars = targetpars,choosop = 2,retop = retop)
  res.ls <- lapply(1:ncol(low.noise),function(i){
    paste(round(low.noise[,i],2),round(high.noise[,i],2),sep = '-')
  })
  res.df <- convertNestedListsDF(res.ls)
  rownames(res.df) <- rownames(low.noise)
  colnames(res.df) <- colnames(low.noise)
  switch(op,res.df,xtable::xtable(res.df) )
}

#converts the noise notation to readable format, e.g.,'1.0.2' to '20'
cirConvNoiseStr <- function(nstr,op=1){
  tmp <- strsplit(nstr,'\\.') #need to escape .
  #cat('\n',as.numeric(tmp[[1]][3]),length(tmp[[1]]),':',tmp[[1]] )
  #the second part holds the 100th place, and the decimal places hold percentages less than 100.
  if(length(tmp[[1]])>2) decpart <- as.numeric(tmp[[1]][3])/10^(1+getpow10(as.numeric(tmp[[1]][3])) )
  else decpart <- 0
  as.numeric(tmp[[1]][2])*100 + decpart*100
}


#function that explores the noise results for a range of values of noise for a certain noise parameter(s) as
#specified by op
#the output gives you the average number of reliablae and unreliable cells for both per odor average and per trials average
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
#notrials: the number of trials, each trial generate noise anew from the original matrix
#resop: option for the how the results output should be formatted, 1 = mbkc_net firing rates, 2 - mbkc, 3 - glomerular 
#4 - list(glomerular,mbkcnet) , 10 - all the circuits
#noiseset: specifies the test or training test to be tested for noise
#noisparams: the noise parameters. If only one element, the same noise params for everything. else different noise params
#for each component
#noiserange: noise range of the components: if multiple components can be a list of noise ranges for each component
#if only one element for multiple component, this will duplicate that for the other components
#self is the original circuit
#outcomp: the outputcomponent we want to see, 1 - glom, 2 - mbkc, 3 - mbkcnet, 4 - mbons
##verbose: F - nothing, T - keeps track of how far along you are in the simulations 
exploreNoiseSeq <- function(self,notrials,noiserange=list(seq(0,0.3,0.1)),noiseparams=list(1,c(0,0.1),1),
                            noiseset=1,outcomp=1,verbose=F,op=1,resop=4){
  #load the noise range and parameters according to the number of components to change
  if((op>=10 && op<=11) || (op>=14 && op<30) || (op==30) || (op==100)) {#determine #components and generate the parameter combination
    nocomp <- computeNoNoiseComp(op) #get number of components 
    #cat('\nno of components',nocomp)
    if(isDataType(noiserange)==1) noiseran <- repList(noiserange,nocomp)
    else if(isDataType(noiserange)==4 && length(noiserange)>1) noiseran <- noiserange
         else noiseran <- repList(noiserange[[1]],nocomp)
    if(identical(list(1,c(0,0.1),1),noiseparams)) noisepar <- repList(noiseparams,nocomp)
    else noisepar <- noiseparams
    #cat('\nnocomp',nocomp,'len noiseran',str(noiseran),str(noisepar),'fi')
  } 
  if(op<10) {#not a multi-component run
    if(isDataType(noiserange)==1) noiseran <- list(noiserange)
    else noiseran <- noiserange   
    noisepar <- noiseparams
  }
  #cat('\nrange',noiserange,str(noiseran),'\n')
  noisevals <- permAllGroups(noiseran,retop = 2) #generating the parameter combinations
  if(verbose) cat('\t',nrow(noisevals),' simulations: ')
  #for the names
  noisevalstr <- sapply(1:nrow(noisevals),function(i) do.call(paste,c(as.list(as.character(noisevals[i,])),sep=',') ) )
  #for each paramter combo, generates the noise params (testparams) and passes them to exploreNoise
  res.ls <- lapply(1:nrow(noisevals),function(i){#go through all the noise parameter combindation 
    if(verbose && (i %% 400 ==0)) cat('\t',i,' runs done.')
    if(op>=10){#different combinations of parameters
      #apply the noise parameter to each component  
      testparams <- lapply(1:ncol(noisevals[i,]),function(j) {
        testpars <- noisepar[[j]] #for the jth component, assign the jth noiseval
        #cat('\nn2',i,str(noisevals[i,j]))
        testpars[[2]][2] <- noisevals[i,j]
        testpars
      })
    } else {#single components
      testparams <- noiseparams
      testparams[[2]][2] <- noisevals[i,1] #replace the noise parameter with the noise we are testing
    }
    res.noise <- exploreNoise(self,notrials = notrials,noiseparams = testparams,noiseset = noiseset,resop = resop,op=op)
    reskc.noise <- res.noise[[1]];respn.noise <- res.noise[[2]]
    #get overall averages
    odor.avg <- apply(getCampDataRel(reskc.noise,op=2),2,mean)
    #averages by trial
    trial.res <- campClassifyTrialCells(reskc.noise,trials = notrials) #lian, check this
    trial.cellclass <- joinListDFs(trial.res[[4]])
    trial.avg <- apply(trial.cellclass,2,mean)
    #get the correlations in similarity between pns and kcs
    odor.cor <- computeMatCor(computeListMatCorr(respn.noise),computeListMatCorr(reskc.noise),op=1)
    odor.cor.fit <- computeMatCor(computeListMatCorr(respn.noise),computeListMatCorr(reskc.noise),op=2)
    res <- list(odor.avg,trial.avg,c(odor.cor,odor.cor.fit),respn.noise,reskc.noise)
    names(res) <- c('per odor','per trial','pn-kc similarity','pn','kc')
    res
  })
  names(res.ls) <- noisevalstr
  res.ls
}

#gives the number of noise components specified by noiseop
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
computeNoNoiseComp <- function(noiseop=1,op=1){
  res <- 0 #everything else, no noise
  if(noiseop<10 && noiseop>0) return(1)
  if(noiseop>=10 && noiseop<30) res <- switch(noiseop-9,2,3,2,10,#10,11,12,13
                                2,2,2,2, #14,15,16,17
                                2,2,3,3,4) #18,19,20,21,22
  if(noiseop==30) res <- 5
  if(noiseop==100) res <- 10
  res
}

#function that explores how the KC output similarity changes with noise in the circuit

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

#processes the results of explorenoiseseq in a format conducive to plotting
#op=1, get ratio of rel to unrel cells per trial, =2, no of rel cells/trial
#op=3, no of unrel cells/trial, 4 - # rel cells/odor, 5 - # unrel/odor
processNoiseSeqResults <- function(res.ls,op=1){
  rel <- unlist(lapply(lapply(res.ls,'[[',2),'[[',1))
  unrel <- unlist(lapply(lapply(res.ls,'[[',2),'[[',2))
  odor.unrel <- unlist(lapply(lapply(res.ls,'[[',1),'[[',3))
  odor.rel <- unlist(lapply(lapply(res.ls,'[[',1),'[[',2))
  res <- (rel/unrel)
  res[res==Inf] <- 0
  switch(op,res,rel,unrel,odor.rel,odor.unrel)
}



#generic method for exploring noise in a circuit
exploreNoise <-function(self,...){
  UseMethod("exploreNoise",self)
}

#new todo: have to add noise to just APL and the combinations of different thigns and APL.
  
#does noise exploration on the fly circuit given by FlyMBg across multiple trials
#noiseset: the set that we would like at explore 0: all odors, 1 - orig set, 2 - training set, 3 - test set
##noiseparams: of the form list(dist,distparams,type of noise:multiplicative or addiditev)
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
#notrials: the number of trials, each trial generate noise anew from the original matrix
#resop: option for the how the results output should be formatted, 1 = mbkc_net firing rates, 2 - mbkc, 3 - glomerular 
#4 - glomerular and mbkc firing rates, 10 - all the circuits
#outcomp: the outputcomponent we want to see, 1 - glom, 2 - mbkc, 3 - mbkcnet, 4 - mbons
exploreNoise.FlyMBg <-function(self,notrials,noiseparams=list(1,c(0,0.1),1),noiseset=1,outcomp=1,op=1,resop=1){
  #add noise to the requisite elements
  #cat('\nexplorenoiseFly',str(noiseparams),'op ',op)
  circuit <- setNoiseParams(self,noiseparams = noiseparams,op=op) 
  odornos <- chooseStimuli(circuit,stimset = noiseset)
  #cat('\n APL firing in explorenoise before trials',unlist(circuit$apl$neurons[1:2]))
  if(notrials>1) {
    res <- lapply(1:(notrials-1),function(i){#circuit is the first trial, now do notirals -1 more trials
      #now compute the firing rate for all, for every trial. And noise will be generated afresh
      res.circuit <- computeFiring(circuit,op=1,odor=odornos,noiseop=op) #noise is added through noiseop
      #cat('\n',i,' APL firing in explorenoise newcircuit',unlist(res.circuit$apl$neurons[7:8]))
      res.circuit
    })
    res <- c(list(computeFiring(circuit,op=1,odor=odornos,noiseop=op)),res) #the first trial is no noise
  } 
  else res <- c(list(computeFiring(circuit,op=1,odor=odornos,noiseop=op))) #only one trial
  switch(resop,processNoiseResults(res,stim = odornos,comp=3),#1
         processNoiseResults(res,stim = odornos,comp=2),#2
         processNoiseResults(res,stim = odornos,comp=1),#3
         list(processNoiseResults(res,stim = odornos,comp=3),processNoiseResults(res,stim = odornos,comp=1)),#4
         res)#5
}

#sets the noise parameters for the required components
#noiseparams: of the form list(dist,distparams,type of noise:multiplicative or addiditev)
#if noise has to be added to multiple components, it has to be a list of lists
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
setNoiseParams <- function(self,noiseparams=list(1,c(0,0.1),1),op=1){
  circuit <- self
  #choose the components to which we will be adding noise
  compids <- setNoiseCompIds(op)
  #now based on what is returned, add noise
  #if only a single number is returned
  if(length(compids)==1) {
    if(compids==0) return(circuit) #no noise to add
    circuit <- setNoiseComp(compid = compids,circuit = circuit,noiseparams = noiseparams)
    return(circuit)
  }
  else {
    if(op < 100) {#multiple compids are chosen, learning and nonlearning structures are separate options
      if(isNestedList(noiseparams)==F && length(noiseparams)==3) noiseparams <- repList(noiseparams,length(compids))
      #cat('\nsetnoise',length(noiseparams))
      for(i in seq_wrap(1,length(compids))){
        circuit <- setNoiseComp(compid = compids[i],circuit = circuit,noiseparams = noiseparams[[i]])
      }
    } else {#op=100, so set noise to everything
      #non-learning structures
      for(i in seq_wrap(1,length(compids[[1]]))){
        circuit <- setNoiseComp(compid = compids[[1]][i],circuit = circuit,noiseparams = noiseparams)
      }
      #now, set the valence related (learning) bits
      #got through each component, then their valences, and set noise
      for(j in seq_wrap(1,length(compids[[2]]))){
        circuit <- setNoiseComp(compid = compids[[2]][j],circuit = circuit,noiseparams = noiseparams)
      }    
    }
  }
  circuit
}

#give the noise option, sets the component ids
#noiseop: components to noise to; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
setNoiseCompIds <- function(op,resop=1){
  #logic go through all the options and return the required result
  if((op>=1 && op<=4) || (op==8) || (op==9)){#only one component or structure for noise addition
    if(op>=1 && op<=4) compid <- switch(op,const.FM$glom,const.FM$glommb_conn,const.FM$kcapl_conn,const.FM$aplkc_conn,
                                        const.FM$dans,const.FM$kcmbon_conn)
    if(op>=8 && op<=9) compid <- switch(op-7,const.FM$apl,const.FM$mbkc) 
    return(compid)
  }
  if(op>=5 && op<=7){#noise to structures that are valence specific
    #now, set the valence related bits
    compid <- switch((op-4),const.FM$kcmbon_conn,const.FM$mbons,const.FM$dans)
    return(compid)
  }
  if((op>=10 && op <= 11) || (op>=14 && op <= 30)) {#multiple structures are being changed
    if(op>=10 && op <= 11)  compids <- switch((op-9),c(const.FM$glom,const.FM$glommb_conn),
                                              c(const.FM$glom,const.FM$glommb_conn,const.FM$aplkc_conn))
    #apl and mbkc structures: 2 structures
    if(op>=14 && op <= 19)  compids <- switch((op-13),c(const.FM$apl,const.FM$aplkc_conn),#14
                                              c(const.FM$apl,const.FM$glommb_conn),c(const.FM$apl,const.FM$glom),#15,16
                                              c(const.FM$apl,const.FM$mbkc),c(const.FM$mbkc,const.FM$glommb_conn),#17,18
                                              c(const.FM$mbkc,const.FM$aplkc_conn))
    #3, 4, and 5 structures
    if(op>=20)  compids <- switch((op-19),c(const.FM$apl,const.FM$aplkc_conn,const.FM$glommb_conn), #20
                                  c(const.FM$apl,const.FM$aplkc_conn,const.FM$mbkc),#21
                                  c(const.FM$apl,const.FM$aplkc_conn,const.FM$glommb_conn,const.FM$glom))#22
    if(op==30)  compids <- c(const.FM$apl,const.FM$aplkc_conn,const.FM$glommb_conn,const.FM$glom,const.FM$mbkc) #30
    return(compids)
  }
  if(op==12) {#the DANs or valence network
    compids <- c(const.FM$kcmbon_conn,const.FM$mbons,const.FM$dans)
    return(compids)
  }
  if(op==100){#set noise to all structures
    #set the non-learning structures
    compids1 <- c(const.FM$glom,const.FM$glommb_conn,const.FM$aplkc_conn,const.FM$apl,const.FM$kcapl_conn,
                 const.FM$mbkc,const.FM$mbkc_net)
    #now, set the valence related bits
    compids2 <- c(const.FM$kcmbon_conn,const.FM$mbons,const.FM$dans)
    return(list(compids1,compids2))
  }
  return(0) #if we are here, none of the options worked, no nosie
}



#set noise params for the given component in the circuit
#op=1, neurons, 2 - connection matrix
#compid: id of the component in the 'circuit'
setNoiseComp <- function(compid,circuit,noiseparams,op=1){
  comp <- getData(circuit,compid) #get the component
  comptype <- getCompType(comp,op=2) #get component type so we can use the appropriate function params
  #the comp types: 1 - neuron, 2 - connection matrix, 3 - Circuit, 4 - learning matrix, 5 - list of neurons, 8 - list of learning matrices
  #cat('\nsetnoisetype',comptype,'compid',compid)
  ##neurons: set params here and basenoise=T with addnosie, then add noise within computefiring.neurons
  if(comptype==1) comp <- setData(comp,val=noiseparams,const.N$noise) 
  if(comptype==2) comp <- setData(comp,val=noiseparams,const.CM$noise)
  if(comptype==5 || comptype==8){#for the learning components
    if(comptype==5) noisetype <- const.N$noise
    else noisetype <- const.CM$noise
    #cat('\nstnoisecomp',length(comp),comptype)
    for(i in 1:length(comp)){
      comp[[i]] <- setData(comp[[i]],val=noiseparams,noisetype)
    }
  } 
  circuit <- setData(circuit,val=comp,compid) #set the component again
}


#processes the results in the format you would like based on op
#component results that we want: 0 - all, 1 - glom, 2 - mbkc, 3 - mbkcnet, 4 - mbons
#stim: specifies the odornos that should be retrieved
#op=1, in the form of list of odors
processNoiseResults <-function(results,stim,comp=1,op=1){
  notrials <- length(results)
  opid <- comp
  #cat('\npricess',opid,notrials)
  compid <- switch(opid,const.FM$glom,const.FM$mbkc,const.FM$mbkc_net,const.FM$mbons)
  notrials <- length(results)
  res <- lapply(1:notrials,function(i){
    flycomp <- getData(results[[i]],compid)
    neurons <- flycomp$neurons
    #cat('\nProcessNoise, ',i,', ',compid,class(results[[i]]),class(flycomp))
    #for each response get the binary responses and actual rate responses
    responses <- lapply(neurons[stim],function(x) {
      tmp <- list(setVecThresh(x,op=2),x)
      #print(tmp)
      names(tmp) <- c('bin','sig') #label the response type
      tmp
    }) 
    names(responses) <- stim #label the stimulus no.
    responses
  })
  res <- flattenLists(res,level = 1) #we want them like t1 od1, t1 od2, t2 od1,...
  res
}


#function takes in the circuit and then looks at all the correlations
analyzeReliability <- function(circuit,topn=25,odors=c(),op=1){
  #get the training data df
  train.df <- getData(circuit,const.FM$trainingval)
  mix <- as.character(train.df[,3])
  #cat('\nAR',str(mix))
  mixu.strings <- sapply(1:nrow(train.df),function(i){
    mixstring <- strsplit(mix[i],split = ',')
    mix.num <- as.numeric(mixstring[[1]])
    as.numeric(c(mix.num[1],train.df[i,2],mix.num[3]) )
  })

  mixl.strings <- sapply(1:nrow(train.df),function(i){
    mixstring <- strsplit(mix[i],split = ',')
    mix.num <- as.numeric(mixstring[[1]])
    as.numeric(c(mix.num[2],train.df[i,2],mix.num[4]) )
  })
  # cat('\nmix',ncol(mixl.strings),ncol(mixu.strings),'\n' )
  # print(table(mixu.strings[3,]))
  # print(table(mixl.strings[3,]))
  mix.strings <- t(cbind(mixl.strings,mixu.strings) )
  mix.cor <- sapply(1:nrow(mix.strings),function(i) {
    tmp.res <- analyzeRelUnRelCor(circuit = circuit,topn = topn,odors = c(mix.strings[i,1],mix.strings[i,2])) 
    c(mix.strings[i,],tmp.res[[1]][c(1,4,7)],tmp.res[[2]][c(1,4,7,8:11,12:15)]) 
  })
  mix.cor <- t(mix.cor)
  rownames(mix.cor) <- mix.strings[,c(2)]
  #cat('\nncol',ncol(mix.cor))
  colnames(mix.cor) <- c('od1','od2','overlap','ob.cor','ob.cor.rat','ob.diff.rat','mb.cor','mb.cor.rat','mb.diff.rat',
                         'mb.diff.top','mb.diff.bot','mb.sum.top','mb.sum.bot','m1.t','m2.t','m1.b','m2.b')
  mix.cor
}

#peocess the results of analyzereliability function
#op=1
processReliabilityData <- function(circuit,topn=25,op=1){
  reldata <- analyzeReliability(circuit,topn = topn)
  #res <- lapply(c(0.2,0.3,0.4), function(i) {
  res <- sapply(seq(1,9)/10, function(i) {
    #cat('\n',i,which(reldata[,3]==i) )
    #print(reldata[which(reldata[,3]==i),])
    apply(reldata[which(reldata[,3]==i),],2,mean) 
  })
  res <- t(res)[,c(3,9,10:17)]
  res
  #table(reldata[,3])
}


#this function deteccts the correlation between pairs of AL and MB responses, then checks to see if the 
#the highest firing neurons in MB carry more of the inforamtion than the lowest firing neurons
analyzeRelUnRelCor <- function(circuit,topn=25,odors=c(),op=1){
  #first get the al and mb responses
  glom.resp <- getData(getData(circuit,const.FM$glom),const.N$neurons)
  mbkc.resp <- getData(getData(circuit,const.FM$mbkc),const.N$neurons)
  #now get the required odor pairs given by odors
  glom1 <- glom.resp[[odors[1]]]
  glom2 <- glom.resp[[odors[2]]]
  mbkc1 <- mbkc.resp[[odors[1]]]
  mbkc2 <- mbkc.resp[[odors[2]]]
  #get the correlations
  glom.cor <- cor(glom1,glom2)
  mbkc.cor <- cor(mbkc1,mbkc2)
  #get correlations of top25 and bottom25
  glompos1 <- union(getTopNPercentile(glom1,topn = topn,op=3),getTopNPercentile(glom2,topn = topn,op=3))
  mbkcpos1 <- union(getTopNPercentile(mbkc1,topn = topn,op=3),getTopNPercentile(mbkc2,topn = topn,op=3))
  glomtop.cor <- cor(glom1[glompos1],glom2[glompos1])
  mbkctop.cor <- cor(mbkc1[mbkcpos1],mbkc2[mbkcpos1])
  #difference too: 
  gdiff1 <- mean(abs(glom1[glompos1]-glom2[glompos1])/(glom1[glompos1]+glom2[glompos1]) )
  mdiff1 <- mean(abs(mbkc1[mbkcpos1]-mbkc2[mbkcpos1])/(mbkc1[mbkcpos1]+mbkc2[mbkcpos1]) )
  mdiff1.abs <- mean(abs(mbkc1[mbkcpos1]-mbkc2[mbkcpos1]) )
  msum1.abs <- mean((mbkc1[mbkcpos1]+mbkc2[mbkcpos1]) )
  mt1 <- c(mean(abs(mbkc1[mbkcpos1])),mean(abs(mbkc2[mbkcpos1])) )
  #bottom posn
  glompos1 <- union(getBottomNPercentile(glom1,bottomn = topn,op=3),getBottomNPercentile(glom2,bottomn = topn,op=3))
  mbkcpos1 <- union(getBottomNPercentile(mbkc1,bottomn = topn,op=3),getBottomNPercentile(mbkc2,bottomn = topn,op=3))
  glombot.cor <- cor(glom1[glompos1],glom2[glompos1])
  mbkcbot.cor <- cor(mbkc1[mbkcpos1],mbkc2[mbkcpos1])
  #difference
  gdiff2 <- mean(abs(glom1[glompos1]-glom2[glompos1])/(glom1[glompos1]+glom2[glompos1]) )
  mdiff2 <- mean(abs(mbkc1[mbkcpos1]-mbkc2[mbkcpos1])/(mbkc1[mbkcpos1]+mbkc2[mbkcpos1]) )
  mdiff2.abs <- mean(abs(mbkc1[mbkcpos1]-mbkc2[mbkcpos1]) )
  msum2.abs <- mean((mbkc1[mbkcpos1]+mbkc2[mbkcpos1]) )
  mt2 <- c(mean(abs(mbkc1[mbkcpos1])),mean(abs(mbkc2[mbkcpos1])) )
  #cat('\n',mdiff,'\n',gdiff)
  resg <- c(glom.cor,glomtop.cor,glombot.cor,glomtop.cor/glombot.cor,gdiff1,gdiff2,gdiff1/gdiff2)
  resmb <- c(mbkc.cor,mbkctop.cor,mbkcbot.cor,mbkctop.cor/mbkcbot.cor,mdiff1,mdiff2,mdiff1/mdiff2,mdiff1.abs,mdiff2.abs,
             msum1.abs,msum2.abs,mt1,mt2)
  list(resg,resmb)
}

#given two vectors, compares the overlap between them for different deciles
#vecs: list of vectors
#deciles: the deciles for which we need to do the calculation
getVecOverlapDeciles <- function(vecs,deciles=seq(25,100,25),decsize=25,op=1){
  #cat('\n',str(vecs))
  #get the overlap for each of the deciles
  res.ls <- lapply(deciles,function(x){
    cell.lst <- lapply(vecs,function(y){
      #cat('\ntopdec',x-1,x-1+decsize)
      topdec <- getTopNDecile(y,topn = x,decile = x,zero = 2,op = 3)
      #topdec <- setdiff(getTopNPercentile(y,topn = x,zero = 2,op = 3),getTopNPercentile(y,topn = x-decsize,zero = 2,op = 3) )
      #cat('\n',x,getTopNDecile(y,topn = x,decile = decsize,zero = 2,op = 3),':',topdec)
      topdec
    })
    common <- Reduce(intersect,cell.lst)
  })
  sapply(res.ls,length)
}

#this function looks at the different stimuli and determines their correlation and then the overlap across
#the different deciles
#type: circuit component type, 1 = mbkc_net, 2 = glomeruli, 3 = mbkc
getStimListOverlaps <- function(circuit,deciles=seq(25,100,25),decsize=25,type=1,op=1){
  #get the training or test stimuli odors, and their supposed overlap
  type.df <- extractStimCor(circuit = circuit,stimtype = 1,type = 1)
  #now, extract the responses of the appropriate circuit component
  glom.lst <- getData(getData(circuit,const.FM$glom),const.N$neurons) #get the glomeruli just to make sure
  resp.lst <- switch(type,getData(circuit,const.FM$mbkc_net),getData(circuit,const.FM$glom),getData(circuit,const.FM$mbkc)) 
  resp.lst <- getData(resp.lst,const.N$neurons)
  #get the training stimulus responses
  stim.lst <- resp.lst[as.numeric(type.df[,1])]
  names(stim.lst) <- type.df[,1]
  #get the original stimulus responses
  orig.lst <- resp.lst[as.numeric(unique(type.df[,2]))]
  names(orig.lst) <- unique(type.df[,2])
  #now get the correlations
  cor.vec <- sapply(1:nrow(type.df),function(i) {
    #cat('\n',type.df[i,1],type.df[i,2],':')#,str(stim.lst[[as.numeric(type.df[i,1]) ]]),str(orig.lst[[as.numeric(type.df[i,2]) ]]))
    stim1 <- stim.lst[[as.character(type.df[i,1]) ]]
    stim2 <- orig.lst[[as.character(type.df[i,2]) ]]
    res.cor <- cor(stim1,stim2) 
    res.over <- getVecOverlapDeciles(list(stim1,stim2))
    #cat('\n',str(as.numeric(type.df[i,1:2])))
    #glomerular correlation to make sure correlations are kosher
    glom.cor <- cor(glom.lst[[as.numeric(type.df[i,1])]],glom.lst[[as.numeric(type.df[i,2])]])
    c(res.cor,glom.cor,res.over)
  })
  cor.vec <- t(cor.vec)
  type.df <- cbind(type.df,cor.vec)
  colnames(type.df)[4:9] <- c('mb.cor','pn.cor','25','50','75','100') 
  type.df
}

#given a circuit, extracts the correlations between the original stimuli and the training or testing stimuli
#circuit: the circuit
#stimtype: 1 - training, 2 - test
#type: circuit component type, 1 = mbkc_net, 2 = glomeruli, 3 = mbkc
extractStimCor <- function(circuit,stimtype=1,type=1,op=1){
  #get the odors in training stimuli
  stim <- switch(stimtype,getData(circuit,const.FM$trainingval),getData(circuit,const.FM$teststim) )
  type.lst <- strsplit(as.character(stim[,3]),split = ',')  
  #cat('\n',str(train.lst),str(lapply(train.lst,'[[',1)) )
  type.df <- cbind.data.frame(as.character(stim[,2]),c(unlist(lapply(type.lst,'[[',1))),
                               c(unlist(lapply(type.lst,'[[',3))),stringsAsFactors=F )
  names(type.df) <- c('trainstim','orig.stim','orig.prop')
  type.df
}

#function makes multiple copies of neuronal vectors
#the generic method
makeMultipleCopies <- function(self,...){
  UseMethod("makeMultipleCopies",self)
}


#n: number of copies
#op: 2 - repeat the first vector n times, then the second one n times and so on
#1: repeat all the neuronal vectors as a unit n times
makeMultipleCopies.Neurons <- function(self,n,op=1){
  #algo: create copies of the neuronal vectors, and then return the modified neuron
  res <- self
  neurons <- joinLists(repList(self$neurons,n=n,op=op))
  res$neurons <- neurons
  names(res$neurons) <- 1:length(neurons)
  res
}

#first the generic method
#computeFiring <-function(x,origin,connmat,op){
computeFiring <-function(self,...){
  UseMethod("computeFiring",self)
}



#origin: the origin neuron population object, the output size is automatically fixed to the size of the origin
#connmat: the connection matrix objext
#self this population of neurons object
#dest: is the destination or target neuron structure
#stimno: odor or list of odor nos for which this should be coomputed
#op=1, do all stimuli, 2 - do stimuli specigied by odor no
#noiseop: specifies whether we should add noise (T) or not (F) to the calculation: 
computeFiring.Neurons <- function(self,connMat,dest=c(),stimno=c(),noiseop=F,op=1){
  #do a check to confirm that the matrices are correct. So, the orign and dest neurons desc should match
  #the same fields in the matrix
  origin <- self
  #cat('\nCF.neurons',self$id,stimno)
  if ( size(connMat)[1] != size(origin,op=2) ){
    cat('\n Conn Matrix and neuron populations do not match',size(connMat)[1],',',size(origin,op=2))
    return(F)
  }
  #addnoise to the neuron if specified
  stim <- getData(origin,const.N$neurons)
  cmat <- getData(connMat,const.CM$connmat)
  if(length(stimno) > 0 && op==2) nostim <- stimno #check if specific or all stimuli
  else  nostim <- 1:length(stim) 
  #now compute the firing rates for the two types of connection matrices
  res.neurons <- lapply(nostim, function(i){
    if(isDataType(cmat)==const.DataType$list){#a lsit of matrices, each matrix is a target neuron
      res <-  lapply(cmat,function(x){
        #calculate the contributions from every mb onto the apl, and then sum them up
        syn.rate <- sapply(1:length(stim[[i]]), function(j) stim[[i]][j]*sum(x[,j]))
        neuron.rate <- sum(syn.rate)
      })
    }
    else {
      neuron.rate <- cmat %*% stim[[i]]
      res <- as.vector(neuron.rate)
    }
    names(res) <- c() #we dont want an indexed list number
    unlist(res)
  })
  #set the target and then assign the computed firing rates to the appropriate stimuli in the target, and threshold appropriately
  if(length(dest)>0) { 
    target <- dest
    #you have to threshold the destination neurons not, the sourc neurons
    #res.neurons <- threshStimuli(res.neurons,getData(dest,const.N$thresh)) #threshold the neurons
  }  
  else target <- origin #if there is no destination, neurons should be the base here
  tar.neurons <- getData(target,const.N$neurons)
  tar.neurons[nostim] <- res.neurons
  #make sure the labels are appropriate
  if( length(which(names(tar.neurons[nostim])=="")) > 0) names(tar.neurons)[nostim] <- nostim
  target <- setData(target,val=tar.neurons,index=const.N$neurons)
  target
}

#calculate the firing of the different neuronal elements of the circuits
#odor: the odor no to calculate the firing for OR a user-specified odor when op=2
#noiseop: components to which we want to add noise;#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 
#12 - DAN-KC-MBON connections, 13 - noise to everything, 0 no noise, 8 - noise to APL, 9 - noise 
#to mbkc, 14 - noise to APL and APL-KC, 15 - noise to APL, PN-MB, 16 - apl, apl-kc, pn-mb, 
#17 - apl, apl-kc, pn-mb, pn, 18 - pn and apl, 19 - pn, apl-kc, apl, 19 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
#op:1, all odors 2, - odor no specified in odor,
computeFiring.FlyMBg <- function(self,odor=1,noiseop=0,op=2){
  #cat('\nCF.N')
  #noise works the following way, based on the noise option chosen in add noise, we will add noise to those 
  #components of the circuit, after we compute firing. It wont change self's firing rates, but will change, newcircuit's
  #firing rates
  newcircuit <- addNoise(self,op=noiseop) #chooses the noise option
  glom <- getData(newcircuit,const.FM$glom)#get and set the glomerular firing rate to the odor
  if(glom$basenoise==T) {#first element in the layers, Set its noise levels here, if appropriate
    glom <- addNoise(glom,stim=odor)
    newcircuit <- setData(newcircuit,val=glom,const.FM$glom)
  }
  mbkc <- computeFiring(glom,getData(newcircuit,const.FM$glommb_conn),getData(newcircuit,const.FM$mbkc),stimno=odor,op=op)
  if(mbkc$basenoise==T) mbkc <- addNoise(mbkc,stim=odor)
  newcircuit <- setData(newcircuit,val = mbkc,const.FM$mbkc) #set the mb firing rate, automatically inherits the parameters of the existing mbkc
  #get apl firing
  apl <- computeFiring(mbkc,getData(newcircuit,const.FM$kcapl_conn),getData(newcircuit,const.FM$apl),stimno=odor,op=op)
  if(apl$basenoise==T) apl <- addNoise(apl,stim=odor)
  newcircuit <- setData(newcircuit,val = apl,const.FM$apl) #set the apl firing rate
  #subtract the apl feedback
  mbkc_net <- twoSetOperations(getData(newcircuit,const.FM$mbkc),computeFiring(apl,getData(newcircuit,const.FM$aplkc_conn),
                               dest=getData(newcircuit,const.FM$mbkc_net)),dest=getData(self,const.FM$mbkc_net))
  if(mbkc_net$basenoise==T) mbkc_net <- addNoise(mbkc_net,stim=odor)
  mbkc_net <- setData(mbkc_net,val=getData(mbkc_net,const.N$neurons),index=const.N$neurons) #set the new firing rates for mbkc_net
  newcircuit <- setData(newcircuit,val = mbkc_net,const.FM$mbkc_net) #set the mb firing rate
  kcmbon_conn.lst <- getData(newcircuit,const.FM$kcmbon_conn)
  mbon.lst <- getData(newcircuit,const.FM$mbons)
  for(i in 1:length(kcmbon_conn.lst) ){
    mbon.lst[[i]] <- computeFiring(mbkc_net,kcmbon_conn.lst[[i]],mbon.lst[[i]],stimno=odor,op=op,noiseop=neuron.noise.op)
    if(mbon.lst[[i]]$basenoise==T) mbon.lst[[i]] <- addNoise(mbon.lst[[i]],stim=odor)
  }
  newcircuit <- setData(newcircuit,val = mbon.lst,const.FM$mbons) #set the mbon firing rates

  newcircuit
} 


#compute the firing rate for each layer of the circuit, basically from that layer onwards to the 
#output layer
#first the generic method
#computeFiring <-function(x,origin,connmat,op){
computeFiringLayer <-function(self,...){
  UseMethod("computeFiringLayer",self)
}


#calculate the firing of the neuronal elements of each successive layer of the circuit
#basically from that layer to the output layer
#layer: The layer from which you want to computer the firing of the neurons, 
#1 - glom, 2- mbkc, 3 - apl, 4 - mbkc_net, 5 - mbon; 3 is not an option.
#odor: the odor no to calculate the firing for OR a user-specified odor when op=2
#op:1, all odors 2, - odor no specified in odor, 
computeFiringLayer.FlyMBg <- function(self,layer=1,odor=1,op=2){
  #do we need to handle noise here?
  cat('\nCFL.N; layer ',layer)
  if(layer==5) return(self) #we are at the end, nothing more to compute, so return
  newcircuit <- self
  #returns components: origin * connMat => destination
  components <- switch(layer,list(getData(self,const.FM$glom),getData(self,const.FM$glommb_conn),getData(self,const.FM$mbkc),const.FM$mbkc),#glom
                       list(getData(self,const.FM$mbkc),getData(self,const.FM$kcapl_conn),getData(self,const.FM$apl),const.FM$apl),#mbkc
                       list(getData(self,const.FM$apl),getData(self,const.FM$aplkc_conn),getData(self,const.FM$mbkc_net),const.FM$mbkc_net
                            ,getData(self,const.FM$mbkc)),#apl, feedback from apl to mbkc
                       list(getData(self,const.FM$mbkc_net),getData(self,const.FM$kcmbon_conn),getData(self,const.FM$mbons),const.FM$mbons)#mbkc_net
  )
  if(layer==3 || layer==4){#special processing for calculating mbkc_net or mbons
    if(layer==3){#mbkc_net apl feedback
      target <- twoSetOperations(components[[5]],computeFiring(components[[1]],components[[2]]),dest=components[[3]])
      target <- setData(components[[3]],val=getData(components[[3]],const.N$neurons),index=const.N$neurons) #set the new firing rates for mbkc_net
    }
    if(layer==4){#mbons, multiple valences
      target <- components[[3]]
      for(i in 1:length(components[[2]]) ){
        target[[i]] <- computeFiring(components[[1]],components[[2]][[i]],target[[i]],stimno=odor,op=op)
      }
    }
  }
  else  target <- computeFiring(components[[1]],components[[2]],dest=components[[3]],stimno=odor,op=op)
  newcircuit <- setData.Circuit(self,val=target,components[[4]])
  newcircuit <- computeFiringLayer(newcircuit,layer=layer+1,odor=odor,op=op)
} 


# To do: havve to add code, so that theflyMBg calculates mbon firing starting from a certain layer of neurons
# 1, from glom, 2 from mbkc, 3 - from mbkc_net



#Takes an object and creates a copy and assigns the specified values from the second to the first
assignObject <-function(self,...){
  UseMethod("assignObject",self)
}

#assigns the values specified by index from the target to self.
#index: the values to be assigned, default: neurons
assignObject.Neurons <-function(self,target,index=const.N$neurons,op=1){
  temp <- self
  if (length(index)>0) temp <- setData(temp,val=getData(target,index),index)
  temp
}

#assigns the values specified by index from the target to self.
#index: the values to be assigned, default: connMat
assignObject.ConnMatrix <-function(self,target,index=const.CM$connmat,op=1){
  temp <- self
  if (length(index)>0) temp <- setData(temp,val=getData(target,index),index)
  temp
}

#assigns the values specified by index from the target to self.
#index: the values to be assigned, default: mbkc
assignObject.FlyMBg <-function(self,target,index=const.FM$mbkc,op=1){
  temp <- self
  if (length(index)>0) temp <- setData(temp,val=getData(target,index),index)
  temp
}

#Takes two sets of neurons and combines them in various ways
#first the generic method
twoSetOperations <-function(self,...){
  UseMethod("twoSetOperations",self)
}

#second, the second neurons population
#dest: the dest set of neurons
#op=1 - subtract the second population from self,
#2- divide the second from the first population
twoSetOperations.Neurons <- function(self,second,dest=c(),noiseop=0,op=1){
  if(length(dest) > 0) target <- dest
  else target <- self
  if ( size(self) != size(second) ){
    cat('\n The 2 neuron populations do not match',size(self))
    return(F)
  }
  neurons <- switch(op,lapply(1:size(self), function(i) threshVal(self$neurons[[i]] - second$neurons[[i]]) ),
                    lapply(1:size(self), function(i) threshVal(self$neurons[[i]]/second$neurons[[i]]) ) )
  names(neurons) <- 1:length(neurons)
  target$neurons <- neurons
  target #return the result
}



#get a particular variable of an object
size <- function(self,...){
  UseMethod("size",self)
}

#gets the number of stimuli encoded by the neuron population or the number of neurons
#op=1, no of stimuli, 2 - no of neurons
size.Neurons <- function(self,op=1){
  switch(op,length(self$neurons),length(self$neurons[[1]]) )
}

#the dimnesions of the matrix
size.ConnMatrix <- function(self,op=1){
  #dim(self$connMat)
  self$size
}

#the dimnesions of the FlyMBg circuit
size.FlyMBg <- function(self,op=1){
  #dim(self$connMat)
  #cat('\nSize FlyMBg')
  res <- c(size(self$glom),size(self$glom,op=2),nrow(self$trainingval),length(self$teststim),size(self$mbkc_net,op=2),
    length(self$dans),size(self$dans$`1`,op=2),size(self$mbons$`1`),length(self$mbons) )
  names(res) <- c('#stim','#glom','#train stim','#test stim','#kcs','#dan-types','#dans','#mbon stim','#mbon-types')
  res
}

#get a particular variable of an object
getData <- function(self,...){
  UseMethod("getData",self)
}


#gets the requested data type for connection matrix
#names(connMatrix) <- c('connmat','mattype','origin','target','params','type','size','learningrule')
#datatype: 1 - connmat, 2 - mattype, 3 - origin, 4 - target, 5 - params, 6 - type, 7 - size, 8 - learningrule
#op=1, get the vals, 2 - get the descriptions, 
getData.ConnMatrix <- function(self,datatype=1,op=1){
  #cat('\ngetdataConnMatrix',datatype)
  switch(op,switch(datatype,
                   self$connMat,self$type,self$source,self$target,self$params,self$size,self$syndist,self$noise,
                   self$learningRule,self$gain,
                   self$dopNo,self$dopParams,self$dopDist,self$dopType,self$dopConn,self$curTime,self$synActFn),
         switch(datatype,
                'connmat','type','source','target','params','size','syndist','noise','learningRule','gain',
                'dopNo','dopParams','dopDist','dopType','dopConn','curTime','synActFn'))
}

#const.CM <- list(connmat=1,type=2,source=3,target=4,params=5,size=6,syndist=7,noise=8,
#                 learningrule=9,dopno=10,dopparams=11,dopdist=12,doptype=13,dopconn=14,curtime=15)



#gets the requested data type for Neurons
#names(connMatrix) <- c('connmat','mattype','origin','target','params','type','size','learningrule')
#datatype: 1 - neurons, 2 - neurondesc, 3 - dist, 4 - basenoise, 5 - neurontype, 6 - params
#op=1, get the vals, 2 - get the descriptions, 3 - get the neurons for the stimulus no specified
getData.Neurons <- function(self,datatype=1,stimno=c(),op=1){
  switch(op,switch(datatype,
                   self$neurons,self$neurondesc,self$dist,self$basenoise,self$neurontype,self$params,self$posns,self$noise,self$id,
                   self$thresh), #actual data types
         switch(datatype,'neurons','neurondesc','dist','basenoise','neurontype','params','noise','id','thresh'),#datatype descriptions
         self$neurons[[stimno]]
  )
}

#gets the rquiested data for the FlyMBg circuit
#valence: the connection matrix for the valence specified in valence
#op=1, get the vals, 2 - get the descriptions, 3 - the required connection matrix, 4 - the mbons for a valence
#5 - the mbon firing rate for a particular odor and valence
getData.FlyMBg <- function(self,datatype=1,valence=1,odorno=1,op=1){
  #if (op==5) cat('\ngetData FlyMBg',datatype,op,odorno,class(self$mbons[[valence]]))
  switch(op,switch(datatype,self$glom,self$mbkc,self$mbkc_net,self$apl,self$dans,self$mbons,
                   self$glommb_conn,self$kcapl_conn,self$aplkc_conn,self$kcmbon_conn,self$noodors,
                   self$valences,self$trainingval,self$teststim,self$tatparams,self$origstim),
         switch(datatype,'glom','mbkc','mbkc_net','apl','dans','mbons',
                'glommb_conn','kcapl_conn','aplkc_conn','kcmbon_conn','noodors','valences','trainingval',
                'teststim','tatparams','origstim'),
         self$kcmbon_conn[[valence]]$connMat,
         self$mbons[[valence]],
         getData(self$mbons[[valence]],op=3,stimno=odorno)
  )
}

#gets the rquiested data for a circuit
getData.Circuit <- function(self,datatype=1,op=1){
  cat('\nCircuit getData Nothing for now')
}


#setup generic method to add noise
setData <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("setData",self)
  }
}


#will update the x object with whatever val has to be updated
#val: value to be updated 
#index: index of the value to be updated. for indices look at getData
setData.Neurons <- function(self,val,index,op=1){
  #maybe introduce a type check
  #cat('\nBefore Setting neuron value',getData(self,const.N$neurons))
  temp <- self
  temp[[index]] <- val #updates the value of x$index
  #cat('\nSetting neuron value',getData(self,const.N$neurons))
  temp
}


#will update the x object with whatever val has to be updated
#val: value to be updated 
#index: index of the value to be updated. for indices look at getData
setData.ConnMatrix <- function(self,val,index,op=1){
  #maybe introduce a type check
  temp <- self
  temp[[index]] <- val #updates the value of x$index
  temp
}

#will update the x object with whatever val has to be updated
#val: value to be updated 
#index: index of the value to be updated. for indices look at getData
setData.FlyMBg <- function(self,val,index,op=1){
  #maybe introduce a type check'
  #cat('\nsetDATA FlyMBg',index)
  temp <- self
  temp[[index]] <- val #updates the value of x$index
  temp
}


#will update the x object with whatever val has to be updated
#val: value to be updated 
#index: index of the value to be updated. for indices look at getData
setData.Circuit <- function(self,val,index,op=1){
  #maybe introduce a type check'
  #cat('\nsetDATA circuit')
  temp <- self
  temp[[index]] <- val #updates the value of x$index
  temp
}


#generic method to add stuff to thee objects 
addData <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("addData",self)
  }
}

#this adds extra stimuli to the Neurons object, If there is more stuff to be addded, can make it generic later
#nostimuli: extra stimuli to be added
#params: either new params or if its empty just use the exisiting ones
#stimuli: the stimuli to be added, can either be glom$neurons or glom
#posn at which the stimuli should be inserted. If 0, add at the end
#op=1: generate the stimuli, 2 - use the stimuli passed in as an argument in stim
addData.Neurons <- function(self,nostimuli=1,stimuli=c(),params=c(),posn=0,op=1){
  #cat('\naddData Neurons')
  size <- size(self,op=2)
  if(length(params) >0 ) tparams <- params
  else tparams <- self$params
  #generate the neurons
  if(length(stimuli)>0) {
    if(c('Neurons') %in% class(stimuli)) neurons <- getData(stimuli,const.N$neurons) 
    else neurons <- stimuli #new stimuli to add 
  }
  else  neurons <- genNeurons(number = size,basenoise = self$basenoise,dist = self$dist,params = tparams,noStimuli = nostimuli)
  #add to the oldvals and update
  oldvals <- getData(self,const.N$neurons)
  if(posn==0) pos <- length(oldvals)+1
  else pos <- posn
  newvals <- insertElemList(lst = oldvals,elem = neurons,posn = pos) #c(oldvals,neurons)
  names(newvals) <- 1:length(newvals)
  newneuron <- setData(self,val=newvals,const.N$neurons)
}

#adds new stimuli to the circuit, and it can either become part of the training set or testing set
#stimtype: the stimulus type to which these stimuli should be added. Default is testing stimuli
#1 - irgininal stimuli, 2 - training stimuli, 3 - testing stimuli, default:3
#stimuli: the stimuli to be added, in the form of a list of neuron firing rates like self$glom$neurons
#can either be self$glom or the circuit itself
#addtype: 1-add, 2 - replace
#op=1, 
addData.FlyMBg <- function(self,nostimuli=1,stimuli=c(),stimtype=3,valences=c(),op=1){
  #cat('\naddData flymbg')
  glom <- getData(self,const.FM$glom)
  if(c('FlyMBg') %in% class(stimuli)) glomstimuli <- getData(getData(stimuli,const.FM$glom),const.N$neurons)
  else if(c('Neurons') %in% class(stimuli)) glomstimuli <- getData(stimuli,const.N$neurons)
       else glomstimuli <- stimuli
  #based on stimtype determine the posn where this should be added
  stim.df <- switch(stimtype,getData(self,const.FM$origstim),getData(self,const.FM$trainingval),
                    getData(self,const.FM$teststim)) 
  posn <- stim.df[nrow(stim.df),2] + 1 
  newglom <- addData(glom,nostimuli=nostimuli,stimuli=glomstimuli,posn=posn)
  circuit <- setData(self,val=newglom,const.FM$glom)
  #we should update all the other components of the circuit too, mbkc, mbkcnet, and mbons
  circuit <- computeFiring(circuit,op=1)
  #update training and testing odors
  circuit <- updateStimParams(circuit,stimtype=stimtype,stimuli=length(glomstimuli),valences=valences,op=1)
}

#generic method to add stuff to thee objects 
setStimuli <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("setStimuli",self)
  }
}

#function sets the stimuli of stimtype and updates the stim dataframe accordingly
#nostimuli: extra stimuli to be added
#params: either new params or if its empty just use the exisiting ones
#stimuli: the stimuli to be added.
#posn at which the stimuli should be inserted. If 0, add at the end
#op=1: generate the stimuli, 2 - use the stimuli passed in as an argument in stim
setStimuli.Neurons <- function(self,nostimuli=1,stimuli=c(),params=c(),posn,op=1){
  size <- size(self,op=2)
  if(length(params) >0 ) tparams <- params
  else tparams <- self$params
  #generate the neurons
  if(length(stimuli)>0) {
    if(c('Neurons') %in% class(stimuli)) neurons <- getData(stimuli,const.N$neurons) 
    else neurons <- stimuli #new stimuli to add 
  }
  else  neurons <- genNeurons(number = size,basenoise = self$basenoise,dist = self$dist,params = tparams,noStimuli = nostimuli)
  #add to the oldvals and update
  oldvals <- getData(self,const.N$neurons)
  newvals <- replaceElemList(lst = oldvals,elem = neurons,stpos = posn[1],endpos = posn[2]) #c(oldvals,neurons)
  #cat('\n old new',length(oldvals),length(newvals),posn,length(neurons))
  names(newvals) <- 1:length(newvals)

  newneuron <- setData(self,val=newvals,const.N$neurons)
}

#adds new stimuli to the circuit, and it can either be designated as the training set or testing set or the original set
#stimtype: the stimulus type to which these stimuli should be added. Default is testing stimuli
#1 - original stimuli, 2 - training stimuli, 3 - testing stimuli, default:3
#stimuli: the stimuli to be added, in the form of a list of neuron firing rates like self$glom$neurons
#can either be self$glom or the circuit itself
#valences: if mix stimuli, the valences is needed to determine valence
#op=1, normal stimuli,
setStimuli.FlyMBg <- function(self,nostimuli=1,stimuli=c(),stimtype=3,valences=c(),op=1){
  #cat('\nsetSTimuli flymbg',nostimuli)
  glom <- getData(self,const.FM$glom)
  if(c('FlyMBg') %in% class(stimuli)) glomstimuli <- getData(getData(stimuli,const.FM$glom),const.N$neurons)
  else if(c('Neurons') %in% class(stimuli)) glomstimuli <- getData(stimuli,const.N$neurons)
       else glomstimuli <- stimuli
  #based on stimtype determine the posn where this should be added
  stim.df <- switch(stimtype,getData(self,const.FM$origstim),getData(self,const.FM$trainingval),
                    getData(self,const.FM$teststim)) 
  stpos <- switch(stimtype,1,stim.df[1,2],stim.df[1,2])
  endpos <- stim.df[nrow(stim.df),2]
  newglom <- setStimuli(glom,nostimuli=nostimuli,stimuli=glomstimuli,posn=c(stpos,endpos))
  circuit <- setData(self,val=newglom,const.FM$glom)
  #get the new odor numbers for this stimtype now
  if(length(stimuli)>0) newstimno <- length(glomstimuli)
  else newstimno <- nostimuli
  #we should update all the other components of the circuit too, mbkc, mbkcnet, and mbons
  circuit <- computeFiring(circuit,op=1)
  #update training and testing odors
  circuit <- updateStimParams(circuit,stimtype=stimtype,stimuli=newstimno,valences=valences,op=2)#
  #cat('\nsetstim',size(circuit$glom,op=1))
  circuit  
}



#generic method to update the original, training and test params 
updateStimParams <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("updateStimParams",self)
  }
}


#FlyMBg method to update the training and test params, when the circuit's stimuli and training stimuli are different
#stimuli: the no of new stimuli that are added
#valences: are the valences for the test stimulations. If not present then just use the training vals or protocol
#op=1,stimuli that are added, 2 - stimuli that are set; needed for determining valence type if valences is NULL
updateStimParams.FlyMBg <-function(self,stimtype=1,stimuli=c(),valences=c(),op=1){
  #get the existing stim dfs, and put them in a list
  origstim <- getData(self,const.FM$origstim)
  trainingval <- getData(self,const.FM$trainingval)
  teststim <- getData(self,const.FM$teststim)
  newstim.lst <- list(origstim,trainingval,teststim)
  newstim <- 1:stimuli
  if(length(valences)==0){#generate new valences according to exisiting rules
    stim.df <- cbind.data.frame(genSequence(length(getData(self,const.FM$valences)),
                                            size = length(newstim),op = op[1]),newstim,newstim ) 
  }
  else {#valence is specified
    stim.df <- cbind.data.frame(valences[,1],newstim,valences[,2])
  }
  colnames(stim.df) <- colnames(newstim.lst[[stimtype]])
  #now add these to target stim
  if(op==1) newstim.lst[[stimtype]] <- joinDF(newstim.lst[[stimtype]],stim.df)
  else newstim.lst[[stimtype]] <- stim.df
  #update the odor numbers for all types of stim based on type of stimuli
  newstim.lst[[1]][,2] <- 1:nrow(newstim.lst[[1]])
  newstim.lst[[2]][,2] <- 1:nrow(newstim.lst[[2]]) + nrow(newstim.lst[[1]])
  newstim.lst[[3]][,2] <- 1:nrow(newstim.lst[[3]]) + nrow(newstim.lst[[1]]) + nrow(newstim.lst[[2]])
  #print(newstim.lst)
  #sanity check: no of glomerli = sum of all three dfs && they should be sequential
  circuit <- setData(self,val=newstim.lst[[1]],const.FM$origstim) 
  circuit <- setData(circuit,val=newstim.lst[[2]],const.FM$trainingval) 
  circuit <- setData(circuit,val=newstim.lst[[3]],const.FM$teststim) 
}



#this function sets the stimuli that are to be used for training and testing 
#testodors: 1 - becomes part of the training and testing set, 2 - becomes part of the testing set only, 
#3 - change it so that this is the new testing set 
setTATParams <- function(self,...){
  if(is.object(self)){#it is a class, so do it
    #cat('\nseetTAT generic',class(self))
    UseMethod("setTATParams",self)
  }
}

#tatparams: two choices each for training and testing c(trainparams,testparams)
#trainparams: 1 - alternate, 2 - sequential, 3 - random
#testparams: 1 - is valence training confirmed, 2 - pairwise discrimination
#train: the odors to be designated as training odors, a vector
#test: the odors to be designated as testing odors, a vector
#valences: are the valences for the test stimulations. If not present then just use the training vals or protocol
#op:1 - set it up according to input params, 2 - upgrade the circuit, i.e., dont change existing train params
#op=3, just change tatparams
setTATParams.FlyMBg <- function(self,train=c(),test=c(),tatparams=c(1,1),valences=c(),op=1){
  if(op==3) return(setData(self,val=tatparams,const.FM$tatparams))
  #cat('\nsetTAT')
  circuit.size <- size(self)
  if(length(train) == 0) trainstim <- 1:circuit.size[1]
  else trainstim <- train
  if(length(test) == 0) teststim <- 1:circuit.size[1]
  else teststim <- test
  #ok, set training and testing params according to the inputs
  trainingval <- cbind(genSequence(length(getData(self,const.FM$valences)),
                                   size = length(trainstim),op = tatparams[1]),trainstim )
  #if valences is specified and its length matches test, those are the test vals otherwise just return the training vals
  if(length(valences)>0 && length(valences)==length(test)) testval <- cbind.data.frame(valences,test,test)
  else testval <- cbind.data.frame(trainingval[teststim,1],teststim,teststim)
  #should do teststimulations
  circuit <- setData(self,val=trainingval,const.FM$trainingval)
  circuit <- setData(circuit,val=testval,const.FM$teststim)
  circuit <- setData(circuit,val=tatparams,const.FM$tatparams)
}
  
  
  
#newNeurons <- function(number=1,basenoise=c(0,0),neurontype=1,neurondesc='excitatory',posns=0,dist=1,params=c(10,1),noStimuli=1){
#generate the neurons
#  neurons <- genNeurons(number = number,basenoise = basenoise,dist = dist,params = params,noStimuli = noStimuli


#returns the type of structure this component is, 1 - neuron, 2 - connection matrix, 3 - Circuit, 4 - learning matrix, 
#5 - list of neurons, 8 - list of learning matrices
#comp: the component to be tested
#op=1, return string, 2 - return number as shown abovve
getCompType <- function(comp,op=1){
  compclass <- class(comp)
  if('Neurons' %in% compclass) return(switch(op,'Neurons',1))
  if('Learning' %in% compclass) return(switch(op,'Learning',4))
  if('ConnMatrix' %in% compclass) return(switch(op,'ConnMatrix',2))
  if('Circuit' %in% compclass) return(switch(op,'Circuit',3))
  if('list' %in% compclass) return(switch(op,paste('list: ',getCompType(comp[[1]],op=op),sep = ''),
                                          getCompType(comp[[1]],op=op)+4 ) ) #returns the type of the first list element +4 designates list
  0 #nada, return 0
}


#setup generic method to add noise
addNoise <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("addNoise",self)
  }
}


#adding noise to neurons, returns a noised version of the neurons
#noiseparams: the noise params (distribution,dist. params)
#noise distribution: 1 - gaussian, 2 - uniform
#odornos: vector of odors nos to which noise should be added, otherwise, add to all odors
#op: 1, neurons that are based on stimuli, 2 - based on learninig like dop., whose numbers are different from stim. nos
addNoise.Neurons <-function(self,noiseparams=c(),stim=c(),op=1){
  #cat('\nadding noise to neurons',stim,op)
  #check where the noiseparameter is
  if(length(noiseparams)>0) noise <- noiseparams
  else noise <- self$noise
  stim.neurons <- self$neurons
  if(op==1 && length(stim) > 0) stimnos <- stim
  else stimnos <- 1:length(stim.neurons)
  #go through eaach odor and to each neuron add 0-mean noise
  #cat('\naddnoise',stimnos,'stim',unlist(stim.neurons[1:2]),str(noise))
  res.neurons <- lapply(stimnos,function(i){
    noisefactor <- genNoise(n=length(stim.neurons[[i]]),noisedist = noise[[1]],noise = noise[[2]][2]) #the noise vector
    #cat('\t nf:',noisefactor,';')
    neuron <- switch(noise[[3]],stim.neurons[[i]]*(1 + noisefactor),stim.neurons[[i]] + noisefactor ) #add based on type
    neuron
  })
  stim.neurons[stimnos] <- res.neurons
  #cat('\nwith noise',unlist(stim.neurons[1:2]))
  #names(res.neurons) <- names(self$neurons)
  neuron <- setData(self,val=stim.neurons,const.N$neurons)
}

#adding noise to the connection matrix, and returns a noised version of the original connection matrix
addNoise.ConnMatrix <-function(self,stimnos=c(),noiseparams=c()){
  #cat('adding noise to conn mat')
  #check where the noiseparameter is
  if(length(noiseparams)>0) noise <- noiseparams
  else noise <- self$noise
  if(noise[[3]]==0) return(self) #no noise, so nothing to do
  #add 0-mean noie again
  origmat <- getData(self,const.CM$connmat)
  #check for the two types of matrices
  if(isDataType(origmat)==const.DataType$list)  
    connmat <- genNoiseMatList(origmat,noisedist = noise[[1]],noise = noise[[2]][2],noisetype = noise[[3]]) 
  else {
    connsize <- dim(origmat)
    matnoise <- matrix(genNoise(n=connsize[1]*connsize[2],noisedist = noise[[1]],noise=noise[[2]][2]),nrow = connsize[1])
    connmat <- switch(noise[[3]],self$connMat * (1+matnoise),self$connMat + matnoise )
    if("dopType" %in% names(self)){#it is a learning matrix, so strengths shouldn't go outside range
      connmat[connmat < 0] <- 0
      connmat[connmat > 0.99] <- 0.99
    }
  }
  # matnoise <- matrix(rnorm(connsize[1]*connsize[2],mean=noiseparams[1],sd=noiseparams[2]),nrow = connsize[2])
  # connmat <- self$connMat * (1+matnoise)
  newconn <- setData(self,val=connmat,const.CM$connmat)
}




#adding noise to the fly circuit
#odors: vector of odors nos to which noise should be added, otherwise, add to all odors
#noiseset: 0, all odors or sets, 1 - orig stim,  2 - training set only, 3 - testing set only
#op: components to which we want to add noise; 0 - no noise, 1, glom/pns, 2 - PN-MB, 3 - KC-APL, 4 - APL-KC
#5- mb-mbon connection matrix, 6 - mbons, 7 - dans, 8 - noise to APL, 9 - noise to mbkc, 
#10 - PN & PN-MB, 11 - PN, PN-MB, APL-KC, 12 - DAN-KC-MBON connections, 13 - empty for not, 0 no noise, 
#14 - APL and APL-KC, 15 - APL, PN-MB, 16 - apl and pn, 17 - apl & mbkc, 18 - mbkc & pn-mb, 
#19 - mbkc & apl-kc , 20 - apl, apl-kc, pn-mb, 21 - apl & apl-kc & mbkc, 22 - apl, apl-kc, pn-mb, pn, 
#30 - apl, apl-kc, pn-mb, pn, mbkc
#100 - add noise according to noise params to every structure
addNoise.FlyMBg <-function(self,odors=c(),noiseset=0,op=1){
  #cat('adding noise to FlyMBg')
  if(op==0) return(self) # no noise to add, just return self
  circuit <- self
  #designate odors for noise, default is all odors, if noiseset > 0 and odors is empty either testing or training sets
  odornos <- chooseStimuli(circuit,stimset = noiseset,stim = odors)
  #choose the components to which we will be adding noise
  compids <- setNoiseCompIds(op)
  #if only a single number is returned
  if(length(compids)==1) {
    if(compids==0) return(circuit) #no noise to add
    circuit <- addCompNoise(circuit,compid = compids,stimnos = odornos)
    return(circuit)
  }
  else {
    if(op < 100) {#multiple compids are chosen, learning and nonlearning structures are separate options
      #cat('\nsetnoise',length(noiseparams))
      for(i in seq_wrap(1,length(compids))){
        circuit <- addCompNoise(circuit,compid = compids[i],stimnos = odornos)
      }
    } else {#op=100, so set noise to everything
      #non-learning structures
      for(i in seq_wrap(1,length(compids[[1]]))){
        circuit <- addCompNoise(circuit,compid = compids[[1]][i],stimnos = odornos)
      }
      #now, set the valence related (learning) bits
      #got through each component, then their valences, and set noise
      for(j in seq_wrap(1,length(compids[[2]]))){
        circuit <- addCompNoise(circuit,compid = compids[[2]][j],stimnos = odornos)
      }    
    }
  }
  circuit
}


#change the noise of the component. Normal neurons/matrix and the learning ones
#for connmection matrices, just add noise. For neurons, set the basenoise=T, and then add noise
#within computeFiring.neurons
#self: the circuit in question; compid: id of the component so that we can retireive and modify it.
#op
addCompNoise <- function(self,compid,stimnos,op=1){
  circuit <- self
  comp <- getData(circuit,compid) #get the component
  comptype <- getCompType(comp,op=2) #type of component: neuron, matrix, or learning component
  #cat('\naddcompnoise',compid,comptype)
  #neurons
  if(comptype==1) {
    #old code: circuit <- setData(circuit,val=addNoise(comp,stim=stimnos),compid)
    comp <- setData(comp,val=T,const.N$basenoise) #set basenoise=T and then put it back in the circuit
    circuit <- setData(circuit,comp,compid)
  }
  if(comptype==2) circuit <- setData(circuit,val=addNoise(comp),compid)
  if(comptype==5 || comptype==8){#do lists of neurons or matrices
    comp <- getData(circuit,compid) #get the component
    for(j in 1:length(comp)){ #iterate over the number of components
      comptype <- getCompType(comp[[j]],op=2)
      if(comptype==1){#liam check the new code vs oldcode
        #we might want to set noise selectively to the DANs based on stimuli
        comp <- setData(comp,val=T,const.N$basenoise) #set basenoise=T and then put it back in the circuit
        #the setting of specific stimnos will have to be done through computeFiring.neurons
        # if(comptype==1 && compid!=const.FM$dans) comp[[j]] <- addNoise(comp[[j]],stim=stimnos,op=1) # normal neurons
        # if(comptype==1 && compid==const.FM$dans) comp[[j]] <- addNoise(comp[[j]],stim=stimnos,op=2) # dop neurons
      }
      if(comptype==2 || comptype == 4) comp[[j]] <- addNoise(comp[[j]]) #learning or normal matrix
    }
    circuit <- setData(circuit,val=comp,compid)
  } 
  circuit
}


#choose or designate which stimuli are chosen for modification
#self: the circuit in question
#stim: vector of stimulus nos to which noise should be added, otherwise, add to all odors
#stimset: 0, all stimuli or sets, 1 - orig stim,  2 - training set only, 3 - testing set only
#op=1
chooseStimuli <- function(self,stim=c(),stimset=1,op=1){
  #designate odors for noise, default is all odors, if noiseset > 0 and odors is empty either testing or training sets
  if(length(stim)==0){
    if(stimset > 0){ #just the set, usually training set
      stimnos <- switch(stimset,getData(self,const.FM$origstim)[,2],getData(self,const.FM$trainingval)[,2],
                        getData(self,const.FM$teststim)[,2])
    }
    else stimnos <- c(getData(self,const.FM$origstim)[,2],getData(self,const.FM$trainingval)[,2],
                      getData(self,const.FM$teststim)[,2]) #basically all odors
  }
  else stimnos <- stim
  stimnos
}

#calculates the firing rate of the target neurons given source neurons and the connection
#matrix, which can either be a matrix or a list of matrices (each matrix denotes a target)
#origin: are of the class neurons and the source neurosn that fuel the target
#connmatix: can be of the form a matrix or a list of matrices, all specified within the matrix class
##dest: is the destination or target neuron structure
#stimno: odor or list of odor nos for which this should be coomputed
multMatNeurons <- function(origin,connmat,dest=c(),stimno=c(),op=1){
  stim <- getData(origin,const.N$neurons)
  cmat <- getData(connmat,const.CM$connmat)
  if(length(stimno) > 0) nostim <- stimno
  else  nostim <- 1:length(stim) #no of stimuli
  res.neurons <- lapply(nostim, function(i){
    if(isDataType(cmat)==4){#a lsit of matrices, each matrix is a target neuron
      res <-  lapply(cmat,function(x){
        #calculate the contributions from every mb onto the apl, and then sum them up
        syn.rate <- sapply(1:length(stim[[i]]), function(j) stim[[i]][j]*sum(x[,j]))
        neuron.rate <- sum(syn.rate)
      })
    }
    else {
      neuron.rate <- cmat %*% stim[[i]]
      res <- as.vector(neuron.rate)
    }
    names(res) <- c() #we dont want an indexed list number
    unlist(res)
  })
  #set the target and then assign the computed firing rates to the appropriate stimuli in the target
  if(length(dest)>0) target <- dest
  else target <- origin #if there is no destination, neurons should be the base here
  tar.neurons <- getData(target,const.N$neurons)
  tar.neurons[nostim] <- res.neurons
  #make sure the labels are appropriate
  if(names(tar.neurons[nostim[1]])=="") names(tar.neurons)[nostim] <- nostim
  target <- setData(target,val=tar.neurons,index=const.N$neurons)
  target 
}


#setup generic method to print object
printObj <-function(self,...){
  if(is.object(self)){#it is a class, so do it
    UseMethod("printObj",self)
  }
}

#op=1, dont print out the neurons, just the population size, 2 - print out a list of neurons
#we dont need a method dispatcher since this is a generic method already
printObj.Neurons <- function(self,op=1){
  #cat('Neuron: ')
  cat('#stimuli: ',size(self),', pop. size: ',size(self,op=2),'\n')
}

#make printing neurons into a generic ffunction
print.Neurons <- function(self,...){
  printObj(self,...)
}

#op=1, dont print out the neurons, just the population size, 2 - print out a list of neurons
#we dont need a method dispatcher since this is a generic method already
printObj.ConnMatrix <- function(self,op=1){
  cat('ConnMatrix: ')
  cat('\t',self$type)
  cat('\t',self$size)
  if(is.null(self$dopNo)){ #not a learning matrix
    cat(' no learning ')
  }
  else{ #learning matrix details
    cat('\nlearning matrix: dop no ',self$dopNo,', size of n/w: ',dim(self$dopConn))
  }
  cat('\n')
}

print.ConnMatrix <- function(self,...){
  printObj(self,...)
}

#op=1, print everything,
#2: print the noise component for all of its components
printObj.FlyMBg <- function(self,op=1){
  if(op==1){
    cat('FlyMBg circuit \n')
    cat('glom: ')
    printObj(self$glom)
    cat('kc: ')
    printObj(self$mbkc)
    cat('apl: ')
    printObj(self$apl)
    cat('mbons: \n')
    lapply(self$mbons, function(i) printObj(i))
    cat('dans: \n')
    lapply(self$dans, function(i) printObj(i))
    cat('glom-kc n/w')
    printObj(self$glommb_conn)
    cat('kcmbon_conn:\n')
    lapply(self$kcmbon_conn, function(i) {
      printObj(i)
    })
  }
  if(op==2){
    cat('FlyMBg circuit noise params\n')
    cat('glom: ')
    cat(unlist(self$glom$noise))
    cat('kc: ')
    cat(unlist(self$mbkc$noise))
    cat('apl: ')
    cat(unlist(self$apl$noise))
    cat('mbons: ')
    lapply(self$mbons, function(i) cat(unlist(i$noise)))
    cat('dans: ')
    lapply(self$dans, function(i) cat(unlist(i$noise)))
    cat('glom-kc n/w ')
    cat(unlist((self$glommb_conn$noise)))
    cat('apl-kc n/w ')
    cat(unlist((self$aplkc_conn$noise)))
    cat('kcmbon_conn:\n')
    lapply(self$kcmbon_conn, function(i) {
      cat(unlist(i$noise))
    })
    
  }  
}

#--------------------------------------------------------------------------------------------------------------
#helper olfactory and general circuit functions

#generates connection matrix, default gaussiam, whose params are specified by params
#the matrix is size rows x cols: rows: target neurons, cols: source neurons
#rows - the rows, e.g., each one might specify a pcx neuron, the target neurons
#cols - the cols e.g., each one might specify a glomerulus, the source neurons
#dop: dopaminergic neurons, for type 5, 
#type - 1- agaussian matrix, 2 -gamma distribution, 3 - compound poisson Gamma distribution, 
#4 - uniform matrix, 5 - compound poisson Gamma distribution, 6 - constant matrix all synapses same strength 
#specified by params, 
#7: generates connections between source neuron(s) and target neuron(s) given the number of synapses "n". The n
# used in a poisson to generate connection number for each source neuron, and then each synapse is generated from a distribution
# and stored in a vector. For more than one source neuron, this would become a matrix  
#8: same as 8, excpet the number of synapses are prespecified as a vector so (vector of synapses, distribution of synapses) 
#10 - fly almb ccircuit
#20: connections are based on distribution and spatial params, e.g., between pcx2e and pcx2i cells, the index
#in the vector indicates the spatial position
#nosyn - not needed: average number of synapses between a glomerulus and neuron 
##noise: list(noisedist,cnoise amplitude,type of noise additive or multiplicative which is default)
#noise: the noise to be added, it is gaussian mean 0 with stdev given by noise
#noise, e.g. noise = .1 is 10 % noise
#if compound poisson gamma params = (poisson mean,alpha,beta)
#params - , no of synapses, and paramss of the distribution
#if type=4, it is a list of list(glomnopars,glompars,synpars)
#for type = 5, it is c(avg no of synapses made with each KC-MBON synapse,strength of the synapse = gaussian(2,0.5))
#for type = 6, specifies the strenght of the synapses
#for type = 7, it is c(avg. no of synapses by each source neuron onto a target,connection distribution)
#for type = 8, no of synapses is specified in vector specified in params for # synapses.
#for type = 10, list(c(glomerular no params,c(# claws,poisson success param) ),c(glomerular params,i.e., the parameters 
#of the hypergeometric distribution used to choose glomeruli),c(the parameters for the synaptic distribution))
#for type = 20, params = list(#connections each source neuron makes,number of neurons per spatial position for the source, 
#for thetarget, l - the width of 95 % of the spread, L - the total width ofthe region, synparams - 
#the distribution for synaptic strengths)
#gain: gain of the connection strength. Multiply whatever matrix you get with the gain.
genConnMat <- function(rows,cols,params=c(1,1.15,.16),type=3,noise=list(0,c(0,0),0),gain=1,op=1){
  if(type==1){#gaussian matrix
    mat <- matrix(rnorm(rows*cols,mean = params[1],sd = params[2]),nrow = rows)
  }
  if(type==2){#gamma matrix
    mat <- matrix(rgamma(rows*cols,shape = params[1],scale = params[2]),nrow = rows)
  }
  if(type==3 || type == 5){#compund poisson gamma/uniform matrix, could be made faster by breaking all 
    #cat('\ngenconnmat.x',str(params))
    mat <- genLongRangeConnMatrix(rows = rows,cols,params = params,type = type) #fills the matrix by columnns
  }
  if(type==4){#uniform matrix
    #cat('\nuniform',params)
    mat <- matrix(runif(rows*cols,min = params[1],max = params[2]),nrow = rows)
  }
  #aplkc_conn <- newConnMatrix(c(1,kcno),syndist = aplkc_conn_params[[1]],dop=c(),params = aplkc_conn_params[-1])
  
  if(type==6){#constant matrix, all synapses are the same strength
    mat <- matrix(rep(params,rows*cols),nrow = rows)
  }
  if(type==7 || type==8){#each source neurons makes "n_i" connections with each target, 
    #where n is derived from a Poisson distribution with mean "n" 
    #7: generate the synapses, 8 - use the vector specified in params for # synapses.
    #cat('\noption',type,' : ',str(params))
    mat <- generateIndSynMatrices(target = rows,source = cols,params,op=type-6)
  }
  if(type==10){#fly al mb circuit
    mat <- genFlyPNMBconn(rows = rows,cols = cols,glomnopars = params[[1]],glompars = params[[2]],synpars = params[[3]])
  }
  if(type==20){#spatial connection matrix
    mat <- genSpatialConnMatrix(rows = rows,cols = cols,params = params)
  }
  #print(mat)
  #adding noise if needed
  if (type==7 || type==8) resmat <- genNoiseMatList(mat,noisedist = noise[[1]],noise = noise[[2]][2],
                                         noisetype = noise[[3]]) #add noise to each individual matrix in the list
  else  resmat <- mat*(1 + genNoiseMat(rows = rows,cols = cols,noisedist = noise[[1]],noise = noise[[2]][2],
                                          noisetype = noise[[3]]))
  #sapply(1:cols,function(x) rnorm(rows,mean = 0,sd = noise)))
  resmat
}



#this generates a connection matrix or a set of connection matrices between a bunch of source neurons 
#and target neurons where the each synapse between the source and target have their own entry. Thus, if there is only 
#one target like the APL, the output is one matrix, but if there are multiple neurons then it is a list of matrices
#returns a matrix or a list of matrices, where each column contains the synapses made by that source neurons with the target
#rows: number of targets
#cols: number of source neurons
#params: params[1]: gives the nujmber of connections between source and target neurons, one number means we use that as the mean of
#a poisson distribution. a vector would mean that each number gives the mean of each source to all targets, and matrix would 
#be all source to all targets; params[2:]: the distribution for the synapses
#c(synapse distribution type,params of synapse distribution: c(no of synapse,distribution type
#distribution type=1 - gamma distibuteion of strngth, 2 - unform distibution of strength
#, distrinbution params)), only 1 synapse -> c(6,c(1)) or list(7,c(6,1,1.15,0.12))
#learningrule: the rule for updating the strength of the synapses.(alpha,beta,ec50,T:sat fn,F:linear function) additive decrease s-alpha, multiplicative inc: s*beta 
#op1: params are a vector, 2 - params are a list list(distribution parameters,synapse vector from
#which we choose the number of synapses: it gives the synapses from each source to the targets) 
generateIndSynMatrices <-function(source,target,params=c(6,1,1.15,0.12),op=1){
  targets <- 1:target
  #break up the params if needed
  if(op==2) synpars <- params[[1]]
  else synpars <- params
  mats.ls <- lapply(targets, function(i){
    #determine the number of synapses, and then we have to build matrics 
    #2 options: either constrain all target neurons to have the same number of source synapses
    #or each target neuron has a potentially diff. no depending on rpois
    if(op==1 || isDataType(params)==const.DataType$vector){
      nosyn <- rpois(source,synpars[1])
      maxsyn <- max(nosyn)
    } else {
      #cat('par',params[[2]])
      maxsyn <- params[[2]][i];nosyn <- maxsyn
    }
    #each source neuron makes nosyn synapses with the target neuron
    #so matrix with nosyn rows and source number of cols, params[-1] leaves out the first synapse
    #cat('\tsy:',maxsyn,source,op)
    mat <- matrix(genSynapses(maxsyn*source,distop=synpars[2],params = synpars[-(1:2)]),ncol = source)
    #fix it to have 0s where there are no synapses.
    mat_update <- updateSynMat(mat,nosyn)
  })
  names(mats.ls) <- targets
  mats.ls
}

#generate the number of synapses specified by no, according the distribution distop, and params
#distop: 1 - gamma distibuteion of strngth, 2 - unform distibution of strength, 3 - all of them are the same strength
#params: distribution parameters
genSynapses<-function(no,distop=1,params=c(1.15,0.16),op=1){
  #cat('\ngenSynapses no',no,' : ',distop)
  if(no==0) return(0) #no synapses, so return 0.
  if(distop==1) strvec <- rgamma(no,shape=params[1],scale=params[2]) #gamma
  if(distop==2) strvec <- runif(no,min = params[1],max = params[2]) #uniform
  if(distop==3) strvec <- rep(params[1],no) #the same strenth no times
  strvec 
}

#updates the connection matrix according to the synapses that each neuron should actually have
#as in if some neurons have less than the max. no of connections, you should update that.
#mat: tjhe connection matrix
#nosyn: an array specifying the number of synapses that every source neuron makes
updateSynMat <- function(mat,nosyn,op=1){
  res <- matrix(rep(0,length(mat)),nrow = nrow(mat))
  norows <- nrow(mat)
  for (i in 1:length(nosyn)) {
    #for each column if norow > nosyn[i], turn the rest to 0
    #cat('\n',nosyn[i])
    res[1:nosyn[i],i] <- mat[1:nosyn[i],i]
  }
  #fix the 1:0 still returning the first element issue
  zeropos <- which(nosyn==0)
  res[1,zeropos] <- 0
  res
}

#gets the number of synapses in this connection matrix from source to target. If there are 
#multiple sources, gets a list
#connmat: the connection matrix
getConnSynapseNo <- function(connMat,op=1){
  nosources <- length(connMat)
  #figure out the type of matrix: is it a list or a matrix
  if(isDataType(connMat)==const.DataType$list){#its a list
    notargets <- ncol(connMat[[1]])
    syn.ls <- lapply(connMat, function(x){
      sapply(1:notargets,function(i) length(which(x[,i]>0)))
    })
  }
  else syn.ls <- nrow(connMat) #each column is a source so its length gives # synapses
  syn.ls
}



#populateConnMatrixIndSyn <-function(source,target,params=list(4,c(1,1.15,0.16)),op=1){

#generates a matrix with connections across regions. e.g., between OB and PCx
#rows - the rows, for eg each one might specify a pcx neuron, the target neurons
#cols - the cols e.g each one might specify a glomerulus, the source neurons
#params - , no of synapses, and paramss of the distribution
#for type = 3 and 5, it is c(avg no of synapses made with each KC-MBON synapse,strength of the synapse = gaussian(2,0.5))
#type = 3 - gamma, 5 - uniform
genLongRangeConnMatrix <- function(rows,cols,params,type=3,op=1){
  #cat('\nparams',rows,cols,params[1])
  synno <- rpois(rows*cols,params[1]) #no of synapses for all connections
  #generate 1.2 times total synapse no of gamma distributions, and then use them
  #although since the no of synapses is from the mean, we shouldnt need them
  #old, mistakes when it generates fewer synapses: strvec <- rgamma(rows*cols*params[1]*1.2,shape=params[2],scale=params[3])
  if(type == 3) strvec <- rgamma(sum(synno)*1.2,shape=params[2],scale=params[3]) #gamma
  if(type == 5) strvec <- runif(sum(synno)*1.2,min = params[2],max = params[3]) #uniform
  index <- cumsum(synno) #the cumulative sum which acts as an index
  zero <- synno #multiplier to turn 0 synapses to 0
  zero[zero > 0] <- 1
  #cat(index,'\n',synno,'\nzero',zero,'\n')
  mat <- sapply(1:(rows*cols),function(i) {
    #strength <- sapply(0:syn,function(x) (x>0)*rgamma(1,shape=params[2],scale=params[3]))
    #start and end give the indices of the strvec that need to be used
    start <- (index[i]-synno[i]+1)*zero[i] 
    end <- index[i]*zero[i]
    #cat(i,'st ',start,',',' end',end,'\t',sum(start,end),' ;')
    strength <- sum(strvec[start:end])
    strength
  })
  #cat('\n',mat)
  mat <- matrix(mat,nrow = rows,byrow = F) #fills the matrix by columnns
  
  
}

#generates a matrix where the connections between source and target depend on their spatial positions. Currently, the dependence 
#on spatial position depends on a gaussian but it can be changed to be linear for instance. 
#rows - the rows, for eg each one might specify a pcx neuron, the target neurons
#cols - the cols e.g each one might specify a glomerulus, the source neurons
#for type = 20, params = list(#connections each source neuron makes,number of neurons per spatial position for the source, 
#for thetarget, l - the width of 95 % of the spread, L - the total width ofthe region, synparams - 
#the distribution for synaptic strengths)
genSpatialConnMatrix <- function(rows,cols,params,op=1){
  no.conn <- params[[1]]    
  nu.s <- params[[2]]
  nu.t <- params[[3]]
  l <- params[[4]]
  L <- params[[5]]
  synparams <- params[[6]]
  #the assumption is that the index of a neuron gives its spatial location from 1...L
  #src and tar neurons are spread the region, the granularity is determined by nu
  src <- cols
  tar <- rows
  #now, calculate the unit of distance for src and target neurons, i.e., the granularity
  src.posns <- src/nu.s
  tar.posns <- tar/nu.t
  L.s <- L/src.posns
  L.t <- L/tar.posns
  #cat('\nconnMat20:L.s',L.s,L,l)
  res.tar <- sapply(1:src,function(i){#each src neuron
    src.posn <- ceiling(i/nu.s)
    src.mn <- (src.posn-0.5)*L.s # the middle point or mean of where the neuron is 
    res.src <- sapply(1:tar,function(j){#go through each target neuron
      tar.pos <- ceiling(j/nu.t)
      tar.interval <- c((tar.pos-1)*L.t,tar.pos*L.t)
      pos.prob <- discreteNormal(x=tar.interval,mn = src.mn,width = l,op=1)
      #cat('\ngenConn20: mean',src.mn,'prob',pos.prob,'target pos',tar.pos)
      pos.prob*no.conn #the probability has to weighted by the number of connections that the cell makes
    })
  })
  colnames(res.tar) <- round(rep(c(1:src.posns)*L.s,each=nu.s),2)
  rownames(res.tar) <- round(rep(c(1:tar.posns)*L.t,each=nu.t),2)
  mat <- populateConnMatrix(connmat = res.tar,synparams.ls = synparams,op=1)
  mat
}

#given a mtrix of synaptic connections, with each entry specifying the number of average connections, it turns it into a synaptic weight
#matrix
#connmat: the connection matrix, with each entry specifying the mean no of connections paramter for the poisson
#synparams: synaptic params matrix for generating the weight for each synapse
#op=1 - gamma distibuteion of strngth, 2 - unform distibution of strength
populateConnMatrix <- function(connmat,synparams.ls=list(1,1.15,.16),op=1){
  if(isDataType(synparams.ls)>2) synparams <- unlist(synparams.ls)
  else synparams <- synparams.ls
  rows <- nrow(mat)
  cols <- ncol(mat)
  #get the number of synapses made by src and target 
  synno <- sapply(connmat,function(i) rpois(1,i) )
  posns <- which(synno>0) #find out all posns with non-zero synapses
  if(op==1) strvec <- rgamma(sum(synno)*1.2,shape=synparams[2],scale=synparams[3]) #gamma
  if(op==2) strvec <- runif(sum(synno)*1.2,min = params[2],max = params[3]) #uniform
  index <- cumsum(synno[posns]) #the cumulative sum which acts as an index
  zero <- synno[posns] #multiplier to turn 0 synapses to 0
  zero[zero > 0] <- 1
  #cat(index,'\n',synno,'\nzero',zero,'\n')
  synvec <- sapply(1:length(posns),function(i) {
    #strength <- sapply(0:syn,function(x) (x>0)*rgamma(1,shape=params[2],scale=params[3]))
    #start and end give the indices of the strvec that need to be used
    start <- (index[i]-synno[posns[i]]+1)*zero[i] 
    end <- index[i]*zero[i]
    #cat(i,'st ',start,',',' end',end,'\t',sum(start,end),' ;')
    strength <- sum(strvec[start:end])
    strength
  })
  #cat('\npCM',synvec,':',posns)
  mat <- connmat
  mat[mat>=0] <- 0
  mat[posns] <- synvec #fills the matrix by columnns
  mat
}



#function for generating topographic circuit connectivity
#source: the structure from which the connections are coming in, e.g., OB
#target: the structure to which the connections are going, e.g., PCx
#syndist - the synaptic distribution; 1- agaussian matrix, 2 -gamma distribution, 3 - compound poisson Gamma distribution
# 4 - uniform distribution
#if compound poisson gamma params = (poisson mean,alpha,beta)
#params - params of the distribution
#noise: the noise to be added, it is gaussian mean 0 with stdev given by noise
#noise, e.g. noise = .1 is 10 % noise
genConnMatTop <- function(source,target,syndist=1,params=c(0,0),noise=0,op=1){
  #calculate the number of targets that a source neuron contacts
  source_target <- size[2]/size[1]
  #generate matrix, with rows as the target and cols are the source
  #lets just do the source as rows first and then transpose in the end
  connmat.ls <- lapply(1:size[1], function(x){
    #figure out which neurons this should connect to
    if (source_target > 1) target_neurons <- c(rep(0,(x-1)*source_target),rep(1,source_target),rep(0,(size[2]/source_target-x)*source_target))
    else {
      #cat(x,source_target,ceiling(x*source_target),x*source_target,'\n')
      target_neurons <- c(rep(0,size[2]))
      target_neurons[ceiling(x*source_target)] <- 1
    }
    target_neurons
  })
  connmat <- convertNestedListsDF(connmat.ls)
  connmat
}

#this function draws the connections from PN to MB according to Stevens and 
#glomnopars: gives the number of glomeruli that each KC makes a connection with.
#glompars: gives the 
#synpats: The parameters for the synaptic strengths 
#synpars[1]: synaptic option, 1 - uniform distribution based on gruntman, 2 - gamma distribution based on chuck pnas supplement fit
#synpars[2:3]: the parameters of the distribution
genFlyPNMBconn <- function(rows,cols,glomnopars=c(8,0.715),glompars=c(0.07),synpars=c(2,4,4),noise=0,op=1){
  #first choose the number of glomeruli for each row of KCs, it follows a binomial
  #distribution with n=8 for gamma cells and 11 for non-gamma cells and p of adding a claw=0.715
  #cat('\bgenFly',rows,cols)
  noglom <- sapply(1:rows, function(x) rbinom(1,size=glomnopars[1],prob = glomnopars[2]))
  #now, to assign these to glomeruli
  geom_cum <- genCumulativeVec(1:cols,dist = 6,params = glompars[1]) #generate the geom cumulative list
  #cat('\nnoglom',noglom,'\t',synpars,'\n',geom_cum,'\n')
  
  kc <- lapply(noglom, function(x){
    no <- runif(x,min = 0,max = 1) #generate a random number b/w 0 and 1, and use cumulative list to pick a glomerulus
    glom <- sapply(1:x, function(y) computeValueIndex(no[y],geom_cum))
    #now, assign synaptic strengths to these
    #glom_strength <- sapply(1:length(glom), function(y) rgamma(1,shape = synpars[2],scale = synpars[3]))
    glom_strength <- switch(synpars[1],sapply(1:length(glom), function(y) runif(1,min = synpars[2],max = synpars[3])),
                            sapply(1:length(glom), function(y) rgamma(1,shape = synpars[2],scale = synpars[3])) )
    #cat('\n',no,glom,', g.strength ',glom_strength)
    res <- createVecPosns(cols,posns = glom,vals = glom_strength,op=2)
  })
  #add noise, maybe
  res <- as.matrix(transposeDF(convertNestedListsDF(kc)))
  colnames(res) <- 1:cols
  rownames(res) <- 1:rows
  res
}

#dop=c() - no dopamine, not a learning matrix, list(no of dopamine neurons, params for the network connections,dopamine type=reward,punish)
# 1 for reward, -1 for punish


#--------------------------------------------------------------------------------------------------------------
#now we can do specific circuits: all of them like classes




newPNMBcircuit <- function(pnno=50,mbno=2000,op=1){
  pn <- newNeurons(number = 50,dist = 2) #projection neuron firing
  pnmb_conn <- newConnMatrix(c*c(pnno,mbno))
  mbkc <- newNeurons(number = 50,dist = 0,params = c(0))
  kcapl_conn <- newConnMatrix(c(mbno,1),type = 1)
  apl <- newNeurons(number = 1,dist = 0,params = c(0))
  aplkc_conn <- newConnMatrix(c(1,mbno),type = 1)
  dans <- newNeurons(number = 2,dist = 0,params = c(1))
  mbons <- newNeurons(number = 2,dist = 0,params = c(0))
  
  pnmb <- list(pn,mbkc,apl,dans,mbons,pnmb_conn,kcapl_conn,aplkc_conn)
  names(pnmb) <- c('pn','mbkc','apl','dans','mbons','pnmb_conn','kcapl_conn','aplkc_conn')
  res <- structure(pnmb,class = c('FlyMBp','circuit')) #make it a class and return it
  res
}



#------------------------------------------
#function to be transferred later

#statsfn
#for a distribution, any distribution, brek the distribution into n intervals and return the 
#the points within the distribution that form these intervals
#distvec: the vector
#n: the values of the intervals
#
computeCumulativeList <-function(distvec,n,op=1){
  cumlist <- RelCumFreq(vec = distvec)
  len <- nrow(cumlist) #get the smaller of the two, length or number of n
  if(n > len) no_indices <- len
  else no_indices <- n
  lastelem <- cumlist[len] #get the last elem
  indices <- seq(0,lastelem,lastelem/no_indices) #compute the indices
  res_indices <- sapply(indices,function(x) computeValueIndex(x,cumlist[,1]))
  res <- cumlist[res_indices,2] 
  res
}

#statsfn
#for a distribution, any distribution, get the probabilities of the points given by 
#vec and then normalize to 1 so that you get a proportionate or cumulative distribution
#of these probabilities
#dist: the distributuion for which we should get the intervals
#1 - gaussian, 2 - gamma, 3 -lognormal, 4 - exponential, 5 - uniform, 6 - geometric
#n: the values of the intervals
#op=1, normalize, 2 - do not normalize
genCumulativeVec <-function(vec,dist=1,params=c(),op=1){
  #choose a distribution
  fnx <- function(x) switch(dist,do.call(pnorm,as.list(c(q=x,c(params)))),
                            do.call(pgamma,as.list(c(q=x,shape=params[1],scale=params[2]))),
                            do.call(plnorm,as.list(c(q=x,c(params)))),
                            do.call(pexp,as.list(c(q=x,params[1]))),
                            do.call(punif,as.list(c(q=x,min=params[1],max=params[2]))),
                            do.call(pgeom,as.list(c(q=x,prob=params[1]))))
  
  #get the cumulative probabilities for the numbers in vec
  res <- sapply(vec,fnx) 
  res <- switch(op,res/fnx(vec[length(vec)]),res) #normalize them to 1 or not
  res
}  

#misc_fns
#gets the highest index <= value in the vector
computeValueIndex <- function(x,vec,op=1){
  res <- length(which(vec<=x))
  res+1
}



#----------------------------------------------------------------------------------
circuitfntesting <-function(){
  #creating a fly circuit
  tmp <- genFlyPNMBconn(2000,50) #generating a fly circuit
  
  
}

#result: with the fly circuit the connectivity matrix is a mix of two uniform distributions. One that includes one synapse, and 
#a second smaller one that includes 2 synapses. 
#
#OLD functions
#----------------------------------------------------------------------------------
#origin: the origin neuron population object, the output size is automatically fixed to the size of the origin
#connmat: the connection matrix objext
#self this population of neurons object
#dest: is the destination or target neuron structure
#stimno: odor or list of odor nos for which this should be coomputed
#op=1, do all stimuli, 2 - do stimuli specigied by odor no
#noiseop: specifies whether we should add noise (T) or not (F) to the calculation: 
computeFiring.Neurons.old <- function(self,connMat,dest=c(),stimno=1,op=1){
  #do a check to confirm that the matrices are correct. So, the orign and dest neurons desc should match
  #the same fields in the matrix
  origin <- self
  #cat('\nCF.neurons',str(origin))
  if ( size(connMat)[1] != size(origin,op=2) ){
    cat('\n Conn Matrix and neuron populations do not match',size(connMat)[1],',',size(origin,op=2))
    return(F)
  }
  #compute 2 things: the inputs and connection matrix, and add to the neuron firing to factor in the basal noise.
  #but, you should also take care because if you keep adding, you'll get a gigantic number
  conn <- getData(connMat,const.CM$connmat)
  neurons <- getData(origin,const.N$neurons)
  if(op==1){#do all stimuli
    #cat('\nCFN',str(conn),str(neurons))
    res.neurons <- lapply(1:size(origin), function(i){
      neurons <- as.vector(conn %*% neurons[[i]]) #generate the destination vector
    })
  }
  if(op==2){#just do one particular stimulus
    if(length(dest)>0) res.neurons <- getData(dest,const.N$neurons)
    else res.neurons <- neurons #if there is a destination, that should be the base here
    res.neurons[stimno] <- lapply(stimno, function(i){
      res <- as.vector(conn %*% neurons[[i]]) #generate the destination vector
    })
  }
  names(res.neurons) <- 1:size(origin)
  #now update the this current neuron population
  if (length(dest)>0) target <- setData(dest,val=res.neurons,index=const.N$neurons)
  else target <- setData(origin,val=res.neurons,index=const.N$neurons)
  target
}

addNoise.FlyMBg.OLD <-function(self,odors=c(),noiseset=0,op=1){
  #cat('adding noise to FlyMBg')
  if(op==0) return(self) # no noise to add, just return self
  circuit <- self
  #designate odors for noise, default is all odors, if noiseset > 0 and odors is empty either testing or training sets
  if(length(odors)==0){
    if(noiseset > 0){ #just the set, usually training set
      odornos <- switch(noiseset,getData(self,const.FM$origstim)[,2],getData(self,const.FM$trainingval)[,2],
                        getData(self,const.FM$teststim)[,2])
    }
    else odornos <- c(getData(self,const.FM$origstim)[,2],getData(self,const.FM$trainingval)[,2],
                      getData(self,const.FM$teststim)[,2]) #basically all odors
  }
  else odornos <- odors
  #choose the components to which we will be adding noise
  if(op>=1 && op<=4){#only one component or structure for noise addition
    compid <- switch(op,const.FM$glom,const.FM$glommb_conn,const.FM$kcapl_conn,const.FM$aplkc_conn,
                     const.FM$dans,const.FM$kcmbon_conn)
    flycomp <- getData(circuit,compid)
    flycomp.type <- getCompType(flycomp,op=2)
    if(flycomp.type == 2) circuit <- setData(circuit,val=addNoise(flycomp),compids[i]) #matrix
    else circuit <- setData(circuit,val=addNoise(flycomp,odors=odornos),compids[i]) # other things
    
    #circuit <- setData(circuit,val=addNoise(flycomp,odors=odornos),compid)
  }
  if(op>=5 && op<=7){#includes the valence based structures
    compid <- switch((op-4),const.FM$kcmbon_conn,const.FM$mbons,const.FM$dans)
    #got through each component, then their valences, and set noise
    flycomp <- getData(circuit,compid)
    for(i in getData(circuit,const.FM$valences)){
      flycomp.type <- getCompType(flycomp[[i]],op=2)
      if(flycomp.type == 2) flycomp[[i]] <- addNoise(flycomp[[i]]) #matrix
      else flycomp[[i]] <- addNoise(flycomp[[i]],odors=odornos) # other things
    }
    circuit <- setData(circuit,val=flycomp,compid)
  }
  if(op>=10 && op <= 11) {#multiple structures are being changed
    compids <- switch((op-9),c(const.FM$glom,const.FM$glommb_conn),
                      c(const.FM$glom,const.FM$glommb_conn,const.FM$aplkc_conn))
    for(i in seq_wrap(1,length(compids))){
      flycomp <- getData(circuit,compids[i])
      flycomp.type <- getCompType(flycomp,op=2)
      if(flycomp.type == 2) circuit <- setData(circuit,val=addNoise(flycomp),compids[i]) #matrix
      else circuit <- setData(circuit,val=addNoise(flycomp,odors=odornos),compids[i]) # other things
    }
  }
  if(op==100){#add noise to all structures
    compids <- c(const.FM$glom,const.FM$glommb_conn,const.FM$aplkc_conn,const.FM$apl,const.FM$kcapl_conn,
                 const.FM$mbkc,const.FM$mbkc_net)
    for(i in seq_wrap(1,length(compids))){
      flycomp <- getData(circuit,compids[i])
      flycomp.type <- getCompType(flycomp,op=2)
      if(flycomp.type == 2) circuit <- setData(circuit,val=addNoise(flycomp),compids[i]) #matrix
      else circuit <- setData(circuit,val=addNoise(flycomp,odors=odornos),compids[i]) # other things
    }
    #now, set the valence related bits
    compids <- c(const.FM$kcmbon_conn,const.FM$mbons,const.FM$dans)
    #got through each component, then their valences, and set noise
    for(j in seq_wrap(1,length(compids))){
      flycomp <- getData(circuit,compids[j])
      for(i in circuit$valences){
        flycomp.type <- getCompType(flycomp[[i]],op=2)
        if(flycomp.type == 2) flycomp[[i]] <- addNoise(flycomp[[i]]) #matrix
        else flycomp[[i]] <- addNoise(flycomp[[i]],odors=odornos) # other things
      }
      circuit <- setData(circuit,val=flycomp,compids[j])
    }    
  }
  circuit
}


