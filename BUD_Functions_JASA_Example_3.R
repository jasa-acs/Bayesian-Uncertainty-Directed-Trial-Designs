EXAMPLE=F
###############################################################################
####        create parameter matrix THETA            ##########################
####        parameter: list all parameter            ##########################
####        num.biom: number of biomarker            ##########################
####        num.drug: number of drugs                ##########################
###############################################################################
create.THETA=function(p0,TE_pos,TE_neg,target=NA){
 
  num.drug = length(TE_pos)+1
  if(is.na(target) || length(target)!=(num.drug-1)) target=seq(1, num.biom)
  num.biom=length(unique(target))
  
  list.sub.group=function(num.biom){
    XN=matrix(rep(0:1,num.biom),nrow=2,ncol=num.biom)
    s=do.call('expand.grid',as.data.frame(XN))
    colnames(s)=paste0('bm',1:num.biom)
    return(s) }
  
  sub.group = list.sub.group(num.biom) #list of subgroups
  parameter = matrix(p0,nrow(sub.group),length(TE_pos)+1)
  
  for(g in 2:ncol(parameter)){ 
      parameter[which(sub.group[,target[g-1]]==1),g]=parameter[which(sub.group[,target[g-1]]==1),g]+TE_pos[g-1]
      parameter[which(sub.group[,target[g-1]]==0),g]=parameter[which(sub.group[,target[g-1]]==0),g]+TE_neg[g-1]
  }
  
  THETA           = matrix(parameter,nrow(sub.group),num.drug)
  rownames(THETA) = (sapply(1:nrow(sub.group),function(i) paste(sub.group[i,],collapse = "")))
  colnames(THETA) = c('control',paste0('drug ',1:(num.drug-1)))
  
  return(THETA)
}

###############################################################################
###            generate the matrix of conterfactual data       ################
###      corresponding to the particular patient's profile     ################
###############################################################################
Gen.Biomarker.data=function(THETA,BIOM_prob,N){
  b=sapply(1:length(BIOM_prob),function(i) rbinom(n = N,size = 1,prob = BIOM_prob[i]))
  b=sapply(1:nrow(b),function(i) paste(b[i,],collapse = ""))
  match_sub_gr=sapply(1:length(b),function(i) which(b[i]==rownames(THETA)))
  OUTCOME=t(sapply(1:N,function(n){
    sapply(1:ncol(THETA),function(t) rbinom(1,1, THETA[match_sub_gr[n],t]))
  }))
  colnames(OUTCOME)=colnames(THETA)
  DATA=data.frame(profile=b,randomization=rep(NA,N),OUTCOME)
  return(DATA)
}

###############################################################################
###           Summarize sufficient statistics       ###########################
###############################################################################
Suff.stat=function(DATA){
  Events=array(0,dim=c(ncol(DATA)-2,2,1))
  dimnames(Events)=list(colnames(DATA)[-(1:2)],c('success','failure'))
  Events[,1,1]=sapply(1:nrow(Events),function(a) sum(DATA[which(DATA[,2]==a),2+a]))
  Events[,2,1]=sapply(1:nrow(Events),function(a) length(DATA[which(DATA[,2]==a),2+a]))-Events[,1,1]
  return(Events)
}

###############################################################################
###           Summarize sufficient statistics by    ###########################
###               subgroup                          ###########################
###############################################################################
Suff.stat.subgroup=function(DATA,THETA){
  Events=array(0 ,dim = c(2,ncol(DATA)-2,nrow(THETA)), dimnames=list(c('success','failure'),colnames(THETA),row.names(THETA)))
  for(i in dimnames(Events)[[3]]){
    Events[1,,i]=sapply(1:ncol(THETA),function(a) sum(DATA[which(DATA[,2]==a & DATA[,1]==i),2+a]))
    Events[2,,i]=sapply(1:ncol(THETA),function(a) length(DATA[which(DATA[,2]==a & DATA[,1]==i),2+a]))-Events[1,,i]
  }
  return(Events)
}

###############################################################################
###           Summarize sufficient statistics by          #####################
###            positive and negative biomarker for        #####################
###            target populations                         #####################
###         drug 1 biomarker 1, drug 2 biomarker 2,...    #####################
###############################################################################
Suff.stat.target=function(EVENTS,target){
  events_pos_neg=array(0,dim=c(2,4,dim(EVENTS)[[2]]-1))
  dimnames(events_pos_neg)=list(c('experimental','control'),
                                c('success +','failure +','success -','failure -'),
                                dimnames(EVENTS)[[2]][-1])
  for(b in 1:(dim(EVENTS)[[2]]-1))
  {
    events_pos_neg[2:1,1,b]=rowSums(EVENTS[1,c(1,b+1),which(substring(dimnames(EVENTS)[[3]],target[b],target[b])==1)])
    events_pos_neg[2:1,2,b]=rowSums(EVENTS[2,c(1,b+1),which(substring(dimnames(EVENTS)[[3]],target[b],target[b])==1)])
    events_pos_neg[2:1,3,b]=rowSums(EVENTS[1,c(1,b+1),which(substring(dimnames(EVENTS)[[3]],target[b],target[b])==0)])
    events_pos_neg[2:1,4,b]=rowSums(EVENTS[2,c(1,b+1),which(substring(dimnames(EVENTS)[[3]],target[b],target[b])==0)])
  }
  return(events_pos_neg)
}

###############################################################################
###          Computation of the probability of a             ##################
###            Beta(a,b) random variable to be larger        ##################
###    than another independent Beta(a’,b’) random variable. ##################
###############################################################################
betacomparison = function(E,ab.C=c(1,1),ab.E=c(1,1)){
  C_a = ab.C[1]+E[1,1]
  C_b = ab.C[2]+E[1,2]
  E_a = ab.E[1]+E[2,1]
  E_b = ab.E[2]+E[2,2]
  f = function(p.Y)  (1-pbeta(p.Y,E_a,E_b))*dbeta(p.Y,C_a,C_b)
  integrate(f,0,1)$value
}

###############################################################################
###          Posterior probabilities with the described              ##########
### model of treatment effects (i) in biomarker positive and         ##########
### (ii) in the biomarker negative subpopulations, assuming  a=b=1.  ##########
###             prior ={ p(E+),p(E- | E+),p(E- | I+)) }.             ##########
###############################################################################
comparisonPositiveNegative<-function(E,ab.C=c(1,1),ab.E=c(1,1),prior=c(.5,.5,.01)){
  a1=betacomparison(E[2:1,1:2],ab.C,ab.E)
  a2=betacomparison(E[2:1,3:4],ab.C,ab.E)
  aa=c(a1*a2, 
  		 a1*(1-a2), 
  		 (1-a1)*a2,(1-a1)*(1-a2))*c(prior[1]*prior[2], 
  		 prior[1]*(1-prior[2]),(1-prior[1])*prior[3],(1-prior[1])*(1-prior[3]))
  aa=aa/sum(aa)
  return(c(pos=aa[1]+aa[2], neg=aa[1]+aa[3], posneg=1-aa[4]))
}


###############################################################################
####  Predictive probabilities with the described model of            #########
###   positive outcome for the next patient, assuming the patient     #########
###   will be assigned to the experimental arm,                       #########
###   (i) in biomarker positive group, and                            #########
###   (ii) in the biomarker negative group, assuming a=b=1.           #########
###############################################################################
prediction<-function(E,prior=c(.5,.5,.01)){
  EE=E
E=EE[,1:2]
xx=rep(0,sum(E[2,])+2)
for(i in 1:(sum(E[2,])+2)){
  x=log(1/(sum(E[2,])+2))
  if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
  if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
  if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
  xx[i]=x}
if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
xx=xx+log(sum(E[1,])+1)
xx1=exp(xx)


E=EE[,3:4]
xx=rep(0,sum(E[2,])+2)
for(i in 1:(sum(E[2,])+2)){
  x=log(1/(sum(E[2,])+2))
  if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
  if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
  if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
  xx[i]=x}
if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
xx=xx+log(sum(E[1,])+1)
xx2=exp(xx)


v1=c(1:(sum(EE[2,1:2])+2))%*%t(rep(1,(sum(EE[2,3:4])+2)))
v2=t(c(1:(sum(EE[2,3:4])+2))%*%t(rep(1,(sum(EE[2,1:2])+2))))
x=xx1%*%t(xx2)

 x=x*(v1>(EE[2,1]+1))*(v2>(EE[2,3]+1))*prior[1]*prior[2]+
   x*(v1>(EE[2,1]+1))*(v2<=(EE[2,3]+1))*prior[1]*(1-prior[2])+
   x*(v1<=(EE[2,1]+1))*(v2>(EE[2,3]+1))*(1-prior[1])*prior[3]+
   x*(v1<=(EE[2,1]+1))*(v2<=(EE[2,3]+1))*(1-prior[1])*(1-prior[3])
x=x/sum(x)
a1=sum(x*(c((1+EE[1,1]):(1+EE[1,1]+EE[2,1]+EE[2,2]+1))/((EE[1,1]+EE[1,2]+EE[2,1]+EE[2,2]+3))))
a2=sum(t(x)*(c((1+EE[1,1+2]):(1+EE[1,1+2]+EE[2,1+2]+EE[2,2+2]+1))/((EE[1,1+2]+EE[1,2+2]+EE[2,1+2]+EE[2,2+2]+3))))

return(c(a1,a2))}
 
 #########

prediction_control<-function(E,prior=c(.5,.5,.01)){EE=E; EE[1,]=E[2,];EE[2,]=E[1,];E=EE

E=EE[,1:2]
xx=rep(0,sum(E[2,])+2)
for(i in 1:(sum(E[2,])+2)){
  x=log(1/(sum(E[2,])+2))
  if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
  if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
  if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
  xx[i]=x}
if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
xx=xx+log(sum(E[1,])+1)
xx1=exp(xx)


E=EE[,3:4]
xx=rep(0,sum(E[2,])+2)
for(i in 1:(sum(E[2,])+2)){
  x=log(1/(sum(E[2,])+2))
  if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
  if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
  if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
  xx[i]=x}
if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
xx=xx+log(sum(E[1,])+1)
xx2=exp(xx)

v1=c(1:(sum(EE[2,1:2])+2))%*%t(rep(1,(sum(EE[2,3:4])+2)))
v2=t(c(1:(sum(EE[2,3:4])+2))%*%t(rep(1,(sum(EE[2,1:2])+2))))
x=xx1%*%t(xx2)
x=x*(v1>(EE[2,1]+1))*(v2>(EE[2,3]+1))*prior[1]*(1-prior[3])+
   x*(v1>(EE[2,1]+1))*(v2<=(EE[2,3]+1))*prior[1]*(prior[3])+
   x*(v1<=(EE[2,1]+1))*(v2>(EE[2,3]+1))*(1-prior[1])*(1-prior[2])+
   x*(v1<=(EE[2,1]+1))*(v2<=(EE[2,3]+1))*(1-prior[1])*(1-prior[2])

x=x/sum(x)

a1=sum(x*(c((1+EE[1,1]):(1+EE[1,1]+EE[2,1]+EE[2,2]+1))/((EE[1,1]+EE[1,2]+EE[2,1]+EE[2,2]+3))))
a2=sum(t(x)*(c((1+EE[1,1+2]):(1+EE[1,1+2]+EE[2,1+2]+EE[2,2+2]+1))/((EE[1,1+2]+EE[1,2+2]+EE[2,1+2]+EE[2,2+2]+3))))

return(c(a1,a2))
}

###############################################################################
####            calculate the utility function                       #########
###############################################################################
Utility=function(prob, l=5,  c=5.7){
	
  if(is.matrix(prob)){   return( (l*(prob[3,]-prob[3,]^c)  + (prob[1,]-prob[1,]^c) + (prob[2,]-prob[2,]^c)))
  }else{                 return( (l*(prob[3]-prob[3]^c)    + (prob[1]-prob[1]^c)   + (prob[2]-prob[2]^c)))         }
}

###############################################################################
###############################################################################
####                Simulate a clinical trial                         #########
###############################################################################
###############################################################################

trial_Biomarker=function(DATA,THETA, l=5, c = 5.7, pow=2, target=NA){
	
                         U      = c()
  if(any(is.na(target))) target = seq(1,ncol(THETA)-1)
  
  for(n in 1:nrow(DATA)){
  	
    events_sub=Suff.stat.target(Suff.stat.subgroup(DATA,THETA),target)
    prob=sapply(1:dim(events_sub)[3],function(a) comparisonPositiveNegative(events_sub[,,a]))
    U=rbind(U,Utility(prob,l = l,c = c))
    v2=c()
    RAND=array(0,dim =dim(events_sub)[3]+1)
    
    for(f in 1:dim(events_sub)[3]){
      eve_succ_sub=eve_fail_sub=events_sub
      mark=as.numeric(substring(DATA[n,1],target[f],target[f]))
      eve_succ_sub[1,2*(1-mark)+1,f]=eve_succ_sub[1,2*(1-mark)+1,f]+1
      eve_fail_sub[1,2*(1-mark)+2,f]=eve_fail_sub[1,2*(1-mark)+2,f]+1
      pr_pred=prediction(events_sub[,,f])[1+(1-mark)]
      v2=c(v2,pr_pred)
      prob_s=sapply(1:dim(events_sub)[3],function(a) comparisonPositiveNegative(eve_succ_sub[,,a]))
      prob_f=sapply(1:dim(events_sub)[3],function(a) comparisonPositiveNegative(eve_fail_sub[,,a]))
      RAND[f+1]=(U[n,f]-(pr_pred*Utility(prob_s,l = l,c = c)[f]+(1-pr_pred)*Utility(prob_f,l = l,c = c)[f]))
    }
    
    for(f in 1:dim(events_sub)[3]){
      eve_succ_sub=eve_fail_sub=events_sub
      mark=as.numeric(substring(DATA[n,1],target[f],target[f]))
      eve_succ_sub[2,2*(1-mark)+1,f]=eve_succ_sub[2,2*(1-mark)+1,f]+1
      eve_fail_sub[2,2*(1-mark)+2,f]=eve_fail_sub[2,2*(1-mark)+2,f]+1
      pr_pred=prediction_control(events_sub[,,f])[1+(1-mark)]
      prob_s=sapply(1:dim(events_sub)[3],function(a) comparisonPositiveNegative(eve_succ_sub[,,a]))
      prob_f=sapply(1:dim(events_sub)[3],function(a) comparisonPositiveNegative(eve_fail_sub[,,a]))
      RAND[1]=RAND[1]+(U[n,f]-(pr_pred*Utility(prob_s,l = l,c = c)[f]+(1-pr_pred)*Utility(prob_f,l = l,c = c)[f]))
    }
    
                     RAND[which(RAND<=0)] = 0
    if(all(RAND==0)) RAND                 = array(1,dim =dim(events_sub)[3]+1)
                     RAND                 = (RAND^pow)/sum(RAND^pow)
                     DATA[n,2]            = sample(1:(dim(events_sub)[3]+1),1, prob=RAND)
  }
  return(list(DATA=DATA,Utility=U))
}


ztest=function(Y){
  p=Y[,1]/rowSums(Y)
  pp=(Y[1,1]+Y[2,1])/sum(Y)
  return(z=(p[1]-p[2])/sqrt(pp*(1-pp)*sum(1/rowSums(Y))))
}



#PSEUDO BOOTSTRAP
PSEUDO_BOOTSTRAP=function(scen, test.arm, boot_sim, test_z1, test_z2, N, pow){
  
  #test.arm= arm yoo want to be test (arm 1 is the control)
  #boot_sim= number of bootstrap simulation
  
  pos_boot           = pos[scen,]
  pos_boot[test.arm] = pos_boot[1]
  neg_boot           = neg[scen,]
  neg_boot[test.arm] = neg_boot[1]
  
  # simulate trials under the null
  boot_param_pos = create.THETA(p0=p.control,TE_pos=pos_boot,TE_neg=neg[scen,],target = BMK_targ[scen,])
  boot_data_neg  = lapply(1:boot_sim,function(j) Gen.Biomarker.data(THETA=boot_param_pos,BIOM_prob=BMK_prevalence[[scen]],N = N))
  boot_trial_pos = lapply(1:boot_sim,function(j) trial_Biomarker(DATA=boot_data_neg[[j]],THETA = boot_param_pos, target = BMK_targ[scen,], pow=pow))
  
  boot_z1=sapply(1:length(boot_trial_pos), function(i){
    AA=Suff.stat.target(target = BMK_targ[scen,],EVENTS = Suff.stat.subgroup(boot_trial_pos[[i]]$DATA, boot_param_pos))
    Z=sapply(1:4,function(a) ztest(Y=rbind(AA[1,1:2,a],AA[2,1:2,a])))
    return(Z)
  })
  z1_ho=apply(boot_z1,1,quantile,probs=0.9,na.rm = T)
  
  # simulate trials under the null  
  boot_param_neg = create.THETA(p0=p.control,TE_pos=pos[scen,],TE_neg=neg_boot,target = BMK_targ[scen,])
  boot_data_neg  = lapply(1:boot_sim,function(j) Gen.Biomarker.data(THETA=boot_param_neg,BIOM_prob=BMK_prevalence[[scen]],N = N))
  boot_trial_neg = lapply(1:boot_sim,function(j) trial_Biomarker(DATA=boot_data_neg[[j]],THETA = boot_param_neg, target = BMK_targ[scen,], pow = pow))
  
  boot_z2=sapply(1:length(boot_trial_neg), function(i){
    AA=Suff.stat.target(target = BMK_targ[scen,],EVENTS = Suff.stat.subgroup(boot_trial_neg[[i]]$DATA, boot_param_neg))
    Z=sapply(1:4,function(a) ztest(Y=rbind(AA[1,1:2,a],AA[2,1:2,a])))
    return(Z)
  })
  z2_ho=apply(boot_z2,1,quantile,probs=0.9,na.rm = T)
  

#POWER
  Power_z1 = colMeans(sapply(1:nrow(test_z1),function(i) test_z1[i,]>z1_ho[i]))
  Power_z2 = colMeans(sapply(1:nrow(test_z2),function(i) test_z2[i,]>z2_ho[i]))
  
  return(rbind(TEST_z1=Power_z1, TEST_z2=Power_z2)[,test.arm])
}


