######################################################################
#### Example 2 (Sec 3.2)
######################################################################

################## BUD design ##################
# interant of the entropy of the posterior of the response rate of the most effective treatment
diff_ent_best_arm = function(t, xs, ns){
  # xs = vector of treatment sucessed by treatment arm
  # ns = sample sizes by treatment arm 
 A = length(xs)
 
  sapply(t,function(x){
    ss  = sum(sapply(1:A,function(a) dbeta(x,1+xs[a],1+ns[a]-xs[a])* prod(sapply(1:A,function(i) pbeta(x,1+xs[i],1+ns[i]-xs[i]))[-a]) ))  
    ss1 = ss*log(ss)
  
    if(is.na(ss1))  
      ss1 = 0
  
  return(-ss1)})
}  

# Simulate single BUD trial using entropy "diff_ent_best_arm()"
BUD_trial=function(pp, N, pow=1){
 
  A = length(pp)
  pv=hh=xx=ss=nn1=xx1=rr1=PPP=HF=c()
  nn=xr=rep(0,A)
  
  for(n in 1:N){
   
    # compute utility
    Pr  = prob_best_arm(xr,nn)
    H   = integrate(diff_ent_best_arm, 0,1, xs=xr, ns=nn)$value
    PPP = rbind(PPP,Pr)
    HF   = c(HF,H)
    
    # compute expected utility
    x1   = matrix(rep(xr,A),A,A,byrow = T)+diag(1,A,A)
    nn12 = matrix(rep(nn,A),A,A,byrow = T)+diag(1,A,A)
    HF1  = c()
    
    for(o in 1:A){  
      H1   = integrate(diff_ent_best_arm,0,1,xs=x1[o,],ns=nn12[o,])$value
      H2   = integrate(diff_ent_best_arm,0,1,xs=xr,ns=nn12[o,])$value
      succ = (1+xr[o])/(2+nn[o])
      HH   = succ*H1+(1-succ)*H2
      HF1  = c(HF1,HH)
    }
    
    # compute randomization probability
                     rand                = H-HF1
    if(any(rand<0))  rand[which(rand<0)] = 0
    if(all(rand==0)) rand                = rep(1,ncol(DATA))
                     RAND                = (rand^pow)/sum(rand^pow)  
                     rr1                 = rbind(rr1,RAND)
                     
    # randomize patient                     
    d     = sample(A,1,prob = RAND)
    
    # update data
    nn[d] = nn[d]+1
    xr[d] = xr[d]+rbinom(1,1,pp[d])
    nn1   = rbind(nn1,nn)
    xx1   = rbind(xx1,xr)
  }
  return(list(Responses=xx1,SampleSize=nn1,Posterior=PPP,Utility=-HF,rand=rr1))
}



################## BAR design ##################
### Pr(theta.a >= theta.j  for all j)
prob_best_arm=function(x, n, alpha=1, beta=1){
  
  k   = length(x)
  ans = numeric(k)
  
  for(i in 1:k){
    indx   = (1:k)[-i]
    f      = function(z){ r =     dbeta(z, x[i] + alpha, n[i] - x[i] + beta)
    for(j in indx)  r = r * pbeta(z, x[j] + alpha, n[j] - x[j] + beta)
    r }
    ans[i] = integrate(f, 0, 1)$value
  }
  return(ans)
}


# Simulate single trial using Thomson sampling
BAR_trial=function(pp,N,pow=1){
 
  A = length(pp)
  nn1=xx1=rr1=PPP=HF=c()
  nn=xr=rep(0,A)
  
  for(n in 1:N){
    # compute posterior probabiliy
    prob = prob_best_arm(xr, nn)
    rand = prob^(pow)/sum(prob^(pow))
    PPP  = rbind(PPP,prob)
    rr1  = rbind(rr1,rand)
    
    # utility
    H    = integrate(diff_ent_best_arm,0,1,xs=xr,ns=nn)$value
    HF   = c(HF,H)
    
    # patient assignment
    d=sample(length(pp),1,prob = rand)
    
    # update data
    nn[d]=nn[d]+1
    xr[d]=xr[d]+rbinom(1,1,pp[d])
    nn1=rbind(nn1,nn)
    xx1=rbind(xx1,xr)
  }
  return(list(Responses=xx1,SampleSize=nn1,Posterior=PPP,Utility=-HF,rand=rr1))
}



################## BR design ##################
# Simulate single trial using balanced randomization
BR_trial=function(pp,N){
 
  A = length(pp)
  nn1=xx1=HF=c()
  nn=xr=rep(0,A)

 for(n in 1:N){
     
    # compute utility
    HF   = c(HF,integrate(diff_ent_best_arm,0,1,xs=xr,ns=nn)$value)
    
    # assignment
    d      = sample(1:4, 1, prob = rep(1/A, A))
    
    # update data
    nn[d] = nn[d]+1
    xr[d] = xr[d]+rbinom(1,1,pp[d])
    nn1   = rbind(nn1,nn)
    xx1   = rbind(xx1,xr)
 }
  
  return(list(Responses=xx1,SampleSize=nn1,Utility=-HF))
}



################## RPW design ##################
# Simulate single trial using randomized play the winner rule
RPW_trial=function(pp, N, nball=1){
  
  A = length(pp)
  nn1=xx1=PPP=HF=c()
  nn= xr = rep(0,A)
  URN=table(seq(1,A))
  
  for(n in 1:N){
    
    # cumpute utility
    prob = prob_best_arm(xr,nn)
    PPP  = rbind(PPP,prob)
    H    = integrate(diff_ent_best_arm,0,1,xs=xr,ns=nn)$value
    HF   = c(HF,H)
    
    # assignment
    rand = URN/sum(URN)
    d    = sample(1:A,1,prob=rand)
    
    # update data
    nn[d] = nn[d]+1
    ris   = rbinom(1,1,pp[d])
    xr[d] = xr[d]+ris
    nn1   = rbind(nn1,nn)
    xx1   = rbind(xx1,xr)
    
    if(ris==1){ URN[d]=URN[d]+nball
    }else      URN[-d]=URN[-d]+nball
    
  }
  return(list(Responses=xx1, SampleSize=nn1, Posterior=PPP, Utility=-HF))
}


