################################################
#### Example 1: Multi-Arms trials 
################################################

################## BUD design ##################
# posterior variances of treatment effect 
varB = function(x, n, x1, n1){
  a  = x+1
  b  = 1+n-x
  a1 = x1+1
  b1 = 1+n1-x1
  
 return(a*b/((a+b)^2 *(a+b+1))+a1*b1/((a1+b1)^2 *(a1+b1+1)))
}

# Simulate single BUD trial using posterior variance `varB()``
BUD_trial=function(pp, N, pow=1){
  
  A = length(pp)
  pv=hh=xx=ss=nn1=xx1=rr1=PPP=HF=c()
  nn=xr=rep(0,A)
  HF= rep(0, n)
  
  for(n in 1:N){
    
    # posterior probability of positive treatment effect 
    Pr  = sapply(2:A,function(a) p.post.sup0(xr[1],nn[1],xr[a],nn[a]))
    PPP = rbind(PPP,Pr)
    
    # utility
    H   = sapply(2:A,function(a) (varB(xr[1],nn[1],xr[a],nn[a])))
    HF[n]  = sum(H)
    
    # expected utility
    x1   = matrix(rep(xr,A),A, A,byrow = T) +diag(1,A,A)
    nn12 = matrix(rep(nn,A),A, A,byrow = T) +diag(1,A,A)
    HF1  = c()
    
    for(i in 1:A){
      
      # all arm specific variances when next patient in arm i is responder (H1) or failuer (H2)
      H1   = sapply(2:A,function(a) varB(x1[i,1], nn12[i,1], x1[i,a], nn12[i,a]))
      H2   = sapply(2:A,function(a) varB(xr[1],   nn12[i,1], xr[a],   nn12[i,a]))
      
      succ = (1+xr[i])/(2+nn[i])               # predicted probability of sucess
      HH   = succ*sum(H1)+(1-succ)*sum(H2)     # expected unformation
      HF1  = c(HF1,HH)
    }
    
    # randomization 
    rand=(HF[n]-HF1)/sum(HF[n]-HF1)  # normalized increment in utility
    rand=(rand^pow)/sum(rand^pow)    
    rr1   = rbind(rr1,rand)
    d     = sample(length(pp),1,prob = rand)
    
    # update data
    nn[d] = nn[d]+1
    xr[d] = xr[d]+rbinom(1,1,pp[d])
    nn1   = rbind(nn1,nn)
    xx1   = rbind(xx1,xr)
  }
  
  return(list(Responses=xx1,SampleSize=nn1,Posterior=PPP,Utility=HF,rand=rr1))
}



################## BR design #######################
# Simulate single trial using balanced randomization
BR_trial=function(pp,N){

  A    = length(pp)
  rand = rep(1/A,A)
  nn1=xx1=HF=c()
  nn=xr=rep(0,A)
  
  for(n in 1:N){
   # compute utility
    H  = sapply(2:A,function(a) (varB(xr[1],nn[1],xr[a],nn[a])))
    HF = c(HF,sum(H))

   # treatment assignment
    d=sample(length(pp),1,prob = rand)
    
    # update data
    nn[d]=nn[d]+1
    xr[d]=xr[d]+rbinom(1,1,pp[d])
    nn1=rbind(nn1,nn)
    xx1=rbind(xx1,xr)
  }
  return(list(Responses=xx1, SampleSize=nn1, Variance=HF))
}



################## BAR design - Trippa 2013 ########
# Pr(theta.a >= theta.0)
p.post.sup0 = function(C_S,C_N,E_S,E_N,ab.C=c(1,1),ab.E=c(1,1)){
  
  # posterior parameter
  C_a = ab.C[1]+C_S      # beta parameter control arm
  C_b = ab.C[2]+C_N-C_S  # beta parameter control arm
  E_a = ab.E[1]+E_S      # beta parameter experimental arm
  E_b = ab.E[2]+E_N-E_S  # beta parameter experimental arm
  
  # probability that experimental rate is larger than control rate
  f   = function(p.Y)  (1-pbeta(p.Y,E_a,E_b))*dbeta(p.Y,C_a,C_b)   
        integrate(f,0,1)$value
}

# Simulate single trial using  Bayesian adaptive randomization (Trippa 2013)  
BAR_trial=function(pp,N){
  
  A   = length(pp)
  nn1 = xx1=HF=c()
  nn  = xr=rep(0,A)

  for(n in 1:N){
    
    # compute utility
    HF  = c(HF,sum(sapply(2:A,function(a) (varB(xr[1],nn[1],xr[a],nn[a])))))
    
    # randomization
                Pr  = sapply(2:A,function(a) p.post.sup0(xr[1],nn[1],xr[a],nn[a]))
    if(n<300){  rand=Pr^(n/300)/sum(Pr^(n/300)) 
    }else{      rand=Pr/sum(Pr) }
                rand=c(exp(.1*( max(nn[-1])-nn[1] ))/length(rand),rand)
                d   =sample(1:A, 1 ,prob = rand)
                
    # update data              
    nn[d]=nn[d]+1
    xr[d]=xr[d]+rbinom(1,1,pp[d])
    nn1=rbind(nn1,nn)
    xx1=rbind(xx1,xr)
  }
  return(list(Responses=xx1,SampleSize=nn1, Variance=HF))
}



################## BAR design - Thall 2006 ########
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

# Simulate single trial using  Bayesian adaptive randomization (Thall 2006)  
BAR_trial.Thall=function(pp,N){
  
  A=length(pp)
  nn1=xx1=HF=c()
  nn=xr=rep(0,A)
  
  for(n in 1:N){

    # compute posterior and randomization probability     
    Pr   = prob_best_arm(xr,nn) ^(0.5*n/N)
    rand = Pr / sum(Pr)
    
    # compute utility
    HF   = c(HF, sum(sapply(2:A,function(a) (varB(xr[1],nn[1],xr[a],nn[a])))) )

    # treatment assignment
    d    = sample(length(pp),1,prob = rand)
    
    # update data
    nn[d]= nn[d]+1
    xr[d]= xr[d]+rbinom(1,1,pp[d])
    nn1  = rbind(nn1,nn)
    xx1  = rbind(xx1,xr)
  }
  return(list(Responses=xx1,SampleSize=nn1, Variance=HF))
}



################## DBCD ##########################
# Simulate single trial using DBCD (Neyman allocation)
DBCD1_trial=function(pp,N,g=2){
  
  A   = length(pp)
  HF  = c()
  nn  = rep(1,A)
  nn1 = nn*diag(1,A,A)
  nn1 = apply(nn1,2,cumsum)
  xr  = rep(sapply(1:A,function(a) rbinom(1,1,pp[a])))
  xx1 = xr*diag(1,A,A)
  xx1 = apply(xx1,2,cumsum)
  
  for(n in A:N){

    # utility
    HF   = c(HF,sum(sapply(2:A,function(a) (varB(xr[1],nn[1],xr[a],nn[a])))))

    # randomization
    est  = (xr+1)/(nn+2)
    tar  = sqrt(est*(1-est))/sum(sqrt(est*(1-est))) 
    curr = nn/n
    d    = sample(length(pp),1,prob = tar*(tar/curr)^g)
    
    # update data
    nn[d]=nn[d]+1
    xr[d]=xr[d]+rbinom(1,1,pp[d])
    nn1=rbind(nn1,nn)
    xx1=rbind(xx1,xr)
  }
  return(list(Responses=xx1,SampleSize=nn1,Variance=HF))
  
  
}

# Simulate single trial using DBCD (target theta.k^.5 )
DBCD2_trial=function(pp,N,g=2){

  A   = length(pp)
  HF  = c()
  nn  = rep(1,A)
  nn1 = nn*diag(1,A,A)
  nn1 = apply(nn1,2,cumsum)
  xr  = rep(sapply(1:A,function(a) rbinom(1,1,pp[a])))
  xx1 = xr*diag(1,A,A)
  xx1 = apply(xx1,2,cumsum)
  
  for(n in length(pp):N){
    
    # compute information for comparision to BUD
    HF  = c(HF,sum( sapply(2:length(xr),function(a) (varB(xr[1],nn[1],xr[a],nn[a])))))

    
    # compute randomization probability
    est  = (xr+1)/(nn+2)
    tar  = sqrt(est)/sum(sqrt(est)) 
    curr = nn/n
    d    = sample(length(pp),1,prob = tar*(tar/curr)^g)   # randomize

    # update data
    nn[d] = nn[d]+1
    xr[d] = xr[d]+rbinom(1,1,pp[d])
    nn1   = rbind(nn1,nn)
    xx1   = rbind(xx1,xr)
  }
  
  return(list(Responses=xx1,SampleSize=nn1,Variance=HF))
  
}


