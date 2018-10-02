# source all functions
path = "Path to BUD_Functions_JASA_Example_1.R here"
setwd(path)

############################################################
#### Example 1: Estimating treatment effect in Multi-Arm trial 
############################################################
rm(list=ls())
source("BUD_Functions_JASA_Example_1.R")

######################## parameter / scenarios
{
n        = 336   # sample size
s        = 1     # scenario
Nsim     = 100   # number of simulation
A        = 4
scenario = rbind( scen1=rep(.4,4),
                  scen2=c(0.4,0.6,0.4,0.4),
                  scen3=c(0.4,0.6,0.4,0.2),
                  scen4=c(0.4,0.6,0.65,0.7))
p        = scenario[s,]
Gamma    = p[-1]-p[1]
}

########################  computing time ##############################
{
msim = 100
system.time({ lapply(1:msim,function(j) BUD_trial(pp=p,N=n,pow=3)) })
system.time({ lapply(1:msim,function(j) BR_trial(pp=p,N=n)) })
system.time({ lapply(1:msim,function(j) BAR_trial(pp=p,N=n)) })
system.time({ lapply(1:msim,function(j) BAR_trial.Thall(pp=p,N=n)) })
system.time({ lapply(1:msim,function(j) DBCD1_trial(pp=p,N=n)) })
system.time({ lapply(1:msim,function(j) DBCD1_trial(pp=p,N=n)) })
}

########################  BUD  ############################## 
{
BUD       = lapply(1:Nsim,function(j) BUD_trial(pp=p,N=n,pow=3))
n.BUD     = sapply(BUD, function(x) x$SampleSize[n,] )
ESS.BUD   = rowMeans(n.BUD)
SD.BUD    = apply(n.BUD, 1, sd)
p.BUD     = sapply(BUD, function(x)  sapply(2:A, function(a){ 
              M = rbind( c(x$Responses[n,a], x$Responses[n,1]), 
                         c(x$SampleSize[n,a]-x$Responses[n,a], x$SampleSize[n,1]-x$Responses[n,1]))
              fisher.test(M, alternative='greater')$p.value
              }))
Po.BUD    = rowMeans(p.BUD<=0.05)
MSE.BUD   = rowMeans(sapply(BUD, function(x){ 
                  p.a = x$Responses[n,-1]/x$SampleSize[n,-1]
                  p.0 = x$Responses[n,1]/x$SampleSize[n,1]
                  ((p.a-p.0)-Gamma)^2 }), na.rm = T)
}

########################  BR  ############################## 
{
BR       = lapply(1:Nsim,function(j) BR_trial(pp=p,N=n))
n.BR     = sapply(BR, function(x) x$SampleSize[n,] )
ESS.BR   = rowMeans(n.BR)
SD.BR    = apply(n.BR, 1, sd)
p.BR     = sapply(BR, function(x)  sapply(2:A, function(a){ 
              M = rbind( c(x$Responses[n,a], x$Responses[n,1]), 
                         c(x$SampleSize[n,a]-x$Responses[n,a], x$SampleSize[n,1]-x$Responses[n,1]))
              fisher.test(M, alternative='greater')$p.value
              }))
Po.BR    = rowMeans(p.BR<=0.05)
MSE.BR   = rowMeans(sapply(BR, function(x){ 
              p.a = x$Responses[n,-1]/x$SampleSize[n,-1]
              p.0 = x$Responses[n,1]/x$SampleSize[n,1]
              ((p.a-p.0)-Gamma)^2 }), na.rm = T)
}

########################  BAR according to Trippa's method # 
{
BAR       = lapply(1:Nsim,function(j) BAR_trial(pp=p,N=n))
n.BAR     = sapply(BAR, function(x) x$SampleSize[n,] )
ESS.BAR   = rowMeans(n.BAR)
SD.BAR    = apply(n.BAR, 1, sd)
p.BAR     = sapply(BAR, function(x)  sapply(2:A, function(a){ 
              M = rbind( c(x$Responses[n,a], x$Responses[n,1]), 
                         c(x$SampleSize[n,a]-x$Responses[n,a], x$SampleSize[n,1]-x$Responses[n,1]))
              fisher.test(M, alternative='greater')$p.value
              }))
Po.BAR    = rowMeans(p.BAR<=0.05)
MSE.BAR   = rowMeans(sapply(BAR, function(x){ 
              p.a = x$Responses[n,-1]/x$SampleSize[n,-1]
              p.0 = x$Responses[n,1]/x$SampleSize[n,1]
              ((p.a-p.0)-Gamma)^2 }), na.rm = T)
}

########################  BAR according to Thall's method  # 
{
BAR2       = lapply(1:Nsim,function(j) BAR_trial.Thall(pp=p,N=n))
n.BAR2     = sapply(BAR2, function(x) x$SampleSize[n,] )
ESS.BAR2   = rowMeans(n.BAR2)
SD.BAR2    = apply(n.BAR2, 1, sd)
p.BAR2     = sapply(BAR2, function(x)  sapply(2:A, function(a){ 
              M = rbind( c(x$Responses[n,a], x$Responses[n,1]), 
                         c(x$SampleSize[n,a]-x$Responses[n,a], x$SampleSize[n,1]-x$Responses[n,1]))
              fisher.test(M, alternative='greater')$p.value
              }))
Po.BAR2    = rowMeans(p.BAR2<=0.05)
MSE.BAR2   = rowMeans(sapply(BAR2, function(x){ 
                  p.a = x$Responses[n,-1]/x$SampleSize[n,-1]
                  p.0 = x$Responses[n,1]/x$SampleSize[n,1]
                  ((p.a-p.0)-Gamma)^2 }), na.rm = T)
}

########################  DBCD1 ############################
{
DBCD1       = lapply(1:Nsim,function(j) DBCD1_trial(pp=p,N=n))
n.DBCD1     = sapply(DBCD1, function(x) x$SampleSize[n,] )
ESS.DBCD1   = rowMeans(n.DBCD1)
SD.DBCD1    = apply(n.DBCD1, 1, sd)
p.DBCD1     = sapply(DBCD1, function(x)  sapply(2:A, function(a){ 
              M = rbind( c(x$Responses[n,a], x$Responses[n,1]), 
                         c(x$SampleSize[n,a]-x$Responses[n,a], x$SampleSize[n,1]-x$Responses[n,1]))
              fisher.test(M, alternative='greater')$p.value
              }))
Po.DBCD1    = rowMeans(p.DBCD1<=0.05)
MSE.DBCD1   = rowMeans(sapply(DBCD1, function(x){ 
                   p.a = x$Responses[n,-1]/x$SampleSize[n,-1]
                   p.0 = x$Responses[n,1]/x$SampleSize[n,1]
                   ((p.a-p.0)-Gamma)^2 }), na.rm = T)
}

########################  DBCD2 ############################
{
DBCD2       = lapply(1:Nsim,function(j) DBCD2_trial(pp=p,N=n))
n.DBCD2     = sapply(DBCD2, function(x) x$SampleSize[n,] )
ESS.DBCD2   = rowMeans(n.DBCD2)
SD.DBCD2    = apply(n.DBCD2, 1, sd)
p.DBCD2     = sapply(DBCD2, function(x)  sapply(2:A, function(a){ 
              M = rbind( c(x$Responses[n,a], x$Responses[n,1]), 
                         c(x$SampleSize[n,a]-x$Responses[n,a], x$SampleSize[n,1]-x$Responses[n,1]))
              fisher.test(M, alternative='greater')$p.value
              }))
Po.DBCD2    = rowMeans(p.DBCD2<=0.05)
MSE.DBCD2   = rowMeans(sapply(DBCD2, function(x){ 
               p.a = x$Responses[n,-1]/x$SampleSize[n,-1]
               p.0 = x$Responses[n,1]/x$SampleSize[n,1]
               ((p.a-p.0)-Gamma)^2 }), na.rm = T)
}

########################  summary for one scenario #########
{
Tab = rbind( c(ESS=ESS.BUD, SD=SD.BUD, Po=Po.BUD, MSE=MSE.BUD),
             c(ESS=ESS.BR,  SD=SD.BR,  Po=Po.BR, MSE=MSE.BR),
             c(ESS=ESS.BAR, SD=SD.BAR, Po=Po.BAR, MSE=MSE.BAR),
             c(ESS=ESS.BAR2, SD=SD.BAR2, Po=Po.BAR2, MSE=MSE.BAR2),
             c(ESS=ESS.DBCD1, SD=SD.DBCD1, Po=Po.DBCD1, MSE=MSE.DBCD1),
             c(ESS=ESS.DBCD2, SD=SD.DBCD2, Po=Po.DBCD2, MSE=MSE.DBCD2))

rownames(Tab) = c("BUD", "BR", "BAR", "BAR2", "DBCD1", "DBCD2")             
Tab

Utility = cbind(BUD.U  = sapply(BUD, function(x) x$Utility[n]),
                BR.U   = sapply(BR, function(x)  x$Variance[n]),
                BAR.U  = sapply(BAR, function(x) x$Variance[n]),
                BAR2.U = sapply(BAR2, function(x) x$Variance[n]),
                DBCD1.U = sapply(DBCD1, function(x) x$Variance[length(x$Variance)]),
                DBCD2.U = sapply(DBCD2, function(x) x$Variance[length(x$Variance)]) )
                
colnames(Utility) = c("BUD", "BR", "BAR", "BAR2", "DBCD1", "DBCD2")  
boxplot(Utility)
}



############################################################
#### Example 2: Selecting the best arm in a Multi-Arm trial 
############################################################
rm(list=ls())
source("BUD_Functions_JASA_Example_2.R")

######################## parameter / scenarios
{
n        = 10                            # total sample size
Nsim     = 20                            # number of simulations
s        = 1                             # selected scenario
scenario = rbind(scen1=c(.3,.4,.5,.6),
                 scen3=c(.4,.4,.4,.8),
                 scen4=c(.35,.45,.7,.8))
p        = scenario[s,]
}

########################  computing time ##############################
{
msim = 10
system.time({ lapply(1:msim,function(j) BUD_trial(pp=p,N=n,pow=2)) })
system.time({ lapply(1:msim,function(j) BR_trial(pp=p,N=n)) })
system.time({ lapply(1:msim,function(j) BAR_trial(pp=p,N=n)) })
system.time({ lapply(1:msim,function(j) RPW_trial(pp=p, N=n,nball=10)) })
}

########################  BUD  ############################## 
{
BUD     = lapply(1:Nsim,function(j) BUD_trial(pp=p, N=n,pow=2))
n_BUD   = sapply(BUD, function(x) x$SampleSize[n,])
ESS_BUD = rowMeans(n_BUD, na.rm = T)
SD_BUD  = apply(n_BUD,1, sd, na.rm = T)
p_BUD   = sapply(BUD, function(x) (x$Responses[n,]+1)/ (x$SampleSize[n,]+2) )
MSE.BUD = mean((apply(p_BUD,2,max)-max(scenario[s,]))^2)*1000
pr_BUD  = sapply(BUD, function(x) x$Posterior[n,])
}

########################  BR  ############################## 
{
BR      = lapply(1:Nsim,function(j) BR_trial(pp=p, N=n))
n_BR    = sapply(BR,  function(x) x$SampleSize[n,])
ESS_BR  = rowMeans(n_BR,  na.rm = T)
SD_BR   = apply(n_BR, 1, sd, na.rm = T)
p_BR    = sapply(BR,  function(x) (x$Responses[n,]+1)/ (x$SampleSize[n,]+2) )
MSE.BR  = mean((apply(p_BR,2,max)-max(scenario[s,]))^2)*1000
pr_BR   = sapply(BR, function(x)  x$Posterior[n,])
}

########################  BAR  ############################## 
{
BAR     = lapply(1:Nsim,function(j) BAR_trial(pp=p, N=n,pow=2))
n_BAR   = sapply(BAR, function(x) x$SampleSize[n,])
ESS_BAR = rowMeans(n_BAR, na.rm = T)
SD_BAR  = apply(n_BAR,1, sd, na.rm = T)
p_BAR   = sapply(BAR,  function(x) (x$Responses[n,]+1)/ (x$SampleSize[n,]+2) )
MSE.BAR = mean((apply(p_BAR,2,max)-max(scenario[s,]))^2)*1000
pr_BAR  = sapply(BAR, function(x) x$Posterior[n,])
}

########################  RPW  ############################## 
{
RPW     = lapply(1:Nsim,function(j) RPW_trial(pp=p, N=n,nball=10))
n_RPW   = sapply(RPW, function(x) x$SampleSize[n,])
ESS_RPW = rowMeans(n_RPW, na.rm = T)
SD_RPW  = apply(n_RPW,1, sd, na.rm = T)
p_RPW   = sapply(RPW,  function(x) (x$Responses[n,]+1)/ (x$SampleSize[n,]+2) )
MSE.RPW = mean((apply(p_RPW,2,max)-max(scenario[s,]))^2)*1000
pr_RPW  = sapply(RPW, function(x) x$Posterior[n,])
}

######################## summary for one scenario / sample size 
{
ESS     = rbind(ESS_BUD, ESS_BR, ESS_BAR, ESS_RPW)
SD      = rbind(SD_BUD, SD_BR, SD_BAR, SD_RPW)
MSE     = rbind(MSE.BUD, MSE.BR, MSE.BAR, MSE.RPW)
POWER   = rbind(BUD = tabulate(apply(pr_BUD, 2, which.max))/ncol(pr_BUD),
              BR  = tabulate(apply(pr_BR,  2, which.max))/ncol(pr_BR),
              BAR = tabulate(apply(pr_BAR, 2, which.max))/ncol(pr_BAR),
              RPW = tabulate(apply(pr_RPW, 2, which.max))/ncol(pr_RPW))

# summary of operating characteristics for the 4 designs
M           = cbind(ESS, SD, POWER, MSE)
colnames(M) = c(paste0('ESS.',1:4), paste0('SD(ESS).',1:4), paste0('P.',1:4), "MSE")
rownames(M) = c("BUD", "BR", "BAR", "RPW")
M

}

############################################################
#### Example 3: Biomarker Example 
############################################################
rm(list=ls())
source("BUD_Functions_JASA_Example_3.R")


######################## parameter / scenarios
{
p.control = 0.35 # response rate of control

# response rates in positive and negative groups 
pos=rbind(scen1=c(0,0,0,0),
          scen2=c(.2,0,0,0),
          scen3=c(.2,.2,0,0),
          scen4=c(.2,.3,0,0),
          scen5=c(.2,.15,.2,.3))

neg=rbind(scen1=c(0,0,0,0),
          scen2=c(0,0,0,0),
          scen3=c(0,0,0,0),
          scen4=c(0,0,0,0),
          scen5=c(0,.15,0,0))

# which biomarker does each arm target
BMK_targ=rbind(scen1=c(1,2,3,4),
               scen2=c(1,2,3,4),
               scen3=c(1,2,3,4),
               scen4=c(1,1,2,3),
               scen5=c(1,1,2,2))

# biomarker prevalance
BMK_prevalence=list(scen1=c(.5,.5,.5,.5),
                    scen2=c(.5,.5,.5,.5),
                    scen3=c(.7,.3,.5,.5),
                    scen4=c(.5,.5,.5),
                    scen5=c(.5,.7))
scen   = 1   # selected scenario
n.sim  = 10  # number of simulation
N      = 500
}

# parameter (response rates), conterfactual outcomes, BUD trial
PARAMETER = create.THETA(p0=p.control, TE_pos=pos[scen,],TE_neg=neg[scen,],target = BMK_targ[scen,])
DATA      = lapply(1:n.sim,function(j) Gen.Biomarker.data(THETA=PARAMETER,BIOM_prob=BMK_prevalence[[scen]],N=N))

BUD         = lapply(DATA,function(x) trial_Biomarker(pow=2, DATA=x, THETA = PARAMETER, target = BMK_targ[scen,]))
# test statistics marker positive population
test_z1.BUD = sapply(BUD, function(x){
  EVENTS = Suff.stat.subgroup(x$DATA, PARAMETER)
  AA     = Suff.stat.target(target = BMK_targ[scen,], EVENTS=EVENTS )
           sapply(1:4,function(a) ztest(Y=rbind(AA[1,1:2,a],AA[2,1:2,a])))  })
# test statistics marker negative population
test_z2.BUD = sapply(BUD, function(x){
  EVENTS = Suff.stat.subgroup(x$DATA, PARAMETER)
  AA     = Suff.stat.target(target = BMK_targ[scen,], EVENTS=EVENTS )
           sapply(1:4,function(a) ztest(Y=rbind(AA[1,3:4,a],AA[2,3:4,a])))  })

# parameter (response rates), conterfactual outcomes for BR trial
BR          = lapply(DATA,function(x) trial_Biomarker(pow=0, DATA=x, THETA = PARAMETER, target = BMK_targ[scen,]))
# test statistics marker positive population
test_z1.BR  = sapply(BUD, function(x){
  EVENTS = Suff.stat.subgroup(x$DATA, PARAMETER)
  AA     = Suff.stat.target(target = BMK_targ[scen,], EVENTS=EVENTS )
           sapply(1:4,function(a) ztest(Y=rbind(AA[1,1:2,a],AA[2,1:2,a])))  })
# test statistics marker negative population
test_z2.BR  = sapply(BUD, function(x){
  EVENTS = Suff.stat.subgroup(x$DATA, PARAMETER)
  AA     = Suff.stat.target(target = BMK_targ[scen,], EVENTS=EVENTS )
           sapply(1:4,function(a) ztest(Y=rbind(AA[1,3:4,a],AA[2,3:4,a])))  })


#scen=same= scenario used for simulations
#test.arm= which arm you want to test
#boot_sim= number of bootstrap simulations

 # average number of patients per arm that with biomarker tareted for each arm
 N.BUD.BMK.p  = sapply(BUD, function(x) colSums(Suff.stat.target(target = BMK_targ[scen,], EVENTS=Suff.stat.subgroup(x$DATA, PARAMETER))[1,1:2,] ) ) 
ESS.BUD.BMK.p = rowMeans(N.BUD.BMK.p)
 SD.BUD.BMK.p = apply(N.BUD.BMK.p, 1, sd, na.rm = T)
 POWER.BUD    = PSEUDO_BOOTSTRAP(scen=scen, test.arm=2, boot_sim=10, test_z1=test_z1.BUD, test_z2=test_z2.BUD)

 # average number of patients per arm that with biomarker tareted for each arm
 N.BUD.BMK.p  = sapply(BUD, function(x) colSums(Suff.stat.target(target = BMK_targ[scen,], EVENTS=Suff.stat.subgroup(x$DATA, PARAMETER))[1,1:2,] ) ) 
ESS.BUD.BMK.p = rowMeans(N.BUD.BMK.p)
 SD.BUD.BMK.p = apply(N.BUD.BMK.p, 1, sd, na.rm = T)
 POWER.BUD    = PSEUDO_BOOTSTRAP(scen=scen, test.arm=2, boot_sim=10, test_z1=test_z1.BUD, test_z2=test_z2.BUD)




############################################################
#### Example 4: Multi-Arms trials with co-primary outcomes 
############################################################
rm(list=ls())
source("BUD_Functions_JASA_Example_4.R")

######################## parameter / scenarios
{
# response probabilities
theta        = rbind(a1=c(0.15, 0.25, 0.4, 0.2),
                     a2=c(0.15, 0.25, 0.4, 0.2),
                     a3=c(0.15, 0.25, 0.4, 0.2),
                     a4=c(0.15, 0.25, 0.4, 0.2))
theta.s1     = theta.s2 = theta.s3 = theta.s4 = theta    # response rates in scenarios 1 to 4 
theta.s2[2,] = c(0.45, 0.15, 0.3, 0.1)
theta.s3[2,] = c(0.30, 0.30, 0.3, 0.1)
theta.s4[2,] = c(0.30, 0.15, 0.45, 0.1)
True.G1      = cbind(c(0, .2, .2, .05), matrix(0, 4, 2))
True.G2      = cbind(c(0, .2, .05, .2), matrix(0, 4, 2))

A            = 4   # number of arms
N            = 348 # number of patients 
n.mc         = 10  # number of simulations - n.mc=10000 in our simulations 
}

########################  computing time ##############################
{
msim = 10
system.time({ BR.MO.all(n.mc=msim, N=N, theta = theta.s1) })
system.time({ BUD.MO.Test.all(n.mc=msim, N=N, theta = theta.s1, h = 4, w12 = 2)   })
system.time({ BUD.MO.Est.all(n.mc=n.mc, N=N, theta = theta.s1, h = 4, w12 = 2)     })
}

########################  BR  ############################## 
{
## simulate trials using balanced randomization 
Results.BR    = NULL
Results.BR$S1 = BR.MO.all(n.mc=n.mc, N=N, theta = theta.s1)	 # scenario 1
Results.BR$S2 = BR.MO.all(n.mc=n.mc, N=N, theta = theta.s2)	 # scenario 2
Results.BR$S3 = BR.MO.all(n.mc=n.mc, N=N, theta = theta.s3)	 # scenario 3
Results.BR$S4 = BR.MO.all(n.mc=n.mc, N=N, theta = theta.s4)	 # scenario 4

# estimate response rates for each simulated scenario
BR.G1   = lapply(Results.BR, function(S) sapply(S, function(x){ p=(x[1,]+x[2,]) / x[5,]; p[-1]-p[1] }) )
BR.G2   = lapply(Results.BR, function(S) sapply(S, function(x){ p=(x[1,]+x[3,]) / x[5,]; p[-1]-p[1] }) )

# estimate MSE for each simulated scenario
MSE1.BR = 1000*sapply(1:4, function(s)  rowMeans( (BR.G1[[s]] -matrix(True.G1[s,], A-1,  ncol(BR.G1[[s]])) )^2, na.rm = T) )
MSE2.BR = 1000*sapply(1:4, function(s)  rowMeans( (BR.G2[[s]] -matrix(True.G2[s,], A-1,  ncol(BR.G2[[s]])) )^2, na.rm = T) )

# power under for each simulated scenario
Test.BR.1  = lapply(Results.BR, function(S) sapply(S, function(x) sapply(2:A, function(a) 1-phyper(q=sum(x[1:2,a])-1,k=sum(x[5,a]),m=sum(x[1:2,c(1,a)]),n=sum(x[3:4,c(1,a)]),lower.tail=TRUE)) ))
Test.BR.2  = lapply(Results.BR, function(S) sapply(S, function(x) sapply(2:A, function(a) 1-phyper(q=sum(x[c(1,3),a])-1,k=sum(x[5,a]),m=sum(x[c(1,3),c(1,a)]),n=sum(x[c(2,4),c(1,a)]),lower.tail=TRUE)) ))
Power.BR.1 = sapply(Test.BR.1, function(x) rowMeans(x<=0.05))
Power.BR.2 = sapply(Test.BR.2, function(x) rowMeans(x<=0.05))
Power.BR   = sapply(1:4, function(s) rowMeans(Test.BR.1[[s]] <= 0.05 & Test.BR.2[[s]] <= 0.05)  )
}

########################  BUD for testing effiacy  #########
{
## simulate trials using 
Results.BUD.1    = NULL
Results.BUD.1$S1 = BUD.MO.Test.all(n.mc=n.mc, N=N, theta = theta.s1, h = 4, w12 = 2)  
Results.BUD.1$S2 = BUD.MO.Test.all(n.mc=n.mc, N=N, theta = theta.s2, h = 4, w12 = 2)  
Results.BUD.1$S3 = BUD.MO.Test.all(n.mc=n.mc, N=N, theta = theta.s3, h = 4, w12 = 2)  
Results.BUD.1$S4 = BUD.MO.Test.all(n.mc=n.mc, N=N, theta = theta.s4, h = 4, w12 = 2)  

# estimate response rates for each simulated scenario
BUD.1.G1        = lapply(Results.BUD.1,    function(S) sapply(S, function(x){ p=(x[1,]+x[2,]) / x[5,]; p[-1]-p[1] }) )
BUD.1.G2        = lapply(Results.BUD.1,    function(S) sapply(S, function(x){ p=(x[1,]+x[3,]) / x[5,]; p[-1]-p[1] }) )

# estimate MSE for each simulated scenario
BUD.1.G1.MSE    = 1000*sapply(1:4, function(s)  rowMeans( (BUD.1.G1[[s]] -matrix(True.G1[s,], A-1,  ncol(BUD.1.G1[[s]])) )^2, na.rm = T) )
BUD.1.G2.MSE    = 1000*sapply(1:4, function(s)  rowMeans( (BUD.1.G2[[s]] -matrix(True.G2[s,], A-1,  ncol(BUD.1.G2[[s]])) )^2, na.rm = T) )

# power under for each simulated scenario
Test.BUD.1.G1  = lapply(Results.BUD.1, function(S) sapply(S, function(x) sapply(2:A, function(a) 1-phyper(q=sum(x[1:2,a])-1,k=sum(x[5,a]),m=sum(x[1:2,c(1,a)]),n=sum(x[3:4,c(1,a)]),lower.tail=TRUE)) ))
Test.BUD.1.G2  = lapply(Results.BUD.1, function(S) sapply(S, function(x) sapply(2:A, function(a) 1-phyper(q=sum(x[c(1,3),a])-1,k=sum(x[5,a]),m=sum(x[c(1,3),c(1,a)]),n=sum(x[c(2,4),c(1,a)]),lower.tail=TRUE)) ))
Power.BUD.1.G1 = sapply(Test.BUD.1.G1, function(x) rowMeans(x<=0.05))
Power.BUD.1.G2 = sapply(Test.BUD.1.G2, function(x) rowMeans(x<=0.05))
Power.BUD.1    = sapply(1:4, function(s) rowMeans(Test.BUD.1.G1[[s]] <= 0.05 & Test.BUD.1.G2[[s]] <= 0.05)  )
}

########################  BUD for estimation  ##############
{
## simulate trials using 
Results.BUD.2    = NULL
Results.BUD.2$S1 = BUD.MO.Est.all(n.mc=n.mc, N=N, theta = theta.s1, h = 4, w12 = 2)  
Results.BUD.2$S2 = BUD.MO.Est.all(n.mc=n.mc, N=N, theta = theta.s2, h = 4, w12 = 2)  
Results.BUD.2$S3 = BUD.MO.Est.all(n.mc=n.mc, N=N, theta = theta.s3, h = 4, w12 = 2)  
Results.BUD.2$S4 = BUD.MO.Est.all(n.mc=n.mc, N=N, theta = theta.s4, h = 4, w12 = 2)  

# estimate response rates for each simulated scenario
BUD.2.G1      = lapply(Results.BUD.2,    function(S) sapply(S, function(x){ p=(x[1,]+x[2,]) / x[5,]; p[-1]-p[1] }) )
BUD.2.G2      = lapply(Results.BUD.2,    function(S) sapply(S, function(x){ p=(x[1,]+x[3,]) / x[5,]; p[-1]-p[1] }) )

# estimate MSE for each simulated scenario
BUD.2.G1.MSE  = 1000*sapply(1:4, function(s)  rowMeans( (BUD.2.G1[[s]] -matrix(True.G1[s,], A-1,  ncol(BUD.2.G1[[s]])) )^2, na.rm = T) )
BUD.2.G2.MSE  = 1000*sapply(1:4, function(s)  rowMeans( (BUD.2.G2[[s]] -matrix(True.G2[s,], A-1,  ncol(BUD.2.G1[[s]])) )^2, na.rm = T) )

# power under for each simulated scenario
Test.BUD.2.G1  = lapply(Results.BUD.2, function(S) sapply(S, function(x) sapply(2:A, function(a) 1-phyper(q=sum(x[1:2,a])-1,k=sum(x[5,a]),m=sum(x[1:2,c(1,a)]),n=sum(x[3:4,c(1,a)]),lower.tail=TRUE)) ))
Test.BUD.2.G2  = lapply(Results.BUD.2, function(S) sapply(S, function(x) sapply(2:A, function(a) 1-phyper(q=sum(x[c(1,3),a])-1,k=sum(x[5,a]),m=sum(x[c(1,3),c(1,a)]),n=sum(x[c(2,4),c(1,a)]),lower.tail=TRUE)) ))
Power.BUD.2.G1 = sapply(Test.BUD.2.G1, function(x) rowMeans(x<=0.05))
Power.BUD.2.G2 = sapply(Test.BUD.2.G2, function(x) rowMeans(x<=0.05))
Power.BUD.2    = sapply(1:4, function(s) rowMeans(Test.BUD.2.G1[[s]] <= 0.05 & Test.BUD.2.G2[[s]] <= 0.05)  )
}

########################  summary for Example 4 ############ 
{
name1 = rep(paste0("Scenario.", 1:4), 5)
name2 = c(rep(c("MSE1", "MSE2") , each=4), rep(c("Power1", "Power2", "P") , each=4)) 
name3 = paste0(name1, "-", name2)

# summarize power and MSE for all three designs 
BR    = rbind(t(MSE1.BR), t(MSE2.BR), 100*t(Power.BR.1), 100*t(Power.BR.2), 100*t(Power.BR))
BUD1  = rbind(t(BUD.1.G1.MSE), t(BUD.1.G2.MSE), 100*t(Power.BUD.1.G1), 100*t(Power.BUD.1.G2), 100*t(Power.BUD.1))
BUD2  = rbind(t(BUD.2.G1.MSE), t(BUD.2.G2.MSE), 100*t(Power.BUD.2.G1), 100*t(Power.BUD.2.G2), 100*t(Power.BUD.2))
rownames(BR) = rownames(BUD1) = rownames(BUD2) = name3

# table S7
Tab = rbind( 
          cbind(BR[c(1,5, 9,13,17),1:2], BUD1[c(1,5, 9,13,17),1:2], BUD2[c(1,5, 9,13,17),1:2]), 
      			 cbind(BR[c(2,6,10,14,18),1:2], BUD1[c(2,6,10,14,18),1:2], BUD2[c(2,6,10,14,18),1:2]),
      			 cbind(BR[c(3,7,11,15,19),1:2], BUD1[c(3,7,11,15,19),1:2], BUD2[c(3,7,11,15,19),1:2]),
			       cbind(BR[c(4,8,12,16,20),1:2], BUD1[c(4,8,12,16,20),1:2], BUD2[c(4,8,12,16,20),1:2]))

Tab = xtable::xtable(Tab, align=rep("r", 7) )
}


