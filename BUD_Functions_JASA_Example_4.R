######################################################################
#### Example 4 (Sec 5): Multi-Arms trials with co-primary outcomes 
######################################################################

################## balanced randomization ##################
## P( Beta(nu.a) > Beta(nu.0)+Delta)
P.Beta.a.0 = function(nu.a, nu.0, Delta=0)
  integrate(f     =function(x) dbeta(x, nu.0[1], nu.0[2])*pbeta(x+Delta, nu.a[1], nu.a[2], lower.tail=F), lower=0, upper=1)$value

## simulate a single trial using balanced randomization
BR.MO       = function(N, theta){
 
	# nr or arms
	A    = nrow(theta)                                     
	# treatment assignments
	Av   = as.double(rmultinom(1, size = N, rep(1/A, A)))  
	# observed data
	Data = sapply(1:A, function(a)  as.double(rmultinom(1, size = Av[a], theta[a,]))   ) 
	
 rbind(Data, colSums(Data))
}

## simulate n.mc single trials using balanced randomization 
BR.MO.all = function(n.mc, N, theta)  
 replicate(n=n.mc,  BR.MO(N=N, theta=theta), simplify = F)


################## BUD for testing #########################
## compute BUD randomization probability and draws random number
Pr.PTE.1.a = function(Data, h, w12){

A     = ncol(Data)
p.hat = p.hat = (Data[1:4,]+1)/( matrix(Data[5,], 4, A, byrow = T)+4)

# compute utility
nu    = rbind(Data[1,]+Data[2,], Data[3,]+Data[4,], Data[1,]+Data[3,], Data[2,]+Data[4,])+2
p.E.1 = apply(nu[1:2,-1], 2, function(v) P.Beta.a.0(nu.a=v, nu.0=nu[1:2,1], Delta=0) )
p.E.2 = apply(nu[3:4,-1], 2, function(v) P.Beta.a.0(nu.a=v, nu.0=nu[3:4,1], Delta=0) )
H1    = p.E.1 - p.E.1^6
H2    = p.E.2 - p.E.2^6
U     = sum(H1)+sum(H2)
	
# expected compute utility
U.new = sapply(1:A, function(a){
    U.a = sapply(1:4, function(y) {
  D           = Data
  D[c(y,5),a] = D[c(y,5),a]+1
  nu          = rbind(D[1,]+D[2,], D[3,]+D[4,], D[1,]+D[3,], D[2,]+D[4,])+2
  p.E.1 = apply(nu[1:2,-1], 2, function(v) P.Beta.a.0(nu.a=v, nu.0=nu[1:2,1], Delta=0) )
  p.E.2 = apply(nu[3:4,-1], 2, function(v) P.Beta.a.0(nu.a=v, nu.0=nu[3:4,1], Delta=0) )
  H1    = p.E.1 - p.E.1^6
  H2    = p.E.2 - p.E.2^6
  U     = sum(H1)+sum(H2)
})
    sum(U.a * p.hat[,a])
})

# compute randimization probability and draw treatment assignment
Delta.a =  U-U.new 
sample(x=1:A, size = 1, prob = Delta.a^h/sum(Delta.a^h)) }

## simulate a single trial using BUD for testing efficacy
BUD.MO.Test = function(N, theta, h, w12){
	
 # generate potential outcome data 
 Outcomes           = cbind(Ai=rep(0,N), apply(theta, 1, function(p) sample(1:4, size = N,  replace = T, p) )) 
 A                  = ncol(Outcomes)-1
 
 # get summary matrix to store trial results
 Data               = matrix(0,5,A) 
 colnames(Data)     = paste0("arm", 1:A)

for(i in 1:N){
  # BUD assignment
  Ai                 = Outcomes[i,1] = Pr.PTE.1.a(Data=Data, h=h, w12) 
  # update data
  Yi                 = Outcomes[i,1+Ai]           
  Data[ c(Yi,5), Ai] = Data[ c(Yi,5), Ai] + 1 
} 
 
return( list(Stat=Data, Data=Outcomes)) }

## simulate n.mc trials using BUD for testing efficacy
BUD.MO.Test.all = function(n.mc, N, theta, h, w12) 
 lapply(1:n.mc, function(i) BUD.MO.Test(N=N, theta=theta, h=h, w12=w12)$Stat)


################## BUD for estimation ######################

## compute BUD randomization probability and draws random number
Pr.PTE.2 = function(Data, w12, h){

 A     = ncol(Data)
p.hat = (Data[1:4,]+1)/( matrix(Data[5,], 4, A, byrow = T)+4)
p1    = p.hat[1,]+p.hat[2,]
p2    = p.hat[1,]+p.hat[3,]

# compute utility
V1    = p1*(1-p1)/( Data[5,] + 4 +1 )
V2    = p2*(1-p2)/( Data[5,] + 4 +1 )
V11   = p.hat[1,]*(1-p.hat[1,])/( Data[5,] + 4 +1 )
V.G1  = V1[-1]+V1[1]
V.G2  = V2[-1]+V2[1]
V.G11 = V11[-1]+V11[1]
U     = sum(V.G11) + w12 * (sum(V.G1)+sum(V.G2))

# expected compute utility
U.new = sapply(1:A, function(a){

    U.a = sapply(1:4, function(y) {
  Data.new           = Data
  Data.new[c(y,5),a] = Data.new[c(y,5),a]+1
  p.hat = (Data.new[1:4,]+1)/( matrix(Data.new[5,], 4, A, byrow = T)+4)
p1    = p.hat[1,]+p.hat[2,]
p2    = p.hat[1,]+p.hat[3,]

V1    = p1*(1-p1)/( Data.new[5,] + 4 +1 )
V2    = p2*(1-p2)/( Data.new[5,] + 4 +1 )
V11   = p.hat[1,]*(1-p.hat[1,])/( Data.new[5,] + 4 +1 )

V.G1  = V1[-1]+V1[1]
V.G2  = V2[-1]+V2[1]
V.G11 = V11[-1]+V11[1]
U.new = sum(V.G11) + w12 * (sum(V.G1)+sum(V.G2))
})
    
    sum(U.a * p.hat[,a])
})

# compute randimization probability and draw treatment assignment
Delta.a =  U-U.new 
sample(x=1:A, size = 1, prob = Delta.a^h/sum(Delta.a^h) )

}

## simulate a single trial using BUD for estimating treatment effects
BUD.MO.Est  = function(N, theta, h, w12){
	
 # generate potential outcome data 
 Outcomes           = cbind(Ai=rep(0,N), apply(theta, 1, function(p) sample(1:4, size = N,  replace = T, p) )) 
 A                  = ncol(Outcomes)-1
 # get summary matrix to store trial results
 Data               = matrix(0,5,A) 
 colnames(Data)     = paste0("arm", 1:A)

for(i in 1:N){
  # select action
  Ai                 = Outcomes[i,1] = Pr.PTE.2(Data, w12, h)
  # update data
  Yi                 = Outcomes[i,1+Ai]
  Data[ c(Yi,5), Ai] = Data[ c(Yi,5), Ai] + 1 } 
 
return( list(Stat=Data, Data=Outcomes)) }

## simulate n.mc trials using BUD for estimating treatment effects
BUD.MO.Est.all  = function(n.mc, N, theta, h, w12) 
 lapply(1:n.mc, function(i) BUD.MO.Est(N=N, theta=theta, h=h, w12=w12)$Stat)
