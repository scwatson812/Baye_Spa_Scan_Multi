###Specify the file path to the working directory
###This directory must contain DataFile.csv
path = 'Path/To/DataFile/'
setwd(path)

###Specify the number of multinomial categories
K = 4

###Specify the Prior Distributions
alpha.in = rep(1,K)
alpha.out = rep(1,K)
alpha.all = rep(1,K)
###Specify prior probability of the null hypothesis being true
prob.H0 = 0.5

###Read in the data
df = read.csv('DataFile.csv')
n = dim(df)

###Define some helper functions
log.B<-function(alpha){
  res = sum(lgamma(alpha)) -lgamma(sum(alpha))
  return(res)
}

log.int.f.y.given.H1 <- function(Y,zone.id){
 Y.in = Y[,zone.index.list[[zone.id]]]
 Y.out = Y[,-zone.index.list[[zone.id]]]
 in.vec = rowSums(as.matrix(Y.in)) + alpha.in
 out.vec = rowSums(as.matrix(Y.out)) + alpha.out
 res = log.B(in.vec) + log.B(out.vec)- log.B(alpha.in) - log.B(alpha.out)
 return(res)
}

log.int.f.y.given.H0 <- function(Y){
  res = log.B(rowSums(Y) + alpha.all) - log.B(alpha.all)
  return(res)
}


zone.index.list <- list()
for(i in 1:(dim(df)[2]-(K+1))){
  zone.index.list[[i]] <- which(df[,(i+K+1)] == 1)
}
n.zones = length(zone.index.list)

Y = t(df[,(2:(K+1))])
Y = as.matrix(Y)

lambda.vec.BMSS = rep(NA,n.zones)
for(s in 1:n.zones){
  if(length(zone.index.list[[s]])>0){
    lambda.vec.BMSS[s] = log.int.f.y.given.H1(Y,s)
  }
}
most.likely.cluster = which(lambda.vec.BMSS == max(lambda.vec.BMSS,na.rm = TRUE))
Bayes.Factor = exp(max(lambda.vec.BMSS,na.rm = TRUE) + log((1-prob.H0)) - log.int.f.y.given.H0(Y) -log(prob.H0))
most.likely.cluster
Bayes.Factor


