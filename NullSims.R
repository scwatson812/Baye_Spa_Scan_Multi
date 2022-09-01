###Specify the file path to the working directory
###This directory must contain DataFile.csv
path = 'Path/to/DataFile/'
setwd(path)

###Number of simulations
sims = 500

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

###Calculate Observed Probabilities
p.all = apply(df[,(2:(K+1))],2,mean)

###Store the output
lambda.mat.BMSS = matrix(NA,sims,n.zones)
most.likely.cluster.list.BMSS95 = list()
K.vec.BMSS95 = rep(NA,sims)


###Store the start time
for(i in 1:sims){
  ###Generate the data
  Y = rmultinom(n = n, size = 1, prob = p.all)
  lambda.vec.BMSS = rep(NA,n.zones)
  for(s in 1:n.zones){
    if(length(zone.index.list[[s]])>0){
      lambda.vec.BMSS[s] = log.int.f.y.given.H1(Y,s)
    }
  }
  most.likely.cluster = which(lambda.vec.BMSS == max(lambda.vec.BMSS,na.rm = TRUE))
  most.likely.cluster.list.BMSS95[[i]] = most.likely.cluster
  K.vec.BMSS95[i] = exp(max(lambda.vec.BMSS,na.rm = TRUE) + log((1-prob.H0)) - log.int.f.y.given.H0(Y) -log(prob.H0))
  print(i)
}

plot.vec = seq(0,1000,by = 0.1)
rejection.rate <- c()
for(i in 1:length(plot.vec)){
  rejection.rate = c(rejection.rate,mean(K.vec.BMSS95 >= plot.vec[i]))
}
###Plot the rejection rate as a function of the rejection threshold
plot(plot.vec,rejection.rate,type = 'l')
###Find the smallest rejection threshold so that the rejection rate is at most 0.05
min(plot.vec[rejection.rate <= 0.05])

