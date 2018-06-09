rm(list=ls())
source("ABF.R")

### Figure 4 example

# Hull defined by orders 1<=2<=3 and 3<=2<=1 
hull <- make.hull(rbind(1:3,3:1))
# Get facet defining inequalities 
ieqFileObject <- traf (new("poiFile", convex_hull=hull))
# Test points “Data1” and “Data2"
test <- cbind(c(.32,.66,.36),c(.8,.3,.8))
# vector indicating if point is in convex hull
is.hull.vector(ieqFileObject,test[,1])
# fast version for large test test, limited to <=
is.hull(ieqFileObject,test)
# Larger example
ntest=1e6; npar=3
test <- matrix(runif(ntest*npar),nrow=npar)
system.time(print(mean(is.hull(ieqFileObject,test)))) 
# [1] 0.666855


### Minimial (2x2) trace model example

hull <- make.hull(make.trace(ntrace=2,ndim=2,trace.increasing=TRUE),mono=TRUE)
ieq <- traf (new("poiFile", convex_hull=hull))
# Test ntest points
ntest <- 1e6; test <- matrix(runif(ntest*8),ncol=ntest)

# Trace model 
mean(is.hull(ieq,test)) 
# [1] 0.056636, p(trace) is 6/24^2 = 0.01041667 = 18.4%

# Removing overlap
hull <- make.hull(get.lap(make.trace(ntrace=2,ndim=2,trace.increasing=TRUE)),
                  mono=TRUE)
ieq <- traf (new("poiFile", convex_hull=hull))
mean(is.hull(ieq,test)) 
# [1] 0.027179, p(trace.mono.lap) is 4/24^2 = 0.006944444 = 25.5%

# Force dim order
hull <- make.hull(
  get.dim(get.lap(make.trace(ntrace=2,ndim=2,trace.increasing=TRUE)))
  ,mono=TRUE)
ieq <- traf (new("poiFile", convex_hull=hull))
v=mean(is.hull(ieq,test)) 
# [1] 0.0017269, p(trace.mono.lap.dim) is (1/(24^2))/v = 1 (witin error) 


### WM data

trace <- make.trace(ntrace=3,ndim=2,trace.increasing=FALSE)
trace.lap <- get.lap(trace)
trace.lap.dim <- get.dim(trace.lap)

trace.hull <- make.hull(trace,mono=TRUE)
trace.lap.hull <- make.hull(trace.lap,mono=TRUE)
trace.lap.dim.hull <- make.hull(trace.lap.dim,mono=TRUE)

trace.ieq <- traf (new("poiFile", convex_hull=trace.hull))
trace.lap.ieq <- traf (new("poiFile", convex_hull=trace.lap.hull))
trace.lap.dim.ieq <- traf (new("poiFile", convex_hull=trace.lap.dim.hull))

# Appendix table
head(trace.lap.ieq@inequalities@num)
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
# [1,]    0    0   -1    0    0    0    0    0    0     0     0     0     0
# [2,]    0    0    0    0    0   -1    0    0    0     0     0     0     0
# [3,]    0    0    0    0    0    0    0    0   -1     0     0     0     0
# [4,]    0    0    0    0    0    0    0    0    0     0     0    -1     0
# [5,]   -1    0    0    0    0    1    0    0    0     0     0     0     0
# [6,]   -1    1    0    0    0    0    0    0    0     0     0     0     0
# ...
tail(trace.lap.ieq@inequalities@num)
# [29,]    0    0    1    0    0   -1    0    0   -1     0     0     1     1
# [30,]    0    1    0   -1    0    0    0   -1    0     1     0     0     1
# [31,]    0    1    0    0   -1    0    0   -1    0     0     1     0     1
# [32,]    0    1    0    0    0   -1    0   -1    0     0     0     1     1
# [33,]    1    0    0   -1    0    0   -1    0    0     1     0     0     1
# [34,]    1    0    0    0   -1    0   -1    0    0     0     1     0     1

# Get data
load("WMdata.RData") 
names(data)[names(data)=="subjNumber"] <- "s"  # subjects
names(data)[names(data)=="setSize"]    <- "T"  # Trace
names(data)[names(data)=="silent"]     <- "D"  # Dimension
names(data)[names(data)=="sequential"] <- "S"  # State
names(data)[names(data)=="change"]     <- "R"  # Response
names(data)[names(data)=="CResp"]      <- "C"  # Response score
data <- data[,c("s","T","D","S","R","C")]

# Calcualte average probability
correct <- tapply(data$C,data[,c("T","D","S")],sum) # number correct
n <- tapply(data$C,data[,c("T","D","S")],length)    # number of trials

## RICH: I AM USING THE ML ESTIMATES CORRECT/N IS THAT OK?
##       ALMOST IDENTICAL BUT COULD USE BAYES p <- (correct+1)/(n+2) 
p <- correct/n 

# Sample from uniform binomial prior
ntest <- 1e6
test <- matrix(rbinom(ntest*length(p),n,.5),ncol=ntest)/as.vector(n)

trace.mono.prior <- mean(is.hull(trace.ieq,test))
# [1] 0.000698, p(trace.mono) is (20/720^2) = 5.5%
trace.mono.lap.prior <- mean(is.hull(trace.lap.ieq,test))
# [1] [1] 0.000598, p(trace.mono.lap) is (18/720^2) = 5.8%
trace.mono.lap.dim.prior <- mean(is.hull(trace.lap.dim.ieq,test))
# [1] 0.000036, p(trace.mono.lap.dim) is (4/720^2) = 21%

# Sample using average probabilities
ntest <- 1e6
test <- matrix(rbinom(ntest*length(p),n,p),ncol=ntest)/as.vector(n)

trace.mono.post <- mean(is.hull(trace.ieq,test))
# [1] 1
trace.mono.lap.post <- mean(is.hull(trace.lap.ieq,test))
# [1] 1
trace.mono.lap.dim.post <- mean(is.hull(trace.lap.dim.ieq,test))
# [1] 0.060109

BF.trace.mono <- trace.mono.post/trace.mono.prior
BF.trace.mono.lap <- trace.mono.lap.post/trace.mono.lap.prior
BF.trace.mono.lap.dim <- trace.mono.lap.dim.post/trace.mono.lap.dim.prior


# Repeatedly sample from prior and posterior until BF = posterior/prior settles
# down into a consistent estimate (by default < 0.1% change on two interations)
# or 100 iterations, each taking 1e6 samples for prior an dposterior
BF.trace.mono <- get.BF(p,n,trace.ieq,verbose=TRUE)
# [1] "1e+06 : 1342.28187919463"
# [1] "2e+06 : 1291.48833413826"
# ...
# [1] "2.6e+07 : 1293.92245153706"
# [1] "2.7e+07 : 1292.94114307052"
BF.trace.mono
# [1] 1292.788
BF.trace.mono.lap <- get.BF(p,n,trace.lap.ieq,verbose=TRUE)
# [1] "1e+06 : 1555.2099533437"
# [1] "2e+06 : 1651.73085079772"
# ...
# [1] "9e+06 : 1583.5068282843"
# [1] "1e+07 : 1583.88630418603"
BF.trace.mono.lap
# [1] 1584.887
BF.trace.mono.lap.dim <- get.BF(p,n,trace.lap.dim.ieq,verbose=TRUE)
# [1] "1e+06 : 2086.55172413793"
# [1] "2e+06 : 1885.38697318008"
# ...
# [1] "1.3e+07 : 2032.61307418777"
# [1] "1.4e+07 : 2031.06690222198"
BF.trace.mono.lap.dim
# [1] 2030.171

