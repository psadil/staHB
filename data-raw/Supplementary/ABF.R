library (rPorta)
library(gtools)

extremal <- function(o,pnames=letters[1:length(o)]) 
# generates a matrix whose rows are the extremal points for an order (o,
# specified as as a set of numbers with no equalites allowed) on a vector 
# of parameters with names pnames.
# e.g., o=1:3, pnames=letters(1:length(o)), a <= b <= c 
#       o=3:1, a >= b >= c  etc.)
{
  o <- rank(o)
  if (length(o)!=length(unique(o)))
    stop("No equalites allowed in the order specificaiton")
  n <- length(o)
  out <- matrix(as.numeric(upper.tri(matrix(nrow=n+1,ncol=n+1))),nrow=n+1)[,-1]
  colnames(out) <- pnames[o]
  out[,pnames]
}

extremal.and <-function(o1,o2=NULL,
  pnames1=letters[1:length(o1)],pnames2=NULL) 
# And between two sets of orders (by default the same twice=monotonic)
{
  if (is.null(o2)) o2 <- o1
  if (is.null(pnames2)) pnames2 <- LETTERS[1:length(o2)]
  o1 <- extremal(o1,pnames1)
  o2 <- extremal(o2,pnames2)
  o <-cbind(matrix(rep(o2,each=dim(o2)[1]),ncol=dim(o2)[2]),
    do.call("rbind", replicate(dim(o2)[1], o2, simplify=F))
	)
  colnames(o) <- c(pnames1,pnames2)
  o
}

  
make.hull <- function(omat,pnames=NULL,mono=FALSE) 
# Makes poi object, unique vertices of convex hull on orders specified in
# the rows of omat. If makes the joint order extremals
{
  if ( is.null(dim(omat)) )
    stop("Orders must be specified as rows of a matrix")
  if ( is.null(pnames) ) {
    if ( dim(omat)[2]>26 )
      stop("Automatic variable naming only works up to 26 variables")
    pnames <- letters[1:dim(omat)[2]]
  }
  if (!mono) out <- extremal(omat[1,],pnames) else
    out <- extremal.and(omat[1,],omat[1,],pnames)
  nrow <- nrow(omat)
  if (nrow>1) for (i in 2:nrow) 
    if (!mono) out <- rbind(out,extremal(omat[i,],pnames)) else
      out <- rbind(out,extremal.and(omat[i,],omat[i,],pnames)) 
  rownames(out) <- apply(out,1,paste,collapse=",")
  out <- out[!duplicated(rownames(out)),]
  as.poi(out[sort(rownames(out)),]) 
}

is.hull.vector <- function(ieqFileObject,test) 
# Tests if a point specified in vector test is in a convex hull
{
  signs <- as.matrix(ieqFileObject@inequalities@sign)
  coef <- as.matrix(ieqFileObject@inequalities@num)/
          as.matrix(ieqFileObject@inequalities@den)
  lhs <- coef[,-(length(test)+1)] %*% test 
  ok <- logical(length(lhs))
  ok[signs==-1] <- lhs <= coef[signs==-1,length(test)+1]
  ok[signs==0] <- lhs == coef[signs==0,length(test)+1]
  ok[signs==1] <- lhs >= coef[signs==1,length(test)+1]
  all(ok)
}

is.hull <- function(ieqFileObject,test) 
# Tests if multiple points specified in matrix test
# (one point per column) are in a convex hull
{
  signs <- as.matrix(ieqFileObject@inequalities@sign)
  if (!all(signs==-1))
    stop("A sign is not -1")
  coef <- as.matrix(ieqFileObject@inequalities@num)/
          as.matrix(ieqFileObject@inequalities@den)
  mult <- coef[,-(dim(test)[1]+1)]
  rhs <- coef[,dim(test)[1]+1]
  nz <- mult !=0
  ok <- !logical(dim(test)[2])
  for (i in 1:nrow(coef)) {
    ok[ok] <- mult[i,nz[i,],drop=FALSE]%*%test[nz[i,],ok,drop=FALSE] <= rhs[i]
  }
  ok
}

make.trace <- function(ntrace,ndim,trace.increasing=TRUE) 
# Create trace model orders given ntrace levels and ndim levels  
{
  npoint <- ntrace*ndim
  trace <- permutations(npoint,npoint)
  nvec <- dim(trace)[1]
  tracei <- matrix(1:npoint,ncol=ndim)
  ok <- !logical(nvec)
  if (trace.increasing) d <- 1 else d <- -1
  for (i in 1:ndim) {
    ok[ok] <- apply(trace[ok,],1,function(x){ 
      all(diff(x[x%in%tracei[,i]])==d) })
  } 
  trace <- trace[ok,]
  attr(trace,"ndim") = ndim
  trace
}

get.lap <- function(trace) 
# Filter output of make.trace to remove non-overlapping orders  
{
  ndim <- attr(trace,"ndim")
  tracei <- matrix(1:dim(trace)[2],ncol=ndim)
  laps <- matrix(!logical(dim(trace)[1]*ndim),ncol=ndim)
  for (i in 1:ndim) {
    laps[,i] <- apply(trace,1,function(x){ 
      all(abs(diff(c(1:length(x))[x%in%tracei[,i]]))==1) })
  }   
  trace <- trace[!apply(laps,1,all),,drop=FALSE]
  attr(trace,"ndim") = ndim
  trace
}

get.dim <- function(trace,dim.ord=1:attr(trace,"dim")) 
# Filter output of make.trace to remove orders not respecting dim.ord  
{
  ndim <- attr(trace,"ndim")
  tracei <- matrix(1:dim(trace)[2],ncol=ndim)
  ok <- !logical(dim(trace)[1])
  for (i in 2:ndim) {
    ok[ok] <- apply(trace[ok,],1,function(x){   
      all(c(1:length(x))[x%in%tracei[,i-1]] < c(1:length(x))[x%in%tracei[,i]])
    })
  } 
  trace <- trace[ok,,drop=FALSE]
  attr(trace,"ndim") = ndim
  trace
}

ci.p <- function(p,S,percent=95)
# percent credible intrval for p obtained from S samples
{
  c(qbeta(percent/200,p*S+1,S-p*S+1),
    qbeta(1-percent/200,p*S+1,S-p*S+1))
}


get.BF <- function(p,n,ieq,stopp=.1,ntest=1e6,maxrep=100,minrep=2,verbose=FALSE)
# Sample nrep*ntest sequentially to get BF if is.na(stopp) otherwise can pull
# out early when absolute % change on an interation is < stopp for two in a row
{
  BF <- 
    mean(is.hull(ieq,matrix(rbinom(ntest*length(p),n,p),ncol=ntest)/as.vector(n)))/
    mean(is.hull(ieq,matrix(rbinom(ntest*length(p),n,0.5),ncol=ntest)/as.vector(n)))
  if (verbose) print(paste(ntest,":",BF))
  dBFr <- Inf
  for (i in 2:maxrep)
  {
    BFi <-   
      mean(is.hull(ieq,matrix(rbinom(ntest*length(p),n,p),ncol=ntest)/as.vector(n)))/
      mean(is.hull(ieq,matrix(rbinom(ntest*length(p),n,0.5),ncol=ntest)/as.vector(n))) 
    oldBF <- BF
    olddBFr <- dBFr
    BF <- (i-1)*BF/i + BFi/i
    dBFr <- abs((BF-oldBF)/oldBF)
    if (!is.na(stopp) && ((dBFr < stopp/100) & (olddBFr < stopp/100)) ) 
      return(BF)
    if (verbose) print(paste(ntest*i,":",BF))
  }
  BF
}










