# Written by Richard Morey and Andrew Heathcote
# Tested on R version 3.1.3, OpenMx Version 2.0.1, gtools 3.4.1, copula 0.999-13
# To get OpenMx type source('http://openmx.psyc.virginia.edu/getOpenMx.R') 
# at the R prompt

# NOTE: This code has not been thoroughly tested on designs with greater than 
#       two dimension-factor levels, Caveat Emptor! 

require('OpenMx') # Multivariate normal integration
require('gtools') # Permutations

post.prob.order = function(order=NULL,mus,sig2,lower=NULL,upper=NULL) 
# Use OpenMX to comptue Laplace approximation of p(order)  
{
  if(is.null(order)) {
    order=1:length(mus)
  }else{
    if(length(mus)!=length(order)) {
      stop("length of order vector must be the same as length of mus vector")
    }
  }
  mus = mus[order]
  sig2 = sig2[order]
  
  # We will determine the Pr(x1>x2>x3>...>xn)
  
  # normal approx

  lowerDiff = 0
  
  if(!is.null(lower) & !is.null(upper)){
    if(upper<=lower) stop("upper must be greater than lower.")
    upperDiff = upper - lower    
  }else{
    upperDiff = Inf
  }
  
  lbounds = rep(lowerDiff, length(mus)-1)
  ubounds = rep(upperDiff, length(mus)-1)
  
  covs = diag(as.vector(sig2))
  rotMat = matrix(c(rep(c(-1,1,rep(0,length(mus)-1)),length(mus)-2),-1,1),
                  byrow=TRUE,ncol=length(mus))
  if(!is.null(lower)){
    rotMat <- rbind(c(1,rep(0,ncol(rotMat)-1)),rotMat)
    lbounds = c(lower,lbounds)
    ubounds = c(ifelse(is.null(upper),Inf,upper),ubounds)
  }
  if(!is.null(upper)){
    rotMat <- rbind(rotMat,c(rep(0,ncol(rotMat)-1),1))
    ubounds = c(ubounds,upper)
    lbounds = c(lbounds,ifelse(is.null(lower),-Inf,lower)) 
  }
  
  newCov = rotMat%*%covs%*%t(rotMat)
  newMu = rotMat%*%mus
  
  # OpenMX integrator  
  omxMnor(newCov,newMu,lbound=lbounds,ubound=ubounds)
}


do.sampler<-function(means,sds,M=10000, lower=-Inf, upper=Inf) 
# Take samples
{
  Ul = pnorm(lower,means,sds)  
  Uu = pnorm(upper,means,sds)
  U = replicate(M,{ runif(length(means),Ul,Uu) })
  qnorm(U,means,sds)
}

restrict.samples<-function(means,sds,
  M=10000,probit=TRUE,
  trace.increasing=NULL,D.order=NULL,ndim=NA,ntrace=NA, 
  lower = -Inf, upper = Inf) 
# Take samples and filter out those not following restrictions  
{
  
  samples = do.sampler(means,sds,M, lower, upper)
  if ( !is.null(trace.increasing) ) {
    if ( is.na(ndim) )
      stop("Must supply dim for trace sampling") 
    ords = apply(samples,2,order)  
    
    is.in <- apply(ords,2,function(x) {
      d <- sign(matrix(apply(matrix(order(x),ncol=ndim),2,diff),ncol=ndim))
      all(d==ifelse(trace.increasing,1,-1))
    })   

    if (!is.null(D.order) & any(is.in)) {
      is.D <- apply(ords[,is.in,drop=F],2,function(x) {
        d <- sign(matrix(apply(matrix(order(x),ncol=ndim)[,D.order],
                               1,diff),nrow=ntrace))
        all(d==ifelse(dim.increasing <- D.order[2]>D.order[1],1,-1))
      })
      is.in[is.in] <- is.D      
    }
    samples <- samples[,is.in]
  }
  if (dim(samples)[2]==0)
    stop("No samples after filtering, consider increasing M")
  if ( probit ) samples <- pnorm(samples)
  samples
}


getPPP <- function(dat,probit=TRUE, # Laplace on probit scale?
  use.sampling=FALSE,          # Very slow!
  M=1e6,nrep=100,verbose=TRUE, # M samples nrep time, verbose: dot/M samples
  trace.increasing=NULL, # does accuracy increase or decrase with trace levels?
  D.order=NULL, # vector of order constaints (integers or D factor levels)
  lower=NULL,   # Lower bound to filter samples (p for binomial, Z for probit)
  upper=NULL)   # Upper bound to filter samples (p for binomial, Z for probit)
# Get prior and posterior probabilities by numerical integration
{

  sumOffDiag <- function(M){
    diag(M) <- 0
    sum(M)
  }
  
  probs.sim <-function(probs,m,s,ndim,ntrace,M=1e6,nrep=10,verbose=TRUE) 
  # Probabilities for each entry in probs by sampling  
  {      
    for (i in 1:nrep) {
      # Compute probabilities of orderings for x and y axis
      sim <- table(apply(restrict.samples(m,s,
        M=M,probit=probit,ndim=ndim,ntrace=ntrace,
        trace.increasing=trace.increasing,D.order=D.order,lower=lower),
        2,function(x){paste(order(x),collapse=",")}))
      probs[names(sim)] <- probs[names(sim)]+sim 
      if (verbose) cat(".")
    }
    if (verbose) cat("\n")
    probs/(nrep*M)
  }

  # Check D.order is properly specified and convert to integer form
  if (!is.null(D.order)) {
    if (length(D.order)<2)
      stop("D.order must have at least two elements")
    if (!is.numeric(D.order)) {
      if ( !all(D.order %in% levels(dat$D)) )
        stop("D.order can only contain names in levels(dat$D)")
      tmp <- 1:length(levels(dat$D))
      names(tmp) <- levels(dat$D)
      D.order <- tmp[D.order]
    }
    if ( !is.integer(D.order) | min(D.order)<1 | 
      max(D.order) > length(levels(dat$D)) )
      stop("D.order must be integers in the range 1:length(levels(dat$D))")
    dim.increasing <- D.order[2]>D.order[1]
  }

  # Reserve object for prior & posterior probabilities
  p <- rep(1,length(levels(dat$s)))
  names(p) <- unique(dat$s)  
  # No constraints
  post.p <- cbind(nonmono=p,mono=p)
  prior.p <- post.p
  # Trace constraint
  post.p.trace <- cbind(prior.p,nontrace=p,nolap=p)
  prior.p.trace <- post.p.trace
  # Trace and dimension constraint
  post.p.trace.dim <- cbind(prior.p.trace,mono.nodim=p)
  prior.p.trace.dim <- post.p.trace.dim
  
  # Design counts
  ntrace <- length(levels(dat$T))
  ndim <- length(levels(dat$D))
  npoints <- ntrace*ndim
  
  # prior counts multiplied by these to get prior probabilities
  n.orders <- factorial(npoints)
  prior.p <- prior.p/n.orders^2
  prior.p.trace <- prior.p.trace/n.orders^2
  prior.p.trace.dim <- prior.p.trace.dim/n.orders^2

  cat("Processing participant: ")
  # Loop over subjects
  for( sub in levels(dat$s) ) {
    
    cat(paste(sub," "))
    
    datSub <- dat[dat$s==sub,]
  
    # set up for normal approximation to posterior
  
    if (probit) { # probit
      YMu <- datSub$probitMean[1:npoints]
      XMu <- datSub$probitMean[1:npoints + npoints]
      YStdErr <- datSub$probitStdErr[1:npoints]
      XStdErr <- datSub$probitStdErr[1:npoints + npoints]
    } else {    # binomial
      YMu <- datSub$phat[1:npoints]
      XMu <- datSub$phat[1:npoints + npoints]
      YStdErr <- datSub$stdErr[1:npoints]
      XStdErr <- datSub$stdErr[1:npoints + npoints]
    }
    
    # Determine all permutations, uses gtools
    
    perms <- permutations(npoints,npoints)
    
    # Compute probabilities of orderings for x and y axis
    
    if (use.sampling) { # Use repeated sampling, e.g., 4674s,for M=1e6, nrep=100
      tmp <- apply(perms,1,function(x){paste(x,collapse=",")})
      probsX <- numeric(length(tmp))
      names(probsX) <- tmp
      probsY <- probs.sim(probsX,my,sy,ndim,ntrace,M,nrep,verbose)
      probsX <- probs.sim(probsX,mx,sx,ndim,ntrace,M,nrep,verbose)
    } else { # Use numerical integration, 1.103s
      probsY <- apply(perms,1,post.prob.order,mus=YMu,sig2=YStdErr^2,
                      lower=lower,upper=upper)
      probsX <- apply(perms,1,post.prob.order,mus=XMu,sig2=XStdErr^2,
                      lower=lower,upper=upper)    
      # Renormalize to fix integration error
      probsY <- probsY / sum(probsY)
      probsX <- probsX / sum(probsX)
      names(probsY) <- names(probsX) <- apply(perms,1,paste,collapse=',')
    }
    
    # Assume independence
    jointOrderProbs <- outer(probsY,probsX)
    
    # No constraint results
    
    # monotonic
    post.p[sub,"mono"] <- sum(diag(jointOrderProbs))
    prior.p[sub,"mono"] <- dim(perms)[1]*prior.p[sub,"mono"]
    # non-monotonic
    post.p[sub,"nonmono"] <- 1-post.p[sub,"mono"]
    prior.p[sub,"nonmono"] <- 1-prior.p[sub,"mono"]

    # Apply trace constraint
    if (!is.null(trace.increasing)) { 
      # subset with trace constraint
      is.trace <- apply(perms,1,function(x) {
        d <- sign(matrix(apply(matrix(order(x),ncol=ndim),2,diff),ncol=ndim))
        all(d==ifelse(trace.increasing,1,-1))
      })   
      n.trace <- sum(is.trace)

      # no overlap within trace
      is.nolap <- is.trace & apply(perms,1,function(x) {
        d <- matrix(apply(matrix(x,ncol=ndim),2,diff),ncol=ndim)
        all(abs(d)==1) & all(apply(d[,-1,drop=F],2,function(y){all(y==d[,1])}))
      })   
     
      # non-monotonic and trace
      nonmono <- jointOrderProbs[is.trace,is.trace]; diag(nonmono) <- 0
      post.p.trace[sub,"nonmono"] <- sum(nonmono)
      prior.p.trace[sub,"nonmono"] <- 
        (n.trace^2-n.trace)*prior.p.trace[sub,"nonmono"]
      # monotonic and trace and overlap
      post.p.trace[sub,"mono"] <- sum(diag(
        jointOrderProbs[is.trace & !is.nolap,is.trace & !is.nolap,drop=FALSE]))
      prior.p.trace[sub,"mono"] <- 
        sum(is.trace & !is.nolap)*prior.p.trace[sub,"mono"]       
      # monotnoinc and trace and no overlap
      post.p.trace[sub,"nolap"] <- 
        sum(diag(jointOrderProbs[is.nolap,is.nolap,drop=FALSE]))
      prior.p.trace[sub,"nolap"] <- 
        sum(is.nolap)*prior.p.trace[sub,"nolap"]            
      # non-trace
      post.p.trace[sub,"nontrace"] <-   
        1-sum(post.p.trace[sub,dimnames(post.p.trace)[[2]]!="nontrace"])
      prior.p.trace[sub,"nontrace"] <-   
        1-sum(prior.p.trace[sub,dimnames(prior.p.trace)[[2]]!="nontrace"])
      prior.trace <- data.frame(prior.p.trace)
      post.trace <- data.frame(post.p.trace)
      
      # subset with trace and dimension constraint
      if (!is.null(D.order)) { 
        is.D <- apply(perms[is.trace,],1,function(x) {
          d <- sign(matrix(apply(matrix(order(x),ncol=ndim)[,D.order],
                                 1,diff),nrow=ntrace))
          all(d==ifelse(dim.increasing,1,-1))
        })
        is.trace.dim <- is.trace
        is.trace.dim[is.trace] <- is.D
      
        # trace but not dimension constraint
        is.trace.not.dim <- is.trace & !is.trace.dim 
        # no overlap within trace and dimension
        is.nolap.dim <- is.nolap & is.trace.dim          
        # non-trace (not effected by dim constraint)
        post.p.trace.dim[sub,"nontrace"] <- post.p.trace[sub,"nontrace"]  
        prior.p.trace.dim[sub,"nontrace"] <- prior.p.trace[sub,"nontrace"]   
      
        # non-monotonic
        post.p.trace.dim[sub,"nonmono"] <- sumOffDiag(jointOrderProbs[
          is.trace.dim,is.trace.dim,drop=F])
        prior.p.trace.dim[sub,"nonmono"] <- sumOffDiag(jointOrderProbs[
          is.trace.dim,is.trace.dim,drop=F]*0+1) / n.orders^2

        # monotonic and trace and overlap and dimension
        post.p.trace.dim[sub,"mono"] <- sum(diag(jointOrderProbs[
          is.trace.dim & !is.nolap.dim,is.trace.dim & !is.nolap.dim,drop=F]))
        prior.p.trace.dim[sub,"mono"] <- 
          sum(is.trace.dim & !is.nolap.dim)*prior.p.trace.dim[sub,"mono"]
        # monotnoinc and trace and no overlap and dimension
        post.p.trace.dim[sub,"mono.nodim"] <- sum(diag(jointOrderProbs[
          is.trace.not.dim & !is.nolap.dim,is.trace.not.dim & !is.nolap.dim,drop=F]))
        prior.p.trace.dim[sub,"mono.nodim"] <- 
          sum(is.trace.not.dim & !is.nolap.dim)*prior.p.trace.dim[sub,"mono.nodim"]
        post.p.trace.dim[sub,"nolap"] <- 
          sum(diag(jointOrderProbs[is.nolap.dim,is.nolap.dim,drop=FALSE]))
        prior.p.trace.dim[sub,"nolap"] <- 
          sum(is.nolap.dim)*prior.p.trace.dim[sub,"nolap"] 
        prior.p.trace.dim <- data.frame(prior.p.trace.dim)
        post.p.trace.dim <- data.frame(post.p.trace.dim)
      } else {
        prior.p.trace.dim <- NA
        post.p.trace.dim <- NA
      }
    } else {
      prior.trace <- NA
      post.trace <- NA
    }
    if (use.sampling & verbose) cat("\n")
  }
  cat("\n")

  out <- list(prior=data.frame(prior.p),
              post=data.frame(post.p),
              prior.trace=prior.trace,
              post.trace=post.trace,
              prior.trace.dim=prior.p.trace.dim,
              post.trace.dim=post.p.trace.dim)
  attr(out,"ords") <- list(trace=perms[is.trace,],trace.dim=perms[is.trace.dim,])
  attr(out,"labs") <- paste(datSub$T[1:npoints],datSub$D[1:npoints],sep=".")
  out
}


show.order <- function(ppp,type="trace.dim",labels=NULL) 
# Print out orders stored in ppp object  
{
  if (is.null(labels)) labels <- attr(ppp,"labs") 
  t(apply(attr(ppp,"ords")[[type]],1,
          function(row,labs) labs[row],labs=labels))
}

getBF <- function(ppp,log10bf=TRUE) 
# Calcualte Bayes Factors from prior and posterior proabilites  
# Value: list with names BFx, where x = number of restrictions
#   BF0: no restricitons, $m.nm = mono/non-mono
#   BF1: trace restriciton, $m.nm (given trace) and $t.nt trace/non-trace
#                           $m.nl (given trace) monotonic vs. non-overlap
#   BF2: trace + dimension. $m.nm (given trace & dimension)
{

  # No restriction
  
  # Mono vs. non-mono
  BF0.m.nm <- apply(ppp$post/ppp$prior,1,function(x){x[2]/x[1]})

  # Trace model

  #mono vs. non-mono
  BF.trace <- data.frame(ppp$post.trace/ppp$prior.trace)
  BF1.m.nm <- apply(BF.trace,1,function(x){x["mono"]/x["nonmono"]})
  # montonic overlap vs. monotonic no-overlap
  BF1.m.nl <- apply(BF.trace,1,function(x){x["mono"]/x["nolap"]})
 
  # trace vs. non-trace
  post.trace.p <- apply(ppp$post.trace[,c("nonmono","mono","nolap")],1,sum)
  prior.trace.p <- apply(ppp$prior.trace[,c("nonmono","mono","nolap")],1,sum)
  BF1.t.nt <- (post.trace.p/prior.trace.p) /
     (ppp$post.trace$nontrace/ppp$prior.trace$nontrace)

  # Trace and Dimension model

  if (!any(is.na(ppp$post.trace.dim))) {
    # Mono vs. non-mono
    BF.trace.dim <- data.frame(ppp$post.trace.dim/ppp$prior.trace.dim)
    BF2.m.nm <- apply(BF.trace.dim,1,function(x){x["mono"]/x["nonmono"]})
  } else BF2.m.nm <- NA
  
  if (log10bf) {
    out <- list(BF0=list(m.nm=log10(BF0.m.nm)),
      BF1=list(m.nm=log10(BF1.m.nm),t.nt=log10(BF1.t.nt),m.nl=log10(BF1.m.nl)),
      BF2=list(m.nm=log10(BF2.m.nm)))
    attr(out,"type") <- "log10(BF)"
  } else {
    out <- list(BF0=list(m.nm=BF0.m.nm),
      BF1=list(m.nm=BF1.m.nm,t.nt=BF1.t.nt,m.nl=BF1.m.nl),
      BF2=list(m.nm=BF2.m.nm))
    attr(out,"type") <- "BF"
  }
  out
}


restrict.samples.mean<-function(means,sds,labs,ndim,ntrace,
  ci=100*pnorm(1),M=10000,probit=TRUE,
  trace.increasing=NULL,D.order=NULL, 
  lower = -Inf, upper = Inf) 
# Get means of filtered samples  
{
  samples <- restrict.samples(means=means,sds=sds,
    trace.increasing=trace.increasing,D.order=D.order,
    ndim=ndim,ntrace=ntrace,lower=lower,upper=upper)
  list(mns=rowMeans(samples),
    lohi=apply(samples,1,quantile,probs=c(100-ci,ci)/100))
}


plot.st <- function(datSub,probit=TRUE,ylim=c(.5,1),xlim=c(.5,1),
                    xlab="Sequential",ylab="Simultaneous",cex=1,
                    col="black",caplen=.05,main="",ci=100*pnorm(1),M=1e4,
                    trace.increasing=NULL,D.order=NULL,lower=-Inf,upper=Inf) 
# Takes M (default 1e4) samples, points with SE (CI = x changes to x% CI) 
# with caplen bars in color col.
{
  
  # Check D.order is properly specified and convert to integer form
  if (!is.null(D.order)) {
    if (length(D.order)<2)
      stop("D.order must have at least two elements")
    if (!is.numeric(D.order)) {
      if ( !all(D.order %in% levels(datSub$D)) )
        stop("D.order can only contain names in levels(datSub$D)")
      tmp <- 1:length(levels(datSub$D))
      names(tmp) <- levels(datSub$D)
      D.order <- tmp[D.order]
    }
    if ( !is.integer(D.order) | min(D.order)<1 | 
      max(D.order) > length(levels(datSub$D)) )
      stop("D.order must be integers in the range 1:length(levels(datSub$D))")
    dim.increasing <- D.order[2]>D.order[1]
  }

  # Design counts
  ntrace <- length(levels(datSub$T))
  ndim <- length(levels(datSub$D))
  npoints <- ntrace*ndim

  if (probit) {
    mname <- "probitMean"; sename <- "probitStdErr"    
  } else {
    mname <- "phat"; sename <- "stdErr"
  }
  # means
  mx <- datSub[[mname]][1:npoints + npoints] 
  my <- datSub[[mname]][1:npoints]           
  sx <- datSub[[sename]][1:npoints + npoints]
  sy <- datSub[[sename]][1:npoints]

  # Simulate
  simx <- restrict.samples.mean(mx,sx,ndim=ndim,ntrace=ntrace,ci=ci,M=M,
    probit=probit,trace.increasing=trace.increasing,D.order=D.order,lower=lower)
  simy <- restrict.samples.mean(my,sy,ndim=ndim,ntrace=ntrace,ci=ci,M=M,
    probit=probit,trace.increasing=trace.increasing,D.order=D.order,lower=lower)
  mx <- simx$mns; my <- simy$mns
  lox <- simx$lohi[1,]; hix <- simx$lohi[2,]
  loy <- simy$lohi[1,]; hiy <- simy$lohi[2,]    
  
  # Plot
  angle=90
  plot(mx,my,pch=as.character(datSub$T[1:npoints]),
       xlab=xlab,ylab=ylab,cex=cex,
       main=main,ylim=ylim,xlim=xlim)
  lines(mx[1:3],my[1:3],lty=1)
  lines(mx[1:3 + 3],my[1:3 + 3],lty=2)
  arrows(lox,my,hix, my,
    code=3,angle=angle,col=col,length=caplen)
  arrows(mx,loy,mx, hiy,
    code=3,angle=angle,col=col,length=caplen)
}

