#----------------------------------
# Benjamin Risk
# brisk@emory.edu
# Functions supporting LCA
# This is beta software. R package to come.
#------------------------------
mlcaFP <- function(xData, n.comp = ncol(xData), W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('tiltedgaussian','logistic','JB'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'), max.comp = FALSE, reinit.max.comp=FALSE, df=0, irlba=FALSE,...) {

    #former option:
    #alg.typ = c('parallel','deflation'),
    #alg.typ = match.arg(alg.typ)
    alg.typ = 'parallel'
    require(steadyICA)

  
    distribution = match.arg(distribution)
    whiten=match.arg(whiten)

    if(restarts.dbyd>0 && whiten!='eigenvec') stop('Use whiten=eigenvec with restarts.dbyd')
    ## whiten:
    
    if(irlba) require(irlba)
    if(max.comp) { #if statement evaluates to true for all max.comp!=0
      s.comp = n.comp
      n.comp = max.comp
    }
    if(max.comp=='TRUE') stop('max.comp should be an integer or FALSE')
    if(reinit.max.comp && max.comp==FALSE) stop('Can not reinitialize from max.comp solution if max.comp==FALSE')
    if(reinit.max.comp && alg.typ=='deflation') stop('reinit.max.comp not yet written for deflation algorithm')
    require(ProDenICA)
    #require(multidcov)
    if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
    if(distribution=='tiltedgaussian' && df==0) stop('df must be greater than 0 for tiltedgaussian')
    if(distribution=='logistic'  && df>0) stop('df should be set to zero when using logistic')
    if(distribution=='logistic') Gfunc = logistic
    if(distribution=='JB'  && df>0) stop('df should be set to zero when using JB')
    if(distribution=='JB') Gfunc = jb.stat
    if(!is.null(W.list) & class(W.list)!='list') stop('W.list must be a list')
    if(length(W.list) && (restarts.pbyd || restarts.dbyd)) stop('restarts.pbyd and restarts.dbyd must be equal to zero when supplying W.list')

    orth.method= match.arg(orth.method)
    p = ncol(xData)
    nRow = nrow(xData)
    d = n.comp
    xData <- scale(xData, center=TRUE, scale=FALSE)
    if (d > p) stop('d must be less than or equal to p')
    if (whiten=='eigenvec') {
      # Use whitener=='eigenvec' so that restarts.dbyd initiates from the
      # span of the first d eigenvectors.
      temp = whitener(X = xData,n.comp = p,irlba=irlba)
      xData = temp$Z
      whitener = temp$whitener
      rm(temp)
      }  else if (whiten=='sqrtprec') {
         est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
         evd.sigma = svd(est.sigma)
         whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
         xData = xData%*%whitener
        }
    else {
      whitener = diag(p)
    }
    # warning('TO DO: Look at whitening methods and check for inconsistent options')
  if (is.null(W.list)) {
    if(restarts.pbyd) W.list = gen.inits(p=p,d=d,runs=restarts.pbyd,orth.method=orth.method)
    if(restarts.dbyd) {
      W.temp = gen.inits(p=d,d=d,runs=restarts.dbyd,orth.method=orth.method)
      #pad with zeros:
      zeros = matrix(0,p-d,d)
      W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros))
      W.list = c(W.list,W.temp)
    }
  }
  ## If restarts.pbyd and restarts.dbyd both equal zero:
  if (is.null(W.list)) W.list = gen.inits(p=p,d=d,runs=1,orth.method=orth.method)
  runs = length(W.list)
  out.list = NULL
  loglik.v = numeric(runs)
  for(k in 1:runs) {
    W0 = as.matrix(W.list[[k]])
    if(alg.typ == 'parallel') {
      out.list[[k]] = lca.par(xData=xData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df, ...)
    out.list[[k]]$df = df
    }
    if(alg.typ == 'deflation') {
      out.list[[k]] = lca.def(xData=xData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df,...)
    out.list[[k]]$df = df
    }
    if(max.comp) {
      flist0 = list()
      for (j in 1:d) flist0[[j]] <- Gfunc(out.list[[k]]$S[, j], ...)
      loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
      orderedLL = order(loglik.S,decreasing=TRUE)
      out.list[[k]]$S = out.list[[k]]$S[,orderedLL[1:s.comp]]
      out.list[[k]]$Ws = out.list[[k]]$Ws[,orderedLL[1:s.comp]]
      out.list[[k]]$loglik = sum(loglik.S[orderedLL[1:s.comp]])
      loglik.v[k] = out.list[[k]]$loglik
      } else {
      loglik.v[k] = out.list[[k]]$loglik
    }
  }
  for(i in 1:runs){
    out.list[[i]]$distribution=distribution
    out.list[[i]]$whitener = whitener
  }
  out = out.list[[which.max(loglik.v)]]
  if(reinit.max.comp) {
    out = lca.par(xData=xData,W0=out$Ws,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=s.comp,df=df, ...)
    out$df = df
    out.list[[k+1]] = out
  }
  if(out.all==TRUE) out.list else out
}
#---------------------------------

#-------------------------
gen.inits <- function(p,d,runs,orth.method=c('svd','givens')) {
  orth.method=match.arg(orth.method)
  W.list = list()
  for(i in 1:runs) {
    if(orth.method=='givens') {
      W.list[[i]] <- as.matrix(theta2W(runif(n=choose(p,2),min=0,max=2*pi)))[,1:d]
      } else {
        temp = matrix(rnorm(p*d),p,d)
        W.list[[i]] <- svd(temp)$u
        }
  }
  W.list
}

#-------------------------------------

logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
  #maximizes likelihood given s then calculates gradient w.r.t. w.hat
  #df is not used
  xData = as.vector(xData)
  # Yuxuan: gs and gps modified for numerical stability reason. Avoid infinity in the numerator.
  # list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = - 2*exp(-xData/scale) / (scale^2*(1+exp(-xData/scale))^2))
  list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), 
       gs = -1/scale + 2/(scale*(1+exp(xData/scale))), 
       gps = - 2 / (scale^2 * (exp(-xData/scale)+2+exp(-xData/scale))) )
}


# jb.stat written by Ze Jin:
jb.stat <- function(x, df=0) {
  n <- length(x)
  s <- sum(x^3)
  k <- sum(x^4)
  Gs <- x^3 * s / n^2 + (x^4 * k / n^2 + 9 / n - 6 * x^4 / n) / 4
  gs <- 6 * x^2 * s / n^2 + (8 * x^3 * (k / n - 3) / n) / 4
  gps <- 6 * (3 * x^4 + 2 * x * s) / n^2 + (24 * x^2 * (k / n - 3) / n + 32 * x^6 / n^2) / 4
  list(Gs = Gs, gs = gs, gps = gps)
}
  
#-----------------------
orthogonalize = function (W) {
  ##For arbitrary W, this is equivalent to (WW^T)^{-1/2} W
  temp <- svd(W)
  tcrossprod(temp$u,temp$v)
}

#------------------------------------------------
# Parallel algorithm:
lca.par <- function(xData,W0,Gfunc,maxit,verbose,density,eps,n.comp,df,...) {
  W0 = as.matrix(W0)
  d = ncol(W0)
  if(n.comp!=d) stop('W0 needs to be p x d')
  p = ncol(xData)
  nRow = nrow(xData)
  s <- xData %*% W0
  flist <- as.list(1:d)
  ##  Densities of likelihood components:
  for (j in 1:d) flist[[j]] <- Gfunc(s[, j], df=df,...)
  flist0 <- flist
  crit0 <- mean(sapply(flist0, "[[", "Gs"))
  nit <- 0
  nw <- 10
  repeat {
    nit <- nit + 1
    gS <- sapply(flist0, "[[", "gs")
    gpS <- sapply(flist0, "[[", "gps")
    #t1 <- t(xData) %*% gS/nRow
    t1 <- crossprod(xData,gS)/nRow
    t2 <- apply(gpS, 2, mean)
    if(d>1) W1 <- t1 - W0%*%diag(t2) else W1 <- t1 - W0*t2
    W1 <- orthogonalize(W1)
    if(d>1) nw <- frobICA(t(W0), t(W1))^2 else nw <- mean((W0-W1)^2) #Uses a measure that works for non-square matrices -- MSE. The measure is defined for M so here we use transpose of W.
    W0 <- W1
    s <- xData %*% W0
    for (j in 1:d) flist0[[j]] <- Gfunc(s[, j], df=df, ...)
    crit0 <- mean(sapply(flist0, "[[", "Gs"))
    if (verbose) cat("Iter", nit, "G", crit0, "Delta", nw, "\n")
    if ((nit >= maxit)) {
      warning('Max iter reached')
      break
    }
    if (nw < eps) break
  }
  out = list(Ws = W0, loglik = d*nRow*crit0, S = s)
  if(density) out$density = lapply(flist0, "[[", "density")
  out
}



#--------------------
tiltedgaussian = function (xData, df = 8, B = 100, ...) {
  #This function is based on ProDenICA::GPois by Trevor Hastie
  #NOTE: Assumes data are zero mean.
  require(gam)
  n <- length(xData)
  sd.x = sd(xData)
  rx <- c(min(xData)-0.1*sd.x, max(xData)+0.1*sd.x) 
  xg <- seq(from = rx[1], to = rx[2], length = B)
  gaps <- diff(rx)/(B - 1)
  xcuts <- c(rx[1] - gaps/2, xg[-B] + gaps/2, rx[2] + gaps/2)
  #NOTE: I use the response variable that corresponds to the LCA paper.
  #This differs from the GPois algorithm in ProDenICA
  ys <- as.vector(table(cut(xData, xcuts)))/(gaps*n)
  pois.fit <- suppressWarnings(gam(ys ~ s(xg, df)+offset(dnorm(xg,log=TRUE)), family = poisson, ...))
  Gs <- predict(pois.fit) #log tilt function predicted at grid locations (note: predict on gam object can not be used to obtain derivatives)
  # the gam object with the predict function can not be used directly to obtain the derivatives
  # of the smoothing spline.
  # Here, we refit another iteration of the IRWLS algorithm used in gam:
  # Note: working residuals = (y - mu0)/mu0
  # weights = mu0
  # fitted(pois.fit) = mu0
  # predict(pois.fit) = eta0 = log(mu0)
  sGs = Gs #+ log(sum(dnorm(xg))/sum(fitted(pois.fit)))
  z0 <- sGs + residuals(pois.fit, type='working')
  pois.refit <- smooth.spline(x=xg, y=z0, w=fitted(pois.fit),df=df) #obtain the log tilt function in an object that can be used to obtain derivatives
  Gs <- predict(pois.refit, xData, deriv = 0)$y
  gs <- predict(pois.refit, xData, deriv = 1)$y
  gps <- predict(pois.refit, xData, deriv = 2)$y
  fGs <- function(x) predict(pois.refit,x,deriv=0)$y
  fgs <- function(x) predict(pois.refit,x,deriv=1)$y
  fgps <- function(x) predict(pois.refit,x,deriv=2)$y
  list(Gs = Gs, gs = gs, gps = gps, fGs = fGs, fgs=fgs, fgps=fgps)
}
#---------------------------------------------

#----------------------------
# estimate mixing matrix from estimates of components:
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) coef(lm(xData~sData))[-1,] else coef(lm(xData~sData-1))
}


#-----------------------------------
# order by likelihood
# option for positive skewness
order.likelihood <- function(S,positive.skew=TRUE,distribution=c('logistic','tiltedgaussian','logcosh','JB'),out.loglik=FALSE,...) {
  distribution = match.arg(distribution)
  nObs = nrow(S)
  d = ncol(S)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='JB') Gfunc = jb.stat
  if(positive.skew) {
    skewness <- function(x, n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    skew = apply(S, 2, skewness)
    sign = -1 * (skew < 0) + 1 * (skew > 0)
    S = S %*% diag(sign)
  }
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
  orderedLL = order(loglik.S,decreasing=TRUE)
  S = S[,orderedLL]
  if(out.loglik) return(list(S=S,loglik=sort(loglik.S,decreasing=TRUE), order=orderedLL)) else list(S=S, order=orderedLL)
}

marginal.likelihoods <- function(S,distribution=c('logistic','tiltedgaussian','logcosh','GPois'),...)
{
  distribution = match.arg(distribution)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='GPois') Gfunc = ProDenICA::GPois
  d = ncol(S)
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  apply(sapply(flist0, "[[", "Gs"),2,sum)
}



