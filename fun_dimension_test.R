BIsearch_test <- function(xData, alpha = 0.05, fun, k.start=NULL, pval.start=NULL, verbose = TRUE, ...){
  # Input : xData(n * p), significance level for each test, function to test H_k
  # Output: 
  # dim_est: estimated dimension 
  # step_search: the test statistics, p-value and k of the whole search process
  # Test true dimension k between 1 and p-1. Error may happen when k = 0 or p
  p = ncol(xData)
  a = 1; b = p-1
  i = 0
  if (is.null(k.start)){
    k = ceiling(p/2) 
    Res_search = data.frame(p.value = 0, k)
  } 
  else{
    k = k.start
    if (is.null(pval.start)){
      #Res_search = data.frame(p.value = 0, statistics = 0, k)
      Res_search = data.frame(p.value = 0, k)
    }  
    else{
      i = i+1
      Res_search[i, ] = data.frame(p.value = pval.start, k)
      if (pval.start > alpha) b = k else a = k
      k = ceiling((a+b)/2)
      if (verbose){
        print(paste('input: starting dimension',k.start,'starting pvalue',pval.start,'first test dimension',k))
      }
    }
  } 
  
  # when b-a>1, we can just take k= ceiling((a+b)/2)
  # but when b-a=1, ceiling((a+b)/2)=b
  
  repeat{
    i = i+1
    Res = fun(xData, k=k, ...)
    # fun should have return 'p.value', 'k', 'statistic'
    # Res_search[i, ] = c(Res$p.value, Res$statistic, Res$k)
    Res_search[i, ] = c(Res$p.value, Res$k)
    if(verbose) print(paste('dimension', k, 'finished, the p.val is', Res$p.value))
    # if p-value larger than significance level, H_k is accepted, search for smaller k
    # if p-value larger than significance level, H_k is rejected, search for bigger k
    if (Res$p.value > alpha) b = k else a = k 
    if ( (b-a)==1 ) break
    k = ceiling((a+b)/2)
  }
  
  # we get the sequence a, a+1=b, either k=a or k=b
  # a should be tested already in the search unless a=1
  # b should be tested already in the search unless b=p-1
  # if both a and b are tested already and b=a+1, it must be true H_a is rejected and H_b is accepted
  # then the true dimension is b
  # find estimated dimension k0
  if (a == 1) {
    Res = fun(xData, k=1, ...)
    #Res_search[i+1, ] = c(Res$p.value, Res$statistics, Res$k)
    Res_search[i+1, ] = c(Res$p.value, Res$k)
    if (Res$p.value > alpha) k0 = 1 else k0 = 2
  } else if (b == p-1) {
    Res = fun(xData, k=p-1, ...)
    #Res_search[i+1, ] = c(Res$p.value, Res$statistics, Res$k)
    Res_search[i+1, ] = c(Res$p.value, Res$k)
    if (Res$p.value > alpha) k0 = p-1 else k0 = p
  } else {
    k0 = b
  }
  
  return(list(dim_est=k0, step_search=Res_search))
}

FOBI_GRF_test <- function(X=NULL, k, noise.dim = c(33,33), FWHM, n_sim=200, parallel = FALSE){
  # NULL H_0: the dimension of non-gaussian subspace is at most k
  # if p.val small, reject null, as least k+1
  # if p.val large, accept null, true k_0 <= k
  # when exactly k, p.val should have uniform distribution 
  p = dim(X)[2]
  require(JADE)
  require(moments)
  require(neuRosim)
  
  Xkurt = sort(apply(FOBI(X)$S, 2, kurtosis), decreasing = TRUE)
  #Xkurt = order.likelihood(icaimax(X, nc = dim(X)[2])$S, distribution = 'logistic', scale = 1, out.loglik = TRUE)$loglik/1089
  #Xkurt = sort(apply(icaimax(X, nc = p)$S, 2, function(x){mean(-x-2*log(1+exp(-x)))}), decreasing = TRUE)
  
  kurt_sim = numeric(n_sim)
  l = p - k
  
  #count = 0 
  if (parallel){
    kurt_sim = foreach (i = 1:n_sim, .combine = cbind)%dopar%{
      noise <- spatialnoise(dim = noise.dim, sigma=1, nscan = l, method = "gaussRF", FWHM = FWHM)
      dim(noise) <- c(prod(noise.dim),l)
      est = FOBI(noise)
      max(apply(est$S, 2, kurtosis))
    }
  }else{
    for (i in 1:n_sim){
      noise <- spatialnoise(dim = noise.dim, sigma=1, nscan = l, method = "gaussRF", FWHM = FWHM)
      dim(noise) <- c(prod(noise.dim),l)
      est = FOBI(noise)
      kurt_sim[i] = max(apply(est$S, 2, kurtosis))
      #est = icaimax(noise, nc = l)
      #kurt_sim[i] = max(apply(est$S, 2, function(x){mean(-x-2*log(1+exp(-x)))}))
      
      #if (Xkurt[k+1] < kurt_sim[i]) count = count + 1
      #if (count > alpha * n_sim){
      #  print(paste('early acceptance at dimension', k))
      #  break
      #}
    }
  }
  
  
  p.val = sum( (Xkurt[k+1] - kurt_sim)<0 )/n_sim
  #p.val = count/n_sim
  
  return(list(p.value = p.val, k = k))
}
