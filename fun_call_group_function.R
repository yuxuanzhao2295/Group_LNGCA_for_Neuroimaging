group_LNGCA = function(data.list, n.individual, n.group,
                      whiten = "eigenvec", restarts.pbyd = 0,  restarts.dbyd = 30, distribution='logistic',
                      parallel = FALSE){
  nsub = length(data.list)
  estSubNGCA = list()
  
  # subject level ngca
  if (parallel){
    estSubNGCA = foreach(i = 1:nsub)%dopar%{
      mlcaFP(xData = data.list[[i]]$X, n.comp = n.individual[i], 
             whiten = whiten, restarts.pbyd = restarts.pbyd,  restarts.dbyd =  restarts.dbyd, 
             distribution=distribution)
    }
  }else{
    for (i in 1:nsub){
      estSubNGCA[[i]] = mlcaFP(xData = data.list[[i]]$X, n.comp = n.individual[i], 
                               whiten = whiten, restarts.pbyd = restarts.pbyd,  restarts.dbyd =  restarts.dbyd, 
                               distribution=distribution)
    }
  }
  
  
  # concatenate data for group NGCA 
  LC.conca = NULL
  for (i in 1:nsub) LC.conca = cbind(LC.conca, estSubNGCA[[i]]$S)
  
  # extract group components from NGCA
  est.group = mlcaFP(xData = whitener(LC.conca, n.comp = n.group)$Z, n.comp = n.group, 
                     whiten = whiten, restarts.pbyd = restarts.pbyd,  restarts.dbyd =  restarts.dbyd, 
                     distribution=distribution)
  
  est.comp = order.likelihood(est.group$S, distribution = distribution)
  
  est.comp
} 

group_ICA = function(data.list, var.ratio.kept = 0.82, n.group,
                     restarts = 30, eps = 1e-6){
  require(steadyICA)
  nsub = length(data.list)
  
  # subject level pca and concatenate data
  ratio.var = sapply(data.list, function(x){a = (svd(x$X)$d)^2; cumsum(a)/sum(a)})
  l = max(which(apply(ratio.var, 1, min) < var.ratio.kept)) + 1
  
  PC.conca = NULL
  #for (i in 1:20) PC.conca = cbind(PC.conca, simData.high[[i]]$whitened.X[,1:l]) # to do
  for (i in 1:nsub) PC.conca = cbind(PC.conca, whitener(data.list[[i]]$X, n.comp = l)$Z) 
  
  # extract group components
  est.group <- infomaxICA(X = whitener(PC.conca, n.comp = n.group)$Z, n.comp = n.group, 
                          whiten = TRUE, restarts=restarts, eps = eps)
  est.group$S
}
