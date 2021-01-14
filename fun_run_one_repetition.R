# one repetition with repetition seed j 
run_comp_est_OneRep = function(nsub = 20, gamma.shape = .02, gamma.rate = .0001, FWHM = 9, var.inactive = 0.001, seed = 1, 
                               n.individual.NGCA, n.group=3, parallel_NGCA = FALSE){
  # n.individual.NGCA should have dimension size nsub*3
  cor.est = array(0, dim = c(n.group,3,2,3), 
                  dimnames = list(NULL,
                                  comp = c('low var comp', 'medium var comp', 'high var comp'),
                                  method = c('group GNCA','group ICA'), 
                                  SVAR = c('High SVAR', 'Medium SVAR', 'Low SVAR')))
  comp.ngca = list()
  comp.gift = list()
  SVAR.setting = c('high','medium','low')
  for (i in 1:3){
    est = run_comp_est_OneSVAR(SVAR.setting = SVAR.setting[i], 
                               n.individual.NGCA = n.individual.NGCA[,i], n.group = n.group,
                               nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, var.inactive = var.inactive,
                               seed = seed, parallel_NGCA = parallel_NGCA)
    
    cor.est[,,,i] = est$corr
    comp.ngca[[i]] = est$comp.ngca
    comp.gift[[i]] = est$comp.gift
  }
  
  names(comp.ngca) = c('high','medium','low')
  names(comp.gift) = c('high','medium','low')
  
  res = list(corr = cor.est, comp.ngca = comp.ngca, comp.gift = comp.gift)
  res
} 

run_comp_est_OneSVAR = function(SVAR.setting, n.individual.NGCA = rep(15, nsub), n.group,
                                nsub = 20, gamma.shape = .02, gamma.rate = .0001, FWHM = 9, seed = 1, var.inactive = 0.001,
                                parallel_NGCA = FALSE){
  cor.est = array(0, dim = c(n.group, 3, 2), 
                  dimnames = list(method = NULL, c('low var comp', 'medium var comp', 'high var comp'), 
                                  c('group GNCA','group ICA')))
  
  # data generation 
  simData = list()
  for (i in 1:nsub){
    set.seed(i + (seed-1)*nsub)
    simData[[i]] = SimFMRI.ngca(snr = SVAR.setting, 
                                gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM,var.inactive = var.inactive)
  }
  
  # group NGCA
  est_group_NGCA = group_LNGCA(simData, n.individual = n.individual.NGCA, n.group = n.group, parallel = parallel_NGCA)
  comp.ngca = order.likelihood(est_group_NGCA$S, distribution = 'logistic')$S
  #temp <- matchICA(S = est_group_NGCA$S, template = simData[[1]]$S[,1:3])
  cor.est[,,1] = cor(comp.ngca, simData[[1]]$S[,1:3], method = "pearson")
  
  # group ICA
  est_group_ICA = group_ICA(simData, n.group = n.group)
  comp.gift = order.likelihood(est_group_ICA, distribution = 'logistic')$S
  #temp <- matchICA(S = est_group_ICA, template = simData[[1]]$S[,1:3])
  cor.est[,,2] = cor(comp.gift, simData[[1]]$S[,1:3], method = "pearson")
  
  res = list(corr = cor.est, comp.ngca = comp.ngca, comp.gift = comp.gift)
}

run_dim_est_OneRep = function(nsub = 20, gamma.rate = 1e-4, gamma.shape = .02, FWHM = 9, var.inactive = 0.001,
                              seed = 1, run.method = 1:3, parallel = FALSE){
  dim.est = array(0, dim = c(nsub ,3, 3), 
                  dimnames = list(NULL, 
                                  SVAR = c('High SVAR', 'Medium SVAR', 'Low SVAR'), 
                                  method = c('FOBI-GRF','FOBIasymp','FOBIboot')))
  
  SVAR.setting = c('high','medium','low')
  for (i in 1:3){
    dim.est[,i,] = run_dim_est_OneSVAR(SVAR.setting = SVAR.setting[i], nsub = nsub, 
                                       gamma.shape = gamma.shape, gamma.rate = gamma.rate, FWHM = FWHM, var.inactive = var.inactive,
                                       seed = seed, run.method = run.method, 
                                       parallel = parallel)
    
    print(paste('finish iteration ',SVAR.setting[i]))
  }
  dim.est
} 

run_dim_est_OneSVAR = function(SVAR.setting, nsub = 20, 
                               gamma.shape = .02, gamma.rate = .0001, FWHM = 9, var.inactive = 0.001,
                               seed = 1, run.method = 1:3, 
                               parallel = FALSE){
  dim.est = array(-1, dim = c(nsub,3), 
                  dimnames = list(NULL, SVAR = c('FOBI-GRF', 'FOBIasymp', 'FOBIboot')))
  
  # data generation 
  simData = list()
  for (i in 1:nsub){
    set.seed(i + (seed-1)*nsub)
    simData[[i]] = SimFMRI.ngca(snr = SVAR.setting, 
                                gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM,
                                var.inactive = var.inactive)
  }
  
  
  if (parallel){
    res = foreach (i = 1:nsub)%dopar%{
      est = numeric(3)
      # FOBI GRF
      if (1 %in% run.method){
        est[1] = BIsearch_test(simData[[i]]$X, fun = FOBI_GRF_test, FWHM = FWHM)$dim_est
      } 
      
      # FOBIasymp
      if (2 %in% run.method){
        require(ICtest)
        est[2] = BIsearch_test(simData[[i]]$X, fun = FOBIasymp)$dim_est
      }
      
      # FOBIboot
      if (3 %in% run.method){
        require(ICtest)
        est[3] = BIsearch_test(simData[[i]]$X, fun = FOBIboot)$dim_est
      }
      
      est
    }
    
    for (i in 1:nsub) dim.est[i,] = res[[i]]
    
  }else{
    for (i in 1:nsub){
      # FOBI GRF
      if (1 %in% run.method){
        dim.est[i,1] = BIsearch_test(simData[[i]]$X, fun = FOBI_GRF_test, FWHM = FWHM)$dim_est
      } 
      
      # FOBIasymp
      if (2 %in% run.method){
        require(ICtest)
        dim.est[i,2] = BIsearch_test(simData[[i]]$X, fun = FOBIasymp)$dim_est
      }
      
      # FOBIboot
      if (3 %in% run.method){
        require(ICtest)
        dim.est[i,3] = BIsearch_test(simData[[i]]$X, fun = FOBIboot)$dim_est
      }
    }
  }
  
  dim.est
}





