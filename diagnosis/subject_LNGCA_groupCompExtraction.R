nsub = 20
gamma.rate = 1e-4
gamma.shape = 0.02
FWHM = 9
nTR = 100
SVAR.setting = 'low'
nNG = 25

seed = 13
simData = list()
for (i in 1:nsub){
  set.seed(i + (seed-1)*nsub)
  simData[[i]] = SimFMRI.ngca(snr = SVAR.setting, 
                              gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM,
                              nNG = nNG, nTR=nTR)
}

resComp = list()
for (i in 1:nsub){
  resComp[[i]] = mlcaFP(xData = simData[[i]]$X, restarts.pbyd = 0,  restarts.dbyd =  30, 
                        distribution='logistic')
}


ranks = array(0, dim = c(20,3))
cors = array(0, dim = c(20,3))
for (i in 1:20){
  for  (j in 1:3){
    out= order.likelihood(resComp[[i]]$S, distribution = 'logistic', out.loglik = TRUE)
    ranks[i,j] = which.max(abs(cor(out$S, simData[[i]]$S[,j])))
    cors[i,j] = max(abs(cor(out$S, simData[[i]]$S[,j])))
  }
}

cors[,1] # all very high, never miss comp 1 if sepcifying 100 comps
ranks[,1] # the component with the highest correlation with comp 1 always has the top 25 non-gaussianity
ranks[,2] - ranks[,1] # for most subjects, the component matched to comp 2 has smaller non-gaussianity than that matched to comp 1

resComp_short = list()
for (i in 1:nsub){
  resComp_short[[i]] = mlcaFP(xData = simData[[i]]$X, n.comp=25, restarts.pbyd = 0,  restarts.dbyd =  30, 
                              distribution='logistic')
}

ranks_short = array(0, dim = c(20,3))
cors_short = array(0, dim = c(20,3))
for (i in 1:20){
  for  (j in 1:3){
    out= order.likelihood(resComp_short[[i]]$S, distribution = 'logistic', out.loglik = TRUE)
    ranks_short[i,j] = which.max(abs(cor(out$S, simData[[i]]$S[,j])))
    cors_short[i,j] = max(abs(cor(out$S, simData[[i]]$S[,j])))
  }
}

# nevertheless, when specifying the number of NG components as 25, we miss comp 1 more ofer 
sum(cors_short[,1]<0.5) # 8
sum(cors_short[,2]<0.5) # 5
