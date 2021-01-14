# generate the simulation data for a single subject
SimFMRI.ngca = function(snr, gamma.rate, gamma.shape,
                        FWHM = 9,
                        noisyICA=FALSE, nNG=25, nTR=50, phi=0.37, dim.data=c(33,33), var.inactive=0.001) {
  if (snr == 'high'){
    snr.ratio = c(0.335, 0.299, 0.366)
  } 
  if (snr == 'low') {
    snr.ratio = c(0.017, 0.460, 0.523)
  }
  if (snr == 'medium'){
    snr.ratio = c(0.176, 0.386, 0.438)
  }
  snr.group = c(0.154,0.298,0.548)
  require(neuRosim)
  require(steadyICA)
  #Latent group components are fixed for each simulation:
  x1 = rep(3,5)
  y1 = c(3:7)
  s1.coords = cbind(x1,y1)
  s1 = specifyregion(dim = dim.data, coord = s1.coords, form = "manual")
  s1[s1!=0] = seq(0.5,1,length=length(x1))
  
  x2 = c(8,8,8,9,10,9,10,10,10,9,8)
  y2 = c(15,14,13,13,13,15,15,16,17,17,17)
  s2.coords = cbind(c(x2,x2+7),c(y2,y2))
  s2 = specifyregion(dim=dim.data, coord = s2.coords, form = 'manual')
  s2[s2!=0] = seq(0.5,1,length=2*length(x2))
  
  x3 = c(13,14,15,15,15,14,13,15,15,14,13)
  y3 = c(19,19,19,20,21,21,21,22,23,23,23)
  s3.coords = cbind(c(x3,x3+7,x3+14),c(y3,y3,y3))
  s3 = specifyregion(dim=dim.data, coord = s3.coords, form = 'manual')
  s3[s3!=0] = seq(0.5,1,length=3*length(x3))
  
  sim.S.group = cbind(as.vector(s1),as.vector(s2),as.vector(s3))
  
  #Latent individual components from random Gamma field 
  sim.S.indi = spatialnoise(dim=dim.data, sigma=1, nscan = nNG-3, method = "gammaRF", FWHM = FWHM,  gamma.shape = gamma.shape, gamma.rate = gamma.rate)
  dim(sim.S.indi) = c(prod(dim.data),nNG-3)
  
  ## Add small amount of Gaussian noise to inactive voxels
  nInactive = sum(sim.S.group == 0)
  baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  sim.S.group[sim.S.group==0] = baseline
  
  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-nNG
  sim.GRF = NULL
  
  sim.GRF <- spatialnoise(dim = dim.data, sigma=1, FWHM = FWHM, nscan = nscan,  method = 'gaussRF')
  dim(sim.GRF) <- c(prod(dim.data),nscan)
  
  ##Mixmat:
  #create timecourses for latent components: 
  sim.Ms = matrix(rnorm(nNG*nTR,0,1),nrow=nNG,ncol=nTR)
  for(t in 2:nTR) sim.Ms[,t] = phi*sim.Ms[,t-1] + sim.Ms[,t]
  # different ratio for three different group components
  # first fix the variance ratio among three group components
  # then let the sum of group variance reach the desirable level
  sim.Xs.group = sim.S.group%*%sim.Ms[1:3,]
  for (i in 1:3){
    sim.Xs.group[,i] = sim.Xs.group[,i]/sd(as.vector(sim.Xs.group[,i])) 
    sim.Xs.group[,i] = sim.Xs.group[,i] * sqrt(snr.group[i])
  } 
  sim.Xs.group = sim.Xs.group/sd(as.vector(sim.Xs.group)) #standardize so we can control SNR
  sim.Xs.group = sim.Xs.group * sqrt(snr.ratio[1])
  
  sim.Xs.indi = sim.S.indi%*%sim.Ms[4:nNG,]
  sim.Xs.indi = sim.Xs.indi/sd(as.vector(sim.Xs.indi)) #standardize so we can control SNR
  sim.Xs.indi = sim.Xs.indi * sqrt(snr.ratio[2])
  
  sim.S = cbind(sim.S.group, sim.S.indi)
  
  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) #standardize so we can control SNR
  sim.Xn = sim.Xn * sqrt(snr.ratio[3])
  
  sim.Xs = sim.Xs.group + sim.Xs.indi
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  
  if(noisyICA) { 
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener, 
                sim.Xs = sim.Xs, sim.Xn = sim.Xn))
  }
}