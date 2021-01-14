# Settings
# set the repository as working directory
setwd("~/Group-Linear-NonGaussian-Component-Analysis-for-Neuroimaging")
library(doParallel)
library(steadyICA)
library(ICtest)
# The simulation in conducted using 20 cores
registerDoParallel(20)
getDoParWorkers()

source('fun_NGCA.R')
source('fun_simNGCA.R')
source('fun_run_one_repetition.R')
source('fun_call_group_function.R')
source('fun_dimension_test.R')

gamma.rate = 1e-4
gamma.shape = 0.02
nsub = 20
nrep = 40
FWHM = 9 

# Dimension estimation 
# Each repetition contains 20 subjects, and thus 40 repetition contain 800 subjects in total
# Therefore, this part of simulation takes long time to finish
res_dimest = array(0, dim = c(nrep,nsub,3,3), dimnames = list(NULL,NULL,SVAR = c('High SVAR', 'Medium SVAR', 'Low SVAR'), method = c('FOBI-GRF (our)','FOBIasymp','FOBIboot')))
for (i in 1:nrep){
  res_dimest[i,,,] = run_dim_est_OneRep(nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, seed = i, run.method = 1:3, parallel = TRUE)
  print(paste('finish rep ', i))
}
#save(res_dimest, file = 'res_dimest.RData')


# Components estimation with 3 group components
result_correst_3comp = foreach(j = 1:nrep) %dopar% {
  run_comp_est_OneRep(nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, seed = j, 
                      n.individual.NGCA = result_dimest[j,,,1], n.group=3) # use 3 group components 
}

result_correst_4comp = foreach(j = 1:nrep) %dopar% {
  run_comp_est_OneRep(nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, seed = j, 
                      n.individual.NGCA = result_dimest[j,,,1], n.group=4) # use 4 group components 
}

result_correst_2comp = foreach(j = 1:nrep) %dopar% {
  run_comp_est_OneRep(nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, seed = j, 
                      n.individual.NGCA = result_dimest[j,,,1], n.group=2) # use 2 group components 
}







