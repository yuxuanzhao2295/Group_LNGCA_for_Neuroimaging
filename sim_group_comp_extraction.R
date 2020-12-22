# Settings
setwd("/Users/yuxuan/Documents/GitHub/Group-Linear-NonGaussian-Component-Analysis-for-Neuroimaging")
library(doParallel)
library(steadyICA)
library(ICtest)
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
nrep = 100
nrep = 40
FWHM = 9 

# 
result_dimest_0715 = array(0, dim = c(nrep,nsub,3,3), dimnames = list(NULL,NULL,SVAR = c('High SVAR', 'Medium SVAR', 'Low SVAR'), 
                                                                  method = c('FOBI-GRF (our)','FOBIasymp','FOBIboot')))
result = foreach(j = 1:nrep)%dopar%{
  run_dim_est_OneRep(nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, seed = j, run.method = 1:3)
}
for (i in 1:nrep) result_dimest_0715[i,,,] = result[[i]]
save(result_dimest_0715, file = 'result_dimest_0715.RData')
apply(result_dimest_0715,c(2,3,4),median)

result_correst_3comp = array(0, dim = c(nrep,2,3,3), 
                                  dimnames = list(NULL, method = c('group NGCA', 'group ICA'), 
                                                  SVAR = c('high SVAR', 'medium SVAR', 'low SVAR'),
                                                  Comp = c('Comp 1', 'Comp 2', 'Comp 3')))
result = foreach(j = 1:nrep) %dopar% {
  run_comp_est_OneRep(nsub = nsub, gamma.rate = gamma.rate, gamma.shape = gamma.shape, FWHM = FWHM, seed = j, 
                      n.individual.NGCA = result_dimest_0715[j,,,1], n.group=3) # use 3 group components 
}
for (i in 1:nrep) result_correst[i,,,] = result[[i]]




