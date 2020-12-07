
rm(list = ls(all=TRUE))

library(AER)
library(MASS)

args <- commandArgs(T)
set.seed(as.numeric(args[1]))
job_id <- (as.numeric(args[1]))
message("job number ", job_id)


source("functions_flexnoX.R")

#Set sample size and number of snps for this simulation
ss = 20000    #sample size
givs = 200 #number of snps
reps = 100

n=ss
L = givs
kx = 2

output <- data.frame()

#simulations with strong instruments 
for(tau_true in c(0,0.5)){
  
  set.seed(10)
  pimax = 2
  delta = 0.7
  mu = c(0,0,0)
  s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
  pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
  pi1 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
  alpha = as.matrix(c(rnorm(L/2,0,tau_true), rep(0,L/2)),nrow=L,ncol=1) 
  
  source("simulations_2ex.R")
  
  
  output <- rbind(output, data.frame(L, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                     tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                     QIVW1_LIML, QIVW1, QIVW2_LIML, QIVW2, Qa, Qa_novar,
                                     Qexact, Qexact_het,
                                     limlest1, limlest2, limlhet1, limlhet2, tau, QIVW1_LIMLnoc, QIVW1_noc, QIVW2_LIMLnoc, QIVW2_noc,
                                     Qexactnoc, Qexact_hetnoc,
                                     limlest1_noc, limlest2_noc, limlhet1_noc, limlhet2_noc, taunoc))
  
}


message("filesave", sprintf("outputstrong%s.Rda", job_id))

save(output, file=sprintf("outputstrong%s.Rda", job_id))
