
rm(list = ls(all=TRUE))

ptm <- proc.time()

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
reps = 50


output <- data.frame()

for(givs in c(50,100,200)){

##Simulations with weak instruments 

pimax = 0.1

for(tau_true in c(0.5)){

set.seed(as.numeric(args[1]))
  
n=ss
L = givs
sim = 1
kx = 2
delta = 0.7
mu = c(0,0,0)
s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
pi1 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
alpha = as.matrix(c(rnorm(L/2,0,tau_true), rep(0,L/2)),nrow=L,ncol=1)  


source("simulations_2ex_jk.R")

output <- rbind(output, cbind(data.frame(sim, givs, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                         tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                        QIVW1, QIVW2, Qa, Qa_novar,
                                      Qexact, Qexact_het,
                                       limlest1, limlest2, limlhet1, limlhet2, tau), jackknife))

}

#simulations with conditionally weak instruments
for(tau_true in c(0.5)){
  
set.seed(as.numeric(args[1]))
pimax = 2
sim = 2
delta = 0.7
L = givs
mu = c(0,0,0)
s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
pi1 = delta*pi2 + c(0,rnorm((L-1),2.8/sqrt(n),0.075))
alpha = as.matrix(c(rnorm(L/2,0,tau_true), rep(0,L/2)),nrow=L,ncol=1)  

source("simulations_2ex_jk.R")


output <- rbind(output, cbind(data.frame(sim, givs, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                   tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                  QIVW1, QIVW2, Qa, Qa_novar,
                                   Qexact, Qexact_het,
                                   limlest1, limlest2, limlhet1, limlhet2, tau), jackknife))

}

}

message("filesave", sprintf("outputjackknife%s.Rda", job_id))
save(output, file=sprintf("outputjackknife%s.Rda", job_id))

time_min <- (proc.time() - ptm)[3]/60
time_min