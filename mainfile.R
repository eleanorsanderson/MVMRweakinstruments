
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
reps = 1000


output <- data.frame()

##Simulations with weak instruments 


pimax = 0.1

for(tau_true in c(0,0.5)){

  set.seed(10)
  
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


source("simulations_2ex.R")


output <- rbind(output, data.frame(sim, L, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                   tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                   QIVW1_LIML, QIVW1, QIVW2_LIML, QIVW2, Qa, Qa_novar,
                                   Qexact, Qexact_het,
                                   limlest1, limlest2, limlhet1, limlhet2, tau, QIVW1_LIMLnoc, QIVW1_noc, QIVW2_LIMLnoc, QIVW2_noc,
                                   Qexactnoc, Qexact_hetnoc,
                                   limlest1_noc, limlest2_noc, limlhet1_noc, limlhet2_noc, taunoc))
}

#simulations with conditionally weak instruments
for(tau_true in c(0,0.5)){
  
  set.seed(10)
pimax = 2
sim = 2
delta = 0.7
mu = c(0,0,0)
s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
pi1 = delta*pi2 + c(0,rnorm((L-1),2.8/sqrt(n),0.075))
alpha = as.matrix(c(rnorm(L/2,0,tau_true), rep(0,L/2)),nrow=L,ncol=1)  

source("simulations_2ex.R")


output <- rbind(output, data.frame(sim, L, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                   tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                   QIVW1_LIML, QIVW1, QIVW2_LIML, QIVW2, Qa, Qa_novar,
                                   Qexact, Qexact_het,
                                   limlest1, limlest2, limlhet1, limlhet2, tau, QIVW1_LIMLnoc, QIVW1_noc, QIVW2_LIMLnoc, QIVW2_noc,
                                   Qexactnoc, Qexact_hetnoc,
                                   limlest1_noc, limlest2_noc, limlhet1_noc, limlhet2_noc, taunoc))

}


#simulations with strong instruments 
for(tau_true in c(0,0.5)){
  
  set.seed(10)
  pimax = 2
  sim = 2.2
  delta = 0.7
  mu = c(0,0,0)
  s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
  pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
  pi1 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
  alpha = as.matrix(c(rnorm(L/2,0,tau_true), rep(0,L/2)),nrow=L,ncol=1) 
  
  source("simulations_2ex.R")
  
  
  output <- rbind(output, data.frame(sim, L, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                     tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                     QIVW1_LIML, QIVW1, QIVW2_LIML, QIVW2, Qa, Qa_novar,
                                     Qexact, Qexact_het,
                                     limlest1, limlest2, limlhet1, limlhet2, tau, QIVW1_LIMLnoc, QIVW1_noc, QIVW2_LIMLnoc, QIVW2_noc,
                                     Qexactnoc, Qexact_hetnoc,
                                     limlest1_noc, limlest2_noc, limlhet1_noc, limlhet2_noc, taunoc))
  
}


message("filesave", sprintf("output2ex%s.csv", job_id))

write.csv(output, file=sprintf("output2ex%s.csv", job_id))

#simulations for table 3
output3 <- data.frame()

pimax = 1
n=ss
L = givs
sim = 3
kx = 3
mu = c(0,0,0)
s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
pi1 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
pi3 = 0.1*pi1+0.9*as.matrix(c(runif((L),0,pimax/5)),nrow=L,ncol=n)


source("simulations_3ex.R")


output3 <- rbind(output3, data.frame(sim, L, tau_true, pimax, ivest1, ivest2, olsest1, olsest2,
                                   tsest1, tsest2, F1, F2, F1_iv, F2_iv, 
                                   QIVW1_LIML, QIVW1, QIVW2_LIML, QIVW2, Qa, Qa_novar,
                                   Qexact, Qexact_het,
                                   limlest1, limlest2, limlhet1, limlhet2, tau, 
                                   ivest1_3, ivest2_3, ivest3_3, olsest1_3, olsest2_3, olsest3_3,
                                   tsest1_3, tsest2_3, tsest3_3, F1_3, F2_3, F3_3, F1_iv_3, F2_iv_3, F3_iv_3, 
                                   QIVW1_3, QIVW2_3, QIVW3_3, Qa3, Qa_novar3,
                                   Qexact_3, Qexact_het_3,
                                   limlest1_3, limlest2_3, limlest3_3, limlhet1_3, limlhet2_3, limlhet3_3, tau_3))





message("filesave", sprintf("output3ex%s.csv", job_id))

write.csv(output3, file=sprintf("output3ex%s.csv", job_id))


