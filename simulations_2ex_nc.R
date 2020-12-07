
##empty parameters for simulation results
ivest1 = ivest2 = NULL
olsest1 = olsest2 = NULL
tsest1 = tsest2 = NULL
F1 = F2 = F1_iv = F2_iv = NULL
correlation = NULL
QIVW1_LIML = QIVW1 = QIVW1_LIMLnoc = QIVW1_noc = NULL
QIVW2_LIML = QIVW2 =  QIVW2_LIMLnoc = QIVW2_noc =  NULL
xx11 = xx12 = NULL
xx21 = xx22 = NULL
limlest1 =  limlest2 = limlest1_noc =  limlest2_noc =NULL
limlhet1 = limlhet2 = limlhet1_noc = limlhet2_noc =NULL
tau = taunoc = NULL
Qexact = Qexact_het = Qexactnoc = Qexact_hetnoc = NULL
Qa = Qa_novar = NULL

##SIMULATIONS

for(j in 1:reps){
  
  G = matrix(nrow = n, ncol=L)
  G[,]= rbinom((n*L),2,0.5)
  c1 = rnorm(n,1,2)
  c2 = rnorm(n,1,2)
  
  V = mvrnorm(n, mu, s)
  
  x1 = G%*%pi1 + c1  + V[,2]
  x2 = G%*%pi2 + c2  + V[,3]
  x3 = G%*%alpha  #new
  y = 0.5*x1 - 0.3*x2 + x3 + 0.5*c1 + 0.5*c2+  V[,1]
  
  x = cbind(x1, x2)
  
  ##Second model for two sample estimation
  G2 = matrix(nrow = n, ncol=L)
  G2[,]= rbinom((n*L),2,0.5)
  c22 = rnorm(n,1,2)
  c21 = rnorm(n,1,2)
  
  V2 = mvrnorm(n, mu, s)
  
  x12 = G2%*%pi1  + c21  + V2[,2]
  x22 = G2%*%pi2  + c22  + V2[,3]
  x32 = G2%*%alpha#new
  y2 = 0.5*x12 - 0.3*x22 +  x32 + 0.5*c21 + 0.5*c22 + V2[,1]
  
  ################
  ###Fit the model
  ################
  ##IV estimation 
  ################
  
  #one sample
  
  ivest = summary(ivreg(y~ x1 + x2 | G))
  
  ivest1[j] = ivest$coefficients[2,1]
  ivest2[j] = ivest$coefficients[3,1]
  
  
  F1[j] = summary(lm(x1~ G))$fstatistic[1]
  F2[j] = summary(lm(x2~ G))$fstatistic[1]
  
  
  reg1 = ivreg(x1~ x2 |G)
  uhat = reg1$residuals
  F1_iv[j] = summary(lm(uhat~ G))$fstatistic[1]
  
  reg2 = ivreg(x2~ x1 |G)
  uhat2 = reg2$residuals
  F2_iv[j] = summary(lm(uhat2~ G))$fstatistic[1]
  
  xx = (1/n)*(t(cbind(x1,x2))%*%cbind(x1,x2))
  xx11[j] = xx[1,1]
  xx12[j] = xx[1,2]
  xx21[j] = xx[2,1]
  xx22[j] = xx[2,2]
  
  olsest = summary(lm(y~ x1 + x2))
  olsest1[j] = olsest$coefficients[2,1]
  olsest2[j] = olsest$coefficients[3,1]
  
  
  
  #two sample estimation
#################################
  
  #indepdendently estimated effects of each SNP
  
  
 pihat = stderr = matrix(nrow = L, ncol = kx)
  gammahat = segamma = NULL
 
  for(k in 1:L){
    
    gammahat[k] = summary(lm(y2~ G2[,k]))$coefficient[2,1]
    segamma[k] = summary(lm(y2~ G2[,k]))$coefficient[2,2]
    
   for(ex in 1:kx){
    
   pihat[k,ex] = summary(lm(x[,ex]~ G[,k]))$coefficient[2,1]
   stderr[k,ex] = summary(lm(x[,ex]~ G[,k]))$coefficient[2,2]
    
  }
  } 
  
  correlation = matrix(c(1, cor(x1,x2), cor(x1,x2), 1), nrow = 2, ncol = 2)
  cov <- as.vector(rep(cor(x1,x2),L))*stderr[,1]*stderr[,2]
  
  tsest_snp = summary(lm(gammahat~  -1 + pihat, weights=(1/segamma^2)))
  tsest1[j] = tsest_snp$coefficients[1,1]
  tsest2[j] = tsest_snp$coefficients[2,1]
  
  
  fitx1 = lm(pihat[,1] ~ -1 + pihat[,2])
  delta12 = fitx1$coef[1]
  VX1 = (stderr[,1]^2)- 2*delta12*(cov) + (delta12^2)*(stderr[,2]^2)
  QIVW1[j] = sum((1/VX1)*(pihat[,1] - (delta12*pihat[,2]))^2)
  
 
   fitx2 = lm(pihat[,2] ~ -1 + pihat[,1])
  delta22 = fitx2$coef[1]
  VX2 = (stderr[,2]^2)- 2*delta22*(cov) + (delta22^2)*(stderr[,1]^2)
  QIVW2[j] = sum((1/VX2)*(pihat[,2] - (delta22*pihat[,1]))^2)
  
  
  ###F_TS estimated using LIML method
  QIVW1_LIML[j] <- optimize(Q_ind, interval = c(-100,100))$objective
  QIVW2_LIML[j] <- optimize(Q_ind2, interval = c(-100,100))$objective
  
  
  ##estimation of Qa
  
  betaest = rbind(tsest1[j], tsest2[j])
  w <- segamma^2
  Qa_novar[j] =  sum((1/w)*((gammahat - pihat%*%betaest)^2))
 
  
   varcov = matrix(nrow = kx, ncol = kx)
  w=NULL
  for(l in 1:L){
    for(pp in 1:kx){
      for(p2 in 1:kx){
        varcov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
      }
    }
    
    w[l] <- segamma[l]^2+t(betaest)%*%varcov%*%betaest     
  }
  
  Qa[j] =  sum((1/w)*((gammahat - pihat%*%betaest)^2))
  
  
  
  
  #LIML estimation of the betas.
  
  liml_beta<- optim(c(0,0), LIML_b)
  
  limlest1[j] <- liml_beta$par[1]
  limlest2[j] <- liml_beta$par[2]
  
  Qexact[j] <- liml_beta$value
  
  #estimate the parameters with LIML - with heterogeneity statistic
  
  limltauest  = optimize(PL_MVMR,interval=c(-10,10))
  tau_i = limltauest$minimum
  
  tau[j] = tau_i
  
  liml_het2<- optim(c(0,0), PL2_MVMR)
  limlhet1[j] <- liml_het2$par[1]
  limlhet2[j] <- liml_het2$par[2]
  Qexact_het[j] <- liml_het2$value
  
  
  #no covariances
  correlation = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  cov <- as.vector(rep(0,L))*stderr[,1]*stderr[,2]
  
  
  #no covariances:
  VX1_noc = (stderr[,1]^2) + (delta12^2)*(stderr[,2]^2)
  QIVW1_noc[j] = sum((1/VX1_noc)*(pihat[,1] - (delta12*pihat[,2]))^2)
  
  VX2_noc = (stderr[,2]^2) + (delta22^2)*(stderr[,1]^2)
  QIVW2_noc[j] = sum((1/VX2_noc)*(pihat[,2] - (delta22*pihat[,1]))^2)
  
  ###F_TS estimated using LIML method
  QIVW1_LIMLnoc[j] <- optimize(Q_ind, interval = c(-100,100))$objective
  QIVW2_LIMLnoc[j] <- optimize(Q_ind2, interval = c(-100,100))$objective
  

  liml_betanoc<- optim(c(0,0), LIML_b)
  
  limlest1_noc[j] <- liml_betanoc$par[1]
  limlest2_noc[j] <- liml_betanoc$par[2]
  
  Qexactnoc[j] <- liml_betanoc$value
  
  limltauest  = optimize(PL_MVMR,interval=c(-10,10))
  tau_i = limltauest$minimum
  
  taunoc[j] = tau_i
  
  liml_het2<- optim(c(0,0), PL2_MVMR)
  limlhet1_noc[j] <- liml_het2$par[1]
  limlhet2_noc[j] <- liml_het2$par[2]
  Qexact_hetnoc[j] <- liml_het2$value
  
  
  
  
  
 
}
