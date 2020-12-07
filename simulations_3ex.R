
##empty parameters for simulation results
ivest1 = ivest2 = NULL
olsest1 = olsest2 = NULL
tsest1 = tsest2 =  NULL
F1 = F2 =  F1_iv = F2_iv = NULL
correlation = NULL
QIVW1_LIML = QIVW1 = QIVW1_LIMLnoc = QIVW1_noc = NULL
QIVW2_LIML = QIVW2 =  QIVW2_LIMLnoc = QIVW2_noc =  NULL
limlest1 =  limlest2 =limlest1_noc =  limlest2_noc = NULL
limlhet1 = limlhet2 = limlhet1_noc = limlhet2_noc = NULL
tau = taunoc = NULL
Qexact = Qexact_het = Qexactnoc = Qexact_hetnoc = NULL
Qa = Qa_novar = NULL

ivest1_3 = ivest2_3 = ivest3_3 = NULL
olsest1_3 = olsest2_3 = olsest3_3 = NULL
tsest1_3 = tsest2_3 = tsest3_3 = NULL
F1_3 = F2_3 = F3_3 =  F1_iv_3 = F2_iv_3 = F3_iv_3 = NULL
correlation = NULL
QIVW1_3 = QIVW2_3 = QIVW3_3 = NULL

limlest1_3 =  limlest2_3 =limlest3_3 = NULL
limlhet1_3 = limlhet2_3 =limlhet3_3= NULL
tau_3 = NULL
Qexact_3 = Qexact_het_3 = NULL
Qa3 = Qa_novar3 = NULL

##SIMULATIONS

for(j in 1:reps){
  
  G = matrix(nrow = n, ncol=L)
  G[,]= rbinom((n*L),2,0.5)
  c1 = rnorm(n,2,4)
  c2 = rnorm(n,2,4)
  c3 = rnorm(n,2,4)
  
  V = mvrnorm(n, mu, s)
  
  x1 = G%*%pi1 + 0.9*c1 + 0.1*c2 + V[,2]
  x2 = G%*%pi2 + 0.9*c2 + 0.1*c1 + V[,3]
  x3 = G%*%pi3 + c3 + rnorm(n,0,1)  
  y = 0.5*x1 - 0.3*x2 + 0.7*x3 + c1 + c2 + c3 +  V[,1]
  
  x = cbind(x1, x2, x3)
  
  ##Second model for two sample estimation
  G2 = matrix(nrow = n, ncol=L)
  G2[,]= rbinom((n*L),2,0.5)
  c22 = rnorm(n,2,4)
  c21 = rnorm(n,2,4)
  c32 = rnorm(n,2,4)
  
  V2 = mvrnorm(n, mu, s)
  
  x12 = G2%*%pi1  + 0.9*c21 + 0.1*c22  + V2[,2]
  x22 = G2%*%pi2  +  0.9*c22 + 0.1*c21 + V2[,3]
  x32 = G2%*%pi3 + c32 + rnorm(n,0,1)  
  y2 = 0.5*x12 - 0.3*x22 + 0.7*x32 + c21 + c22 + c32 + V[,1]
  
  
  ################
  ###Fit the model
  ################
  ##IV estimation - two exposures
  ################
  
  kx=2
  
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
  
  
 pihat = stderr = matrix(nrow = L, ncol = 2)
  gammahat = segamma = NULL
 
  for(k in 1:L){
    
    gammahat[k] = summary(lm(y2~ G2[,k]))$coefficient[2,1]
    segamma[k] = summary(lm(y2~ G2[,k]))$coefficient[2,2]
    
   for(ex in 1:2){
    
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
  
  
  #estimation of Q
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
  
##########################################  
##estimation including all three exposures
##########################################
  kx=3
  
  #one sample
  
  ivest = summary(ivreg(y~ x1 + x2 +x3 | G))
  
  ivest1_3[j] = ivest$coefficients[2,1]
  ivest2_3[j] = ivest$coefficients[3,1]
  ivest3_3[j] = ivest$coefficients[4,1]
  
  F1_3[j] = summary(lm(x1~ G))$fstatistic[1]
  F2_3[j] = summary(lm(x2~ G))$fstatistic[1]
  F3_3[j] = summary(lm(x3~ G))$fstatistic[1]
  
  reg1 = ivreg(x1~ x2 + x3|G)
  uhat = reg1$residuals
  F1_iv_3[j] = summary(lm(uhat~ G))$fstatistic[1]
  
  reg2 = ivreg(x2~ x1 + x3|G)
  uhat2 = reg2$residuals
  F2_iv_3[j] = summary(lm(uhat2~ G))$fstatistic[1]
  
  reg3 = ivreg(x3~ x1 + x2|G)
  uhat3 = reg3$residuals
  F3_iv_3[j] = summary(lm(uhat3~ G))$fstatistic[1]
  
    olsest = summary(lm(y~ x1 + x2 + x3))
  olsest1_3[j] = olsest$coefficients[2,1]
  olsest2_3[j] = olsest$coefficients[3,1]
  olsest3_3[j] = olsest$coefficients[4,1]
  
  
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
  
  
  correlation = matrix(c(1, cor(x1,x2), cor(x1,x3), cor(x1,x2), 1, cor(x2,x3), cor(x1,x3), cor(x2,x3), 1), nrow = 3, ncol = 3)
  cov12 = correlation[1,2]*stderr[,1]*stderr[,2]
  cov13 = correlation[1,3]*stderr[,1]*stderr[,3]
  cov23 = correlation[2,3]*stderr[,2]*stderr[,3]
  
  tsest_snp = summary(lm(gammahat~  -1 + pihat, weights=(1/segamma^2)))
  tsest1_3[j] = tsest_snp$coefficients[1,1]
  tsest2_3[j] = tsest_snp$coefficients[2,1]
  tsest3_3[j] = tsest_snp$coefficients[3,1]
  
  fitx1 = lm(pihat[,1] ~ -1 + pihat[,2] + pihat[,3])
  delta12 = fitx1$coef[1]
  delta13 = fitx1$coef[2]
  VX1 = (stderr[,1]^2)- 2*delta12*(cov12) - 2*delta13*(cov13) + (delta12^2)*(stderr[,2]^2) + (delta13^2)*(stderr[,3]^2) + delta12*delta13*cov23
  QIVW1_3[j] = sum((1/VX1)*(pihat[,1] - (delta12*pihat[,2]+delta13*pihat[,3]))^2)
  
  
  fitx2 = lm(pihat[,2] ~ -1 + pihat[,1] + pihat[,3])
  delta21 = fitx2$coef[1]
  delta23 = fitx2$coef[2]
  VX2 = (stderr[,2]^2)- 2*delta21*(cov12) - 2*delta23*(cov23) + (delta21^2)*(stderr[,1]^2) + (delta23^2)*(stderr[,3]^2) + delta21*delta23*cov13
  QIVW2_3[j] = sum((1/VX2)*(pihat[,2] - (delta21*pihat[,1]+delta23*pihat[,3]))^2)
  
  
  
  fitx3 = lm(pihat[,3] ~ -1 + pihat[,1] + pihat[,2])
  delta31 = fitx3$coef[1]
  delta32 = fitx3$coef[2]
  VX3 = (stderr[,3]^2)- 2*delta31*(cov13) - 2*delta32*(cov23) + (delta31^2)*(stderr[,1]^2) + (delta32^2)*(stderr[,2]^2) + delta31*delta32*cov12
  QIVW3_3[j] = sum((1/VX3)*(pihat[,3] - (delta31*pihat[,1]+delta32*pihat[,2]))^2)
  
  
  #estimation of Q
  betaest = rbind(tsest1_3[j], tsest2_3[j], tsest3_3[j])
  w <- segamma^2
  Qa_novar3[j] =  sum((1/w)*((gammahat - pihat%*%betaest)^2))
  
  
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
  
  Qa3[j] =  sum((1/w)*((gammahat - pihat%*%betaest)^2))
  
  
  
  
  #LIML estimation of the betas.
  
  liml_beta<- optim(c(0,0,0), LIML_b)
  
  limlest1_3[j] <- liml_beta$par[1]
  limlest2_3[j] <- liml_beta$par[2]
  limlest3_3[j] <- liml_beta$par[3]
  
  Qexact_3[j] <- liml_beta$value
  
  #estimate the parameters with LIML - with heterogeneity statistic
  
  limltauest  = optimize(PL_MVMR,interval=c(-10,10))
  tau_i = limltauest$minimum
  
  tau_3[j] = tau_i
  
  liml_het2<- optim(c(0,0,0), PL2_MVMR)
  limlhet1_3[j] <- liml_het2$par[1]
  limlhet2_3[j] <- liml_het2$par[2]
  limlhet3_3[j] <- liml_het2$par[3]
  Qexact_het_3[j] <- liml_het2$value
  
 
}


