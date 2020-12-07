
Fstrong <- function(exposure)
{
  
  F.strong <- data.frame()
  for(j in 1:length(exposure)){
    
    F.analysis <- data.frame(Fstat[,exposure[j]])
    keep <- as.vector(as.numeric(F.analysis> 5))
    F.analysis <- F.analysis*keep
    F.analysis[F.analysis == 0] <- NA
    F.analysis<- na.omit(F.analysis)
    
    F.strong[j,1] <- apply(F.analysis,2, mean)
    F.strong[j,2] <- dim(F.analysis)[1]
    
    
  }
  
  row.names(F.strong) <- exposure
  colnames(F.strong) <- c("F.stat", "no.snps")
  return(F.strong)
  
}



conditionalF <- function(exposures)
{
  
  F.analysis <- data.frame(Fstat[,exposures])
  ex.analysis <- data.frame(exp[,exposures])
  se.analysis <- data.frame(dat_se[,exposures])
  cor.analysis <- data.frame(correlations[exposures,exposures])
  
  maxF_row <- apply(F.analysis,1,function(x) max(as.numeric(x)))
  keep <- as.vector(as.numeric(maxF_row > 5))
  
  n.snp <- sum(keep)
  
  F.analysis <- F.analysis*keep
  F.analysis[F.analysis == 0] <- NA
  F.analysis<- na.omit(F.analysis)
  
  ex.analysis <- ex.analysis[,1:length(exposures)]*keep
  ex.analysis[ex.analysis == 0] <- NA
  ex.analysis<- na.omit(ex.analysis)
  
  se.analysis <- se.analysis[,1:length(exposures)]*keep
  se.analysis[se.analysis == 0] <- NA
  se.analysis<- na.omit(se.analysis)
  
  Q <- data.frame()
  F.conditional <- data.frame()
  F.mean <- data.frame()
  
  for(j in 1:length(exposures)){
    
    F.mean <- apply(F.analysis,2,function(b) mean(as.numeric(b)))
    
    X1 <- ex.analysis[,exposures[j]]
    X2 <- data.frame(ex.analysis[,-which(names(ex.analysis) == exposures[j])])
    X2 <- (matrix(unlist(X2), ncol = length(exposures)-1))
    
    X1se <- se.analysis[,exposures[j]]
    X2se <- se.analysis[,-which(names(ex.analysis) == exposures[j])]
    X2se <- (matrix(unlist(X2se), ncol = length(exposures)-1))
    
    delta <- as.vector(lm(X1~ -1 + X2)$coefficients)
    delta1 <- c(1,delta)
    
    vx <- NULL
    
    for(i in 1:n.snp){
      
      covariances <- data.frame()
      
      for(l in 1:length(exposures)){
        for(m in 1:length(exposures)){
          covariances[l,m] <- correlations[exposures[l],exposures[m]]*se.analysis[i,exposures[l]]*se.analysis[i,exposures[m]]
        } 
      }
      colnames(covariances) <- exposures
      rownames(covariances) <- exposures
      
      covX1a <- data.frame(covariances[exposures[j],])
      colnames(covX1a) <- exposures
      covX1 <- covX1a[,-which(names(covX1a) == exposures[j])]
      
      covX2 <- as.matrix(covariances[-which(names(covariances) == exposures[j]),-which(names(covariances) == exposures[j])])
      
      vx[i] <- X1se[i]^2 - 2*delta%*%t(covX1) + t(delta)%*%covX2%*%delta
      
      
    }
    
    Q[1,j] <- sum((1/vx)*((lm(X1~ -1 + X2)$residuals)^2))
    F.conditional[1,j]<-Q[1,j]/length(X1) 
    no.snps <- length(X1)
    
  }
  F.conditional <- t(F.conditional)
  rownames(F.conditional)<-exposures
  colnames(F.conditional)<-"Conditional F"

  results <- data.frame(F.conditional, F.mean, no.snps)
  return(results)

}




MRfunction_jk <- function(exposures)
{
  
  
  PL_MVMR = function(a){
    tau2   = a[1]
    
    PL2_MVMR = function(ab){
      b=ab
      
      cov = data.frame()
      w=NULL
      for(i in 1:n.snp){
        for(l in 1:length(exposures)){
          for(m in 1:length(exposures)){
            cov[l,m] <- correlations[exposures[l],exposures[m]]*X_se[i,subexp_se[l]]*X_se[i,subexp_se[m]]
          }
        }
        cov <- as.matrix(cov)
        w[i] <- Y1_se[i]^2+t(b)%*%cov%*%b + tau2    
      }
      
      q =  sum((1/w)*((Y1 - X1%*%b)^2))
      
    }
    
    
    st_PL2  =  cbind(rep(1,kx))
    bc    = optim(st_PL2, PL2_MVMR)
    
    bcresults <- bc$par
    
    cov = data.frame()
    w=NULL
    for(i in 1:n.snp){
      for(l in 1:length(exposures)){
        for(m in 1:length(exposures)){
          cov[l,m] <- correlations[exposures[l],exposures[m]]*X_se[i,subexp_se[l]]*X_se[i,subexp_se[m]]
        }
      }
      cov <- as.matrix(cov)
      w[i] <- Y1_se[i]^2+t(bcresults)%*%cov%*%bcresults + tau2     
    }
    
    q = (sum((1/w)*((Y1 - X1%*%bcresults)^2))-(n.snp-kx))^2
  }
  
  
  
  PL2_MVMR = function(ab){
    b=ab
    
    w=NULL
    cov = data.frame()
    for(i in 1:n.snp){
      for(l in 1:length(exposures)){
        for(m in 1:length(exposures)){
          cov[l,m] <- correlations[exposures[l],exposures[m]]*X_se[i,subexp_se[l]]*X_se[i,subexp_se[m]]
        }
      }
      cov <- as.matrix(cov)
      w[i] <- Y1_se[i]^2+t(b)%*%cov%*%b + tau_i 
    }
    
    q =  sum((1/w)*((Y1 - X1%*%b)^2))
    
  }
  
  
  
  
  ex.analysis <- analysis.dat[,c("amd",subexp)]
  se.analysis <- analysis.dat[,c("amd_se",subexp_se)]
  F.analysis <- analysis.dat[,c(subexp_f)]
  cor.analysis <- data.frame(correlations[exposures,exposures])
  
  kx <- length(exposures)
  
  n.snp <- length(analysis.dat[,"amd"])
  
  X <- ex.analysis[,exposures]
  X1 <- (matrix(unlist(X), ncol = length(exposures)))
  X_se <- se.analysis[,subexp_se]
  X1_se <- (matrix(unlist(X_se), ncol = length(exposures)))
  Y <- ex.analysis[,"amd"]
  Y1 <- (matrix(unlist(Y), ncol = 1))
  Y_se <- se.analysis[,"amd_se"]
  Y1_se <- (matrix(unlist(Y_se), ncol = 1))
  
  ##simple IVW
  
  MR.IVW <- summary(lm(Y1 ~ -1 + X1, weights = (Y1_se)^-2))$coefficients
  
  beta <- as.vector(lm(Y1 ~ -1 + X1, weights = (Y1_se)^-2))$coefficients
  
  
  
  for(k in 1:10){
    vx <- NULL
    for(i in 1:n.snp){
      
      covariances <- data.frame()
      
      for(l in 1:length(exposures)){
        for(m in 1:length(exposures)){
          covariances[l,m] <- correlations[exposures[l],exposures[m]]*X_se[i,subexp_se[l]]*X_se[i,subexp_se[m]]
        } 
      }
      colnames(covariances) <- exposures
      rownames(covariances) <- exposures
      
      covar <- as.matrix(covariances)
      
      vx[i] <- Y1_se[i]^2 + t(beta)%*%covar%*%beta
      
      
    }
    
    beta <- as.vector(lm(Y1 ~ -1 + X1, weights = (vx)^-1))$coefficients
    
  }
  
  MR.IVW_update <- summary(lm(Y1 ~ -1 + X1, weights = (vx)^-2))$coefficients
  
  
  #liml estimation with the heterogeneity statistic 
  
  tau_i = 0 
  liml <- optim(rep(0,kx), PL2_MVMR)
  Q_cov <- liml$value
  MR.liml <- liml$par
  MR.liml <- cbind(MR.liml)
  
  tau_i = optimize(PL_MVMR,interval=c(-10,10))$minimum
  MR.liml_results<- optim(rep(0,kx), PL2_MVMR, hessian=TRUE)
  
  MR.liml_het <- cbind(MR.liml_results$par)
  Hess <- MR.liml_results$hessian
  var.cov <- solve(Hess)
  MR.liml_het.se <- as.matrix(sqrt(diag(var.cov)),ncol=1)
  
  rownames(MR.IVW)<-exposures
  rownames(MR.IVW_update)<-exposures
  rownames(MR.liml)<-exposures
  rownames(MR.liml_het)<-exposures
  rownames(MR.liml_het.se)<-exposures
  
  results <- data.frame(n.snp, MR.IVW, MR.IVW_update, MR.liml, Q_cov, MR.liml_het, MR.liml_het.se, tau_i)
  return(results)
  
}

