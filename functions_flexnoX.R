
#estimation of the weak instrument F through minimisation of Q to estimate delta. 

Q_ind = function(a){
  

  d1 = a[1]
  
  w = 1/(stderr[,1]^2- 2*d1*cov + (d1^2)*stderr[,2]^2)
  q =  sum(w*(pihat[,1] - (d1*pihat[,2]))^2)
  
}



Q_ind2 = function(a){
  
  
  d2 = a[1]
  
  w = 1/(stderr[,2]^2- 2*d2*cov + (d2^2)*stderr[,1]^2)
  q =  sum(w*(pihat[,2] - (d2*pihat[,1]))^2)
  
}





#Estimation of beta through minimisation of the Q statistic - no heterogeneity
LIML_b = function(ab){
  
  b=ab
  
  cov = matrix(nrow = kx, ncol = kx)
w=NULL
  for(l in 1:L){
    for(pp in 1:kx){
      for(p2 in 1:kx){
        cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
      }
    }
    
  w[l] <- segamma[l]^2+t(b)%*%cov%*%b     
  }

  q =  sum((1/w)*((gammahat - pihat%*%b)^2))
  
}





#estimation with the heterogeneity statistic



PL_MVMR = function(a){
  tau2   = a[1]
  
  PL2_MVMR = function(ab){
    b=ab
    
    cov = matrix(nrow = kx, ncol = kx)
    w=NULL
    for(l in 1:L){
      for(pp in 1:kx){
        for(p2 in 1:kx){
          cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
        }
      }
      
      w[l] <- segamma[l]^2+t(b)%*%cov%*%b + tau2    
    }
    
    q =  sum((1/w)*((gammahat - pihat%*%b)^2))
    
  }
  
  
  st_PL2  =  rep(0,kx)
  bc    = optim(st_PL2, PL2_MVMR)
  
  bcresults <- bc$par
  
  cov = matrix(nrow = kx, ncol = kx)
  w=NULL
  for(l in 1:L){
    for(pp in 1:kx){
      for(p2 in 1:kx){
        cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
      }
    }
    
    w[l] <- segamma[l]^2+t(bcresults)%*%cov%*%bcresults + tau2    
  }
  
  q = (sum((1/w)*((gammahat - pihat%*%bcresults)^2))-(L-kx))^2
}



PL2_MVMR = function(ab){
  b=ab
  
  w=NULL
  cov = matrix(nrow = kx, ncol = kx)
  for(l in 1:L){
    for(pp in 1:kx){
      for(p2 in 1:kx){
        cov[pp,p2] <- correlation[pp,p2]*stderr[l,pp]*stderr[l,p2]
      }
    }
    
    w[l] <- segamma[l]^2+t(b)%*%cov%*%b + tau_i    
  }
  
  q =  sum((1/w)*((gammahat - pihat%*%b)^2))
  
}


