smoothprod_self2 <- nimbleCode({
  
  lambda0 ~ dgamma(10*lambda0_base,10)
  theta ~ dgamma(mean = mean_theta,sd = sd_theta)
  phi_S ~ dgamma(mean = mean_phi_S,sd = sd_phi_S)
  phi_T ~ dgamma(mean = mean_phi_T,sd = sd_phi_T)
  
  for (i in 1:nbasis){
    alpha[i] ~ dnorm(0,sd = 1)
  }
  
  for (i in 1:nbasisT){
    alphaT[i] ~ dnorm(0,sd = 1)
  }
  
  for (s in 1:S){
    logit(b_S[s]) <- inprod(alpha[1:nbasis],X_splines[s,1:nbasis])
  }
  
  for (i in 1:1) {
    
    b_T_i[i] <- b_T[time_all[i] + 1]
    lambda_event[i] <- lambda0
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in 2:N_events) {
    
    b_T_i[i] <- b_T[time_all[i] + 1]
    lambda_event[i] <- lambda0+b_T_i[i]*b_S[IdS_events[i]]*theta*inprod(exp(-matrix_spatial_distances_events[i,]/phi_S),exp(-matrix_distances[i,]/phi_T))
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N_0) { 
    
    b_T_i[i] <- b_T[time_all[i] + 1]
    lambda_int[1:S,i-N_events] <- rep(lambda0,S)
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1 # 1 is the length of the temporal intervals 
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_0+1):N_1) { 
    
    b_T_i[i] <- b_T[time_all[i] + 1]
    aux[1:S,i-N_events] <- theta*(exp(-matrix_spatial_distances_events_grid[1,]/phi_S))*exp(-matrix_distances[2,1]/phi_T)
    lambda_int[1:S,i-N_events] <- lambda0+b_T_i[i]*b_S[1:S]*aux[1:S,i-N_events]
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1 # 1 is the length of the temporal intervals 
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_1+1):N) { 
    
    b_T_i[i] <- b_T[time_all[i] + 1]
    aux[1:S,i-N_events] <- theta*t(exp(-matrix_spatial_distances_events_grid[1:n_previous[i],]/phi_S))%*%t(t(exp(-matrix_distances_events_time_points[1:n_previous[i],time_all[i]+1]/phi_T)))
    lambda_int[1:S,i-N_events] <- lambda0+b_T_i[i]*b_S[1:S]*aux[1:S,i-N_events]
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1 # 1 is the length of the temporal intervals 
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in 1:(Tmax+1)){
    logit(b_T[i]) <- inprod(alphaT[1:nbasisT],X_splinesT[i,1:nbasisT])
  }
  
})