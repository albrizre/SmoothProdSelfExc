standard_self <- nimbleCode({
  
  lambda0 ~ dgamma(10*lambda0_base,10)
  theta ~ dunif(0,1)
  phi_S ~ dunif(30,10000)
  phi_T ~ dunif(1,100)

  for (i in 1:1) {
    
    lambda_event[i] <- lambda0
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in 2:N_events) {
    
    lambda_event[i] <- lambda0+theta*inprod(exp(-matrix_spatial_distances_events[i,]/phi_S),exp(-matrix_distances[i,]/phi_T))
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N_0) { 
    
    lambda_int[1:S,i-N_events] <- rep(lambda0,S)
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1 # 1 is the length of the temporal intervals 
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_0+1):N) { 
    
    lambda_int[1:S,i-N_events] <- lambda0+theta*t(exp(-matrix_spatial_distances_events_grid[1:n_previous[i],]/phi_S))%*%t(t(exp(-matrix_distances_events_time_points[1:n_previous[i],time_all[i]+1]/phi_T)))
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1 # 1 is the length of the temporal intervals 
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})