log_Prior_Chloe_model<- function(param_pop,Delta_mat)
{
  log_denom = log(((2*pi)^3/2)*10^-2)
  Inv_Sigma_Var = diag(c(0.01,0.01,1))
  res_log=-0.5*t(param_pop[1:3]- c(0,0,-4.1))%*%Inv_Sigma_Var%*%(param_pop[1:3]- c(0,0,-4.1))-log_denom 
  return(res_log)
  
}
