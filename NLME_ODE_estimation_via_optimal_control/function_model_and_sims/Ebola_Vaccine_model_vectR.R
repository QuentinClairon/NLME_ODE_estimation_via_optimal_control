Ebola_Vaccine_model_vectR<- function(t,param_sub,param_pop,exo_par)
{
  log_delta_L = exo_par[[2]]
  coeff_mult  = exo_par[[3]]
  
  
  delta_S = exp(param_pop[1])
  delta_L = exp(log_delta_L)
  phi_S = exp(param_pop[2] +param_sub[1] )
  phi_L = exp(param_pop[3] +param_sub[2] )
  # delta_Ab = exp(param_pop[4]  +param_sub[3])
  
  S_cell_part =  phi_S*exp(-delta_S*t)
  L_cell_part =  phi_L*exp(-delta_L*t)
  
  res_t = matrix(0,1,1)
  
  res_t[1,1] =  coeff_mult*(S_cell_part+L_cell_part)
  
  #print(res_t)
  return(res_t)
  
}
