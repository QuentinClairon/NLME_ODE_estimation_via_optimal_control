Ebola_Vaccine_model_matA_inhomogen<- function(t,param_sub,param_pop,exo_par)
{
 # coeff_mult  = exo_par[[3]]
  delta_Ab = exp(param_pop[3]  +param_sub[3])
  
  res_t = matrix(0,1,1)
  res_t[1,1] =  -delta_Ab
  return(res_t)
  
}
