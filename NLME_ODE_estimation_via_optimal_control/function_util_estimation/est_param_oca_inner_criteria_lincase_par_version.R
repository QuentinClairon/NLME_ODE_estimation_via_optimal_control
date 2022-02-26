est_param_oca_inner_criteria_lincase_par_version <- function(seq_argument){
  ############################################################
  # DESCRIPTION: estimation of the subject specific parameters via the inner criteria optimization as well as the solution of the related optimal control problem  for linear ODE based mixed-effect models
  # 
  # FUNCTION SIGNATURE:  est_param_oca_inner_criteria_lincase_par_version <- function(seq_argument) where seq_argument is a list composed of
  # 
  # INPUT :
  # Times_integ_set           :(1 x nb_time_integ_points matrix) integration time points 
  # pseudo_Y_set              :(dim_obs x nb_time_integ_points sized matrix) extended observations matrix 
  # weight_integ_set          :(1 x nb_time_integ_points matrix)  values of weights 
  # bi_ini                    :(vector) initial guess for the subject specific parameter estimator
  # Delta_mat                 :(square matrix ) current value for Delta
  # mat_U                     :(dim_control square matrix) weighing matrix for control magnitude penalization
  # mat_A_pop                 :(function) matrix valued function representing the matrix A in the linear ODE dX/dt = A(t,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              mat_A<- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                 OUTPUT: (dim_syst x dim_syst  matrix) the value of  A(t,bi,theta)
  # vect_r_pop                 :(function) vector valued function representing the Vector r in the linear ODE dX/dt = A(t,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              vect_r <-<- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                OUTPUT: (dim_syst x 1  matrix) the value of  r(t,bi,theta)
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # x0_known                  :(vector) subject specific known initial condition
  # inner_criteria_type       :(integer) precize the definition of the inner criteria being use (= 2 for gi)
  # opt_traj_est_only         :(integer) equal 0 if the subject specific parameter have to be estimated, 1 if it is just the optimal control and trajectory
  # theta_pop_cur             :(vector) current value for the mean population parameter
  # exo_pop_cur               :(list) list of fixed parameter and covariates
  # type_algorithm_optim      :(integer) type of optimization algorithm (=0 for Nelder-Mead (default), 1 for Nlminb)
  # inner_information crit    :(2 dimensional vector) supplementary details about optimization algorithm
  #                             - First entry (integer) number of iteration (default 200)
  #                             - Second entry(integer) print algorithm information (= 1 for print, 0 no print (default))
  # OUTPUT:
  # Return a list  composed of elements:
  # bi_est                   : (vector) estimated subject specific parameters 
  # ui_opt                   : (dim_control x (nb_time_integ_points-1) matrix) estimated optimal control
  # Xi_opt                   : (dim_syst x nb_time_integ_points matrix) estimated optimal trajectory 
  # Xi_opt_obs               : (dim_syst x nb_obs matrix) estimated optimal trajectory at observation time points
  # val_gi                   : (real) estimated inner criteria value 
  # val_hi                   : (real) estimated  h value (useful for asymptotic variance estimation)
  ############################################################  
  pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
  sapply(pathnames, FUN=source);
  pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
  sapply(pathnames, FUN=source);
  library('optimx')
  
  Times_integ= seq_argument[[1]]
  pseudo_Y = seq_argument[[2]]
  weight_integ_snd =  seq_argument[[3]]
  param_sub_ini = seq_argument[[4]]
  Delta_mat  = seq_argument[[5]]
  mat_U = seq_argument[[6]]
  mat_A = seq_argument[[7]]
  vect_r = seq_argument[[8]]
  mat_B = seq_argument[[9]]
  mat_C = seq_argument[[10]]
  x0_known =  seq_argument[[11]]
  inner_criteria_type  = seq_argument[[12]]
  opt_traj_est_only = seq_argument[[13]]
  theta_pop_cur = seq_argument[[14]]
  exo_pop_cur = seq_argument[[15]]
  
  type_algorithm_optim = 0
  if (length(seq_argument) > 15){
    type_algorithm_optim = seq_argument[[16]]
  }
 
  
  nb_iter_max_nlminb =200
  trace_val = 0
  if (length(seq_argument) > 16){
    nb_iter_max_nlminb = seq_argument[[17]][1]
    trace_val =seq_argument[[17]][2]
  }

  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  nb_time_integ = length(Times_integ)
  nb_unknown = dim_syst - length(x0_known)

  Yn_fin = pseudo_Y[,nb_time_integ]
 
  func_eco_E <-function(param_sub_cur){
 
    mat_Rn  = t(mat_C)%*%mat_C
    vect_Hn = -t(mat_C)%*%Yn_fin 
 
    History_R = list(mat_Rn) 
    History_H = list(vect_Hn) 
    History_G = list()
    mat_Rk  = mat_Rn
    vect_Hk = vect_Hn
    
    val_Sn_prof = t(Yn_fin)%*%Yn_fin
    
    for (nk in 1:(nb_time_integ-1)){
      
      delta_k = Times_integ[nb_time_integ-nk+1] - Times_integ[nb_time_integ-nk] 
      Yk= pseudo_Y[,nb_time_integ-nk]
      
      mat_Ak = mat_A(Times_integ[nb_time_integ-nk],param_sub_cur,theta_pop_cur,exo_pop_cur)
      mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
      
      vect_rk = vect_r(Times_integ[nb_time_integ-nk],param_sub_cur,theta_pop_cur,exo_pop_cur)
      
      weight_cur = weight_integ_snd[nb_time_integ-nk]
      
      #print(mat_U+delta_k*t(mat_B)%*%mat_Rk%*%mat_B)
      mat_G = solve(mat_U+delta_k*t(mat_B)%*%mat_Rk%*%mat_B)
      History_G = append(list(mat_G),History_G)
     
      Sn_prof_first_term = weight_cur*t(Yk)%*%Yk +(2*t(vect_Hk)+delta_k*t(vect_rk)%*%mat_Rk)%*%vect_rk
      Sn_prof_second_term =  t(vect_Hk+delta_k*mat_Rk%*%vect_rk)%*%mat_B%*%mat_G%*%t(mat_B)%*%(vect_Hk+delta_k*mat_Rk%*%vect_rk)
      val_Sn_prof = val_Sn_prof +delta_k*(Sn_prof_first_term-Sn_prof_second_term)    
      
      mat_Rk_m1_1part = mat_Rk + delta_k*weight_cur*t(mat_C)%*%mat_C + delta_k*(mat_Rk%*%mat_Ak+t(mat_Ak)%*%mat_Rk)+(delta_k^2)*t(mat_Ak)%*%mat_Rk%*%mat_Ak
      mat_Rk_m1_2part = -delta_k*t(mat_Ak_delta)%*%mat_Rk%*%mat_B%*%mat_G%*%t(mat_B)%*%mat_Rk%*%mat_Ak_delta
      mat_Rk_m1 = mat_Rk_m1_1part+mat_Rk_m1_2part
      
      vect_Hk_m1 =  vect_Hk - delta_k*weight_cur*t(mat_C)%*%Yk+delta_k*t(mat_Ak)%*%vect_Hk +delta_k*t(mat_Ak_delta)%*%mat_Rk%*%vect_rk
      vect_Hk_m1 =  vect_Hk_m1 -  delta_k*t(mat_Ak_delta)%*%mat_Rk%*%mat_B%*%mat_G%*%t(mat_B)%*%(vect_Hk +delta_k*mat_Rk%*%vect_rk)
      
      mat_Rk = mat_Rk_m1
      vect_Hk = vect_Hk_m1
      History_R = append(list(mat_Rk),History_R)
      History_H = append(list(vect_Hk),History_H)
    }
    
    
    mat_R0  =  History_R[[1]]
    vect_H0 = History_H[[1]]
    
    if (length(x0_known)==0){
      inv_mat_R0 = solve(mat_R0)
      est_x0 = -inv_mat_R0%*%vect_H0
      val_Sn_prof  = val_Sn_prof -t(vect_H0)%*%inv_mat_R0%*%vect_H0
    }else{
      
      if (length(x0_known)== dim_syst){
        est_x0 = x0_known
        val_Sn_prof  = val_Sn_prof+ t(x0_known)%*%mat_R0%*%x0_known+2*t(vect_H0)%*%x0_known
      }else{
      mat_R0_11 = mat_R0[seq(1,nb_unknown),seq(1,nb_unknown)]
      mat_R0_12 = mat_R0[seq(1,nb_unknown),seq(nb_unknown+1,dim_syst)]
      mat_R0_22 = mat_R0[seq(nb_unknown+1,dim_syst),seq(nb_unknown+1,dim_syst)]
      vect_H0_1 = vect_H0[seq(1,nb_unknown)]
      vect_H0_2 = vect_H0[seq(nb_unknown+1,dim_syst)]
      
      inv_mat_R0_11 = solve(mat_R0_11)
      val_inter = (mat_R0_12%*%x0_known+vect_H0_1)
      
      est_x0_unknown = -inv_mat_R0_11%*%val_inter
      est_x0 = c(t(est_x0_unknown),t(x0_known))
      
      val_Sn_prof  = val_Sn_prof -t(val_inter)%*%inv_mat_R0_11%*%val_inter+ t(x0_known)%*%mat_R0_22%*%x0_known+2*t(vect_H0_2)%*%x0_known
      }
    }
    
    
    
    State_i = matrix(data = rep(0,dim_syst*nb_time_integ ), nrow = dim_syst, ncol = nb_time_integ )
    State_i[,1] =  est_x0
    
    control_i =matrix(data = rep(0,dim_control*(nb_time_integ-1)), nrow = dim_control, ncol = nb_time_integ-1 )
    
    val_RSS = 0
    for (nkk in 1:(nb_time_integ-1)){
      delta_k = Times_integ[nkk+1] - Times_integ[nkk] 
      Ykk= pseudo_Y[,nkk];
      mat_Ak = mat_A(Times_integ[nkk],param_sub_cur,theta_pop_cur,exo_pop_cur)
      mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
      
      vect_rk = vect_r(Times_integ[nkk],param_sub_cur,theta_pop_cur,exo_pop_cur)
      
      mat_Rk_p1 = History_R[[nkk+1]]
      vect_Hk_p1 = History_H[[nkk+1]]
      mat_G = History_G[[nkk]]
      
      control_ukk = -mat_G%*%t(mat_B)%*%(mat_Rk_p1%*%(mat_Ak_delta%*%State_i[,nkk]+delta_k*vect_rk)+vect_Hk_p1)
      control_i[,nkk] = control_ukk
      State_i[,nkk+1]  = mat_Ak_delta%*%State_i[,nkk]+delta_k*vect_rk+delta_k*mat_B%*%control_ukk
      
      if (weight_integ_snd[nkk] > 0){
        val_RSS = val_RSS+ sum((mat_C%*%State_i[,nkk] - Ykk)^2);
      }
      
      
    }
    
    val_RSS  = val_RSS+ sum((mat_C%*%State_i[,nb_time_integ] - Yn_fin)^2);
    indice_obs  =weight_integ_snd>0      
    State_i_obs_time = State_i[,indice_obs]
    
    out_Riccati_inner = list(Val_cost_prof = val_Sn_prof ,control_opt =control_i, Ricatti = History_R,Adjoint =  History_H,Sequence_G = History_G,State_opt =State_i,State_obs = State_i_obs_time,val_RSS=val_RSS)
    
   
    val_h_cur = out_Riccati_inner$val_RSS + sum((Delta_mat%*%param_sub_cur)^2)
    
    if (inner_criteria_type == 1){
      fval = val_h_cur 
    }
    else{
      val_Sn_prof=out_Riccati_inner$Val_cost_prof
      fval = val_Sn_prof + sum((Delta_mat%*%param_sub_cur)^2)
    }
    out_eco_E = c(out_Riccati_inner,val_hi =val_h_cur,fval= fval)
    
    return(out_eco_E )                          
  }
  
  func_eco <- function (param_c){
    fval_E = func_eco_E(param_c)
    return(fval_E$fval) 
  }
  
  if (opt_traj_est_only == 1){
    out_ensemble_param = func_eco_E(param_sub_ini)
    bi_est = param_sub_ini
    val_gi = out_ensemble_param$fval
  }else {
    
    if (type_algorithm_optim==0){
      res_estim = optimr(param_sub_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control =  list(trace=trace_val,maxit  =nb_iter_max_nlminb))
      bi_est =res_estim$par
      val_gi = res_estim$value
    }else{
      res_estim = optimr(param_sub_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control =  list(trace=trace_val,maxit  =nb_iter_max_nlminb),method="nlminb")
      bi_est =res_estim$par 
      val_gi = res_estim$value
    }
 
    if (length( bi_est)==1){
      bi_est = matrix( bi_est,1,1)
    }
   
    out_ensemble_param = func_eco_E( bi_est)   
  }
 
  val_hi =  out_ensemble_param$val_hi 
  ui_opt = out_ensemble_param$control_opt
  Xi_opt = out_ensemble_param$State_opt
  Xi_opt_obs = out_ensemble_param$State_obs
  out_func_inner = list(bi_est=bi_est,ui_opt=ui_opt,Xi_opt=Xi_opt,Xi_opt_obs=Xi_opt_obs,val_gi=val_gi, val_hi= val_hi)
  
  return(out_func_inner)
}
