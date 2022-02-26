est_param_oca_variance_estimation_lincase_prof_par_version<- function(Times_integ_set,pseudo_Y_set,param_exogen_set,weight_integ_set,param_pop_est,log_Delta_vect_est,param_sub_est,
                                                                      mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,nb_obs_tot,list_x0_known = list(),est_pop_parameter_only=0,delta_known=0,method_est_jacob =0,cl1= list()){
  
  ############################################################
  # DESCRIPTION: compute the variance for the mean population parameter and log-transformed of the subject specific variance for linear ODE based mixed-effect models
  # 
  # FUNCTION SIGNATURE: est_param_oca_variance_estimation_lincase_prof_par_version<- function(Times_integ_set,pseudo_Y_set,param_exogen_set,weight_integ_set,param_pop_est,log_Delta_vect_est,param_sub_est,
  #                                                     mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,nb_obs_tot,list_x0_known = list(),est_pop_parameter_only=0,delta_known=0,method_est_jacob =0,cl1= list())
  #  
  # INPUT :
  # Times_integ_set           :(list of 1 x nb_time_integ_points matrices) list of integration time points for each subject
  # pseudo_Y_set              :(list of dim_obs x nb_time_integ_points sized matrix) list of observations matrix  for each subject
  # param_exogen_set          :(list of vectors) list of vector of covariates or other fixed parameters for each subject
  # weight_integ_set          :(list of 1 x nb_time_integ_points matrices)  values of weights sequence for each subject (w in the paper)
  # param_pop_est             :(vector) estimated value for theta, the mean population parameter value
  # log_Delta_vect_est        :(vector) estimated value for delta, the log transformation of the subject specific variance delta
  # param_sub_est             :(list of vectors) list of estimated subject specific parameters
  # mat_U                     :(dim_control square matrix) weighing matrix for control magnitude penalization
  # mat_A_pop                 :(function) matrix valued function representing the matrix A in the linear ODE dX/dt = A(t,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              mat_A<- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                 OUTPUT: (dim_syst x dim_syst  matrix) the value of  A(t,bi,theta)
  # vect_r_pop                 :(function) matrix valued function representing the Vector r in the linear ODE dX/dt = A(t,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              vect_r <- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                OUTPUT: (dim_syst x 1  matrix) the value of  r(t,bi,theta)
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # nb_obs_tot                :(integer) total number of observations per subject and state variable
  # list_x0_known             :(list of vector) list of known initial conditions for each subject
  # est_pop_parameter_only    :(integer) equal to 1 if subject specific parameters are known or 0 if they are estimated 
  # delta_known               :(integer) equal to 1 if the variance of subject specific parameters is known or 0 if it is estimated
  # method_est_jacob          :(integer) precise which kind of numerical method is required for jacobian/hessian approximation during variance estimation 1 for Richarsdon based method, 0 for simple finite difference
  # cl1                       :(cluster object or empty list) if a cluster object is provided, parallelization of the inner bi estimation and optimal control problem solving
  
  # OUTPUT:
  # Return a list  composed of:
  # Cov_est : (matrix) the estimated variance/covariance matrix 
  # est_A   : (matrix) corresponding A_hat matrix 
  # est_B   : (matrix) corresponding B_hat matrix 
  ############################################################
  #  print(log_Delta_vect_est)
  nb_subject =length(Times_integ_set);
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  nb_param_pop = length(param_pop_est);
  nb_bm = length(log_Delta_vect_est);
  nb_param_cov =nb_bm;
  
  Delta_mat_est = diag(x=exp(log_Delta_vect_est),nrow = nb_bm, ncol = nb_bm)
  inv_Delta = solve( Delta_mat_est)
  der_det_delta = matrix(0,nb_bm,1)
  der_sec_det_delta = matrix(0,nb_bm,nb_bm+nb_param_pop)
  
  if (delta_known == 1){
    param_ext_est = param_pop_est;
    nb_param_tot = length(param_ext_est)
  }else{
    param_ext_est = c(param_pop_est, log_Delta_vect_est);
    nb_param_tot = length(param_ext_est)
    
    
    for (dk in 1:nb_bm){
      ind_delta_k_mat = matrix(0,nb_bm,nb_bm)
      ind_delta_k_mat[dk,dk] = exp(log_Delta_vect_est[[dk]])
      der_det_delta[dk] =  sum(diag(inv_Delta%*%ind_delta_k_mat))
    }
  }
  
  
  
  func_comp_h<- function(seq_arg){
    pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
    sapply(pathnames, FUN=source);
    pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
    sapply(pathnames, FUN=source);
    
    param_tot_cur = seq_arg[[1]]
    Times_integ_cur = seq_arg[[2]]
    pseudo_Y_cur = seq_arg[[3]]
    weight_integ_cur = seq_arg[[4]]
    param_sub_est_cur = seq_arg[[5]]
    mat_U_cur= seq_arg[[6]]
    mat_A_pop_cur= seq_arg[[7]]
    vect_r_pop_cur= seq_arg[[8]] 
    mat_B_cur= seq_arg[[9]]  
    mat_C_cur= seq_arg[[10]]  
    x0_known_cur= seq_arg[[11]] 
    param_exogen_set_cur= seq_arg[[12]]
    est_pop_parameter_only= seq_arg[[13]]
    delta_known= seq_arg[[14]]
    
    nb_bm_cur = length(param_sub_est_cur)
    
    nb_param_tot_cur = length(param_tot_cur)
    dim_obs_cur = nrow(mat_C_cur)
    
    
    if  (delta_known == 1){
      Delta_mat_cur = Delta_mat_est;
      nb_param_pop_cur = nb_param_tot_cur
      param_pop_cur = param_tot_cur
    }
    else{
      nb_param_pop_cur = nb_param_tot_cur -nb_bm_cur
      param_pop_cur = param_tot_cur[seq(1,nb_param_pop_cur,1)];
      log_Delta_vect = param_tot_cur[seq(nb_param_pop_cur+1,nb_param_pop_cur+nb_bm_cur,1)];
      Delta_mat_cur = diag(x= exp(log_Delta_vect),nrow = nb_bm_cur,ncol =nb_bm_cur)
      
    }
    
    
    seq_argument_ind_sub = list(Times_integ_cur,pseudo_Y_cur,weight_integ_cur ,param_sub_est_cur, Delta_mat_cur,
                                mat_U_cur,mat_A_pop_cur,vect_r_pop_cur,mat_B_cur,mat_C_cur,x0_known_cur,2,est_pop_parameter_only,param_pop_cur,param_exogen_set_cur)
    
    out_ind_sub =  est_param_oca_inner_criteria_lincase_par_version(seq_argument_ind_sub)
    val_h =  out_ind_sub$val_hi
    
    return(val_h)
  }
  
  
  
  func_comp_h_grad_hessian <- function(seq_arg_tr){
    library('numDeriv')
    param_to_est = seq_arg_tr[[1]]
    method_est_jacob_cur= seq_arg_tr[[15]]
    
    func_comp_to_est_grad_hessian <- function(par_c,seq_arg_cur){
      new_seq_arg = seq_arg_cur
      new_seq_arg[[1]] <-par_c
      return(func_comp_h(new_seq_arg))
    }
    
    comp_h = func_comp_h(seq_arg_tr)
    if (method_est_jacob_cur>0){
      if (method_est_jacob_cur==3){
        comp_h_grad= grad(function(par_c){func_comp_to_est_grad_hessian(par_c,seq_arg_tr)},param_to_est,method='simple')
        comp_h_hessian= hessian(function(par_c){func_comp_to_est_grad_hessian(par_c,seq_arg_tr)},param_to_est,method='complex')
      }else{
        comp_h_grad= grad(function(par_c){func_comp_to_est_grad_hessian(par_c,seq_arg_tr)},param_to_est,method='Richardson') 
        comp_h_hessian= hessian(function(par_c){func_comp_to_est_grad_hessian(par_c,seq_arg_tr)},param_to_est,method='Richardson')
      }
    }
    else{
      comp_h_grad= grad(function(par_c){func_comp_to_est_grad_hessian(par_c,seq_arg_tr)},param_to_est,method='simple')
      comp_h_hessian= hessian(function(par_c){func_comp_to_est_grad_hessian(par_c,seq_arg_tr)},param_to_est,method='complex')
    }
    
    return(list(comp_h,comp_h_grad,comp_h_hessian))
  }
  
  
  est_A_hat = matrix(rep(0,nb_param_tot^2),nb_param_tot,nb_param_tot)
  est_B_hat = matrix(rep(0,nb_param_tot^2),nb_param_tot,nb_param_tot)
  val_sum_h_div_nb_obs = 0
  val_sum_grad_h_div_nb_obs = matrix(0,nb_param_tot,1)
  
  
  
  list_seq_arg = list()
  list_seq_arg_tronq = list()
  
  
  for (nbs4 in 1: nb_subject){
    
    seq_argument_sens_nbs = list(param_ext_est,Times_integ_set[[nbs4]],pseudo_Y_set[[nbs4]],weight_integ_set[[nbs4]],param_sub_est[[nbs4]], 
                                 mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known[[nbs4]],
                                 param_exogen_set[[nbs4]],est_pop_parameter_only,delta_known,method_est_jacob)
    
    
    
    list_seq_arg = append(list_seq_arg,list(seq_argument_sens_nbs))
    
  }
  
  if (length(cl1)==0){
    out_parLapply_func_comp_h_grad_hessian = lapply(list_seq_arg,func_comp_h_grad_hessian)
  }else{
    out_parLapply_func_comp_h_grad_hessian = parLapply(cl1,list_seq_arg,func_comp_h_grad_hessian)
  }
  
  
  #  for (ns in 1:nb_subject){
  #  print(out_parLapply_func_comp_h_grad_hessian[[ns]])
  #    h_ns = c(out_parLapply_func_comp_h_grad_hessian[[ns]][[1]])
  #    grad_h_ns = out_parLapply_func_comp_h_grad_hessian[[ns]][[2]]
  
  #   val_sum_h_div_nb_obs =  val_sum_h_div_nb_obs+h_ns
  #    val_sum_grad_h_div_nb_obs = val_sum_grad_h_div_nb_obs+ matrix(grad_h_ns,nb_param_tot,1)
  #  }
  #  val_sum_h_div_nb_obs=  val_sum_h_div_nb_obs/(nb_obs_tot+ nb_bm*nb_subject)
  #  val_sum_grad_h_div_nb_obs=  val_sum_grad_h_div_nb_obs/(nb_obs_tot+ nb_bm*nb_subject)
  
  for (ns in 1:nb_subject){
    h_ns = out_parLapply_func_comp_h_grad_hessian[[ns]][[1]]
    grad_h_ns = matrix(out_parLapply_func_comp_h_grad_hessian[[ns]][[2]],nb_param_tot,1)
    hessian_h_ns = out_parLapply_func_comp_h_grad_hessian[[ns]][[3]]
    
    
    if (delta_known == 1){
      val_J_ns =   grad_h_ns
      der_param_val_J_ns = hessian_h_ns
    }else{
      val_J_ns =   grad_h_ns
      #  print(grad_h_ns)
      val_J_ns[(nb_param_pop+1):(nb_param_pop+nb_param_cov)] =  
        val_J_ns[(nb_param_pop+1):(nb_param_pop+nb_param_cov)]-2*(nb_subject/(nb_obs_tot+ nb_bm*nb_subject))*der_det_delta*h_ns
      
      der_param_val_J_ns = hessian_h_ns
      # print( hessian_h_ns)
      der_param_val_J_ns[(nb_param_pop+1):(nb_param_pop+nb_param_cov),] =  der_param_val_J_ns[(nb_param_pop+1):(nb_param_pop+nb_param_cov),]
      -2*(nb_subject/(nb_obs_tot+ nb_bm*nb_subject))*(der_det_delta%*%t( grad_h_ns) +der_sec_det_delta* h_ns)
    }
    
    est_A_hat = est_A_hat+ der_param_val_J_ns
    est_B_hat = est_B_hat+ val_J_ns%*%t(val_J_ns)
  }
  
  est_A_hat = -est_A_hat/nb_subject
  est_B_hat = est_B_hat/nb_subject
  
  Cov_est = solve(est_A_hat)%*%est_B_hat%*%t(solve(est_A_hat))/nb_subject
  return(list(est_var_cov_Matrix = Cov_est,est_A = est_A_hat,est_B = est_B_hat))
}
