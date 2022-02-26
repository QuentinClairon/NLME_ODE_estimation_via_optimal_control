erreur_prediction_estimation_oca <- function(Times_integ_set,pseudo_Y_set,pseudo_State_ini,weight_integ_set,param_sub_ini_set,param_exogen_set, Delta_mat,param_pop_est, mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known = list()){
  ############################################################
  # DESCRIPTION: compute the summed prediction error corresponding to the current population and subject specific parameter estimations for ODE based mixed-effect models
  # 
  # FUNCTION SIGNATURE:  erreur_prediction_estimation_oca <- function(Times_integ_set,pseudo_Y_set,pseudo_State_ini,weight_integ_set,param_sub_ini_set,param_exogen_set, Delta_mat,param_pop_est, mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known = list())
  # 
  # INPUT :
  # Times_integ_set           :(list of 1 x nb_time_integ_points matrices) list of integration time points for each subject
  # pseudo_Y_set              :(list of dim_obs x nb_time_integ_points sized matrix) list of observations matrix  for each subject
  # pseudo_State_ini          :(list of dim_syst x nb_time_integ_points sized matrix) list of initial state value at integration time points for the optimal control  for each subject
  # weight_integ_set          :(list of 1 x nb_time_integ_points matrices)  values of weights sequence for each subject (w in the paper)
  # param_sub_ini_set         :(list of vectors) list of estimated values for the subject specific parameters
  # param_exogen_set          :(list of vectors) list of vector of covariates or other fixed parameters for each subject
  # Delta_mat                 :(square matrix)   estimator of Delta matrix
  # param_pop_est             :(vector) estimator of theta
  # mat_U                     :(dim_control square matrix) weighing matrix for control magnitude penalization
  # mat_A_pop                 :(function) matrix valued function representing the matrix A in the pseudolinear formulation dX/dt = A(t,x,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              mat_A<- function(t,vx,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : vx (vector) current state value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  # vect_r_pop                 :(function) matrix valued function representing the Vector r in the  pseudolinear formulation dX/dt = A(t,x,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              vect_r <- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                OUTPUT: (dim_syst x 1  matrix) the value of  r(t,bi,theta)
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # list_x0_known             :(list of vector) list of known initial conditions for each subject
  
  # OUTPUT:
  # error_pred:               :(positive real) the prediction error estimator
  ############################################################  
  nb_subject =length(Times_integ_set)
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  
  v_field <- function(Time, State, Pars,nbs_cur){
    dState = mat_A_pop(Time,State,param_sub_ini_set[[nbs_cur]],param_pop_est,param_exogen_set[[nbs_cur]])%*%State+ vect_r_pop(Time,param_sub_ini_set[[nbs_cur]],param_pop_est,param_exogen_set[[nbs_cur]])
    
    return(list(c(dState)))
  }
  
  error_pred = 0
  for (nbs  in 1:nb_subject){
    Times_integ_set_nbs = Times_integ_set[[nbs]]
    pseudo_Y_set_nbs= pseudo_Y_set[[nbs]]
    pseudo_State_ini_nbs = pseudo_State_ini[[nbs]]
    weight_integ_set_nbs= weight_integ_set[[nbs]]
    param_sub_ini_set_nbs = param_sub_ini_set[[nbs]]
    param_exogen_set_nbs = param_exogen_set[[nbs]]
    
    seq_argument_nbs = list(Times_integ_set_nbs,pseudo_Y_set_nbs,pseudo_State_ini_nbs,weight_integ_set_nbs,param_sub_ini_set_nbs, Delta_mat,
                            mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known[[nbs]],2,est_pop_parameter_only=1,param_pop_est,param_exogen_set_nbs)
    
    Times_obs_nbs = Times_integ_set_nbs[weight_integ_set_nbs>0]
    Y_obs_nbs = subset(t(pseudo_Y_set_nbs),weight_integ_set_nbs>0)
    Y_obs_nbs = t(Y_obs_nbs)
    
    out_inner_control_problem =est_param_oca_inner_criteria_par_version(seq_argument_nbs)
    State_opt_nbs = out_inner_control_problem$Xi_opt
    State_obs_nbs = subset(t(State_opt_nbs),weight_integ_set_nbs>0)
    State_obs_nbs = t(State_obs_nbs)
    
    v_field_nbs   =function(Time, State, Pars,nbs_cur){v_field(Time, State, Pars,nbs)}
    
    middle_interval = ceiling(length(Times_obs_nbs)/2)
    
    Times_obs_nbs_first_part = Times_obs_nbs[1:(middle_interval-1)]
    Times_obs_nbs_second_part = Times_obs_nbs[middle_interval:length(Times_obs_nbs)]
    
    Y_obs_nbs_first_part = Y_obs_nbs[,1:(middle_interval-1)]
    Y_obs_nbs_second_part = Y_obs_nbs[,middle_interval:length(Times_obs_nbs)]
    
    X0_opt =  State_obs_nbs[,1]
    Xmiddle_opt = State_obs_nbs[,middle_interval]
    
    
    pred_first_part <- ode(func =v_field_nbs, y = X0_opt, parms =c(), times = Times_obs_nbs_first_part)
    pred_second_part <- ode(func =v_field_nbs, y = Xmiddle_opt, parms =c(), times = Times_obs_nbs_second_part)
    
    error_pred = error_pred+sum((Y_obs_nbs_first_part - mat_C%*%t(pred_first_part[,2:(dim_syst+1)]))^2)
    error_pred = error_pred+sum((Y_obs_nbs_second_part - mat_C%*%t(pred_second_part[,2:(dim_syst+1)]))^2)
  }
  
  
  
  return(error_pred)
}
