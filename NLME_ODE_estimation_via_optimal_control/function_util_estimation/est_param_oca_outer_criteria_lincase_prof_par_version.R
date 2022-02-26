est_param_oca_outer_criteria_lincase_prof_par_version <- function(Subject_data,param_pop_ini,log_Delta_vect_ini,mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C
                                                              ,mesh_iter,set_theta_sub,list_x0_known = list(),est_pop_parameter_only=1,delta_known=0,log_prior=function(param_pop,Delta_mat){return(0)}
                                                              ,type_optim =0,est_var =1,type_optim_inner=0,nb_iter_max = c(1000,200,0),cl=list()){
  
  ############################################################
  # DESCRIPTION: estimation of the mean population parameter and log-transformed of the subject specific variance via outer criteria optimization for linear ODE based mixed-effect models
  # 
  # FUNCTION SIGNATURE:  est_param_oca_outer_criteria_lincase_prof_par_version <- function(Subject_data,param_pop_ini,log_Delta_vect_ini,mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C
  #                                                                               ,mesh_iter,set_theta_sub,list_x0_known = list(),est_pop_parameter_only=1,delta_known=0,log_prior=function(param_pop,Delta_mat){return(0)}
  #                                                                             ,type_optim =0,est_var =1,type_optim_inner=0,nb_iter_max = c(1000,200,0),cl=list())
  # 
  # INPUT :
  # Subject_data              :(list) list of available data composed for each subject of a list with the following three elements: 
  #                                                        - (1 x nb_obs sized matrix) matrix containing the observation times
  #                                                        - (dim_obs x  nb_obs sized matrix) matrix containing the observations
  #                                                        - (list) a list of fixed parameters and covariates values
  # param_pop_ini             :(vector) initial guess for the estimator of theta
  # log_Delta_vect_ini        :(vector) initial guess for the estimator of log transformed of the subject specific variance 
  # mat_U                     :(dim_control square matrix) weighing matrix for control magnitude penalization
  # mat_A_pop                 :(function) matrix valued function representing the matrix A in the linear ODE dX/dt = A(t,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              mat_A<- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                 OUTPUT: (dim_syst x dim_syst  matrix) the value of  A(t,bi,theta)
  # vect_r_pop                :(function) vector valued function representing the vector r in the linear ODE dX/dt = A(t,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              vect_r <-<- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                OUTPUT: (dim_syst x 1  matrix) the value of  r(t,bi,theta)
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # mesh_iter                 :(integer) refinment of the mesh size, number of added dicretization points between two observation points  
  # set_theta_sub             :(list of vectors) list of initial guesses for the subject specific parameters
  # list_x0_known             :(list of vectors) list of known initial conditions for each subject
  # est_pop_parameter_only    :(integer) equal to 1 if subject specific parameters are known or 0 if they are estimated 
  # delta_known               :(integer) equal to 1 if the variance of subject specific parameters is known or 0 if it is estimated
  # log_prior                 :(function) specification of the log prior  distribution of signature log<- function(theta,Delta_mat)
  # type_optim                :(integer) specify the optimization algorithm for the outer criteria, if 0 Nelder-Mead or nlminb otherwise
  # est_var                   :(integer) specify if the variance is estimated, 1 if so, 0 otherwise
  # type_optim_inner          :(integer) type of optimization algorithm (=0 for Nelder-Mead (default), 1 for Nlminb)
  # nb_iter_max               :(3 dimensional vector) supplementary details about optimization algorithms
  #                             - First entry (integer) number of iteration for outer criteria (default 1000)
  #                             - Second entry (integer) number of iteration for inner criteria (default 200)
  #                             - Third  entry (integer) print algorithm information for outer criteria (= 1 for print, 0 no print (default))
  # cl                        :(cluster object or empty list) if a cluster object is provided, parallelization of the inner bi estimation and optimal control problem solving
  
  # OUTPUT:
  # Return a list  composed of elements:
  # theta_est                : (vector) estimated value for theta, the mean population parameter value
  # delta_est                : (vector) estimated value for delta, the log transformation of the subject specific variance delta
  # sigma_square_est         : (positive real) estimated variance of measurement error distribution
  # Times_integ_set          : (list) list of integration times for optimal controls and trajectories estimation
  # list_bi_est              : (list) list of estimated subject specific parameters
  # list_ui_opt              : (list) list of estimated optimal controls
  # list_Xi_opt              : (list) list of estimated optimal trajectories 
  # Variance_component       : (list) Variance-covariance estimator for the population parameter (if estimated, otherwise NULL)
  # erreur_cross_validation  : (positive real) error prediction for the whole population
  ############################################################  
  nb_subject =length(Subject_data);
  
  # close_cl_at_the_end = 0
  # if (length(cl)==0){
  #   cl <- makeSOCKcluster(rep("localhost",nb_subject))
  #   close_cl_at_the_end = 1
  # }
  
  # Sortie (State_i, control_i,Riccati_i,val_Sn_prof, State_i_obs_time]
  
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  
  
  nb_param_pop = length(param_pop_ini);
  nb_bm = length(log_Delta_vect_ini);
  nb_param_cov =nb_bm;
  Delta_mat_ini = diag(x=exp(log_Delta_vect_ini),nrow = nb_bm, ncol = nb_bm)
  
  
  if (delta_known == 1){
    param_ext_ini = param_pop_ini;
  }else{
    param_ext_ini = c(param_pop_ini, log_Delta_vect_ini);
  }
  
  Times_integ_set = list()
  weight_integ_set = list()
  pseudo_Y_set  = list()
  param_exogen_set = list()
  nb_obs_tot = 0
   
  for (nbs in 1:nb_subject){
    
    Times_obs_s =Subject_data[[nbs]][[1]]
    Yn_s =Subject_data[[nbs]][[2]]
    
    nb_obs_tot = nb_obs_tot+length(Times_obs_s)
    nb_time_integ_s = mesh_iter*(length(Times_obs_s)-1)+1;
    
    Times_integ_s =  matrix(data = rep(0,nb_time_integ_s ), nrow = 1, ncol = nb_time_integ_s )
    weight_integ_snd_s = matrix(data = rep(0,nb_time_integ_s ), nrow = 1, ncol = nb_time_integ_s )
    
    pseudo_Y_s =  matrix(data = rep(0,dim_obs*nb_time_integ_s ), nrow = dim_obs, ncol = nb_time_integ_s )
    
    for (obs_i in 1:(length(Times_obs_s)-1)){
      
      ecart_i_s = Times_obs_s[obs_i+1] - Times_obs_s[obs_i] 
      
      refined_interval_i_s =  seq(Times_obs_s[obs_i],Times_obs_s[obs_i+1],by = (ecart_i_s/mesh_iter))
      refined_interval_i_s[length( refined_interval_i_s)] = Times_obs_s[obs_i+1]
      Times_integ_s[seq(mesh_iter*(obs_i-1)+1,(mesh_iter*obs_i+1),by =1)] = refined_interval_i_s
      
      weight_integ_snd_s[mesh_iter*(obs_i-1)+1] =  1/(refined_interval_i_s[2]-Times_obs_s[obs_i])
      
      
      pseudo_Y_s[,mesh_iter*(obs_i-1)+1]= Yn_s[,obs_i]
      pseudo_Y_s[,mesh_iter*obs_i+1] = Yn_s[,obs_i+1]
      
    }
    Times_integ_s = c(Times_integ_s)
    weight_integ_snd_s = c(weight_integ_snd_s)
    
    Times_integ_set = append(Times_integ_set,list(Times_integ_s))
    weight_integ_set= append(weight_integ_set,list( weight_integ_snd_s))
    pseudo_Y_set = append(pseudo_Y_set,list(pseudo_Y_s))
    param_exogen_set = append(param_exogen_set,list(Subject_data[[nbs]][[3]]))
  }
  
  nb_obs_tot = dim_obs*nb_obs_tot;
  param_sub_ini = set_theta_sub
  
  list_arguments_refined = list()
  for (nbs1  in 1:nb_subject){
    seq_argument_nbs1_refined = list(Times_integ_set[[nbs1]],pseudo_Y_set[[nbs1]],weight_integ_set[[nbs1]],param_sub_ini[[nbs1]], diag(exp(log_Delta_vect_ini),length(log_Delta_vect_ini),length(log_Delta_vect_ini)),
                                     mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known[[nbs1]],2,est_pop_parameter_only,param_pop_ini,param_exogen_set[[nbs1]],type_optim_inner)
    list_arguments_refined = append(list_arguments_refined,list(seq_argument_nbs1_refined))
  }
  
  
  if (length(cl)==0){
    out_parLapply_refined =  lapply(list_arguments_refined,est_param_oca_inner_criteria_lincase_par_version)
  }else{
    out_parLapply_refined =  parLapply(cl,list_arguments_refined,est_param_oca_inner_criteria_lincase_par_version) 
  }
  
  param_sub_set_refined = list()
  
  for (nbs2  in 1:nb_subject){
    #   print( out_parLapply[[nbs2]]$val_h_s)
    theta_sub_s = out_parLapply_refined[[nbs2]]$bi_est
    param_sub_set_refined  = append( param_sub_set_refined ,list(theta_sub_s))
  }
  
  
  
  func_eco_E <- function(param_pop,Delta_mat){
    param_sub_set = list()
    val_h =  matrix(0, nrow = nb_subject, ncol =1)
    list_ui_opt = list()
    list_Xi_opt = list()
  
   
    list_arguments = list()
    for (nbs1  in 1:nb_subject){
      seq_argument_nbs1 = list(Times_integ_set[[nbs1]],pseudo_Y_set[[nbs1]],weight_integ_set[[nbs1]],param_sub_set_refined[[nbs1]], Delta_mat,
                               mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known[[nbs1]],2,est_pop_parameter_only,param_pop,param_exogen_set[[nbs1]],type_optim_inner)
      list_arguments = append(list_arguments,list(seq_argument_nbs1))
    }
    
    if (length(cl)==0){
      out_parLapply =  lapply(list_arguments,est_param_oca_inner_criteria_lincase_par_version)
    }else{
      out_parLapply =  parLapply(cl,list_arguments,est_param_oca_inner_criteria_lincase_par_version) 
    }
    
    
    for (nbs2  in 1:nb_subject){
      theta_sub_s = out_parLapply[[nbs2]]$bi_est
      val_h[[nbs2]] =  out_parLapply[[nbs2]]$val_hi
      param_sub_set = append(param_sub_set,list(theta_sub_s))
      list_ui_opt= append(list_ui_opt,list(out_parLapply[[nbs2]]$ui_opt))
      list_Xi_opt = append(list_Xi_opt,list(out_parLapply[[nbs2]]$Xi_opt))
    }
    

    #param_sub_ini = param_sub_set 
    sigma_square = sum(val_h)/(nb_obs_tot+nb_bm*nb_subject)
    
    fval = -(nb_subject*log(det(Delta_mat)) - (nb_obs_tot+nb_bm*nb_subject)*log(sigma_square)/2 +log_prior(param_pop,Delta_mat))
    
    
    return(list(fval=fval,sigma_square =sigma_square,param_sub_est =  param_sub_set,list_ui_opt=list_ui_opt,list_Xi_opt=list_Xi_opt))
  }
  
  func_eco <- function(param_c){
    
    param_pop_cur =  param_c[seq(1,nb_param_pop,1)];
    if  (delta_known == 1){
      Delta_mat_cur = Delta_mat_ini;
    }
    else{
      Delta_vect = param_c[seq(nb_param_pop+1,nb_param_pop+nb_param_cov,1)];
      Delta_mat_cur = diag(x= exp(Delta_vect),nrow = nb_bm,ncol = nb_bm)
      
    }
    
    out_eco_E = func_eco_E(param_pop_cur,Delta_mat_cur)
    return(out_eco_E$fval)
  }
  
  
  func_eco_tot <- function(param_c){
    param_pop_cur =  param_c[seq(1,nb_param_pop,1)];
    if  (delta_known == 1){
      Delta_mat_cur = Delta_mat_ini
    }
    else{
      Delta_vect = param_c[seq(nb_param_pop+1,nb_param_pop+nb_param_cov,1)]
      Delta_mat_cur = diag(x= exp(Delta_vect),nrow = nb_bm,ncol = nb_bm)
    }
    
    out_eco_E = func_eco_E(param_pop_cur,Delta_mat_cur)
    return(list(param_pop_cur,Delta_mat_cur,out_eco_E))
  }
  
  
  if (type_optim ==0){
    res_estim = optimr(param_ext_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control = list(trace=1,maxit = nb_iter_max[1]))  
  }else{
    res_estim = optimr(param_ext_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control = list(trace=1,maxit = nb_iter_max[1]),method='nlminb')  
  }
  
   out_ensemble_param = func_eco_tot(res_estim$par)
  #if (close_cl_at_the_end==1){
  #  stopCluster(cl)
  #}
 
   # print(res_estim)
   # print( out_ensemble_param)
   # print( out_ensemble_param[[2]])
  theta_est = out_ensemble_param[[1]]
  Delta_mat_est = out_ensemble_param[[2]]
  delta_est = log(diag(Delta_mat_est))
  
  list_bi_est = out_ensemble_param[[3]]$param_sub_est
  sigma_square_est = out_ensemble_param[[3]]$sigma_square
  list_ui_opt = out_ensemble_param[[3]]$list_ui_opt
  list_Xi_opt = out_ensemble_param[[3]]$list_Xi_opt
  
  erreur_prediction_estimation_mat_U_cur = erreur_prediction_estimation_oca_lincase(Times_integ_set,pseudo_Y_set,weight_integ_set,list_bi_est,param_exogen_set, 
                                                                            Delta_mat_est,  theta_est, mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,list_x0_known)
 

    if (est_var==1){
    
      out_var = est_param_oca_variance_estimation_lincase_prof_par_version(Times_integ_set,pseudo_Y_set,param_exogen_set,weight_integ_set,theta_est, delta_est,list_bi_est,
                                                                   mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,nb_obs_tot,list_x0_known,est_pop_parameter_only,delta_known,method_est_jacob = 1)
    }else{
      
      out_var = list()
    
  }
  
  return(list(theta_est=theta_est,delta_est=delta_est,sigma_square_est=sigma_square_est, list_bi_est=list_bi_est,
              Times_integ_set=Times_integ_set,list_ui_opt=list_ui_opt,list_Xi_opt=list_Xi_opt,
              Variance_component=out_var,erreur_cross_validation =erreur_prediction_estimation_mat_U_cur))
}
