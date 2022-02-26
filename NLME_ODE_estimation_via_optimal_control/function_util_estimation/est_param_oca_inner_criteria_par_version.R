est_param_oca_inner_criteria_par_version <- function(arg_sequence){
  ############################################################
  # DESCRIPTION: estimation of the subject specific parameters iva the inner criteria optimization as well as the solution of the related optimal control problem  for ODE based mixed-effect models
  # 
  # FUNCTION SIGNATURE: est_param_oca_inner_criteria_par_version <- function(arg_sequence) where seq_argument is a list composed of
  # 
  # INPUT :
  # Times_integ_set           :(1 x nb_time_integ_points matrix) integration time points 
  # pseudo_Y_set              :(dim_obs x nb_time_integ_points sized matrix) extended observations matrix 
  # pseudo_State_ini          :(dim_syst x nb_time_integ_points sized matrix) initial state value at integration time points for the optimal control problem
  # weight_integ              :(1 x nb_time_integ_points matrix)  values of weights sequence 
  # bi_ini                    :(vector) initial guess for the subject specific parameter estimator
  # Delta_mat                 :(square matrix ) current value for Delta
  # mat_U                     :(dim_control square matrix) weighing matrix for control magnitude penalization
  # mat_A                     :(function) matrix valued function representing the matrix A in the pseudolinear formulation dX/dt = A(t,x,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              mat_A<- function(t,vx,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : vx (vector) current state value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                 OUTPUT: (dim_syst x dim_syst  matrix) the value of  A(t,x,bi,theta)
  # vect_r                    :(function) matrix valued function representing the Vector r in the  pseudolinear formulation dX/dt = A(t,x,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              vect_r <- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                OUTPUT: (dim_syst x 1  matrix) the value of  r(t,bi,theta)
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # x0_known                  :(vector) subject specific known initial condition
  # inner_criteria_type       :(integer) precize the definition of the inner criteria being use (now just 2 for Sn)
  # opt_traj_est_only         :(integer) equal 0 if the subject specific parameter have to be estimated, 1 if it is just the optimal control 
  # param_pop_cur             :(vector) current value for the mean population parameter
  # param_exo_cur             :(list) list of fixed parameter and covariates
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
  
  pseudo_Times_obs=arg_sequence[[1]] 
  pseudo_Yn=arg_sequence[[2]] 
  pseudo_State_cur_sub_ini = arg_sequence[[3]]
  weight_integ  = arg_sequence[[4]]
  bm_ini= arg_sequence[[5]]
  Delta_mat= arg_sequence[[6]]
  mat_U= arg_sequence[[7]]
  mat_A= arg_sequence[[8]] 
  vect_r= arg_sequence[[9]] 
  mat_B= arg_sequence[[10]] 
  mat_C =arg_sequence[[11]] 
  x0_known =arg_sequence[[12]]
  inner_criteria_type=arg_sequence[[13]]
  opt_traj_est_only=arg_sequence[[14]]
  param_pop_cur=arg_sequence[[15]]
  param_exo_cur = arg_sequence[[16]]
  
  type_algorithm_optim = 0
  if (length(arg_sequence) > 16){
    type_algorithm_optim = arg_sequence[[17]]
  }
  # Sortie [param_sub_est, fval_opt,val_h,control_final,State_final,Riccati_final] 
  
  nb_iter_max_nlminb =200
  trace_val = 0
  if (length(arg_sequence) > 17){
    nb_iter_max_nlminb = arg_sequence[[18]][1]
    trace_val =arg_sequence[[18]][2]
  }
 
 
  #[fval,val_h_cur, control_i,State_i,Riccati_i]  
  func_eco_E <-function(param_sub_cur,pseudo_state_0){
    val_err_iter_cur = 10^100;
    stop_criteria = 0;
    nb_iter_max = 50;
    nb_iter = 0;
   
    new_State_ini = pseudo_state_0
    
    nb_iter = 0
    previous_Riccati_inner = list()
    out_Riccati_inner = list()
  
    while (stop_criteria ==0){
      previous_Riccati_inner = out_Riccati_inner
      out_Riccati_inner =  compute_Riccati_nonlinear_profiled_ci(pseudo_Times_obs,pseudo_Yn,  new_State_ini,param_sub_cur,mat_U,mat_A,vect_r,mat_B,mat_C,weight_integ,x0_known,param_pop_cur,param_exo_cur)  
      val_err_iter = Mod(sum((out_Riccati_inner$State_opt -   new_State_ini)^2))
     
   
      #print(out_Riccati_inner$State_opt)
      nb_iter = nb_iter+1
     
      # Riccati_seq= append(Riccati_seq,list(out_Riccati_inner$Ricatti))
      # State_seq  = append(State_seq,list(out_Riccati_inner$State_opt))
      # Adjoint_seq  = append(Adjoint_seq,list(out_Riccati_inner$Adjoint))
      # Mat_G_seq = append(Mat_G_seq,list(out_Riccati_inner$Sequence_G ))

      new_State_ini  = out_Riccati_inner$State_opt
      
      #   print(nb_iter)
      #     print(val_err_iter) 
   
      #  print(max(abs(out_Riccati_inner$State_opt)))
      #  print(is.nan(sum(out_Riccati_inner$State_opt)))
      # print(is.nan(sum(out_Riccati_inner$State_opt)))
      
      if (is.nan(sum(out_Riccati_inner$State_opt))){
        stop_criteria =2 
      }else{
        max_ric_val =max(abs(out_Riccati_inner$State_opt))
        #  print(max(abs(out_Riccati_inner$State_opt)))
        if (is.na( max_ric_val)){
          stop_criteria =2
        }else{
          if (max_ric_val >10^3){
            stop_criteria =2
          }
        }
     
       
      }
      
      if (is.na(val_err_iter)){
        stop_criteria =2
        }else{
          if (is.nan(val_err_iter)){
          }else{
            if (val_err_iter <10^-1){
              stop_criteria = 1;
            }
          }
        
        }
      
      if (nb_iter  > 40){
        stop_criteria = 1;
      }
    }
    
    if (  stop_criteria==2){
      Riccati_inner_sel =  previous_Riccati_inner
      val_h_cur = Riccati_inner_sel$val_RSS + sum((Delta_mat%*%param_sub_cur)^2) 
      if (nb_iter<2){
        val_h_cur = sum(pseudo_Yn^2)
      }
    }else{
      Riccati_inner_sel =  out_Riccati_inner
      val_h_cur = Riccati_inner_sel$val_RSS + sum((Delta_mat%*%param_sub_cur)^2) 
    }
    
    State_opt =   Riccati_inner_sel$State_opt
   
    if (inner_criteria_type == 1){
      fval = val_h_cur 
    }
    else{
      val_Sn_prof= Riccati_inner_sel$Val_cost_prof
      fval = val_Sn_prof + sum((Delta_mat%*%param_sub_cur)^2)
    }
    
  #  out_eco_E = c(out_Riccati_inner,val_h_s =val_h_cur,fval= fval, State_opt= State_opt,Ensemble_fun = list(Riccati_seq=Riccati_seq,Adjoint_seq =Adjoint_seq,State_seq=State_seq, Mat_G_seq= Mat_G_seq),nb_iter_req = nb_iter )
    out_eco_E = c( Riccati_inner_sel,val_hi =val_h_cur,fval= fval,nb_iter_req = nb_iter )
    
    return(out_eco_E )                          
  }
  
  
  func_eco <- function (param_c){
    fval_E = func_eco_E(param_c,pseudo_State_cur_sub_ini)
    return(fval_E$fval) 
  }
  
  if (opt_traj_est_only == 1){
    out_ensemble_param = func_eco_E(bm_ini, pseudo_State_cur_sub_ini)
    bi_est = bm_ini
    val_gi = out_ensemble_param$fval
  }else {
  
    if (type_algorithm_optim==0){
      res_estim = optimr(bm_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control =  list(trace=trace_val,maxit  =nb_iter_max_nlminb))
      bi_est =res_estim$par 
      val_gi = res_estim$value
    }else{
     res_estim = optimr(bm_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control =  list(trace=trace_val,maxit  =nb_iter_max_nlminb),method="nlminb")
     bi_est =res_estim$par 
     val_gi = res_estim$value
    }
    
    if (length( bi_est)==1){
      bi_est = matrix( bi_est,1,1)
    }
    out_ensemble_param = func_eco_E(bi_est,pseudo_State_cur_sub_ini)   
  }
  
  val_hi =  out_ensemble_param$val_hi 
  ui_opt = out_ensemble_param$control_opt
  Xi_opt = out_ensemble_param$State_opt
  Xi_opt_obs = out_ensemble_param$State_obs
  out_func_inner = list(bi_est=bi_est,ui_opt=ui_opt,Xi_opt=Xi_opt,Xi_opt_obs=Xi_opt_obs,val_gi=val_gi, val_hi= val_hi)
  
  return(out_func_inner)
}
