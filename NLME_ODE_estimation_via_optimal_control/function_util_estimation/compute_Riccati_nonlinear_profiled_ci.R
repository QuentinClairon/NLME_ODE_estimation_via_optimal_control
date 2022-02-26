compute_Riccati_nonlinear_profiled_ci <- function(Times_integ,pseudo_Y,State_i_m1,param_struct,mat_U,mat_A,vect_r,mat_B,mat_C,weight_integ_snd,x0_known,param_pop_cur,param_exo_cur){
  ############################################################
  # DESCRIPTION: compute the solution of the discrete Riccati equation at iteration l for ODE based mixed-effect models
  # 
  # FUNCTION SIGNATURE: compute_Riccati_nonlinear_profiled_ci2 <- function(Times_integ,pseudo_Y,State_i_m1,param_struct,mat_U,mat_A,vect_r,mat_B,mat_C,weight_integ_snd,x0_known,param_pop_cur,param_exo_cur)
  # 
  # INPUT :
  # Times_integ               :(1 x nb_time_integ_points matrix) integration time points 
  # pseudo_Y                  :(dim_obs x nb_time_integ_points sized matrix) extended observations matrix  
  # State_i_m1                :(dim_syst x nb_time_integ_points sized matrix) value of the optimal control solution at previous iteration
  # param_struct              :(vector) current value for the subject specific parameter
  # mat_U                     :(dim_control square matrix) weighing matrix for control magnitude penalization
  # mat_A                     :(function) matrix valued function representing the matrix A in the pseudolinear formulation dX/dt = A(t,x,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              mat_A<- function(t,vx,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : vx (vector) current state value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  # vect_r                    :(function) matrix valued function representing the Vector r in the  pseudolinear formulation dX/dt = A(t,x,bi,theta)X +r(t,bi,theta) of signature: 
  #                                              vect_r <- function(t,param_sub,param_pop,exo_par) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param_sub (vector) subject specific parameter value
  #                                                     : param_pop (vector) population parameter value
  #                                                     : exo_par (list) list of fixed parameter or covariates
  #                                                OUTPUT: (dim_syst x 1  matrix) the value of  r(t,bi,theta)
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # weight_integ_snd          :(1 x nb_time_integ_points matrix)  values of weights sequence 
  # x0_known                  :(vector) subject specific known initial condition
  # param_pop_cur             :(vector) value of mean population parameters
  # param_exo_cur             :(list) list of fixed parameter and covariates
  
  # OUTPUT:
  # Return a list  composed of:
  # out_ensemble_param: a list composed of: 
  #                Val_cost_prof      : (real) current profiled cost function
  #                control_opt        : (dim_control x nb_time_integ-1) estimated optimal control at iteration l
  #                State_opt          : (dim_syst x nb_time_integ) estimated optimal trajectory at iteration l
  ############################################################  
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  nb_time_integ = length(Times_integ)
  Yn_fin = pseudo_Y[,nb_time_integ];
  nb_unknown = dim_syst - length(x0_known)
  
  mat_Rn  = t(mat_C)%*%mat_C
  vect_Hn = -t(mat_C)%*%Yn_fin 

  
  History_R = list()
  History_H = list()
  History_R = list(mat_Rn) 
  History_H = list(vect_Hn) 
  History_G = list()
  mat_Rk  = mat_Rn
  vect_Hk = vect_Hn
  val_Sn_prof = t(Yn_fin)%*%Yn_fin
  
  for (nk in 1:(nb_time_integ-1)){
    delta_k = Times_integ[nb_time_integ-nk+1] - Times_integ[nb_time_integ-nk] 
    Yk= pseudo_Y[,nb_time_integ-nk]
    
    State_pred_cur = State_i_m1[,nb_time_integ-nk]
    mat_Ak_val = mat_A(Times_integ[nb_time_integ-nk],State_pred_cur ,param_struct,param_pop_cur,param_exo_cur)
    mat_Ak =  mat_Ak_val
   
    vect_rk = vect_r(Times_integ[nb_time_integ-nk],param_struct,param_pop_cur,param_exo_cur)
    mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
    weight_cur = weight_integ_snd[nb_time_integ-nk]
    
 
    mat_G = solve(mat_U+delta_k*t(mat_B)%*%mat_Rk%*%mat_B)
    History_G = append(list(mat_G),History_G)
    
   
   
    val_Sn_prof1 = weight_cur*t(Yk)%*%Yk +t(2*vect_Hk+delta_k*mat_Rk%*%vect_rk)%*%vect_rk
    val_Sn_prof2 = t(vect_Hk+delta_k*mat_Rk%*%vect_rk)%*%mat_B%*%mat_G%*%t(mat_B)%*%(vect_Hk+delta_k*mat_Rk%*%vect_rk)
    
    val_Sn_prof = val_Sn_prof +delta_k*(val_Sn_prof1-val_Sn_prof2)    
    
    
    
    mat_Rk_m1_1part = mat_Rk + delta_k*weight_cur*t(mat_C)%*%mat_C + delta_k*(mat_Rk%*%mat_Ak+t(mat_Ak)%*%mat_Rk)+(delta_k^2)*t(mat_Ak)%*%mat_Rk%*%mat_Ak
    mat_Rk_m1_2part = -delta_k*t(mat_Ak_delta)%*%mat_Rk%*%mat_B%*%mat_G%*%t(mat_B)%*%mat_Rk%*%mat_Ak_delta
    mat_Rk_m1 = mat_Rk_m1_1part+mat_Rk_m1_2part
    
   
    
    vect_Hk_m1 = vect_Hk - delta_k*weight_cur*t(mat_C)%*%Yk+delta_k*t(mat_Ak)%*%vect_Hk+  delta_k*t(mat_Ak_delta)%*%mat_Rk%*%vect_rk
    vect_Hk_m1 = vect_Hk_m1 - delta_k*t(mat_Ak_delta)%*%mat_Rk%*%mat_B%*%mat_G%*%t(mat_B)%*%(vect_Hk +delta_k*mat_Rk%*%vect_rk)
    
    mat_Rk = mat_Rk_m1
    vect_Hk = vect_Hk_m1
    

    History_R = append(list(mat_Rk),History_R)
    History_H = append(list(vect_Hk),History_H)
   
  }
  
  mat_R0  =  History_R[[1]]
  vect_H0 = History_H[[1]]
 
 # print(mat_R0)
  if (nb_unknown== dim_syst){
    
    inv_mat_R0 = matrix(10^8,dim_syst,dim_syst)
    out_try_catch <-tryCatch({
      inv_mat_R0 = solve(mat_R0)
    },error=function(cond){
      print("erreur inversion R0")
      return(1)}
    )
    
    est_x0 = -inv_mat_R0%*%vect_H0
  }else{
    if (nb_unknown>0){
    
    mat_R0_11 = mat_R0[seq(1,nb_unknown),seq(1,nb_unknown)]
    mat_R0_12 = mat_R0[seq(1,nb_unknown),seq(nb_unknown+1,dim_syst)]
    #print(nb_unknown)
    
    
    if (length(x0_known) ==1){
      est_x0_unknown = -solve(mat_R0_11)%*%(mat_R0_12*x0_known+vect_H0[seq(1,nb_unknown)])
    }else{
      est_x0_unknown = -solve(mat_R0_11)%*%(mat_R0_12%*%x0_known+vect_H0[seq(1,nb_unknown)])
    }
    est_x0 = c(t(est_x0_unknown),t(x0_known))
    est_x0 = matrix(est_x0,length(est_x0),1)
    }
    else{
      est_x0 = x0_known
    }
  }
  
  
  val_Sn_prof  = val_Sn_prof + t(est_x0)%*%mat_R0%*%est_x0+2*t(vect_H0)%*%est_x0 
  
  State_i = matrix(data = rep(0,dim_syst*nb_time_integ ), nrow = dim_syst, ncol = nb_time_integ )
  State_i[,1] =  est_x0
  
  control_i =matrix(data = rep(0,dim_control*(nb_time_integ-1)), nrow = dim_control, ncol = nb_time_integ-1 )
  val_RSS = 0
  for (nkk in 1:(nb_time_integ-1)){
    delta_k = Times_integ[nkk+1] - Times_integ[nkk] 
    Ykk= pseudo_Y[,nkk];
    
    State_pred_cur = State_i_m1[,nkk]
    mat_Ak_val = mat_A(Times_integ[nkk],State_pred_cur ,param_struct,param_pop_cur,param_exo_cur)
    mat_Ak =  mat_Ak_val
  
    mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
    vect_rk = vect_r(Times_integ[nkk],param_struct,param_pop_cur,param_exo_cur)
    
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
  indice_obs  =weight_integ_snd>0.2      
  State_i_obs_time = State_i[,indice_obs]
  
  ensemble_est = list(Val_cost_prof = val_Sn_prof ,control_opt =control_i, Ricatti = History_R,Adjoint =  History_H,State_opt =State_i,Sequence_G = History_G,State_obs = State_i_obs_time,val_RSS=val_RSS)
  
  return(ensemble_est)
}
